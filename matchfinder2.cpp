#include "Tree.hpp"
#include "Tree.cpp"
#include "TreeGenerator.hpp"
#include "TreeGenerator.cpp"

#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <unordered_map>

using namespace std;

int q_split = 509;
//int q_split = 7;
const int C_split[] = {0, 376, 438, 229}; 
//const int C_split[] = {0, 3, 5, 6};

int q_lists = 499999993;
const int C_lists[][4] = {{0, 226438576, 285588390, 258698732},
					  {0, 200528554, 240711014, 356710614}, 
					  {0, 259178286, 56395895, 241554533},
					  {0, 65755988, 72153473, 367934924},
					  {0, 261348789, 274137929, 206977495},
					  {0, 86602449, 301910846, 14804088}};
					  
const uint64_t NUM_TREES[] = {1,1,1,1,2,3,6,11,23,47,106,235,551,1301,3159,7741,19320,48629,123867,317955,823065,2144505,
 5623756,14828074,39299897,104636890,279793450,751065460,2023443032,5469566585,14830871802,
 40330829030,109972410221,300628862480, 823779631721,2262366343746,6226306037178};

struct child_lists_t {
	int nchildren[32];
	int children[32][32];
};

void initChildLists(int numVertices, int *parentList, child_lists_t &chlist) {
	memset(chlist.nchildren, 0, numVertices*4);
	// 0 is the root vertex
	for (int i=1; i<numVertices; i++)
		chlist.children[parentList[i]][chlist.nchildren[parentList[i]]++] = i;
}

void unpackTree(int numVertices, uint64_t packedTree, int *parentList) {
	int rootStack[32];
	int stackSize = 0, index = numVertices;
	while (packedTree > 0) {
		if (packedTree & 1)
			parentList[rootStack[--stackSize]] = index;
		else
			rootStack[stackSize++] = index--;
		packedTree >>= 1;
	}
}

uint64_t packTreeFrom(child_lists_t &chlist, int vtx) {
	uint64_t packedSubtrees[32], temp;
	int nkids = chlist.nchildren[vtx];
	uint64_t p = (((1ULL<<nkids)-1)<<1);
	for (int i=0; i<nkids; ++i)
		packedSubtrees[i] = packTreeFrom(chlist, chlist.children[vtx][i]);
	sort(packedSubtrees, packedSubtrees+nkids);
	int length = 1;
	for (int i=0; i<nkids; ++i) {
		while ((1ULL<<length) <= packedSubtrees[i])
			++length;
		p <<= length;
		p |= packedSubtrees[i];
	}
	return p;
}

uint64_t packTree(int numVertices, int *parentList) {
	child_lists_t chlist;
	initChildLists(numVertices, parentList, chlist);
	return packTreeFrom(chlist, 0);
}

void csfSearchPartial(child_lists_t &chlist, int v, int q, const int *C, int64_t *r) {
	int64_t s[4], t[4];
	r[1] = 1; r[2] = 0; r[3] = 0;
	for (int i=0; i<chlist.nchildren[v]; i++) {
		memcpy(t+1, r+1, 24);
		csfSearchPartial(chlist, chlist.children[v][i], q, C, s);
		s[0] = (s[1]*C[1] + s[2]*C[2] + s[3]*C[3]) % q;
		for (int j=1; j<4; ++j) {
			s[j] = q-s[j];
			r[j] = 0;
		}
		for (int j=1; j<4; ++j) {
			for (int p=0; p<4-j; ++p) {
				r[j+p] += t[j]*s[p];
			}
		}
		for (int j=1; j<4; j++) {
			r[j] %= q;
		}
	}
}

// returns phi_{q,C}(T), where T is the tree represented by parentList
int ecsf(int numVertices, int *parentList, int q, const int *C) {
	child_lists_t chlist;
	initChildLists(numVertices, parentList, chlist);
	int64_t r[4];
	csfSearchPartial(chlist, 0, q, C, r);
	return (r[1]*C[1] + r[2]*C[2] + r[3]*C[3]) % q;
}

// returns true if all of the trees in the vector definitely have different CSFs
bool recursiveConjCheck(int numVertices, vector<uint64_t> &trees, int depth) {
	if (depth == 5)
		return false;
	unordered_map<int, vector<uint64_t> > matchingTrees;
	int q = q_lists, csfVal;
	const int *C = C_lists[depth];
	int parentList[32];
	for (int i=0; i<trees.size(); i++) {
		unpackTree(numVertices, trees[i], parentList);
		csfVal = ecsf(numVertices, parentList, q, C);
		matchingTrees[csfVal].push_back(trees[i]);		
	}
	for (auto it = matchingTrees.begin(); it != matchingTrees.end(); ++it) {
		if ((it->second).size() > 1 && !recursiveConjCheck(numVertices, it->second, depth+1)) {
			return false;
		}
	}
	return true;
}

// returns true if all of the trees contained in the tree file have distinct CSFs
bool checkTreeFile(int numVertices, const char *fileName) {
	ifstream fin(fileName, ios::in | ios::binary);
	vector<uint64_t> trees;
	uint64_t treebuf[256];
	while (fin.read((char*)treebuf, 2048)) {
		for (int i=0; i<256; ++i) {
			trees.push_back(treebuf[i]);
		}
	}
	for (int i=0; i<fin.gcount()/8; ++i) {
		trees.push_back(treebuf[i]);
	}
	fin.close();
	return recursiveConjCheck(numVertices, trees, 0); 
}

// need q <= 512
void splitToFiles(int numVertices, int q, const int *C, const char *filePrefix) {
	ofstream outFiles[512];
	char nameBuffer[32];
	for (int i=0; i<strlen(filePrefix); ++i) {
		nameBuffer[i] = filePrefix[i];
	}
	for (int i=0; i<q; ++i) {
		sprintf(nameBuffer+strlen(filePrefix), "%d.bin", i);
		outFiles[i].open(nameBuffer, ios::out | ios::binary);
	}
	TreeGenerator gen(numVertices);
	int parentList[32], csfVal;
	uint64_t packedTree;
	uint64_t treeCtr = 0;
	uint64_t prog = 5;
	while (gen.nextTreeParList(parentList)) {
		++treeCtr;
		if (treeCtr > prog*NUM_TREES[numVertices]/100) {
			cout << "trees " << prog << "% categorized\n";
			prog += 5;
		}
		csfVal = ecsf(numVertices, parentList, q, C);
		packedTree = packTree(numVertices, parentList);
		outFiles[csfVal].write((char*)(&packedTree), 8);
	}
	for (int i=0; i<q; ++i) {
		outFiles[i].close();
	}
}

bool checkConjecture(int n) {
	splitToFiles(n, q_split, C_split, "group");
	char filename[32];
	for (int i=0; i<q_split; ++i) {
		sprintf(filename, "group%d.bin", i);
		if (!checkTreeFile(n, filename)) {
			return false;
		}
		if (i%20 == 19 || i == q_split-1) {
			cout << i+1 << " categories checked\n";
		}
	}
	return true;
}

string adjListsString(int *parentList, int nVertices) {
	vector<int> adjlists[32];
	char arr[4];
	for (int i=1; i<nVertices; ++i) {
		adjlists[parentList[i]].push_back(i);
		adjlists[i].push_back(parentList[i]);
	}
	for (int i=0; i<nVertices; ++i)
		sort(adjlists[i].begin(), adjlists[i].end());
	string s = "[";
	for (int i=0; i<nVertices; ++i) {
		s += "[";
		for (int j=0; j<adjlists[i].size(); ++j) {
			sprintf(arr, "%d", adjlists[i][j]);
			s += string(arr);
			if (j < adjlists[i].size()-1)
				s += ",";
			else
				s += "]";
		}
		if (i < nVertices-1)
			s += ", ";
		else
			s += "]";
	}
	return s;
}

int main() {
	for (int n = 20; n < 30; ++n) {
		cout << "Checking conjecture for n = " << n << "...\n";
		clock_t t = clock();
		if (checkConjecture(n)) {
			cout << "The conjecture holds for n = " << n << ".\n";
		} else {
			cout << "Potential counterexample found!\n";
		}
		t = clock()-t;
		cout << "Total time: " << t/1000.0 << " seconds.\n";
		cout << "Trees per second: " << NUM_TREES[n]*1000LL/t << '\n';
	}
	return 0;
}
		




