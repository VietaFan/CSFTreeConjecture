#include <pthread.h>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include "printers.h"
#include "DegSeqTree.h"

#define MAX_DEPTH_LIMIT 4
#define MODULUS 1000000007
#define CHECK_INTERVAL 1
#define UPDATE_INTERVAL 60
#define CUTOFF_INTERVAL 1000000
#define VERBOSITY_1

#ifdef VERBOSITY_5
#define VERBOSITY_4
#endif
#ifdef VERBOSITY_4
#define VERBOSITY_3
#endif
#ifdef VERBOSITY_3
#define VERBOSITY_2
#endif
#ifdef VERBOSITY_2
#define VERBOSITY_1_5
#endif
#ifdef VERBOSITY_1_5
#define VERBOSITY_1
#endif
#ifdef VERBOSITY_1
#define VERBOSITY_0
#endif
using namespace std;

#define depthLimit 3

const uint64_t NUM_TREES[] = {1,1,1,1,2,3,6,11,23,47,106,235,551,1301,3159,7741,19320,48629,123867,317955,823065,2144505,
 5623756,14828074,39299897,104636890,279793450,751065460,2023443032,5469566585,14830871802,
 40330829030,109972410221,300628862480, 823779631721,2262366343746,6226306037178};

uint64_t treesSearched, degseqsSearched, dsSizeRecord, nextCutoff;
vector<vector<int> > degSeqs;
map<vector<int>, bool> verified;
vector<vector<vector<int> > > counterexamples;
pthread_mutex_t outputMutex, degseqMutex, resultsMutex, ctrexMutex, progressMutex;


struct child_lists_t {
	int nchildren[32];
	int children[32][32];
};


template<typename T>
ostream& operator<<(ostream &out, vector<T> &vec) {
	out << "[";
	if (vec.size() > 0) {		
		for (int i=0; i<vec.size()-1; ++i)
			out << vec[i] << ", ";
		out << vec[vec.size()-1];
	}
	out << "]";
	return out;
}

template<typename T>
ostream& operator<<(ostream &out, vector<T> &&vec) {
	out << "[";
	if (vec.size() > 0) {		
		for (int i=0; i<vec.size()-1; ++i)
			out << vec[i] << ", ";
		out << vec[vec.size()-1];
	}
	out << "]";
	return out;
}

void initChildLists(int numVertices, int *parentList, child_lists_t &chlist) {
	for (int i=0; i<numVertices; i++)
		chlist.nchildren[i] = 0;
	// 0 is the root vertex
	for (int i=1; i<numVertices; i++)
		chlist.children[parentList[i]][chlist.nchildren[parentList[i]]++] = i;
}

int csfSearchPartial(child_lists_t &chlist, const int *p, int v, int64_t *L) {
	int a;
	int64_t Y[MAX_DEPTH_LIMIT], L2[MAX_DEPTH_LIMIT];
	memset(L, 0, ((depthLimit+1)<<3));
	L[0] = 1; a = 1;
	int y, z;
	for (int i=0; i<chlist.nchildren[v]; i++) {
		y = csfSearchPartial(chlist, p, chlist.children[v][i], Y);
		memset(L2, 0, depthLimit<<3);
		for (int i=0; i<a; i++) {
			z = min(y, depthLimit-i);
			for (int j=0; j<z; j++)	
				L2[i+j] += L[i]*Y[j];
		}
		if (a < depthLimit) {
			if (a+y-1 < depthLimit)
				a += y-1;
			else
				a = depthLimit;
		}
		for (int i=0; i<a; i++)
			L[i] = L2[i]%MODULUS;
	}
	for (int i=a; i>0; i--)
		L[i] = -L[i-1];
	int64_t X=0;
	for (int i=1; i<=a; i++) {
		X -= L[i]*p[i];
	}
	X %= MODULUS;
	if (X < 0)
		X += MODULUS;
	L[0] = X;
	return min(a+1, depthLimit);
}

void unpackTree(int n, uint64_t packedTree, int *parentList) {
	int rootStack[32];
	int stackSize = 0, index = n;
	while (packedTree > 0) {
		if (packedTree & 1)
			parentList[rootStack[--stackSize]] = index;
		else
			rootStack[stackSize++] = index--;
		packedTree >>= 1;
	}
}

/* Returns true if all the CSFs of the trees are distinct,
false if there is a possible match*/
bool compareTreeCSFs(int numVertices, vector<uint64_t> &trees, int comparisonLevel, int threadID) {
	if (comparisonLevel == 6)
		return false;
	child_lists_t chlist;		
	int64_t L[32];
	int tree[32];
	const int pvalLists[][4] = {{0, 226438576, 285588390, 258698732},
					  {0, 200528554, 240711014, 356710614}, 
					  {0, 259178286, 56395895, 241554533},
					  {0, 65755988, 72153473, 367934924},
					  {0, 261348789, 274137929, 206977495},
					  {0, 86602449, 301910846, 14804088}};
	vector<int> csfValVec, sortedCsfValVec;
	sortedCsfValVec.reserve(trees.size());
	csfValVec.reserve(trees.size());
	unordered_set<int> collisions;
	#ifdef VERBOSITY_3
		pthread_mutex_lock(&outputMutex);
		cout << "[thread " << threadID << "]: Comparing " << trees.size() << " trees at comparison level " << comparisonLevel << ".\n";
		pthread_mutex_unlock(&outputMutex);
	#else
	#ifdef VERBOSITY_1_5
		if (comparisonLevel > 0) {
			pthread_mutex_lock(&outputMutex);
			cout << "[thread " << threadID << "]: Comparing " << trees.size() << " trees at comparison level " << comparisonLevel << ".\n";
			pthread_mutex_unlock(&outputMutex);
		}
	#endif
	#endif
	/*cout << "comparing trees with comparison level " << comparisonLevel << endl;
	cout << trees.size() << " trees.\n";*/
	for (uint64_t packedTree: trees) {
		unpackTree(numVertices, packedTree, tree);
		initChildLists(numVertices, tree, chlist);
		csfSearchPartial(chlist, pvalLists[comparisonLevel], 0, L);
		#ifdef VERBOSITY_3
			pthread_mutex_lock(&outputMutex);
			vector<int> treeVec;
			for (int i=0; i<numVertices; ++i)
				treeVec.push_back(tree[i]);
			cout << "[thread " << threadID << "]: tree " << treeVec << " ==> ";
			for (int i=0; i<=comparisonLevel; ++i) {
				csfSearchPartial(chlist, pvalLists[i], 0, L);
			 	cout << L[0] << ' ';
			}
			cout << endl;
			pthread_mutex_unlock(&outputMutex);
		#else
		#ifdef VERBOSITY_1_5
			if (comparisonLevel > 0) {
				pthread_mutex_lock(&outputMutex);
				vector<int> treeVec;
				for (int i=0; i<numVertices; ++i)
					treeVec.push_back(tree[i]);
				cout << "[thread " << threadID << "]: tree " << treeVec << " ==> ";
				for (int i=0; i<=comparisonLevel; ++i) {
					csfSearchPartial(chlist, pvalLists[i], 0, L);
					cout << L[0] << ' ';
				}
				cout << endl;
				pthread_mutex_unlock(&outputMutex);
			}
		#endif
		#endif
		csfValVec.push_back(L[0]);
	}
	sortedCsfValVec = csfValVec;
	sort(sortedCsfValVec.begin(), sortedCsfValVec.end());
	for (int i=1; i<trees.size(); ++i) {
		if (sortedCsfValVec[i] == sortedCsfValVec[i-1]) {
			collisions.insert(sortedCsfValVec[i]);
		}
	}
	unordered_map<int, vector<uint64_t> > csfValMap;
	for (int i=0; i<trees.size(); ++i) {
		if (collisions.count(csfValVec[i])) {
			csfValMap[csfValVec[i]].push_back(trees[i]);
		}
	}
	for (int csfVal: collisions) {
		if (!compareTreeCSFs(numVertices, csfValMap[csfVal], comparisonLevel+1, threadID)) {
			return false;
		}
	}
	return true;
}

bool searchTreesDegSeq(vector<int> &degSeq, int threadID, unordered_map<uint64_t, vector<vector<int> > > &subtreesDP) {
	int numVertices = 2;
	for (int i=0; i<degSeq.size(); ++i)
		numVertices += degSeq[i];
	vector<uint64_t> trees;
	/*#ifdef VERBOSITY_2
	pthread_mutex_lock(&outputMutex);
	cout << "[thread " << threadID << "]: Searching degree sequence " << degSeq << ".\n";
	cout << "[thread " << threadID << "]: Generating trees...\n";
	pthread_mutex_unlock(&outputMutex);
	#endif*/
	generateTrees(degSeq, trees, subtreesDP);
	pthread_mutex_lock(&outputMutex);
	if (trees.size() > dsSizeRecord) {
		dsSizeRecord = trees.size();
		cout << dsSizeRecord << " trees: " << degSeq << "\n";
	}
	pthread_mutex_unlock(&outputMutex);
	#ifdef VERBOSITY_2
	pthread_mutex_lock(&outputMutex);
	cout << "[thread " << threadID << "]: " << trees.size() << " trees with degree sequence " << degSeq << " found.\n";
	//cout << "[thread " << threadID << "]: Comparing CSFs of these trees...\n";
	pthread_mutex_unlock(&outputMutex);
	#endif
	bool allDistinct = compareTreeCSFs(numVertices, trees, 0, threadID);
	/*#ifdef VERBOSITY_2
	pthread_mutex_lock(&outputMutex);
	if (allDistinct) {
		cout << "[thread " << threadID << "]: The conjecture holds for trees with degree sequence " << degSeq << ".\n";
	} else {
		cout << "[thread " << threadID << "]: Potential counterexample found with degree sequence " << degSeq << ".\n";
	}
	pthread_mutex_unlock(&outputMutex);
	#endif*/
	#ifdef VERBOSITY_1
	pthread_mutex_lock(&progressMutex);
	treesSearched += trees.size();
	++degseqsSearched;
	if (treesSearched > nextCutoff) {
		while (nextCutoff <= treesSearched)
			nextCutoff += CUTOFF_INTERVAL;
		pthread_mutex_lock(&outputMutex);
		cout << treesSearched << " of " << NUM_TREES[numVertices] << " searched.\n";
		pthread_mutex_unlock(&outputMutex);
	}
	pthread_mutex_unlock(&progressMutex);
	#endif
	return allDistinct;
}

int totalCtr = 0;

void *worker(void *arg) {
	int x = *((int*)arg);
	int id = x%256, nvtx = x>>8;
	#ifdef VERBOSITY_2
	pthread_mutex_lock(&outputMutex);
	cout << "[thread " << id << "]: hi\n";
	pthread_mutex_unlock(&outputMutex);
	#endif
	unordered_map<uint64_t, vector<vector<int> > > subtreesDP;
	if (nvtx > 20)
		initSubtreesDP(18, subtreesDP);
	else
		initSubtreesDP(nvtx-2, subtreesDP);
	vector<vector<int> > sequences;
	vector<bool> results;
	bool isRunning = true;
	while (isRunning) {
		pthread_mutex_lock(&degseqMutex);
		if (degSeqs.empty()) {
			isRunning = false;
		} else {
			sequences.push_back(degSeqs.back());
			degSeqs.pop_back();
		}
		pthread_mutex_unlock(&degseqMutex);
		/*#ifdef VERBOSITY_2
		pthread_mutex_lock(&outputMutex);
		cout << "[thread " << id << "]: Working on sequences " << sequences << ".\n";
		pthread_mutex_unlock(&outputMutex);
		#endif*/
		for (int i=0; i<sequences.size(); ++i) {
			results.push_back(searchTreesDegSeq(sequences[i], id, subtreesDP));
		}
		pthread_mutex_lock(&resultsMutex);
		for (int i=0; i<sequences.size(); ++i) {
			verified[sequences[i]] = results[i];
		}
		pthread_mutex_unlock(&resultsMutex);
		sequences.clear();
		results.clear();
	}
	subtreesDP.~unordered_map();
	pthread_exit(NULL);
}

bool testConjecture(int n, int nthreads) {
	treesSearched = 0;
	degseqsSearched = 0;
	nextCutoff = 1000000;
	degSeqs.clear();
	verified.clear();
	counterexamples.clear();
	#ifdef VERBOSITY_1
	cout << "[main thread]: Testing conjecture for n = " << n << ".\n";
	#endif
	#ifdef VERBOSITY_1_5
	cout << "[main thread]: Generating degree sequences for n = " << n << ".\n";
	#endif
	getDegreeSequences(n, degSeqs);
	int numSequences = degSeqs.size();
	#ifdef VERBOSITY_1_5
	cout << "[main thread]: Degree sequences generated. Creating child threads...\n";
	#endif
	pthread_t* threads = new pthread_t[nthreads];
	int *threadIDs = new int[nthreads];
	pthread_attr_t attributes;
	pthread_mutex_init(&outputMutex, NULL);
	pthread_mutex_init(&degseqMutex, NULL);
	pthread_mutex_init(&resultsMutex, NULL);
	pthread_mutex_init(&ctrexMutex, NULL);
	pthread_mutex_init(&progressMutex, NULL);
	pthread_attr_init(&attributes);
	pthread_attr_setdetachstate(&attributes, PTHREAD_CREATE_JOINABLE);
	
	//clock_t T0 = clock()/1000, checkTime = T0+CHECK_INTERVAL, updateTime = T0+UPDATE_INTERVAL;
	time_t T0;
	time(&T0);
	time_t checkTime = T0+CHECK_INTERVAL, updateTime = T0+UPDATE_INTERVAL;
	for (int i=0; i<nthreads; ++i) {
		threadIDs[i] = (i+1)+(n<<8);
		pthread_create(&threads[i], &attributes, worker, (void*)(threadIDs+i));
	}
	
	#ifdef VERBOSITY_1_5
	pthread_mutex_lock(&outputMutex);
	cout << "[main thread]: All child threads created.\n";
	pthread_mutex_unlock(&outputMutex);
	#endif 
	
	/*
	#ifdef VERBOSITY_1
	bool inProgress = true;
	time_t currentTime;
	while (inProgress) {
		time(&currentTime);
		//currentTime = clock()/1000;
		if (currentTime >= checkTime) {
			pthread_mutex_lock(&progressMutex);
			if (treesSearched == NUM_TREES[n]) {
				inProgress = false;
			}
			pthread_mutex_unlock(&progressMutex);
			checkTime += CHECK_INTERVAL;
		}
		if (currentTime >= updateTime) {
			pthread_mutex_lock(&progressMutex);
			cout << "[main thread]: Time elapsed: " << currentTime-T0 << " s.\n";
			cout << "[main thread]: " << treesSearched << " of " << NUM_TREES[n] << " trees searched.\n";
			cout << "[main thread]: " << degseqsSearched << " of " << numSequences << " degree sequences searched.\n";
			cout << "[main thread]: Current rate: " << treesSearched/(currentTime-T0) << " trees per second.\n";
			pthread_mutex_unlock(&progressMutex);
			updateTime += UPDATE_INTERVAL;
		}
	}
	#endif
	*/
	
	pthread_attr_destroy(&attributes);
	
	void *status;
	for (int i=0; i<nthreads; ++i) {
		pthread_join(threads[i], &status);
	}
	time_t Tf;
	time(&Tf);
	//clock_t Tf = clock()/1000;
	
	pthread_mutex_destroy(&outputMutex);
	pthread_mutex_destroy(&degseqMutex);
	pthread_mutex_destroy(&resultsMutex);
	pthread_mutex_destroy(&ctrexMutex);
	pthread_mutex_destroy(&progressMutex);
	//pthread_exit(NULL);
	delete[] threads;
	delete[] threadIDs;
	getDegreeSequences(n, degSeqs);
	
	bool conjHolds = true;
	for (vector<int> degSeq: degSeqs) {
		if (!verified[degSeq]) {
			cout << "[main thread]: Possible counterexample found for degree sequence " << degSeq << ".\n";
			conjHolds = false;
		}
	}
	for (int i=0; i<counterexamples.size(); ++i) {
		cout << "[main thread]: Potential counterexample found: the CSFs of the following trees matched at the tested points.\n";
		for (vector<int> tree: counterexamples[i]) {
			cout << tree << endl;
		}
	}
	if (conjHolds) {
		cout << "[main thread]: The conjecture holds for n = " << n << ".\n\n";
	}
	
	cout << "[main thread]: " << degseqsSearched << " degree sequences and " << treesSearched << " trees searched.\n";
	cout << "[main thread]: Total computation time: " << Tf-T0 << " s.\n";
	if (Tf > T0) {
		cout << "[main thread]: Average rate: " << NUM_TREES[n]/(Tf-T0) << " trees per second.\n";
	}
	cout << '\n';
	return conjHolds;
}

int main(int argc, char **argv)  {
	int N=10, T=4;
	if (argc >= 2) {
		N = atoi(argv[1]);
	}
	if (argc >= 3) {
		T = atoi(argv[2]);
	}
	if (argc == 4 && strcmp(argv[3], "all") == 0) {
		for (int n=3; n<=N; ++n)
			testConjecture(n, T);
	} else {
		testConjecture(N, T);
	}
	return 0;
}
	
