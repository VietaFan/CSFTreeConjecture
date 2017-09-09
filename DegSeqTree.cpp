#include "DegSeqTree.h"
#include "printers.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <ctime>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cstdint>
using namespace std;

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

template<typename T1, typename T2>
ostream& operator<<(ostream &out, unordered_map<T1, T2> &dict) {
	out << "{";
	int k = dict.size();
	for (auto it=dict.begin(); it!=dict.end(); ++it) {
		out << (it->first) << ": " << (it->second);
		if (k > 1) {
			out << ", ";
		}
		--k;
	}
	out << "}";
	return out;
}

template<typename T>
ostream& operator<<(ostream &out, unordered_set<T> &S) {
	out << "{";
	int k = S.size();
	for (auto it=S.begin(); it!=S.end(); ++it) {
		out << (*it);
		if (k > 1) {
			out << ", ";
		}
		--k;
	}
	out << "}";
	return out;
}


struct child_lists_t {
	int nchildren[32];
	int children[32][32];
};

void changeRootInPlace(vector<int> &tree, int root) {	
	int v = root, top = -1, nodeStack[32];
	while (v != -1) {
		nodeStack[++top] = v;
		v = tree[v];
	}
	for (int i=top; i>0; --i) {
		tree[nodeStack[i]] = nodeStack[i-1];
	}
	tree[root] = -1;
}

void getChildLists2(vector<int> &tree, child_lists_t &chlist, int rootPos) {
	memset(chlist.nchildren, 0, 4*tree.size());
	for (int i=0; i<rootPos; ++i)
		chlist.children[tree[i]][chlist.nchildren[tree[i]]++] = i;
	for (int i=rootPos+1; i<tree.size(); ++i)
		chlist.children[tree[i]][chlist.nchildren[tree[i]]++] = i;
}
uint64_t packTreeFrom2(child_lists_t &chlist, int vtx);
// algorithm from http://courses.csail.mit.edu/6.046/fall01/handouts/ps9sol.pdf
void getCenters(vector<int> &tree, int &cent1, int &cent2) {
	child_lists_t chlist;
	getChildLists2(tree, chlist, 0);
	int distance[32], nodeQueue[32], pos=0, top=0, N=tree.size(), d;
	int maxiDist = -1, farNode, endPoint2;
	nodeQueue[0] = 0;
	distance[0] = 0;
	while (pos < N) {
		d = distance[nodeQueue[pos]];
		if (d > maxiDist) {
			maxiDist = d;
			farNode = nodeQueue[pos];
		}
		for (int i=0; i<chlist.nchildren[nodeQueue[pos]]; ++i) {
			distance[chlist.children[nodeQueue[pos]][i]] = d+1;
			nodeQueue[++top] = chlist.children[nodeQueue[pos]][i];
		}
		++pos;
	}
	changeRootInPlace(tree, farNode);
	getChildLists2(tree, chlist, farNode);
	nodeQueue[0] = farNode;
	distance[farNode] = 0;
	maxiDist = -1;
	pos = 0;
	top = 0;
	while (pos < N) {
		d = distance[nodeQueue[pos]];
		if (d > maxiDist) {
			maxiDist = d;
			endPoint2 = nodeQueue[pos];
		}
		for (int i=0; i<chlist.nchildren[nodeQueue[pos]]; ++i) {
			distance[chlist.children[nodeQueue[pos]][i]] = d+1;
			nodeQueue[++top] = chlist.children[nodeQueue[pos]][i];
		}
		++pos;
	}
	int node = endPoint2;
	for (int i=0; i<maxiDist/2; ++i)
		node = tree[node];
	cent1 = node;
	if (maxiDist%2) {
		cent2 = tree[node];
		if (chlist.nchildren[cent1] > chlist.nchildren[cent2]) {
			cent2 = -1;
		} else if (chlist.nchildren[cent1] < chlist.nchildren[cent2]) {
			cent1 = cent2;
			cent2 = -1;
		}
			
	} else {
		cent2 = -1;
	}
}

uint64_t packTreeFrom2(child_lists_t &chlist, int vtx) {
	uint64_t packedSubtrees[32], temp;
	int nkids = chlist.nchildren[vtx];
	uint64_t p = (((1ULL<<nkids)-1)<<1);
	for (int i=0; i<nkids; ++i)
		packedSubtrees[i] = packTreeFrom2(chlist, chlist.children[vtx][i]);
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

void addTreesTop(vector<vector<vector<int> > > &subtrees, vector<int> &indexVec, vector<int> &preds, int pos, vector<uint64_t > &trees) {
	// we've chosen all of our indices, so we can make the tree
	if (pos == subtrees.size()) {		
		int treepos = 0;
		vector<int> T;
		T.push_back(-1);
		for (int i=0; i<subtrees.size(); ++i) {
			T.push_back(0); // each subtree has a new direct ancestor of the root
			treepos = T.size(); // now add in the subtree relative to it
			for (int par: subtrees[i][indexVec[i]]) {
				T.push_back(par+treepos);
			}
		}
		child_lists_t chlist0, chlist;
		getChildLists2(T, chlist0, 0);
		uint64_t packed0 = packTreeFrom2(chlist0, 0);
		// try all other vertices of maximal degree as roots; if this is isomorphic to another tree one will have a smaller packed result
		int targetSize = chlist0.nchildren[0]-1;
		for (int i=1; i<T.size(); ++i) {
			if (chlist0.nchildren[i] == targetSize) {
				changeRootInPlace(T, i);
				getChildLists2(T, chlist, i);
				if (packTreeFrom2(chlist, i) < packed0)
					return;
			}
		}
		trees.push_back(packed0);
		return;
	}
	// recurse and pick the next index
	int bound;
	if (preds[pos] == -1) {
		bound = subtrees[pos].size();
	} else {
		bound = indexVec[preds[pos]]+1;
	}
	for (int i=0; i<bound; ++i) {
		indexVec[pos] = i;
		addTreesTop(subtrees, indexVec, preds, pos+1, trees);
	}
}

void addTrees(vector<vector<vector<int> > > &subtrees, vector<int> &indexVec, vector<int> &preds, int pos, vector<vector<int> > &trees) {
	// we've chosen all of our indices, so we can make the tree
	if (pos == subtrees.size()) {		
		int treepos = 0;
		// combine the chosen subtrees and create the tree
		vector<int> tree;
		// add in the ith subtree
		for (int i=0; i<subtrees.size(); ++i) {
			vector<int> &subtree = subtrees[i][indexVec[i]];
			tree.push_back(-1); // each subtree has a new direct ancestor of the root
			treepos = tree.size(); // now add in the subtree relative to it
			for (int par: subtree) {
				tree.push_back(par+treepos);
			}
		}
		trees.push_back(tree);
		return;
	}
	// recurse and pick the next index
	int bound;
	if (preds[pos] == -1) {
		bound = subtrees[pos].size();
	} else {
		bound = indexVec[preds[pos]]+1;
	}
	for (int i=0; i<bound; ++i) {
		indexVec[pos] = i;
		addTrees(subtrees, indexVec, preds, pos+1, trees);
	}
}

void getParts(vector<int> &seq, vector<vector<int> > &stack, vector<vector<vector<int> > > &partseqs, int nparts) {
	if (seq.size() == 0) {
		vector<vector<int> > v;
		for (int i=0; i<nparts; ++i)
			v.push_back(vector<int>());
		partseqs.push_back(v);
		return;
	}
	if (nparts == 1) {
		// check to see if placing everything left in the last subtree is ok
		int k = seq.size();
		if (stack.size() > 0) {
			vector<int> &last = stack[stack.size()-1];
			int i = 0;
			while (i < k && seq[i] == last[i]) {
				++i;
			}
			// the last part doesn't fit
			if (i < k && seq[i] > last[i]) {
				return;
			}
		}
		stack.push_back(seq);
		partseqs.push_back(stack);
		stack.pop_back();
		return;
	}
	int k = seq.size();
	vector<int> nextvec;
	for (int i=0; i<k; ++i)
		nextvec.push_back(0);
	if (stack.size() == 0) {
		int pos = 0, oldpos;
		stack.push_back(vector<int>());
		while (pos > -1) {
			if (nextvec[pos] > seq[pos]) {
				nextvec[pos] = 0;
				--pos;
				if (pos > -1) {
					++nextvec[pos];
				}
				continue;
			}
			pos = k-1;
			// process with current choice for the next vector and recurse
			stack[0] = nextvec;
			vector<int> rem_seq;
			for (int i=0; i<k; ++i) {
				rem_seq.push_back(seq[i]-nextvec[i]);
			}
			getParts(rem_seq, stack, partseqs, nparts-1);
			++nextvec[pos];
		}
		stack.pop_back();
		return;
	}
	vector<int> vec = stack[stack.size()-1];
	int depth = stack.size();
	stack.push_back(vector<int>());
	// first s positions match, then we're strictly below the previous one
	int s;
	for (s=0; s<k; ++s) {
		for (int i=0; i<s; ++i) {
			nextvec[i] = vec[i];
		}
		// pick strictly smaller value for sth position
		for (int i=0; i<vec[s] && i<=seq[s]; ++i) {
			nextvec[s] = i;
			// now fill in the last k-s-1 positions
			int pos = s+1, oldpos;
			while (pos > s) {
				if (pos < k && nextvec[pos] > seq[pos]) {
					nextvec[pos] = 0;
					--pos;
					if (pos > s) {
						++nextvec[pos];
					}
					continue;
				}
				pos = k-1;
				// process with current choice for the next vector and recurse
				stack[depth] = nextvec;
				vector<int> rem_seq;
				for (int i=0; i<k; ++i) {
					rem_seq.push_back(seq[i]-nextvec[i]);
				}
				getParts(rem_seq, stack, partseqs, nparts-1);
				++nextvec[pos];
			}
		}
		if (seq[s] < vec[s]) // if we can't get the first s+1 positions to match, leave
			break;
	}
	if (s == k) {
		stack[depth] = vec;
		vector<int> rem_seq;
		for (int i=0; i<k; ++i)
			rem_seq.push_back(seq[i]-vec[i]);
		getParts(rem_seq, stack, partseqs, nparts-1);
	}
	stack.pop_back();
}

uint64_t packDegseq(vector<int> &degvec, vector<int> &degrees) {
	uint64_t ans = 0;
	int k = 0;
	for (int i=0; i<degvec.size(); ++i) {
		ans <<= degvec[i];
		ans |= ((1ULL<<degvec[i])-1);
		ans <<= (degrees[i]-k);
		k = degrees[i];
	}
	return ans;
}



void genSubtrees(vector<int> &degvec, vector<int> &degrees, vector<vector<int> > &trees, unordered_map<uint64_t, vector<vector<int> > > &subtreesDP, int printDepth) {
	uint64_t packed = packDegseq(degvec, degrees);
	if (subtreesDP.count(packed)) {
		trees = subtreesDP[packed];
		return;
	}
	trees.clear();
	if (degvec.empty()) {
		trees.push_back(vector<int>()); // single element corresponding to empty tree
		return;
	}
	vector<vector<vector<int> > > subtrees;
	vector<vector<vector<int> > > parts_lists;
	int dsize = degrees.size(), d2size, numParts;
	for (int i=0; i<dsize; ++i) {
		--degvec[i];
		vector<int> degvec2,degrees2;
		if (degvec[i] == 0) {
			for (int j=0; j<i; ++j) {
				degvec2.push_back(degvec[j]);
				degrees2.push_back(degrees[j]);
			}
			for (int j=i+1; j<dsize; ++j) {				
				degvec2.push_back(degvec[j]);
				degrees2.push_back(degrees[j]);
			}
		} else {
			degvec2 = degvec;
			degrees2 = degrees;
		}
		numParts = degrees[i];
		d2size = degrees2.size();
		parts_lists.clear();	
		vector<vector<int> > v;	
		getParts(degvec2, v, parts_lists, numParts);
		for (auto parts: parts_lists) {
			subtrees.clear();
			int partNum = 1;
			for (auto part: parts) {		
				vector<vector<int> > subtree;
				vector<int> degvec_next, degrees_next;
				for (int j=0; j<d2size; ++j) {
					if (part[j] > 0) {
						degvec_next.push_back(part[j]);
						degrees_next.push_back(degrees2[j]);
					}
				}
				genSubtrees(degvec_next, degrees_next, subtree, subtreesDP, printDepth+1);
				subtrees.push_back(subtree);
				++partNum;
			}
			int start = 0;
			vector<int> indexvec;
			for (int i=0; i<numParts; ++i)
				indexvec.push_back(0);	
			vector<int> preds;
			for (int i=0; i<numParts; ++i) {
				preds.push_back(-1);		
				for (int j=i-1; j>-1; --j) {
					if (parts[i] == parts[j]) {
						preds[i] = j;
						break;
					}
				}
			}
			addTrees(subtrees, indexvec, preds, 0, trees);
		}
		++degvec[i];
	}
}

void genSubtreesTop(vector<int> &degvec, vector<int> &degrees, vector<uint64_t> &trees, unordered_map<uint64_t, vector<vector<int> > > &subtreesDP) {
	vector<int> degvec_orig = degvec, degrees_orig = degrees;
	string spaces = "";
	trees.clear();
	vector<vector<vector<int> > > subtrees;
	vector<vector<vector<int> > > parts_lists;
	int root_deg = degrees.back()+1;
	degvec.back()--;
	if (degvec.back() == 0) {
		degvec.pop_back();
		degrees.pop_back();
	}
	int dsize = degrees.size(), numParts;
	vector<vector<int> > v;
	getParts(degvec, v, parts_lists, root_deg);
	cout << "parts_list size = " << parts_lists.size() << '\n';
	for (auto parts: parts_lists) {
		subtrees.clear();
		int partNum = 1;
		for (auto part: parts) {		
			vector<vector<int> > subtree;
			vector<int> degvec_next, degrees_next;
			for (int j=0; j<dsize; ++j) {
				if (part[j] > 0) {
					degvec_next.push_back(part[j]);
					degrees_next.push_back(degrees[j]);
				}
			}
			genSubtrees(degvec_next, degrees_next, subtree, subtreesDP, 0);
			subtrees.push_back(subtree);
			++partNum;
		}
		int start = 0;
		vector<int> indexvec;
		for (int i=0; i<root_deg; ++i)
			indexvec.push_back(0);	
		vector<int> preds;
		for (int i=0; i<root_deg; ++i) {
			preds.push_back(-1);		
			for (int j=i-1; j>-1; --j) {
				if (parts[i] == parts[j]) {
					preds[i] = j;
					break;
				}
			}
		}
		addTreesTop(subtrees, indexvec, preds, 0, trees);
	}
	degvec = degvec_orig;
	degrees = degrees_orig;
}

void generateTrees(vector<int> &degSeq, vector<uint64_t> &trees, unordered_map<uint64_t, vector<vector<int> > > &subtreesDP) {
	map<int, int> degmap;
	int N = 0;
	for (int deg: degSeq) {
		N += deg;
		if (!degmap.count(deg)) {
			degmap[deg] = 1;
		}
		else {
			++degmap[deg];
		}
	}
	vector<int> degvec, degrees;
	for (auto it = degmap.begin(); it != degmap.end(); ++it) {
		degrees.push_back(it->first);
		degvec.push_back(it->second);
	}
	genSubtreesTop(degvec, degrees, trees, subtreesDP);
}

vector<vector<int> > getPartitions(int n, int k) {
	if (n == 0) {
		vector<vector<int> > v;
		v.push_back(vector<int>());
		return v;
	}
	vector<vector<int> > ans;
	vector<int> v;
	if (k > n)
		k = n;
	for (int t=1; t<=k; ++t) {
		for (vector<int> prefix: getPartitions(n-t, t)) {
			prefix.push_back(t);
			ans.push_back(prefix);
		}
	}
	return ans;
}

void prepareSubtreesDP(int n, unordered_map<uint64_t, vector<vector<int> > > &subtreesDP) {
	vector<vector<int> > vec = getPartitions(n-2,n-2);	
	for (vector<int> degSeq: vec) {
		map<int, int> degmap;
		for (int deg: degSeq) {
			if (!degmap.count(deg)) {
				degmap[deg] = 1;
			}
			else {
				++degmap[deg];
			}
		}
		vector<int> degvec, degrees;
		for (auto it = degmap.begin(); it != degmap.end(); ++it) {
			degrees.push_back(it->first);
			degvec.push_back(it->second);
		}
		vector<vector<int> > trees;
		genSubtrees(degvec, degrees, trees, subtreesDP, 0);
		subtreesDP[packDegseq(degvec, degrees)] = trees;
	}
}

void initSubtreesDP(int n, unordered_map<uint64_t, vector<vector<int> > > &subtreesDP) {
	for (int k=3; k<=n; ++k) {
		prepareSubtreesDP(k, subtreesDP);
	}
}

void getDegreeSequences(int numVertices, vector<vector<int> > &degSeqs) {
	degSeqs = getPartitions(numVertices-2, numVertices-2);
}


