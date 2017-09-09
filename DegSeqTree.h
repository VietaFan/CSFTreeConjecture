#ifndef DEG_SEQ_TREE
#define DEG_SEQ_TREE
#include <vector>
#include <unordered_map>
#include <cstdint>

void getDegreeSequences(int numVertices, std::vector<std::vector<int> > &degSeqs);
void generateTrees(std::vector<int> &degSeq, std::vector<uint64_t> &trees, std::unordered_map<uint64_t, std::vector<std::vector<int> > > &subtreesDP);
void initSubtreesDP(int n, std::unordered_map<uint64_t, std::vector<std::vector<int> > > &subtreesDP);

#endif
