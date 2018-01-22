/*
Tree.hpp - written by Keeler Russell
https://github.com/keeler/csf/blob/master/src/Tree.hpp
*/
#ifndef TREE_HPP
#define TREE_HPP

#include <vector>
#include <string>

struct Edge
{
    unsigned int u;
    unsigned int v;
};

class Tree
{
    public:
        Tree();
        Tree(unsigned int order, int *levelSequence);
        Tree(const Tree & that);
        ~Tree();

        unsigned int getOrder() const { return mOrder; };
        unsigned int getDegree(unsigned int v) const { if(v < mOrder) return mDegrees[v]; return 0; };

        std::vector<Edge> getEdges() const;
        void getDegreeSequence(std::vector<unsigned int> & degreeSeqnence) const;
        void getPathNums(std::vector<unsigned int> & pathNums) const;

        void print() const;

        Tree & operator=(const Tree & rhs);
    
    private:
        void pathNumsDfs(int i, int count, bool *& visited, std::vector<unsigned int> & pathNums) const;

        unsigned int    mOrder;
        unsigned int    *mDegrees;
        unsigned int    *mOffsets;
        unsigned int    *mAdjArray;
};

#endif
