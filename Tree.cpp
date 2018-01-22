/*
Tree.cpp - written by Keeler Russell
https://github.com/keeler/csf/blob/master/src/Tree.cpp
*/
#include "Tree.hpp"
//#include "Debug.hpp"
#include <vector>
#include <queue>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <string>
#include <iostream>
using namespace std;

Tree::Tree() :
    mOrder(0),
    mDegrees(NULL),
    mOffsets(NULL),
    mAdjArray(NULL)
{
}

Tree::Tree(unsigned int order, int *levelSequence) :
    mOrder(order)
{
    mOffsets = new unsigned int[mOrder];
    mDegrees = new unsigned int[mOrder];
    memset(mOffsets, 0, sizeof(unsigned int) * mOrder);
    memset(mDegrees, 0, sizeof(unsigned int) * mOrder);

    vector<Edge> edges;
    Edge e;
    for(unsigned int i = 2; i <= mOrder; i++)
    {
        e.u = i - 1;
        e.v = levelSequence[i - 1] - 1;
        mDegrees[e.u]++;
        mDegrees[e.v]++;
        edges.push_back(e);
    }

    for(unsigned int i = 1; i < mOrder; i++)
    {
        mOffsets[i] = mOffsets[i - 1] + mDegrees[i - 1];
    }

    mAdjArray = new unsigned int[2 * (mOrder - 1)];
    unsigned int *currSizes = new unsigned int[mOrder];
    memset(currSizes, 0, sizeof(unsigned int) * mOrder);
    for(unsigned int i = 0; i < edges.size(); i++)
    {
        unsigned int u = edges[i].u;
        unsigned int v = edges[i].v;
        mAdjArray[mOffsets[u] + currSizes[u]++] = v;
        mAdjArray[mOffsets[v] + currSizes[v]++] = u;
    }
    delete [] currSizes;
}

Tree::Tree(const Tree & that)
{
    mOrder = that.mOrder;
    mDegrees = new unsigned int[mOrder];
    mOffsets = new unsigned int[mOrder];
    mAdjArray = new unsigned int[2 * (mOrder - 1)];

    memcpy(mDegrees, that.mDegrees, sizeof(unsigned int) * mOrder);
    memcpy(mOffsets, that.mOffsets, sizeof(unsigned int) * mOrder);
    memcpy(mAdjArray, that.mAdjArray, sizeof(unsigned int) * (2 * (mOrder - 1)));
}

Tree::~Tree()
{
    if(mDegrees != NULL)
    {
        delete [] mDegrees;
    }
    if(mOffsets != NULL)
    {
        delete [] mOffsets;
    }
    if(mAdjArray != NULL)
    {
        delete [] mAdjArray;
    }
}

std::vector<Edge> Tree::getEdges() const
{
    std::vector<Edge> edges;
    for(unsigned int i = 0; i < mOrder; i++)
    {
        for(unsigned int j = 0; j < mDegrees[i]; j++)
        {
            Edge e;
            e.u = i;
            e.v = mAdjArray[mOffsets[i] + j];
            if(e.u < e.v)
            {
                edges.push_back(e);
            }
        }
    }

    return edges;
}

void Tree::getDegreeSequence(std::vector<unsigned int> & degreeSequence) const
{
    degreeSequence.clear();
    degreeSequence.insert(degreeSequence.end(), &mDegrees[0], &mDegrees[mOrder]);
    std::sort(degreeSequence.begin(), degreeSequence.end(), std::greater<int>());
}

void Tree::getPathNums(std::vector<unsigned int> & pathNums) const
{
    pathNums.clear();
    for(unsigned int i = 0; i < mOrder; i++)
    {
        pathNums.push_back(0);
    }

    bool * visited = new bool[mOrder];
    for(unsigned int i = 0; i < mOrder; i++)
    {
        memset(visited, 0, sizeof(bool) * mOrder); 
        pathNumsDfs(i, 1, visited, pathNums);
    }
    delete [] visited;

    for(unsigned int i = 0; i < pathNums.size(); i++)
    {
        pathNums[i] /= 2;
    }
}

void Tree::print() const
{
    vector<Edge> edgeSet = getEdges();
    for(unsigned int i = 0; i < edgeSet.size(); i++)
    {
        cout << "(" << edgeSet[i].u << ", " << edgeSet[i].v << ") ";
    }
    cout << endl;
}

Tree & Tree::operator=(const Tree & that)
{
    if(this != &that)
    {
        delete [] mDegrees;
        delete [] mOffsets;
        delete [] mAdjArray;

        mOrder = that.mOrder;

        mDegrees = new unsigned int[mOrder];
        mOffsets = new unsigned int[mOrder];
        mAdjArray = new unsigned int[2 * (mOrder - 1)];

        memcpy(mDegrees, that.mDegrees, sizeof(unsigned int) * mOrder);
        memcpy(mOffsets, that.mOffsets, sizeof(unsigned int) * mOrder);
        memcpy(mAdjArray, that.mAdjArray, sizeof(unsigned int) * (2 * (mOrder - 1)));
    }

    return *this;
}

void Tree::pathNumsDfs(int i, int count, bool *& visited, std::vector<unsigned int> & pathNums) const
{
    visited[i] = true;
    for(unsigned int j = 0; j < mDegrees[i]; j++)
    {
        unsigned int adjVert = mAdjArray[mOffsets[i] + j];
        if(!visited[adjVert])
        {
            pathNums[count]++;
            pathNumsDfs(adjVert, count + 1, visited, pathNums);
        }
    }
}
