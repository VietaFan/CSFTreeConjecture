/*
TreeGenerator.cpp - written by Keeler Russell
https://github.com/keeler/csf/blob/master/src/TreeGenerator.cpp
*/
#include "TreeGenerator.hpp"
#include <cstring>
#include <climits>

#include <iostream>
#include <vector>
using namespace std;

TreeGenerator::TreeGenerator(unsigned int numVertices) :
    mNumVertices(numVertices),
    firstTime(true)
{
    L = new int[mNumVertices];
    currentLevelSequence = new int[numVertices];
}

TreeGenerator::~TreeGenerator()
{
    delete [] L;
    delete [] currentLevelSequence;
}

bool TreeGenerator::nextTree(Tree & t)
{
    if(!firstTime && this->q == 0)
    {
        return false;
    }

    if(firstTime)
    {
        generateFirstLevelSequence();
        firstTime = false;
    }
    else
    {
        generateNextLevelSequence();
    }
		
    Tree tree(mNumVertices, currentLevelSequence );
    
    t = tree;
    return true;
}
 // another function I added to get rid of stuff we don't need
bool TreeGenerator::nextTreeParList(int* parList)
{
    if(!firstTime && this->q == 0)
    {
        return false;
    }

    if(firstTime)
    {
        generateFirstLevelSequence();
        firstTime = false;
    }
    else
    {
        generateNextLevelSequence();
    }
    
    for (int i=1; i<mNumVertices; i++)
    	parList[i] = currentLevelSequence[i]-1;
		/*
    Tree tree(mNumVertices, currentLevelSequence );
    
    t = tree;*/
    return true;
}

void TreeGenerator::generateFirstLevelSequence()
{
    int k = (mNumVertices / 2) + 1;

    if(mNumVertices == 4)
    {
        this->p = 3;
    }
    else
    {
        this->p = mNumVertices;
    }

    this->q = mNumVertices - 1;
    this->h1 = k;
    this->h2 = mNumVertices;

    if(mNumVertices % 2 == 0)
    {
        this->c = mNumVertices + 1;
    }
    else
    {
        this->c = INT_MAX;      // i.e. infinity
    }

    this->r = k;

    for(int i = 1; i <= k; i++)
    {
        L[i - 1] = i;
    }
    for(unsigned int i = k + 1; i <= mNumVertices; i++)
    {
        L[i - 1] = i - k + 1;
    }
    for(unsigned int i = 0; i < mNumVertices; i++)
    {
        currentLevelSequence[i] = i;
    }

    if(mNumVertices > 2)
    {
        currentLevelSequence[k] = 1;
    }
    if(mNumVertices <= 3)
    {
        this->q = 0;
    }
}

void TreeGenerator::generateNextLevelSequence()
{
	  int fixit = 0;

    int needr = 0;
    int needc = 0;
    int needh2 = 0;

    int n = mNumVertices;

    if(c == n + 1 ||
        (p == h2 &&
        ((L[h1 - 1] == L[h2 - 1] + 1 &&
          n - h2 > r - h1) ||
          (L[h1 - 1] == L[h2 - 1] &&
          n - h2 + 1 < r - h1))))
    {
        if(L[r - 1] > 3)
        {
            p = r;
            q = currentLevelSequence[r - 1];
            if(h1 == r)
            {
                h1--;
            }
            fixit = 1;
        }
        else
        {
            p = r;
            r--;
            q = 2;
        }
    }

    if(p <= h1)
    {
        h1 = p - 1;
    }
    if(p <= r)
    {
        needr = 1;
    }
    else if(p <= h2)
    {
        needh2 = 1;
    }
    else if(L[h2 - 1] == L[h1 - 1] - 1 && n - h2 == r - h1)
    {
        if(p <= c)
        {
            needc = 1;
        }
    }
    else
    {
        c = INT_MAX;
    }

    int oldp = p;
    int delta = q - p;
    int oldlq = L[q - 1];
    int oldwq = currentLevelSequence[q - 1];
    p = INT_MAX;

    for(int i = oldp; i <= n; i++)
    {
        L[i - 1] = L[i - 1 + delta];
        if(L[i - 1] == 2)
        {
            currentLevelSequence[i - 1] = 1;
        }
        else
        {
            p = i;
            if(L[i - 1] == oldlq)
            {
                q = oldwq;
            }
            else
            {
                q = currentLevelSequence[i - 1 + delta] - delta;
            }
            currentLevelSequence[i - 1] = q;
        }
        if(needr == 1 && L[i - 1] == 2)
        {
            needr = 0;
            needh2 = 1;
            r = i - 1;
        }
        if(needh2 == 1 && L[i - 1] <= L[i - 2] && i > r + 1)
        {
            needh2 = 0;
            h2 = i - 1;
            if(L[h2 - 1] == L[h1 - 1] - 1 && n - h2 == r - h1)
            {
                needc = 1;
            }
            else
            {
                c = INT_MAX;
            }
        }
        if(needc == 1)
        {
            if(L[i - 1] != L[h1 - h2 + i - 1] - 1)
            {
                needc = 0;
                c = i;
            }
            else
            {
                c = i + 1;
            }
        }
    }

    if(fixit == 1)
    {
        r = n - h1 + 1;
        for(int i = r + 1; i <= n; i++)
        {
            L[i - 1] = i - r + 1;
            currentLevelSequence[i - 1] = i - 1;
        }
        currentLevelSequence[r] = 1;
        h2 = n;
        p = n;
        q = p - 1;
        c = INT_MAX;
    }
    else
    {
        if(p == INT_MAX)
        {
            if(L[oldp - 2] != 2)
            {
                p = oldp - 1;
            }
            else
            {
                p = oldp - 2;
            }
            q = currentLevelSequence[p - 1];
        }
        if(needh2 == 1)
        {
            h2 = n;
            if(L[h2 - 1] == L[h1 - 1] - 1 && h1 == r)
            {
                c = n + 1;
            }
            else
            {
                c = INT_MAX;
            }
        }
    }
}
