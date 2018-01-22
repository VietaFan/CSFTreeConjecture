/*
TreeGenerator.hpp - written by Keeler Russell
https://github.com/keeler/csf/blob/master/src/TreeGenerator.hpp
*/
#ifndef TREE_GENERATOR_HPP
#define TREE_GENERATOR_HPP

#include "Tree.hpp"

// Reference:
// Robert Alan Wright, L. Bruce Richmond, Andrew M. Odlyzko, and Brendan D. McKay.
// Constant time generation of free trees.
// SIAM J. Comput., Vol. 15, No. 2, pages 540–548, 1986.

class TreeGenerator
{
    public:
        TreeGenerator(unsigned int numVertices);
        ~TreeGenerator();

        // Get the next tree generated, returns true until there are no more trees
        // with the specified number of vertices. If false, the tree t won't necessarily
        // contain anything meaningful, so take care using it in such a case.
        bool nextTree(Tree & t);
				bool nextTreeParList(int *parList);
				
        
    private:
        unsigned int    mNumVertices;
        int             *L;
        int             *currentLevelSequence;

        int             p;
        int             q;
        int             h1;
        int             h2;
        int             c;
        int             r;
				bool            firstTime;

        void generateFirstLevelSequence();
        void generateNextLevelSequence();
};

#endif
