# CSFTreeConjecture
Contains C++ code to verify Richard Stanley's conjecture that the chromatic symmetric function distinguishes distinct trees, on at most 28 vertices. This project uses the pthread library for multithreading.

This project was inspired by Keeler Russell's verification of the conjecture for all trees with at most 25 vertices, given here: https://github.com/keeler/csf and uses Russell's library for generating non-isomorphic trees, but compares the chromatic symmetric functions of the trees using an alternate algorithm to allow a greater volume of trees to be processed.

The files Tree.h, Tree.cpp, TreeGenerator.h, and TreeGenerator.cpp are slightly edited versions of Keeler Russell's code from https://github.com/keeler/csf. 

The file matchfinder2.cpp contains the main code from our most recent computational verification of Stanley's conjecture on trees with at most 28 vertices. The other files contain code from other, previous attempts.


