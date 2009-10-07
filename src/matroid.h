// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev: 315 $ $Date: 2009-09-22 21:12:42 -0400 (Tue, 22 Sep 2009) $


/// \file matroid.h
#ifndef MATROID_H
#define MATROID_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <set>
#include <list>
#include <cmath>
#include "graph.h"
#include "matrix.h"

using namespace::std;

#define VECTOR_MATROID 1 
#define GRAPH_MATROID 2 
#define ABSTRACT_MATROID 3
#define POINTSET_MATROID 4
#define UNIFORM_MATROID 4

//needed for bases estimation
const int MAXDIM = 300; 
//RAND_MAX declared in cstdlib as max value returned by rand function
const int MAXIMAL_RAND = RAND_MAX;
// I didn't know what to put here
//const int MAXIMAL_RAND = 100000;
const float SOLVER_ACCURACY = 0.000001;
const float PI = 3.14159254;

class Matroid {
public:
    Matroid ();
    Matroid (std::istream &);
    virtual ~Matroid ();
    friend std::ostream& operator<< (std::ostream& o, Matroid &someMatroid);
    friend std::istream& operator>> (std::istream& in, Matroid &someMatroid);
    virtual int rank () {return -1;};
    virtual int calcBases () {return -1;}; //Returns number of bases
    virtual int isBasis (set <unsigned> S) {return -1;}; //Returns 1 if it is a basis, 0 else
    virtual int setRank (set <unsigned> S) {return -1;}; //Returns the rank of the set
    virtual set <unsigned> randomBasis () = 0; //Returns a random basis
    int getNumElements() {return numElements;};
    // This will initialize all data for enumerating
    // Neighbors of initBasis. This is so specialized 
    // techniques can be used for graphical matroids for example
    virtual void initializePivot(set <unsigned> initBasis) = 0;
    // This will set &pivot to the next pivot from initBasis
    // Returns 1 if &pivot is set, i.e. there is a next pivot
    // Returns 0 else. Only returns basis
    virtual int nextPivot(set <unsigned> &pivot) = 0;
    int matroidType ();

    // Given a column vector Weight, this returns the max basis under the weighting.
    set <unsigned> GreedyAlgorithmMax(Matrix Weight);
	
	//functions used in EstimateBases
	float random_weight_logistic(float x);
	float H_fn_logistic(float t);
	float upper_logistic(float avg_GAMMA, int k);
	float lower_logistic(float avg_GAMMA, int k);
	int modified_rand();
	
	//Will estimate the number of bases of a matrix using Barvinok
	//void EstimateBases();

    // This will calculate all bases and project by the matrix Weight
    virtual set <Matrix, ltcolvec> calcAllBasesProj(Matrix &Weight); 
protected:
    virtual void printMatroid (std::ostream &o) { o << "Matroid called\n";};
    virtual void getMatroid (std::istream &in)  { };
    int numElements;
    int matroidRank;
    int thisMatroidType;
    list < set <unsigned> > Bases;
    list < set <unsigned> > Independents;
    list < set <unsigned> > Circuits;
    list < set <unsigned> > Cocircuits;
};

class VectorMatroid : public Matroid {
public:
    VectorMatroid();
    VectorMatroid(std::istream &in);
    VectorMatroid(Matrix M);
    ~VectorMatroid();

    int rank ();
    int calcBases ();
    int isBasis (set <unsigned> S); //Returns 1 if it is a basis, 0 else
    int setRank (set <unsigned> S); //Returns the rank of the set
    set <unsigned> randomBasis (); //Returns a random basis

    // This will initialize all data for enumerating
    // Neighbors of initBasis. This is so specialized 
    // techniques can be used for graphical matroids for example
    // This uses no specialized technique. Simply checks all (i,j) and
    // see if initbasis - i + j is a basis
    void initializePivotGeneric(set <unsigned> initBasis);
    // This will set &pivot to the next pivot from initBasis
    // Returns 1 if &pivot is set, i.e. there is a next pivot
    // Returns 0 else. Only returns basis
    int nextPivotGeneric(set <unsigned> &pivot);

    // This uses a specialized technique for vector matroids
    // If B is the matrix given by initBasis and A_j is
    // some column of matrixRep, then if we solve
    // Bx=A_j, for any x_i that are non-zero (considering some error)
    // A_i can be swaped for A_j.
    // This calls LAPACK dgesv_ function to solve Bx=A_j
    void initializePivotLAPACK(set <unsigned> initBasis);
    int nextPivotLAPACK(set <unsigned> &pivot);

    void initializePivot(set <unsigned> initBasis);
    int nextPivot(set <unsigned> &pivot);

protected:
    void printMatroid (std::ostream &o);
    void getMatroid(std::istream &in);

    set <unsigned> pivotBasis;
    set <unsigned> currentBasis;
    set <unsigned> remSet, addSet; // These hold all elements to try and remove and add
    set <unsigned>::iterator remEl, addEl,si;
    set <unsigned> currentCycle;
    Matrix      currentInitBasis;
    unsigned     remSetCount;
    unsigned     addSetCount;
    
    Matrix matrixRep;
};

class GraphicalMatroid : public Matroid {
public:
    GraphicalMatroid(); 
    GraphicalMatroid(std::istream &in);
    GraphicalMatroid(Graph G);
    ~GraphicalMatroid();

    int rank ();
    //int calcBases ();

    // someElements is indexed by the total elements, not by internal graph index
    int isBasis (set <unsigned> S); //Returns 1 if it is a basis, 0 else
    int setRank (set <unsigned> S); //Returns the rank of the set

    //int setRank (set <unsigned> S); //Returns the rank of the set
    set <unsigned> randomBasis (); //Returns a random basis

    // This will initialize all data for enumerating
    // Neighbors of initBasis. This is so specialized 
    // techniques can be used for graphical matroids for example
    void initializePivot(set <unsigned> initBasis);
    // This will set &pivot to the next pivot from initBasis
    // Returns 1 if &pivot is set, i.e. there is a next pivot
    // Returns 0 else. Only returns basis
    int nextPivot(set <unsigned> &pivot);

    // This functions returns data with respect to the pivot basis
    // This return the cycle formed by adding the edge (n1,n2)
    // Since pivot basis is a spanning forest, there will always
    // be a cycle if (n1,n2) is NOT in pivotBasis.
    set <unsigned> getCycle (unsigned n1, unsigned n2);

    set <Matrix, ltcolvec> calcAllBasesProj(Matrix &Weight); 

protected:
    void printMatroid (std::ostream &o);
    void getMatroid(std::istream &in);

    set <unsigned> pivotBasis;
    set <unsigned> currentBasis;
    set <unsigned>  addSet; // These hold all elements to try and remove and add
    set <unsigned>::iterator remEl, addEl,si;
    set <unsigned> currentCycle;
    unsigned     remSetCount;
    unsigned     addSetCount;

    Graph   graphRep;
    // This holds the predecessor matrix for initialize pivot
    // If want path from node i to j, take path i to predMatrix(i,j) to j
    Matrix  predMatrix;
};

class PointSetMatroid : public Matroid {
public:

protected:
};

class AbstractMatroid : public Matroid {
public:
protected:
};

class UniformMatroid : public Matroid {
public:
    UniformMatroid ();
    UniformMatroid (int rank, int elements);
    int rank () {return matroidRank;};
    //int calcBases (); //Returns number of bases
    int isBasis (set <unsigned> S); //Returns 1 if it is a basis, 0 else
    int setRank (set <unsigned> S); //Returns the rank of the set
    set <unsigned> randomBasis () ; //Returns a random basis

    // This will initialize all data for enumerating
    // Neighbors of initBasis. This is so specialized 
    // techniques can be used for graphical matroids for example
    void initializePivot(set <unsigned> initBasis);
    // This will set &pivot to the next pivot from initBasis
    // Returns 1 if &pivot is set, i.e. there is a next pivot
    // Returns 0 else. Only returns basis
    int nextPivot(set <unsigned> &pivot);
    
protected:
    void printMatroid (std::ostream &o) ;
    void getMatroid (std::istream &in) ;

    set <unsigned> pivotBasis;
    set <unsigned> currentBasis;
    set <unsigned> remSet, addSet; // These hold all elements to try and remove and add
    set <unsigned>::iterator remEl, addEl,si;
    unsigned     remSetCount;
    unsigned     addSetCount;
};

#endif
