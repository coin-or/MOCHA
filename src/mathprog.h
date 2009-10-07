// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev: 315 $ $Date: 2009-09-22 21:12:42 -0400 (Tue, 22 Sep 2009) $

/// \file mathprog.h
#ifndef MATHPROG_H
#define MATHPROG_H

#include <iostream>
#include <set>
#include <list>
#include <ctime>
#include "matrix.h"
#include "matroid.h"

using namespace::std;

class Functional {
public:
    Functional() {};
    ~Functional() {};
    virtual double Evaluate (Matrix M) = 0; 
protected:
};

class FunctionalGeneral : public Functional {
public:
    FunctionalGeneral (double (*tempFct) (Matrix M));
    ~FunctionalGeneral ();
    double Evaluate (Matrix M); 
protected:
    double (*BalFct) (Matrix M) ; 
};

class FunctionalMinVariance : public Functional {
public:
    FunctionalMinVariance (Matrix M);
    ~FunctionalMinVariance ();
    double Evaluate (Matrix M); 
protected:
    Matrix WeightSum;
};

class FunctionalQuadratic : public Functional {
public:
    FunctionalQuadratic (Matrix M);
    ~FunctionalQuadratic ();
    double Evaluate (Matrix M);
protected:
    Matrix Apex;
};

class FunctionalLinear : public Functional {
public:
    FunctionalLinear (Matrix M);
    ~FunctionalLinear ();
protected:
    double Evaluate (Matrix M);
    Matrix C;
};

class MathProg {
public:
    MathProg ();
    MathProg (std::istream &in);
    virtual ~MathProg ();
    friend std::ostream& operator<< (std::ostream& o, MathProg &someMathProg);
    friend std::istream& operator>> (std::istream& in, MathProg &someMathProg);
    
protected:
    virtual void printMathProg (std::ostream &o) { o << "MathProg called\n";};
    virtual void getMathProg (std::istream &in)  { };
};

class MatroidOpt : public MathProg {
public:
    MatroidOpt () {M = 0;};
protected:
    unsigned    MatroidType;
    Matroid *M;
};

class ProjMatroidOpt : public MatroidOpt{
public:
    // This will calculate all projected spanning trees
    // only if the underlying matroid is a graph. Uses
    // Matsui algorithm to find all spanning trees and the
    // weighting to project
    set <Matrix, ltcolvec> calcAllProjBases ();
protected:
    Matrix  Weight;
};

class ProjBalMatroidOpt : public ProjMatroidOpt {
public:
    ProjBalMatroidOpt ();
    ProjBalMatroidOpt (std::istream &in);
    ProjBalMatroidOpt (std::istream &in, double (*tempFct) (Matrix));
    void setBalFct (double (*tempFct) (Matrix));
     
    // This assumes BalFct is a convex function and we wish to minimize
    // Actually, BalFct need not be convex, but if not the algorithms behavior
    // is unexpected.
    // This implements the steepest decent pivoting algorithm
    // Also known as local search
    list <set <unsigned> > LocalSearch (set <unsigned> firstBasis);
    list <set <unsigned> > LocalSearchRandomStart ();

    //list <set <unsigned> > LocalSearch (set <unsigned> firstBasis, double (*tempFct) (Matrix));
    list <set <unsigned> > LocalSearch (set <unsigned> firstBasis, Functional *F);
    //list <set <unsigned> > LocalSearchRandomStart (double (*tempFct) (Matrix));
    list <set <unsigned> > LocalSearchRandomStart (Functional *F);


    // This assumes BalFct is a convex function and we wish to minimize
    // Actually, BalFct need not be convex, but if not the algorithms behavior
    // is unexpected.
    // This implements the steepest decent pivoting algorithm
    // Also known as local search
    list <set <unsigned> > FirstComeFirstServe (set <unsigned> firstBasis);
    list <set <unsigned> > FirstComeFirstServeRandomStart ();

    //list <set <unsigned> > FirstComeFirstServe (set <unsigned> firstBasis, double (*tempFct) (Matrix));
    list <set <unsigned> > FirstComeFirstServe (set <unsigned> firstBasis, Functional *F);
    //list <set <unsigned> > FirstComeFirstServeRandomStart (double (*tempFct) (Matrix));
    list <set <unsigned> > FirstComeFirstServeRandomStart (Functional *F);


    // This implements the Tabu search method. It always pivots to the minimum 
    // neighbor under the BalFct. It will also never return to a previously seen basis.
    // pivotLimit species the threshold of pivots allowed without improvement
    // Returns the bases pivoted to.
    list <set <unsigned> > TabuSearchHeuristic (set <unsigned> firstBasis, unsigned pivotLimit, list <set <unsigned> > &tabuBasesList);
    list <set <unsigned> > TabuSearchHeuristicRandomStart (unsigned pivotLimit, list <set <unsigned> > &tabuBasesList);

    //list <set <unsigned> > TabuSearchHeuristic (set <unsigned> firstBasis, double (*tempFct) (Matrix), unsigned pivotLimit, list <set <unsigned> > &tabuBasesList);
    list <set <unsigned> > TabuSearchHeuristic (set <unsigned> firstBasis, Functional *F, unsigned pivotLimit, list <set <unsigned> > &tabuBasesList);
    //list <set <unsigned> > TabuSearchHeuristicRandomStart (double (*tempFct) (Matrix), unsigned pivotLimit, list <set <unsigned> > &tabuBasesList);
    list <set <unsigned> > TabuSearchHeuristicRandomStart (Functional *F, unsigned pivotLimit, list <set <unsigned> > &tabuBasesList);

    // This function takes in a list of subsets and outputs the
    // min bases under the functional and weighting  for this object 
    set <unsigned> FindMin (list <set <unsigned> > &);

    // This function takes in a list of subsets and outputs the
    // min bases under the given functional and weighting 
    set <unsigned> FindMin (list <set <unsigned> > &, Functional *);

    // This function takes in a list of subsets and outputs the
    // min bases under the functional for this object 
    Matrix FindMin (set <Matrix, ltcolvec> &);

    // This function takes in a list of subsets and outputs the
    // min bases under the given functional 
    Matrix FindMin (set <Matrix, ltcolvec> &, Functional *);

    //list <set <unsigned> > ProjBalMatroidOpt::SimulatedAnnealing(set <unsigned> firstBasis, list < double > temperatures, list < unsigned > times, list <set <unsigned> > &minBases);
    list <set <unsigned> > SimulatedAnnealing(set <unsigned> firstBasis, list < double > temperatures, list < unsigned > times, list <set <unsigned> > &minBases);
    list <set <unsigned> > SimulatedAnnealing(set <unsigned> firstBasis, list < double >, list < unsigned >, list <set <unsigned> > &, Functional *F);

    // This takes in a basis, a temperature T, and a function tempFct. It outputs an adjacent basis
    // according to the Boltzmann distribution using the Metropolis Chain method.
    //set <unsigned> MetropolisBoltzmannUpdateFunction(set <unsigned>, double T, double (*tempFct) (Matrix));
    set <unsigned> MetropolisBoltzmannUpdateFunction(set <unsigned>, double T, Functional *F);
    // This function assumes tempFct is BalFct
    set <unsigned> MetropolisBoltzmannUpdateFunction(set <unsigned>, double T);

    // This will enumerate a superset of the boundary of the convex hull starting with firstBasis
    // Currently only works when the projected dimension is 2
    list <set <unsigned> > Boundary(set <unsigned> firstBasis, set <Matrix, ltcolvec> &CH);
    // Does the same, except will calculate the firstBasis with a linear program
    list <set <unsigned> > Boundary(set <Matrix, ltcolvec> &CH);


    // This function takes in pareto Boundary points (in 2d) and will calculate triangular regions
    // where other pareto optimum may lie.
    set <Matrix, ltcolvec> BoundaryTrianglesTwoDim ( set <Matrix, ltcolvec> &CH);

    set <Matrix, ltcolvec> BFSDifferentFiber(set <unsigned> firstBasis);
    set <Matrix, ltcolvec> BFSDifferentFiberRandomStart();

    // This will attempt up to numTests using LocalSearch to minimize
    // a convex function that zeros only at the weightedSet point.
    // Returns 1 if weighted set is a projected point, 0 if can not determine
    int PivotTestLocalSearch(Matrix weightedSet, int numTests);
    set <Matrix, ltcolvec> PivotTestLocalSearch(set <Matrix, ltcolvec> &, int numTests);
    
    // This will attempt up to numTests using TabuSearch to minimize
    // a convex function that zeros only at the weightedSet point.
    // Returns 1 if weighted set is a projected point, 0 if can not determine
    int PivotTestTabuSearch(Matrix weightedSet, int numTests, int pivotLimit);
    set <Matrix, ltcolvec> PivotTestTabuSearch(set <Matrix, ltcolvec> &, int numTests, int pivotLimit);


    // This will tests all points in the box defined by lowerCorner and upperCorner
    // using PivotTestLocalSearch numTests times. It also will not test points in projPoints. 
    void BoxPivotTestLocalSearch(const Matrix &lowerCorner, const Matrix &upperCorner, set <Matrix, ltcolvec> &projPoints, int numTests, set <Matrix, ltcolvec> &newPoints);

    // This will tests all points in the box defined by lowerCorner and upperCorner
    // using PivotTestTabuSearch numTests times, with pivotLimit. It also will not test points in projPoints. 
    void BoxPivotTestTabuSearch(const Matrix &lowerCorner, const Matrix &upperCorner, set <Matrix, ltcolvec> &projPoints, int numTests, set <Matrix, ltcolvec> &newPoints, int pivotLimit);

    // This will use local search to find upper and lower bounds for each dimension.
    // Then it will call BoxPivotTestLocalSearch
    // using PivotTestLocalSearch numTests times. It also will not test points in projPoints. 
    void AutoBoundsPivotTestLocalSearch( set <Matrix, ltcolvec> &projPoints, int numTests, set <Matrix, ltcolvec> &newPoints);

    // This will use local search to find upper and lower bounds for each dimension.
    // Then it will call BoxPivotTestTabuSearch
    // using PivotTestTabuSearch numTests times, with pivotLimit. It also will not test points in projPoints. 
    void AutoBoundsPivotTestTabuSearch( set <Matrix, ltcolvec> &projPoints, int numTests, set <Matrix, ltcolvec> &newPoints, int pivotLimit);

    // This will perform numSearches many BFS random starts
    // It will always try to find a new random basis to start mod projection
    // from. If it can not find a new random basis before newRandTolerance attempts
    // it returns what it currently has calculated. -1 =: infinity
    // findAllBoundary == 1 means do not find random boundary elements but
    // compute the entire boundary and run bfs on each point
    void MultiBFSRandomStarts(int numSearches, int BFSSearchDepth, int newRandToleranceBoundary, int findAllBoundary, int newRandTolerance, set <Matrix, ltcolvec> &);
    Matrix projectSet(set <unsigned>);

    // Allisons work.
    // This takes in points in R^2 and outputs the Pareto Optimum.
    set <Matrix, ltcolvec> ParetoOptimum( set <Matrix, ltcolvec> &);

	// This takes in points in R^2 and outputs the minmax Optimum.
    set <Matrix, ltcolvec> MinMax(set <Matrix, ltcolvec> &);

    // This will create a random vector and pivot to a terminal basis and return it.
    // By basis we mean a projected basis.
    set <unsigned> randomLinearBasis();

    static void printPivots(const list <set <unsigned> > &);
    static void printBFSList (const set <Matrix, ltcolvec> &);
    static void printBFSListMatlab (const set <Matrix, ltcolvec> &,string);
    void printPivotsMin(const list <set <unsigned> > &);
    //void printPivotsMin(const list <set <unsigned> > &,double (*tempFct) (Matrix));
    void printPivotsMin(const list <set <unsigned> > &, Functional *F);
    void printPivotsMatlab(const list <set <unsigned> > &,string);

    static void writePivots(const list <set <unsigned> > &, std::ostream &);
    static void writeBFSList (const set <Matrix, ltcolvec> &, std::ostream &);
    static void writeBFSListMatlab (const set <Matrix, ltcolvec> &, std::ostream &,string);
    void writePivotsMin(const list <set <unsigned> > &, std::ostream &);
    //void writePivotsMin(const list <set <unsigned> > &, std::ostream &,double (*tempFct) (Matrix));
    void writePivotsMin(const list <set <unsigned> > &, std::ostream &, Functional *F);
    void writePivotsMatlab(const list <set <unsigned> > &, std::ostream &,string);

    double  evalFct(set <unsigned>);
    static unsigned BFSLevel;
    static int BFSTerminateLevel; // -1 signifies infinity
    static unsigned BFSPrinted;
    static unsigned pivotsPrinted;
    unsigned    projDim (); //Returns the number of rows of the weighting matrix
    set <unsigned> randomBasis ();
protected:
    unsigned totalBoxTests;
    time_t   BoxTestsModTime;
    unsigned BoxTestsModValue;
    unsigned BoxTests;
    // This is the function that optimizes.
    // It takes in a list since we do not know the proper dimensions are run time
    //void BFSDifferentFiberInternal(set <unsigned>,set <Matrix, ltcolvec> &,set <Matrix, ltcolvec> &);
    void BFSDifferentFiberInternal(set <unsigned>,set <Matrix, ltcolvec> &);

    // Internal recursive function to go through lowerCorner - upperCorner points and test
    void BoxPivotTestLocalSearchRec(int colIndex, Matrix &currentMatrix, const Matrix &lowerCorner, const Matrix &upperCorner, set <Matrix, ltcolvec> &projPoints, set <Matrix, ltcolvec> &newPoints);
    // Internal recursive function to go through lowerCorner - upperCorner points and test
    void BoxPivotTestTabuSearchRec(int colIndex, Matrix &currentMatrix, const Matrix &lowerCorner, const Matrix &upperCorner, set <Matrix, ltcolvec> &projPoints, set <Matrix, ltcolvec> &newPoints, int pivotLimit);
    double (*BalFct) (Matrix M) ; 
    Functional *BalanceFunction;
    virtual void printMathProg (std::ostream &o);
    void getMathProg (std::istream &in);
};

class MinVarianceBalClustering : public ProjBalMatroidOpt {
public:
    MinVarianceBalClustering ();
    // This constructor reads in the data points. There is an implied uniform matroid structure.
    // If n points are read in, there is a rank n/2 uniform matroid on n elements. 
    MinVarianceBalClustering (std::istream &in);

    // This takes in a subset of columns of Weights (the data points)
    // and returns two matrices; one that is a submatrix of weights using S (M1) and
    // the remaining (M2)
    //void getClusters(set <unsigned> S, Matrix &M1, Matrix &M2);

    // S is the basis that represents half of the points.
    // label is the matlab label to use. It will use label1 and label2 
    void writeClustersMatlab (set <unsigned> S, std::ostream &,string label);
protected:
    virtual void printMathProg(std::ostream &o);
};


class TwoDSetAngle 
{
    public:
    set <unsigned> S; //This holds a set
    double angle; // This is the angle this set with some fixed set makes with the origin
    unsigned dely;
    unsigned delx;
};

struct ltTwoDSetAngle {
    bool operator()(const TwoDSetAngle &S1, const TwoDSetAngle &S2) const
    {
        if (S1.angle < S2.angle)
        {
            return 1;
        }
        return 0;    
    }

};

class TwoDAngleOrderedList
{ 
public:
    //TwoDAngleOrderedList();
    TwoDAngleOrderedList(set <unsigned> initialSet, ProjBalMatroidOpt *PBMO);
    ~TwoDAngleOrderedList();
    // Inserts S into the list using refPBMO for angle calculation
    void insert(set <unsigned> S);
    // Returns the two sets with the largest angle between them and the angle
    // Returns 0 if there are less than 2 sets in the list.
    double largestAngle (set <unsigned> &S1, set <unsigned> &S2);

    void print ();
protected:
    ProjBalMatroidOpt   *refPBMO;
    double x,y; //x and y values for the initial set
    set < TwoDSetAngle, ltTwoDSetAngle> orderedRays; 
};

#endif
