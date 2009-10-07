// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev: 315 $ $Date: 2009-09-22 21:12:42 -0400 (Tue, 22 Sep 2009) $

/// \file mathprog.cpp

#include <iostream>
#include <set>
#include <list>
#include <string>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <ctime>
#include "matrix.h"
#include "matroid.h"
#include "mathprog.h"

unsigned ProjBalMatroidOpt::BFSLevel = 0;
unsigned ProjBalMatroidOpt::pivotsPrinted = 1;
int ProjBalMatroidOpt::BFSTerminateLevel = -1; // -1 signifies infinity
unsigned ProjBalMatroidOpt::BFSPrinted = 0;

FunctionalGeneral::FunctionalGeneral (double (*tempFct) (Matrix M))
{
    BalFct = tempFct;
}

double FunctionalGeneral::Evaluate (Matrix M)
{
    return (*BalFct)(M);
}

FunctionalMinVariance::FunctionalMinVariance (Matrix M)
{
    WeightSum = M;
    cout << "WeightSum:" << endl << WeightSum;
}

double FunctionalMinVariance::Evaluate(Matrix M)
{
    //cout << "       FunctionalMinVariance::Evaulation called with matrix:" << endl << "***" << endl << M << "***" << endl;
    //cout << "FunctionalMinVariance::Evaluate called." << endl;
    //cout << M;
    double  tempD=0;    
    //cout << "tempD =" << tempD << endl;

    tempD = M.twoNormSquared ();
    //cout << "tempD =" << tempD << endl;

    Matrix tempMatrix = M - WeightSum;
    //cout << "Weight sum" << endl << WeightSum;
    //cout << tempMatrix;
    tempD += tempMatrix.twoNormSquared();
    //cout << "tempD =" << tempD << endl;

    return (-1)*tempD;
}

FunctionalQuadratic::FunctionalQuadratic (Matrix M)
{
    if (M.cols != 1)
    {
        cerr << "FunctionalQuadratic::FunctionalQuadratic called with non-column vector Matrix M." << endl;
        exit (0);
    }
    Apex = M;
}

FunctionalQuadratic::~FunctionalQuadratic ()
{
}

double FunctionalQuadratic::Evaluate (Matrix M)
{
    double tempD=0;
    
    if (M.rows != Apex.rows || M.cols != Apex.cols)
    {
        cerr << "FunctionalQuadratic::Evaluate called with mismatching matrices." << endl;
        exit (0);
    }

    for (unsigned i=0;i<M.rows;i++)
    {
        tempD += (Apex(i,0) - M(i,0))*(Apex(i,0) - M(i,0));
    }
    return tempD;
}

FunctionalLinear::FunctionalLinear (Matrix M)
{
    C = M;
}

FunctionalLinear::~FunctionalLinear ()
{
}

double FunctionalLinear::Evaluate (Matrix M)
{
    Matrix tempM = (M.transpose())*C;
    return tempM(0,0);
}


MathProg::MathProg ()
{
};

MathProg::MathProg (std::istream &in)
{
    //getMathProg(in);
};

MathProg::~MathProg ()
{
};

std::ostream& operator<< (std::ostream& o, MathProg &someMathProg)
{
    someMathProg.printMathProg(o); 
    return o;
}

std::istream& operator>> (std::istream& in, MathProg &someMathProg)
{
    someMathProg.getMathProg(in);
    return in;
}

set <Matrix, ltcolvec> ProjMatroidOpt::calcAllProjBases()
{
    return M->calcAllBasesProj(Weight);
}

void ProjBalMatroidOpt::printMathProg (std::ostream &o)
{
    o << "Projected Balanced Matroid Optimization Program" << endl;

    if (M != 0)
    {
        o << "MATROID" << endl;
        o << *M;
        o << "WEIGHTING" << endl;
        Matrix::printPadLength=4;
        o << Weight;
        //o << "Row sum of Weight" << endl;
        //o << Weight.rowSum();
        Matrix::printPadLength=0;
    }
}
void ProjBalMatroidOpt::getMathProg (std::istream &in)
{
    //Determine type of matroid we are reading in
    string S;

    in >> S;
    if (S == "VECTOR")
    {
        MatroidType = VECTOR_MATROID;
        // Read in matrix representing the matroid
        M = new VectorMatroid(in);

        // Read in matrix representing the weights
        in >> Weight;
    }
    else if (S == "GRAPH")
    {
        MatroidType = GRAPH_MATROID; 
        // Read in matrix representing graph
        M = new GraphicalMatroid(in);

        // Read in matrix representing the weights
        in >> Weight;
    }

    
}


ProjBalMatroidOpt::ProjBalMatroidOpt (std::istream &in)
{
    BFSPrinted = 1;
    getMathProg(in);
    BalFct = 0;
    BalanceFunction = 0;
}

ProjBalMatroidOpt::ProjBalMatroidOpt (std::istream &in, double (*tempFct) (Matrix))
{
    BFSPrinted = 1;
    getMathProg(in);
    BalanceFunction = new FunctionalGeneral(tempFct);
    setBalFct(tempFct);
}

ProjBalMatroidOpt::ProjBalMatroidOpt ()
{
    BFSPrinted = 1;
    BalFct = 0;
}

void ProjBalMatroidOpt::setBalFct (double (*tempFct) (Matrix M))
{
    BalFct = tempFct;
}

list <set <unsigned> > ProjBalMatroidOpt::LocalSearchRandomStart ()
{
    return LocalSearchRandomStart(BalanceFunction);
}

list <set <unsigned> > ProjBalMatroidOpt::LocalSearch (set <unsigned> firstBasis)
{
    return LocalSearch(firstBasis, BalanceFunction);
}

//list <set <unsigned> > ProjBalMatroidOpt::LocalSearchRandomStart (double (*tempFct) (Matrix))
list <set <unsigned> > ProjBalMatroidOpt::LocalSearchRandomStart (Functional *F)
{
    set <unsigned> S;

    S = M->randomBasis();
    //cout << "Random basis: ";
    //copy(S.begin(), S.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;

    return LocalSearch (S,F);
}

//list <set <unsigned> > ProjBalMatroidOpt::LocalSearch (set <unsigned> firstBasis, double (*tempFct) (Matrix) )
list <set <unsigned> > ProjBalMatroidOpt::LocalSearch (set <unsigned> firstBasis, Functional *F )
{
    if (M == 0)
    {
        cerr << "LocalSearch called with no matroid." << endl;
        exit (0);
    }
    if (F == 0)
    {
        cerr << "LocalSearch called with no Functional." << endl;
        exit (0);
    }
    if (M->isBasis(firstBasis) != 1)
    {
        cerr << "LocalSearch called with bad first basis." << endl;
        copy(firstBasis.begin(), firstBasis.end(), ostream_iterator<const int>(cout, " "));
        cout << endl;
        exit (0);
    }
    
    double  currentMin;
    int     pivoted=1;
    set <unsigned> minBasis, currentBasis;
    list <set <unsigned> > pivotBases;
    
    minBasis = firstBasis; // in case no neighbors are min 
    // While we get new pivot basis
    while (pivoted == 1)
    {
        //cout << "Pivot " << pivotBases.size () << endl;
        pivoted=0;
        //copy(minBasis.begin(), minBasis.end(), ostream_iterator<const int>(cout, " "));
        //cout << endl;
        //cout << projectSet(minBasis);
        pivotBases.push_back(minBasis);
        M->initializePivot(minBasis);
        currentMin = F->Evaluate( (Weight.subColumns(minBasis)).rowSum() );
        //cout << "Current Min: " << currentMin << endl;
        while(M->nextPivot(currentBasis))
        {
            Matrix evaluateMatrix = Weight.subColumns(currentBasis).rowSum();
            Matrix tempM = projectSet(currentBasis);
            //cout << "Checking pivot" << endl;
            //cout << "   ";
            //cout << projectSet(currentBasis);
            //copy(currentBasis.begin(), currentBasis.end(), ostream_iterator<const int>(cout, " "));
            //cout << endl;
            //Matrix::printPadLength = 4;
            //cout << projectSet(currentBasis); 
            //Matrix::printPadLength = 0;
            //cout << "Input for Evaluate:" << endl;
            //cout << evaluateMatrix;
            //cout << Weight.subColumns(currentBasis).rowSum();
            //cout << Weight.subColumns(currentBasis).rowSum();
            //cout << Weight.subColumns(currentBasis).rowSum();
            //cout << "   " << F->Evaluate( (Weight.subColumns(currentBasis)).rowSum()) << endl;
            //cout << "   " << F->Evaluate( (Weight.subColumns(currentBasis)).rowSum()) << endl;
            //cout << "   " << F->Evaluate( (Weight.subColumns(currentBasis)).rowSum()) << endl;
            //cout << "   Attempt 1: " << F->Evaluate(evaluateMatrix) << endl;
            //cout << "   Attempt 2: " << F->Evaluate(evaluateMatrix) << endl;
            //cout << "   Attempt 3: " << F->Evaluate(evaluateMatrix) << endl;
            //if (F->Evaluate( (Weight.subColumns(currentBasis)).rowSum() ) < currentMin)
            if (F->Evaluate(evaluateMatrix) < currentMin)
            {
                //cout << "       Pivoted!!!" << endl;
                pivoted = 1;
                // Then set minBasis = currentBasis
                minBasis = currentBasis;
                currentMin = F->Evaluate( (Weight.subColumns(currentBasis)).rowSum() );
            }
        }
        if (pivoted == 1)
        {
            // If we pivot and the matroid is graphical, then find all paths
            //pivotBases.push_back(minBasis);
            currentBasis = minBasis;
        }
    } 
    return pivotBases;         
}

list <set <unsigned> > ProjBalMatroidOpt::FirstComeFirstServeRandomStart ()
{
    return FirstComeFirstServeRandomStart(BalanceFunction);
}

list <set <unsigned> > ProjBalMatroidOpt::FirstComeFirstServe (set <unsigned> firstBasis)
{
    return FirstComeFirstServe(firstBasis, BalanceFunction);
}

//list <set <unsigned> > ProjBalMatroidOpt::FirstComeFirstServeRandomStart (double (*tempFct) (Matrix))
list <set <unsigned> > ProjBalMatroidOpt::FirstComeFirstServeRandomStart (Functional *F)
{
    set <unsigned> S;

    S = M->randomBasis();
    //cout << "Random basis: ";
    //copy(S.begin(), S.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;

    return FirstComeFirstServe(S,F);
}

list <set <unsigned> > ProjBalMatroidOpt::FirstComeFirstServe (set <unsigned> firstBasis, Functional *F )
{
    if (M == 0)
    {
        cerr << "FirstComeFirstServe called with no matroid." << endl;
        exit (0);
    }
    if (F == 0)
    {
        cerr << "FirstComeFirstServe called with no Functional." << endl;
        exit (0);
    }
    if (M->isBasis(firstBasis) != 1)
    {
        cerr << "FirstComeFirstServe called with bad first basis." << endl;
        copy(firstBasis.begin(), firstBasis.end(), ostream_iterator<const int>(cout, " "));
        cout << endl;
        exit (0);
    }
    
    double  currentMin;
    int     pivoted=1;
    set <unsigned> minBasis, currentBasis;
    list <set <unsigned> > pivotBases;
    
    minBasis = firstBasis; // in case no neighbors are min 
    // While we get new pivot basis
    while (pivoted == 1)
    {
        pivoted=0;
        pivotBases.push_back(minBasis);
        M->initializePivot(minBasis);
        currentMin = F->Evaluate( (Weight.subColumns(minBasis)).rowSum() );
        while(M->nextPivot(currentBasis) && pivoted == 0) // Exit on the first pivot!
        {
            Matrix evaluateMatrix = Weight.subColumns(currentBasis).rowSum();
            Matrix tempM = projectSet(currentBasis);
            if (F->Evaluate(evaluateMatrix) < currentMin)
            {
                pivoted = 1;
                // Then set minBasis = currentBasis
                minBasis = currentBasis;
                currentMin = F->Evaluate( (Weight.subColumns(currentBasis)).rowSum() );
            }
        }
        if (pivoted == 1)
        {
            // If we pivot and the matroid is graphical, then find all paths
            //pivotBases.push_back(minBasis);
            currentBasis = minBasis;
        }
    } 
    return pivotBases;         
}

// Tabu Search
list <set <unsigned> > ProjBalMatroidOpt::TabuSearchHeuristicRandomStart (unsigned pivotLimit, list <set <unsigned> > &tabuBasesList)
{
    return TabuSearchHeuristicRandomStart(BalanceFunction,pivotLimit,tabuBasesList);
}

list <set <unsigned> > ProjBalMatroidOpt::TabuSearchHeuristic (set <unsigned> firstBasis, unsigned pivotLimit, list <set <unsigned> > &tabuBasesList)
{
    return TabuSearchHeuristic(firstBasis, BalanceFunction,pivotLimit,tabuBasesList);
}

//list <set <unsigned> > ProjBalMatroidOpt::TabuSearchHeuristicRandomStart (double (*tempFct) (Matrix), unsigned pivotLimit, list <set <unsigned> > &tabuBasesList)
list <set <unsigned> > ProjBalMatroidOpt::TabuSearchHeuristicRandomStart (Functional *F, unsigned pivotLimit, list <set <unsigned> > &tabuBasesList)
{
    set <unsigned> S;

    S = M->randomBasis();
    //cout << "Random basis: ";
    //copy(S.begin(), S.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;

    return TabuSearchHeuristic (S,F,pivotLimit,tabuBasesList);
}

//list <set <unsigned> > ProjBalMatroidOpt::TabuSearchHeuristic (set <unsigned> firstBasis, double (*tempFct) (Matrix), unsigned pivotLimit, list <set <unsigned> > &tabuBasesList )
list <set <unsigned> > ProjBalMatroidOpt::TabuSearchHeuristic (set <unsigned> firstBasis, Functional *F, unsigned pivotLimit, list <set <unsigned> > &tabuBasesList )
{
    if (M == 0)
    {
        cerr << "TabuSearchHeuristic called with no matroid." << endl;
        exit (0);
    }
    if (F == 0)
    {
        cerr << "TabuSearchHeuristic called with no Functional" << endl;
        exit (0);
    }
    if (M->isBasis(firstBasis) != 1)
    {
        cerr << "TabuSearchHeuristic called with bad first basis." << endl;
        copy(firstBasis.begin(), firstBasis.end(), ostream_iterator<const int>(cout, " "));
        cout << endl;
        exit (0);
    }
    tabuBasesList.push_back(firstBasis);
        
    double  currentMin;
    double  tempMin;
    unsigned pivotsSinceNoImprovement = 0;
    set <unsigned> minBasis, currentBasis, currentMinBasis;
    list <set <unsigned> > pivotBases;
    set <set <unsigned> > tabuBases;
    set <set <unsigned> >::const_iterator tbi;


    minBasis = firstBasis; // in case no neighbors are min 
    currentMin = F->Evaluate( (Weight.subColumns(minBasis)).rowSum() );
    currentMinBasis = firstBasis;
    tabuBases.insert(minBasis);
    // Keep pivoting until some criteria
    while (pivotsSinceNoImprovement < pivotLimit)
    {
        M->initializePivot(minBasis);
        // Take the first pivot and set as tempmin
        if (M->nextPivot(currentBasis))
        {
            tempMin = F->Evaluate( (Weight.subColumns(currentBasis)).rowSum() );
            minBasis = currentBasis;
        }
        while(M->nextPivot(currentBasis))
        {
            // Pivot if it is the minimum under tempFct and it has not been visited.
            if (F->Evaluate( (Weight.subColumns(currentBasis)).rowSum() ) < tempMin)
            {
                tbi = tabuBases.find(currentBasis);
                if (tbi == tabuBases.end())
                {
                    tempMin = F->Evaluate( (Weight.subColumns(currentBasis)).rowSum() );
                    // Then set minBasis = currentBasis
                    minBasis = currentBasis;
                }
            }
        }
        if (tempMin < currentMin)
        {
            currentMin = tempMin;
            currentMinBasis = minBasis;
            pivotsSinceNoImprovement=0;
            tabuBasesList.push_back(currentMinBasis);
            tabuBases.insert(currentMinBasis);
            //cout << "Pivots " << tabuBasesList.size() << endl;
        }
        pivotBases.push_back(minBasis);
        tabuBases.insert(minBasis);
        pivotsSinceNoImprovement++;
    } 
    return pivotBases;         
}

// Tabu Search End

// This function takes in a list of subsets and outputs the
// min bases under the functional and weighting for this object 
set <unsigned> ProjBalMatroidOpt::FindMin (list <set <unsigned> > &subSets)
{
    return FindMin(subSets, BalanceFunction);

}

// This function takes in a list of subsets and outputs the
// min bases under the given functional and weighting 
set <unsigned> ProjBalMatroidOpt::FindMin (list <set <unsigned> > &subSets, Functional *F)
{
    list <set <unsigned> >::const_iterator lsit = subSets.begin();
    //cout << "Calling Eval" << endl;
    double currentMin = F->Evaluate(projectSet(*lsit));
    //cout << "Done Calling Eval" << endl;
    set <unsigned> minBasis = *lsit;

    while (lsit != subSets.end())
    {
        if (F->Evaluate(projectSet(*lsit)) < currentMin)
        {
            currentMin = F->Evaluate(projectSet(*lsit));
            minBasis = *lsit;
        }
        lsit++;
    }
    return minBasis;
}

// This function takes in a list of subsets and outputs the
// min bases under the functional and weighting for this object 
Matrix ProjBalMatroidOpt::FindMin (set <Matrix, ltcolvec> &points)
{
    return FindMin(points, BalanceFunction);

}

// This function takes in a list of subsets and outputs the
// min bases under the given functional and weighting 
Matrix ProjBalMatroidOpt::FindMin (set <Matrix, ltcolvec> &points, Functional *F)
{
    set <Matrix, ltcolvec>::const_iterator mit = points.begin();
    double currentMin = F->Evaluate(*mit);
    Matrix minPoint = *mit;

    while (mit != points.end())
    {
        if (F->Evaluate(*mit) < currentMin)
        {
            currentMin = F->Evaluate(*mit);
            minPoint = *mit;
        }
    }
    return minPoint;
}

void ProjBalMatroidOpt::printPivotsMin(const list <set <unsigned> > &pivotBases)
{
    writePivotsMin(pivotBases,cout,BalanceFunction);
}

//void ProjBalMatroidOpt::printPivotsMin(const list <set <unsigned> > &pivotBases,double (*tempFct) (Matrix))
void ProjBalMatroidOpt::printPivotsMin(const list <set <unsigned> > &pivotBases, Functional *F)
{
    writePivotsMin(pivotBases,cout,F);
}

void ProjBalMatroidOpt::writePivotsMin(const list <set <unsigned> > &pivotBases, std::ostream &o)
{
    writePivotsMin(pivotBases,o,BalanceFunction);
}    

//void ProjBalMatroidOpt::writePivotsMin(const list <set <unsigned> > &pivotBases, std::ostream &o,double (*tempFct) (Matrix))
void ProjBalMatroidOpt::writePivotsMin(const list <set <unsigned> > &pivotBases, std::ostream &o,Functional *F)
{
    list <set <unsigned> >::const_iterator li;
    set <unsigned> TS1;

    o << "Pivots" << endl;
    li = pivotBases.begin();
    while (li != pivotBases.end())
    {
        TS1 = *li;    

        o << "       ";
        copy(TS1.begin(), TS1.end(), ostream_iterator<const int>(o, " "));
        //o << "       f(.) = " << evalFct(TS1) <<endl;
        o << "       f(.) = " << F->Evaluate ( (Weight.subColumns(TS1)).rowSum() ) <<endl;

        li++;
    }
    o << endl;
}

void ProjBalMatroidOpt::printBFSList (const set <Matrix, ltcolvec> &BFSList)
{
    ProjBalMatroidOpt::writeBFSList(BFSList,cout);
}

void ProjBalMatroidOpt::writeBFSList (const set <Matrix, ltcolvec> &BFSList, std::ostream &o)
{
    set <Matrix, ltcolvec>::const_iterator cmit = BFSList.begin();

    while (cmit != BFSList.end())
    {
        o << (*cmit) << "---" << endl;; 
        cmit++;
    }
}

void ProjBalMatroidOpt::printPivots(const list <set <unsigned> > &pivotBases)
{
    ProjBalMatroidOpt::writePivots(pivotBases,cout);
}

void ProjBalMatroidOpt::writePivots(const list <set <unsigned> > &pivotBases, std::ostream &o)
{
    list <set <unsigned> >::const_iterator li;
    set <unsigned> TS1;

    o << "Pivots" << endl;
    li = pivotBases.begin();
    while (li != pivotBases.end())
    {
        TS1 = *li;    

        o << "       ";
        copy(TS1.begin(), TS1.end(), ostream_iterator<const int>(o, " "));
        o << endl;
        li++;
    }

}

void ProjBalMatroidOpt::printBFSListMatlab (const set <Matrix, ltcolvec> &BFSList, string label)
{
    ProjBalMatroidOpt::writeBFSListMatlab(BFSList,cout,label); 
}

void ProjBalMatroidOpt::writeBFSListMatlab (const set <Matrix, ltcolvec> &BFSList, std::ostream &o, string label)
{
    set <Matrix, ltcolvec>::const_iterator cmit = BFSList.begin();
    
    o << label << "= [" << endl;
    while (cmit != BFSList.end())
    {
        for (unsigned i=0;i<(*cmit).rows;i++)
        {
            o << (*cmit)(i,0) << " ";
        }
        o << ";" << endl;
        cmit++; 
    }
    o << "];" << endl;
    
    BFSPrinted++;
}

void ProjBalMatroidOpt::printPivotsMatlab(const list <set <unsigned> > &pivotBases,string label)
{
    ProjBalMatroidOpt::writePivotsMatlab(pivotBases,cout,label);
}

void ProjBalMatroidOpt::writePivotsMatlab(const list <set <unsigned> > &pivotBases, std::ostream &o,string label)
{
    list <set <unsigned> >::const_iterator li;
    Matrix  M;

    o << label << " = [" << endl;
    li = pivotBases.begin();
    while (li != pivotBases.end())
    {
        M = projectSet(*li);
        for (unsigned i=0;i<M.rows;i++)
        {
            o << M(i,0) << " ";
        }
        o << ";" << endl;
    
        li++;
    }
    o << "];" << endl;
    pivotsPrinted++;
}


double  ProjBalMatroidOpt::evalFct(set <unsigned> someSet)
{
    if (BalFct != 0)
    {
        return (*BalFct) ( (Weight.subColumns(someSet)).rowSum() );
    }
    cerr << "ProjBalMatroidOpt::evalFct called with BalFct not set.\n";
    exit (0);
    return 0;
}

void ProjBalMatroidOpt::MultiBFSRandomStarts(int numSearches, int BFSSearchDepth, int newRandToleranceBoundary, int findAllBoundary, int newRandTolerance, set <Matrix, ltcolvec> &BFSResults)
{
    //set <Matrix, ltcolvec> BFSResults;
    if (numSearches <= 0)
    {
        return ; 
    }

    int oldBFSTerminateLevel = BFSTerminateLevel;
    set <unsigned> randomBasis;

    if (findAllBoundary != 1 || Weight.rows != 2) 
    {
        int newRandToleranceBoundaryTenPercent = (int) newRandToleranceBoundary / 10;
        if (newRandToleranceBoundaryTenPercent == 0)
        {
            newRandToleranceBoundaryTenPercent = 1; 
        }
        BFSTerminateLevel = BFSSearchDepth; 
        cout << "Performing " << numSearches << " BFS searches on boundary of convex hull." << endl;
        for(int i=0;i<numSearches;i++)
        {
            cout << "Finding new pseudo random boundary point." << endl;
            BFSLevel = 0;
            int numRandAttempts = 0;
            //randomBasis = M->randomBasis();
            randomBasis = randomLinearBasis();
            while (randomBasis.size() == 0)
            {
                randomBasis = randomLinearBasis();
            }
            set <Matrix, ltcolvec>::const_iterator mit = BFSResults.find(projectSet(randomBasis));
            numRandAttempts++;
            time_t startTime = time(0);
            while (mit != BFSResults.end())
            {
                if ((numRandAttempts % newRandToleranceBoundaryTenPercent) == 0)
                {

                    cout << numRandAttempts / newRandToleranceBoundaryTenPercent <<  "0" << "%" << endl;
                    cout << "           " << time(0) - startTime << " seconds." << endl;
                    startTime = time(0);
                }
                //randomBasis = M->randomBasis();
                randomBasis = randomLinearBasis();
                while (randomBasis.size() == 0)
                {
                    randomBasis = randomLinearBasis();
                }
                mit = BFSResults.find(projectSet(randomBasis));
                //cout << "@";
                numRandAttempts++;
                if (newRandToleranceBoundary != -1 && numRandAttempts > newRandToleranceBoundary)
                {
                    cout << "Exceeded Random Retry Limit." << endl;
                    cout << "BFSResults.size() = " << BFSResults.size() << endl;
                    //return; // Don't exit. Still more tests to compute
                    i = numSearches;
                    break;
                }
            }
            cout << "(" << i << ") Random Linear Search Basis (" << numRandAttempts << "): " << endl << projectSet(randomBasis) << endl;
            if (i < numSearches)
            {
                BFSResults.insert(projectSet(randomBasis));
                //set <Matrix, ltcolvec>  grey;
                //BFSDifferentFiberInternal(randomBasis, BFSResults,grey);
                BFSDifferentFiberInternal(randomBasis, BFSResults);
                cout << "BFSResults.size() = " << BFSResults.size() << endl;
            }
        }
    }
    else if (findAllBoundary == 1 && Weight.rows == 2)
    {
        cout << "Computing entire boundary." << endl;
        BFSTerminateLevel = BFSSearchDepth; 

        set <Matrix, ltcolvec> boundaryPoints; 
        list < set <unsigned> > boundaryBases = Boundary(boundaryPoints);
        boundaryPoints.clear ();

        list < set <unsigned> >::iterator lsit = boundaryBases.begin();
        // For each boundaryBases, run bfs

        cout << "Testing each new boundary point." << endl;
        int skipCount = 0;
        for (;lsit != boundaryBases.end();lsit++)
        {
            BFSLevel = 0;
            // Check to see if this projected bases has already been found
            set <Matrix, ltcolvec>::const_iterator mit = BFSResults.find(projectSet(*lsit));
            if (mit == BFSResults.end())
            {
                cout << "Skipped " << skipCount << endl;
                skipCount = 0;
                cout << projectSet(*lsit);
                BFSResults.insert(projectSet(*lsit));
                BFSDifferentFiberInternal(*lsit, BFSResults);
                cout << "BFSResults.size() = " << BFSResults.size() << endl;
            }
            else {
                skipCount++;
            }
        }
    }

    int newRandToleranceTenPercent = (int) newRandTolerance / 10;
    if (newRandToleranceTenPercent == 0)
    {
        newRandToleranceTenPercent = 1; 
    }
    cout << "Performing " << numSearches << " BFS searches on random interior." << endl;
    for(int i=0;i<numSearches;i++)
    {
        cout << "Finding new random point." << endl;
        BFSLevel = 0;
        int numRandAttempts = 0;
        randomBasis = M->randomBasis();
        while (randomBasis.size() == 0)
        {
            randomBasis = M->randomBasis();
        }
        if (randomBasis.size () == 0)
        {
            cerr << "Not supposed to be empty randomBasis returned by M.randomBasis()" << endl;
            exit(0);
        }
        //randomBasis = randomLinearBasis();
        set <Matrix, ltcolvec>::const_iterator mit = BFSResults.find(projectSet(randomBasis));
        numRandAttempts++;
        time_t startTime = time(0);
        while (mit != BFSResults.end())
        {
            if ((numRandAttempts % newRandToleranceTenPercent) == 0)
            {
                cout << numRandAttempts / newRandToleranceTenPercent << "0" << "%" << endl;
                cout << "           " << time(0) - startTime << " seconds." << endl;
                startTime = time(0);
            }
            randomBasis = M->randomBasis();
            while (randomBasis.size() == 0)
            {
                randomBasis = M->randomBasis();
            }

            if (randomBasis.size () == 0)
            {
                cerr << "Not supposed to be empty randomBasis returned by M.randomBasis()" << endl;
                exit(0);
            }
            //randomBasis = randomLinearBasis();
            mit = BFSResults.find(projectSet(randomBasis));
            //cout << "@";
            numRandAttempts++;
            if (newRandTolerance != -1 && numRandAttempts > newRandTolerance)
            {
                cout << "Exceeded Random Retry Limit." << endl;
                cout << "BFSResults.size() = " << BFSResults.size() << endl;
                return; 
            }
        }
        cout << "(" << i << ") Random projected point (" << numRandAttempts << "): " << endl << projectSet(randomBasis) << endl;
        BFSResults.insert(projectSet(randomBasis));
        //set <Matrix, ltcolvec>  grey;
        //BFSDifferentFiberInternal(randomBasis, BFSResults, grey);
        BFSDifferentFiberInternal(randomBasis, BFSResults);
        cout << "BFSResults.size() = " << BFSResults.size() << endl;
    }
    cout << "Done performing BFS enumeration." << endl;

    BFSTerminateLevel = oldBFSTerminateLevel;
    return; 
}

set <Matrix, ltcolvec> ProjBalMatroidOpt::BFSDifferentFiberRandomStart()
{
    set <unsigned> S;

    S = M->randomBasis();

    return BFSDifferentFiber (S);
}

set <Matrix, ltcolvec> ProjBalMatroidOpt::BFSDifferentFiber(set <unsigned> firstBasis)
{
    set <Matrix, ltcolvec>  currentLevel;
    //set <Matrix, ltcolvec>  grey;
    
    BFSLevel = 0;
    
    //BFSDifferentFiberInternal(firstBasis, currentLevel, grey);
    BFSDifferentFiberInternal(firstBasis, currentLevel);
    return currentLevel;
}

//void ProjBalMatroidOpt::BFSDifferentFiberInternal(set <unsigned> currentBasis, set <Matrix, ltcolvec> &currentLevel, set <Matrix, ltcolvec> &grey)
void ProjBalMatroidOpt::BFSDifferentFiberInternal(set <unsigned> currentBasis, set <Matrix, ltcolvec> &currentLevel)
{
    //cout << "BFS Level: " << BFSLevel << endl;
    if (BFSTerminateLevel != -1 && (int)BFSLevel > BFSTerminateLevel)
    {
        return;
    }
    // Only insert if BFSDifferentFiberInternal called
    //currentLevel.insert(projectSet(currentBasis));
    // currentBasis is marked black, so unmark it as grey
    //set <Matrix, ltcolvec>::const_iterator greyit = grey.find(projectSet(currentBasis));
    //if (greyit != grey.end())
    //{
    //    grey.erase(greyit);
    //}

    BFSLevel++;
    set <unsigned> pivotBasis;
    double  currentMin;

    //currentMin = (*BalFct) ( (Weight.subColumns(currentBasis)).rowSum() );
    currentMin = evalFct(currentBasis);
    list <set <unsigned> > basesToTest;

    M->initializePivot(currentBasis);

    unsigned totalAdjacent = 0;
    unsigned totalDiffEval = 0;
    while (M->nextPivot(pivotBasis) == 1)
    {
        totalAdjacent++;
        // BFS should be independent of evalFct
        // if project(pivotBasis) != project(currentBasis)
        //if (evalFct(pivotBasis) != currentMin)
        if (!(projectSet(pivotBasis) == projectSet(currentBasis)))
        {
            totalDiffEval++;
            set <Matrix, ltcolvec>::const_iterator mit = currentLevel.find(projectSet(pivotBasis));
            //set <Matrix, ltcolvec>::const_iterator greyit = grey.find(projectSet(pivotBasis));
            // If it is not black (currentLevel) or grey then it is white and we add to
            // list to recursively call bfs
            if ((mit == currentLevel.end()))// && (greyit == grey.end()))
            {
                ///grey.insert(projectSet(pivotBasis));
                currentLevel.insert(projectSet(pivotBasis));
                basesToTest.push_back(pivotBasis);
            }
        }
    }
    //cout << basesToTest.size() << " number of bases to test at level " << BFSLevel << endl;
    if (1 == 0)
    {
        double adjacentToNewBases = (double)basesToTest.size () / (double)totalAdjacent;
        adjacentToNewBases *= 100;

        ios::fmtflags oldFlags;
        int oldStreamSize;
        oldStreamSize = cout.precision(1);
        oldFlags = cout.setf (ios::fixed);

        cout << adjacentToNewBases << "% adjacent bases eval differently and new." << endl;
        cout.setf (oldFlags);
        cout.precision(oldStreamSize);
    }

    //cout << basesToTest.size () << "/" << totalAdjacent << endl;
    list <set <unsigned> >::iterator btti = basesToTest.begin ();
    set <unsigned> tempBasis;

    while (btti != basesToTest.end() )
    {
        tempBasis = *btti;
        //BFSDifferentFiberInternal(tempBasis, currentLevel, grey);
        BFSDifferentFiberInternal(tempBasis, currentLevel);
        btti++;
    }

    BFSLevel--;
}




// This takes in points in R^2 and outputs the Pareto Optimum.
set <Matrix, ltcolvec> ProjBalMatroidOpt::ParetoOptimum( set <Matrix, ltcolvec> &inputPoints)
{
    // This is the set of pareto optimum we will calculate.
    // This is a set structure from the standard library which
    // holds Matrix classes that we defined. ltcolvec is a bool
    // function that compares two Matrix classes with one column.
    
    // Please do not modify inputPoints
    set <Matrix, ltcolvec> popt = inputPoints;

    // Perform some error checking
    if (inputPoints.empty())
    {
        return popt;
    }
  
    // An iterator to go through all inputPoints.
    set <Matrix, ltcolvec>::iterator msetit = popt.begin();
	
	//An iterator to go in front of points so when we erase we can get back to where we started
	set <Matrix, ltcolvec>::iterator msetit2 = popt.begin();
	msetit2++;
	
	//An iterator to compare to
	set <Matrix, ltcolvec>::iterator stepit = popt.begin();
	
	int counter = 0;
	
	for( ; msetit != popt.end(); msetit++)
	{
	  
	  //// Set this temp matrix to the current matrix msetit is refering to
      Matrix tempM = *msetit;
	  for(stepit = popt.begin(); stepit != popt.end(); stepit++)
	  {
	    if (stepit != msetit)
		{
	    counter = 0;
		Matrix tempM2 = *stepit;
		for(int i = 0; i < tempM.rows; i++)
		{
		  if(tempM(i,0) >= tempM2(i,0)) //if all coordinates are larger than another point
		    counter++;
        }
		if(counter == tempM.rows)
		{
		  popt.erase(msetit);
		  msetit = msetit2;
		  msetit--;
	      break;
		}
		}
	  }
	  msetit2++;
	}

    // tempM(0,0); // Example of the (0,0) position in the matrix
    // tempM(1,0);
    // tempM(2,0);
    // tempM.rows // number of rows in this matrix
    // tempM.cols // number of columns in this matrix

    //// Make sure we are dealing with matrices with 2 rows and 1 column
    //if (tempM.rows != 2 || tempM.cols != 1)
    //{
    //    cerr << "ParetoOptimum called with incorrect inputPoints." << endl;
    //    exit (0);
    //}


    // Now go through and find Pareto optumum

    return popt;
}


set <Matrix, ltcolvec> ProjBalMatroidOpt::MinMax(set <Matrix, ltcolvec> &inputPoints)
{
	set <Matrix, ltcolvec> mmpoints;
	mmpoints.clear();

	// Perform some error checking
    if (inputPoints.empty())
    {
        return mmpoints;
    }
	
	// An iterator to go through all inputPoints.
    set <Matrix, ltcolvec>::iterator msetit = inputPoints.begin();
	
	//An array to hold the max coordinate for each point
	int m = inputPoints.size();
	int maxcoor[m];
	
	int mspot = 0;
	
	//put max coordinate for each point in maxcoor array
	for( ; msetit != inputPoints.end(); msetit++)
	{
		Matrix tempM = *msetit;
		
		maxcoor[mspot] = tempM(0,0);
		for(int i = 1; i < tempM.rows; i++)
		{
			if(tempM(i,0) > maxcoor[mspot])
				maxcoor[mspot] = tempM(i,0);
		}
		mspot++;
	}
	
	//find min amongst all max coordinates
	int min = maxcoor[0];
	for(int j = 1; j < m; j++)
	{
		if(maxcoor[j] < min)
			min = maxcoor[j];
	}
	  
	//Add all points with max coor equal to min to the set mmpoints
	msetit = inputPoints.begin();
	for(int k = 0; k < m; k++)
	{
		if(maxcoor[k] == min)
			mmpoints.insert(*msetit);
			
		msetit++;
	} 
	
	return mmpoints;
}


Matrix ProjBalMatroidOpt::projectSet(set <unsigned> someSet)
{
    return (Weight.subColumns(someSet)).rowSum();
}

// Expects to be column vector
Matrix SpecFctMatrix;

double SpecFct (Matrix M)
{
    //cout << "SpecFct called" << endl;
    double returnValue = 0;
    
    for (int i=0;i<SpecFctMatrix.getRows();i++)
    {
        returnValue += (M(i,0)-SpecFctMatrix(i,0))*(M(i,0)-SpecFctMatrix(i,0));
    }
    //cout << "SpecFct done" << endl;
    return returnValue;
}   

// Expects to be column vector
Matrix SpecLinearFctMatrix;

double SpecLinearFct (Matrix M)
{
    //cout << "SpecLinearFct called" << endl;
    double returnValue = 0;
    
    for (int i=0;i<SpecLinearFctMatrix.getRows();i++)
    {
        returnValue += (M(i,0)*SpecLinearFctMatrix(i,0));
    }
    //cout << "SpecLinearFct done" << endl;
    return returnValue;
}   

int ProjBalMatroidOpt::PivotTestLocalSearch(Matrix weightedSet, int numTests)
{
    //cout << "pivottest called" << endl;
    if (numTests < 0)
    {
        return 0;
    }
    //SpecFctMatrix = weightedSet;
    FunctionalQuadratic MyQuad(weightedSet);
    list <set <unsigned> > resultPivots;
    list <set <unsigned> >::iterator pit;
    set <unsigned> lastPivot;

    for (int i=0;i<numTests;i++)
    {
        resultPivots = LocalSearchRandomStart (&MyQuad);
        //printPivotsMin(resultPivots,SpecFct); 
        pit = resultPivots.end();
        pit--;
        //copy((*pit).begin(), (*pit).end(), ostream_iterator<const int>(cout, " "));
        //cout << Weight.subColumns(*pit);
        //cout << "i=" << i << "   Value: " << SpecFct( (Weight.subColumns(*pit)).rowSum()) <<  endl;
        if (MyQuad.Evaluate( (Weight.subColumns(*pit)).rowSum()) == 0 )
        {
            //cout << "pivottest done" << endl;
            return 1;
        }
    }

    //cout << "pivottest done" << endl;
    return 0;
}

set <Matrix, ltcolvec> ProjBalMatroidOpt::PivotTestLocalSearch(set <Matrix, ltcolvec> &testPoints, int numTests)
{
    set <Matrix, ltcolvec> returnPoints;
    set <Matrix, ltcolvec>::const_iterator mit = testPoints.begin();

    cout << testPoints.size() << " points to test." << endl;    
    unsigned tenPercent = testPoints.size() / 10;
    unsigned i=0;
    for (;mit != testPoints.end();mit++)
    {
        //if ((i % 100) == 0)
        //{
        //    cout << i << endl;
        //}
        if ((i % tenPercent) == 0 && tenPercent != 0)
        {
            cout << i/tenPercent << "0%" << endl;
        }
        i++;
        if (PivotTestLocalSearch(*mit,numTests) == 1)
        {
            returnPoints.insert(*mit);
        }
    }
    return returnPoints;
}

int ProjBalMatroidOpt::PivotTestTabuSearch(Matrix weightedSet, int numTests, int pivotLimit)
{
    if (numTests < 0)
    {
        return 0;
    }
    //SpecFctMatrix = weightedSet;
    FunctionalQuadratic MyQuad(weightedSet);
    list <set <unsigned> > resultPivots;
    list <set <unsigned> >::iterator pit;
    set <unsigned> lastPivot;

    for (int i=0;i<numTests;i++)
    {
        TabuSearchHeuristicRandomStart (&MyQuad,pivotLimit,resultPivots);
        if (resultPivots.size() == 0)
        {
            cerr << "PivotTestTabuSearch: resultPivots.size () == 0." << endl;
            exit(0);
        }
        //printPivotsMin(resultPivots,SpecFct); 
        //cout << "resultPivots.size() : " << resultPivots.size() << endl;
        pit = resultPivots.end();
        pit--;
        //copy((*pit).begin(), (*pit).end(), ostream_iterator<const int>(cout, " "));
        //cout << Weight.subColumns(*pit);
        //cout << "i=" << i << "   Value: " << SpecFct( (Weight.subColumns(*pit)).rowSum()) <<  endl;
        if (MyQuad.Evaluate( (Weight.subColumns(*pit)).rowSum()) == 0 )
        {
            //cout << "pivottest done" << endl;
            return 1;
        }
    }

    return 0;
}

set <Matrix, ltcolvec> ProjBalMatroidOpt::PivotTestTabuSearch(set <Matrix, ltcolvec> &testPoints, int numTests, int pivotLimit)
{
    set <Matrix, ltcolvec> returnPoints;
    set <Matrix, ltcolvec>::const_iterator mit = testPoints.begin();
    Matrix tempM;

    cout << testPoints.size() << " points to test." << endl;    
    unsigned tenPercent = testPoints.size() / 10;
    if (tenPercent == 0)
    {
        tenPercent = 1;
    }
    unsigned i=0;
    for (;mit != testPoints.end();mit++)
    {
        if ((i % tenPercent) == 0 && tenPercent != 0)
        {
            cout << i/tenPercent << "0%" << endl;
        }
        //if ((i % 100) == 0)
        //{
        //    cout << i << endl;
        //}
        i++;
        if (PivotTestTabuSearch(*mit,numTests, pivotLimit) == 1)
        {
            returnPoints.insert(*mit);
        }
    }
    return returnPoints;
}

// This will use local search to find upper and lower bounds for each dimension.
// Then it will call BoxPivotTestLocalSearch
// using PivotTestLocalSearch numTests times. It also will not test points in projPoints. 
void ProjBalMatroidOpt::AutoBoundsPivotTestLocalSearch( set <Matrix, ltcolvec> &projPoints, int numTests, set <Matrix, ltcolvec> &newPoints)
{
    Matrix lowBound(Weight.rows,1);
    Matrix highBound(Weight.rows,1);
    Matrix objective(Weight.rows,1);
    FunctionalLinear *LF;
    list <set <unsigned> > LSResult;
    list <set <unsigned> >::iterator lsit;

    // For each dimension, find the high and low bounds
    for (unsigned i=0;i<Weight.rows;i++)
    {
        for (unsigned j=0;j<Weight.rows;j++)
        {
            objective(j,0) = 0;
        }
        objective(i,0) = 1;
        LF = new FunctionalLinear(objective);

        LSResult = LocalSearchRandomStart(LF);
        if (LSResult.size() > 0)
        {
            lsit = LSResult.end();
            lsit--;

            lowBound(i,0) = ((Weight.subColumns(*lsit)).rowSum())(i,0);
        }
        delete LF;

        for (unsigned j=0;j<Weight.rows;j++)
        {
            objective(j,0) = 0;
        }
        objective(i,0) = -1;
        LF = new FunctionalLinear(objective);

        LSResult = LocalSearchRandomStart(LF);
        if (LSResult.size() > 0)
        {
            lsit = LSResult.end();
            lsit--;

            highBound(i,0) = ((Weight.subColumns(*lsit)).rowSum())(i,0);
        }
        delete LF;
    }

    cout << "Low Bounds" << endl;
    cout << lowBound;
    cout << "High Bounds" << endl;
    cout << highBound;
    BoxPivotTestLocalSearch(lowBound, highBound, projPoints, numTests, newPoints);
}

// This will use local search to find upper and lower bounds for each dimension.
// Then it will call BoxPivotTestTabuSearch
// using PivotTestTabuSearch numTests times, with pivotLimit. It also will not test points in projPoints. 
void ProjBalMatroidOpt::AutoBoundsPivotTestTabuSearch( set <Matrix, ltcolvec> &projPoints, int numTests, set <Matrix, ltcolvec> &newPoints, int pivotLimit)
{
    Matrix lowBound(Weight.rows,1);
    Matrix highBound(Weight.rows,1);
    Matrix objective(Weight.rows,1);
    FunctionalLinear *LF;
    list <set <unsigned> > LSResult;
    list <set <unsigned> >::iterator lsit;

    // For each dimension, find the high and low bounds
    for (unsigned i=0;i<Weight.rows;i++)
    {
        for (unsigned j=0;j<Weight.rows;j++)
        {
            objective(j,0) = 0;
        }
        objective(i,0) = 1;
        LF = new FunctionalLinear(objective);

        LSResult = LocalSearchRandomStart(LF);
        if (LSResult.size() > 0)
        {
            lsit = LSResult.end();
            lsit--;

            lowBound(i,0) = ((Weight.subColumns(*lsit)).rowSum())(i,0);
        }
        delete LF;

        for (unsigned j=0;j<Weight.rows;j++)
        {
            objective(j,0) = 0;
        }
        objective(i,0) = -1;
        LF = new FunctionalLinear(objective);

        LSResult = LocalSearchRandomStart(LF);
        if (LSResult.size() > 0)
        {
            lsit = LSResult.end();
            lsit--;

            highBound(i,0) = ((Weight.subColumns(*lsit)).rowSum())(i,0);
        }
        delete LF;
    }

    cout << "Low Bounds" << endl;
    cout << lowBound;
    cout << "High Bounds" << endl;
    cout << highBound;
    BoxPivotTestTabuSearch(lowBound, highBound, projPoints, numTests, newPoints, pivotLimit);
}

void ProjBalMatroidOpt::BoxPivotTestLocalSearch(const Matrix &lowerCorner, const Matrix &upperCorner, set <Matrix, ltcolvec> &projPoints, int numTests, set <Matrix, ltcolvec> &newPoints)
{
    Matrix currentMatrix((lowerCorner.rows),1); 
    BoxTests = 0;
    totalBoxTests = (unsigned)(upperCorner(0,0) + 1 - lowerCorner(0,0));
    for (unsigned i=1;i<lowerCorner.rows;i++)
    {
        totalBoxTests *= (unsigned)(upperCorner(i,0) + 1 - lowerCorner(i,0));
    }
    totalBoxTests *= numTests;
    cout << "Total tests to perform: " << totalBoxTests << endl;
    BoxTestsModValue = (unsigned)fabs((double)totalBoxTests*0.10);
    cout << "BoxTestsModValue = " << BoxTestsModValue  << endl;
    for (int i=0;i<numTests;i++)
    {
        BoxTestsModTime = time(0);
        BoxPivotTestLocalSearchRec(0, currentMatrix, lowerCorner, upperCorner, projPoints, newPoints);
        cout << endl;
    }   
}

void ProjBalMatroidOpt::BoxPivotTestLocalSearchRec(int colIndex, Matrix &currentMatrix, const Matrix &lowerCorner, const Matrix &upperCorner, set <Matrix, ltcolvec> &projPoints, set <Matrix, ltcolvec> &newPoints)
{
    //cout << "   BoxPivotTestLocalSearchRec called" << endl;
    if (colIndex == (int)lowerCorner.rows)
    {
        BoxTests++;
        if ((BoxTests % BoxTestsModValue) == 0)
        {
            cout << (lround((double)BoxTests / (double)BoxTestsModValue) )*10 << "% Box Tests Complete." << endl;
            cout << "           " << time(0) - BoxTestsModTime << " seconds." << endl;
            BoxTestsModTime = time(0);
        }
        //cout << "[" << BoxTests << "] ";
        //cout.flush ();
        //cout << "   Testing Point" << endl;
        //cout << "       Reached colIndex == lowerCorner.rows" << endl;
        // Check if the current matrix is in projPoints. If not try to show it is a proj point
        set <Matrix, ltcolvec>::const_iterator mit;
        mit = projPoints.find(currentMatrix);
        if (mit == projPoints.end() )
        {
            //cout << "            Point not in projPoints" << endl;
            //cout << "------" << endl << currentMatrix;
            if (PivotTestLocalSearch(currentMatrix,1) == 1)
            {
                //cout << "*";
                //cout.flush();
                projPoints.insert(currentMatrix); 
                newPoints.insert(currentMatrix); 
            }    
        }
        return;
    }

    for(int i=(int)lowerCorner(colIndex,0);i<=(int)upperCorner(colIndex,0);i++)
    {
        currentMatrix(colIndex,0) = i;
        BoxPivotTestLocalSearchRec(colIndex+1,currentMatrix, lowerCorner, upperCorner, projPoints, newPoints);
    }
}

void ProjBalMatroidOpt::BoxPivotTestTabuSearch(const Matrix &lowerCorner, const Matrix &upperCorner, set <Matrix, ltcolvec> &projPoints, int numTests, set <Matrix, ltcolvec> &newPoints, int pivotLimit)
{
    Matrix currentMatrix((lowerCorner.rows),1); 
    BoxTests = 0;
    totalBoxTests = (unsigned)(upperCorner(0,0) + 1 - lowerCorner(0,0));
    for (unsigned i=1;i<lowerCorner.rows;i++)
    {
        totalBoxTests *= (unsigned)(upperCorner(i,0) + 1 - lowerCorner(i,0));
    }
    totalBoxTests *= numTests;
    cout << "Total tests to perform: " << totalBoxTests << endl;
    BoxTestsModValue = (unsigned)fabs((double)totalBoxTests*0.10);
    cout << "BoxTestsModValue = " << BoxTestsModValue  << endl;
    for (int i=0;i<numTests;i++)
    {
        BoxTestsModTime = time(0);
        BoxPivotTestTabuSearchRec(0, currentMatrix, lowerCorner, upperCorner, projPoints, newPoints,pivotLimit);
        cout << endl;
    }   
}

void ProjBalMatroidOpt::BoxPivotTestTabuSearchRec(int colIndex, Matrix &currentMatrix, const Matrix &lowerCorner, const Matrix &upperCorner, set <Matrix, ltcolvec> &projPoints, set <Matrix, ltcolvec> &newPoints, int pivotLimit)
{
    //cout << "   BoxPivotTestTabuSearchRec called" << endl;
    if (colIndex == (int)lowerCorner.rows)
    {
        BoxTests++;
        if ((BoxTests % BoxTestsModValue) == 0)
        {
            cout << "           " << time(0) - BoxTestsModTime << " seconds." << endl;
            BoxTestsModTime = time(0);
            cout << (lround((double)BoxTests / (double)BoxTestsModValue) )*10 << "% Box Tests Complete." << endl;
        }
        //cout << "[" << BoxTests << "] ";
        //cout.flush ();
        //cout << "   Testing Point" << endl;
        //cout << "       Reached colIndex == lowerCorner.rows" << endl;
        // Check if the current matrix is in projPoints. If not try to show it is a proj point
        set <Matrix, ltcolvec>::const_iterator mit;
        mit = projPoints.find(currentMatrix);
        if (mit == projPoints.end() )
        {
            //cout << "            Point not in projPoints" << endl;
            //cout << "------" << endl << currentMatrix;
            if (PivotTestTabuSearch(currentMatrix,1,pivotLimit) == 1)
            {
                //cout << "*";
                //cout.flush();
                projPoints.insert(currentMatrix); 
                newPoints.insert(currentMatrix); 
            }    
        }
        return;
    }

    for(int i=(int)lowerCorner(colIndex,0);i<=(int)upperCorner(colIndex,0);i++)
    {
        currentMatrix(colIndex,0) = i;
        BoxPivotTestTabuSearchRec(colIndex+1,currentMatrix, lowerCorner, upperCorner, projPoints, newPoints, pivotLimit);
    }
}


set <unsigned> ProjBalMatroidOpt::randomLinearBasis()
{
    list <set <unsigned> >  linearPivots;
    Matrix  randDirection(Weight.getRows(),1);

    int allZero = 1;
    
    while (allZero == 1)
    {
        allZero = 1;
        for (int i=0;i<Weight.getRows();i++)
        {
            randDirection(i,0) = (rand() % 101) - 50;
            if ((int)randDirection(i,0) != 0)
            {
                allZero = 0;
            }
        }
        //cout << "RandomDirection" << endl;
        //cout << randDirection;
    }

    //SpecLinearFctMatrix = randDirection;
    FunctionalLinear MyLin(randDirection);

    linearPivots = LocalSearchRandomStart (&MyLin);

    list <set <unsigned> >::const_iterator si;

    si = linearPivots.end();
    si--;

    return *si; 
}

//Returns the number of rows of the weighting matrix
unsigned    ProjBalMatroidOpt::projDim ()
{
    return  Weight.getRows(); 
}
// Does the same, except will calculate the firstBasis with a linear program
list <set <unsigned> > ProjBalMatroidOpt::Boundary(set <Matrix, ltcolvec> &CH)
{
    if (Weight.rows != 2)
    {
        cerr << "ProjBalMatroidOpt::Boundary called with number of weightings !=2" << endl;
        exit (0);
    }
    set <unsigned> firstBasis = randomLinearBasis();
    //cout << "first basis: ";
    //copy(firstBasis.begin(), firstBasis.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    return Boundary(firstBasis, CH);
}

// This will enumerate a subset of the boundary of the convex hull starting with firstBasis
// All points are extremal
// Currently only works when the projected dimension is 2
list <set <unsigned> > ProjBalMatroidOpt::Boundary(set <unsigned> firstBasis, set <Matrix, ltcolvec> &CH)
{
    if (Weight.rows != 2)
    {
        cerr << "ProjBalMatroidOpt::Boundary called with number of weightings !=2" << endl;
        exit (0);
    }
    // This assumes that firstBasis is an extremal point of WPM

    list <set <unsigned> > boundaryBases;
    boundaryBases.push_back(firstBasis);
    
    // Add the first basis to the convex hull.
    CH.insert(projectSet(firstBasis));
    Matrix firstBasisProjected = projectSet(firstBasis);
    Matrix tM1, tM2;
    set <unsigned> S1,S2;
    double tangle;
    TwoDAngleOrderedList    *tempTDAOL;

    set <unsigned> currentBasis = firstBasis;
    set <unsigned> adjBasis;
    set <Matrix, ltcolvec>::const_iterator Mit;

    list <set <unsigned> > basesToTest;
    list <set <unsigned> >::iterator bttit;
    basesToTest.push_front(currentBasis);
    while (basesToTest.size() != 0)
    {
        bttit = basesToTest.begin();
        currentBasis = *bttit;

        //Matrix tempM = projectSet(currentBasis);
        //if( tempM(0,0) == 28 && tempM(1,0) == 139)
        //{
        //    cout << "current basis: ";
        //    copy(currentBasis.begin(), currentBasis.end(), ostream_iterator<const int>(cout, " "));
        //    cout << endl;
        //    cout << projectSet(currentBasis);
        //    set <unsigned> remSet,addSet;

        //    remSet = currentBasis;
        //    for (unsigned i=0;i<M->getNumElements();i++)
        //    {
        //        if (remSet.find(i) == remSet.end())
        //        {
        //            addSet.insert(i);
        //        }
        //    }
        //    cout << "remSet: ";
        //    copy(remSet.begin(), remSet.end(), ostream_iterator<const int>(cout, " "));
        //    cout << endl;
        //    cout << "addSet: ";
        //    copy(addSet.begin(), addSet.end(), ostream_iterator<const int>(cout, " "));
        //    cout << endl;
        //    cout << "Printing all possible projected sets." << endl;
        //
        //    Matrix potMatrix;
        //    set <unsigned>::iterator remsuit = remSet.begin();
        //    set <unsigned>::iterator addsuit = addSet.begin();
        //    for (;remsuit!= remSet.end();remsuit++)
        //    {
        //        currentBasis.erase(currentBasis.find(*remsuit));
        //        for (;addsuit!= addSet.end();addsuit++)
        //        {
        //            currentBasis.insert(*addsuit);
        //            cout << "Potential basis: ";
        //            copy(currentBasis.begin(), currentBasis.end(), ostream_iterator<const int>(cout, " "));
        //            cout << endl;

        //            potMatrix = projectSet(currentBasis);
        //            cout << potMatrix;
        //            //if (potMatrix(0,0) == 37 && potMatrix (1,0) == 144)
        //            //{
        //            cout << "Is basis = " << M->isBasis(currentBasis) << endl;
        //            //}
        //            cout << endl;
        //            currentBasis.erase(currentBasis.find(*addsuit));
        //        }
        //        currentBasis.insert(*remsuit);
        //    }

        //    exit (0);
        //}




        basesToTest.pop_front();
        M->initializePivot(currentBasis);
        //cout << "current basis: ";
        //copy(currentBasis.begin(), currentBasis.end(), ostream_iterator<const int>(cout, " "));
        //cout << endl;
        //cout << projectSet(currentBasis);
        tempTDAOL = new TwoDAngleOrderedList(currentBasis,this);
            
        while(M->nextPivot(adjBasis))
        {
            //cout << "Inserted adjBasis: ";
            //copy(adjBasis.begin(), adjBasis.end(), ostream_iterator<const int>(cout, " "));
            //cout << endl;
            tempTDAOL->insert(adjBasis);
        }
        
        //Matrix tempM = projectSet(currentBasis);
        //if (tempM(0,0) == 237 && tempM(1,0) == 465)
        //{
        //    tempTDAOL->print ();
        //    cout << "Max angle: " << tempTDAOL->largestAngle(S1,S2) << endl;
        //    copy(S1.begin(), S1.end(), ostream_iterator<const int>(cout, " "));
        //    cout << endl;
        //    cout << projectSet(S1);
        //    copy(S2.begin(), S2.end(), ostream_iterator<const int>(cout, " "));
        //    cout << endl;
        //    cout << projectSet(S2);
        //}
        //tempTDAOL->print ();
        //cout << "Max angle: " << tempTDAOL->largestAngle(S1,S2) << endl;
        //copy(S1.begin(), S1.end(), ostream_iterator<const int>(cout, " "));
        //cout << endl;
        //cout << projectSet(S1);
        //copy(S2.begin(), S2.end(), ostream_iterator<const int>(cout, " "));
        //cout << endl;
        //cout << projectSet(S2);

        tangle = tempTDAOL->largestAngle(S1,S2);
        //cout << "Angle: " << tangle << "   " << M_PI << endl;
        //Matrix tempM1 = projectSet(currentBasis);
        //Matrix::printPadLength = 4;
        //cout << tempM1;
        //tempTDAOL->print ();
        //Matrix::printPadLength = 0;
        //char c;
        //cin >> c;
        //if (tempM1(0,0) == 323 && tempM1(1,0) == 45)
        //{
        //    cout << "   Angle: " << tangle << endl;
        //    tempTDAOL->print();
        //}
        // Numerical instabilities. :( Maybe try to use gmp rationals
        //if (tangle >= M_PI)
        if (tangle + 0.000001 >= M_PI)
        {
            //if (tempM(0,0) == 237 && tempM(1,0) == 465)
            //{
            //    cout << "       tangle >= " << M_PI << endl;
            //}
            tM1 = projectSet(S1);
            tM2 = projectSet(S2);
            //if (tM1(0,0) == 237 && tM1(1,0) == 465)
            //{
            //    cout << "Trying to insert " << endl << tM1;
            //    cout << projectSet(currentBasis);
            //    cout << "************" << endl;
            //    tempTDAOL->print ();
            //    cout << "Max angle: " << tempTDAOL->largestAngle(S1,S2) << endl;
            //    copy(S1.begin(), S1.end(), ostream_iterator<const int>(cout, " "));
            //    cout << endl;
            //    cout << projectSet(S1);
            //    copy(S2.begin(), S2.end(), ostream_iterator<const int>(cout, " "));
            //    cout << endl;
            //    cout << projectSet(S2);
            //}
            //if (tM2(0,0) == 237 && tM2(1,0) == 465)
            //{
            //    cout << "Trying to insert " << endl << tM2;
            //    cout << projectSet(currentBasis);
            //    cout << "************" << endl;
            //    tempTDAOL->print ();
            //    cout << "Max angle: " << tempTDAOL->largestAngle(S1,S2) << endl;
            //    copy(S1.begin(), S1.end(), ostream_iterator<const int>(cout, " "));
            //    cout << endl;
            //    cout << projectSet(S1);
            //    copy(S2.begin(), S2.end(), ostream_iterator<const int>(cout, " "));
            //    cout << endl;
            //    cout << projectSet(S2);
            //}
            Matrix::printPadLength=4;
            //cout << tM1;
            Mit = CH.find(tM1);
            if (Mit == CH.end())
            {
                //cout << "       Added" << endl;
                CH.insert(tM1);
                basesToTest.push_front(S1);
                boundaryBases.push_back(S2);
                currentBasis = S1;
            }
            //cout << "    ----" << endl;
            //cout << tM2;
            Matrix::printPadLength=0;
            Mit = CH.find(tM2);
            if (Mit == CH.end())
            {
                //cout << "       Added" << endl;
                CH.insert(tM2);
                basesToTest.push_front(S2);
                boundaryBases.push_back(S2);
                currentBasis = S2;
            }
            //cout << "    ----" << endl;
        }
        else 
        {
            // Theoretically we should not get here!
            // Only boundary points should be added to basesToTest 
            cout << "Angle not larger than PI." << endl;
            cout << projectSet(currentBasis);
            cout << "************" << endl;
            tempTDAOL->print ();
            cout << "Max angle: " << tempTDAOL->largestAngle(S1,S2) << endl;
            copy(S1.begin(), S1.end(), ostream_iterator<const int>(cout, " "));
            cout << endl;
            cout << projectSet(S1);
            copy(S2.begin(), S2.end(), ostream_iterator<const int>(cout, " "));
            cout << endl;
            cout << projectSet(S2);
        }

        delete tempTDAOL;
    }
    return boundaryBases;
}

// This will enumerate a superset of the boundary of the convex hull starting with firstBasis
// Currently only works when the projected dimension is 2
set <Matrix, ltcolvec> ProjBalMatroidOpt::BoundaryTrianglesTwoDim ( set <Matrix, ltcolvec> &CH)
{
    //cout << "BoundaryTrianglesTwoDim called." << endl;
    set <Matrix, ltcolvec> triangularRegions;
    if (CH.empty() == 1)
    {
        return triangularRegions;
    }

    set <Matrix, ltcolvec>::const_iterator mit = CH.begin();

    if ((*mit).rows != 2 || (*mit).cols != 1)
    {
        cerr << "BoundaryTrianglesTwoDim called with non 2d points." << endl;
        exit (0);
    }

    // Copy CH into triangularRegions
    triangularRegions = CH;

    int lowFound;
    int highFound;
    set <Matrix, ltcolvec>::const_iterator tmit;
    Matrix tempM;
    double lowDistance;
    double highDistance;
    Matrix lowNeighbor(2,1), highNeighbor(2,1);
    //cout << "lowNeighbor " << &lowNeighbor << endl;
    //cout << "highNeighbor " << &highNeighbor << endl;

    for (;mit != CH.end();mit++)
    {
        //cout << "mit = " << &(*mit) << endl;
        // For each point on the pareto-boundary, find its closest neighbors
        // and add the triangles to triangularRegions.
        lowFound = 0;
        highFound = 0;

        for (tmit=CH.begin();tmit != CH.end();tmit++)
        {
            //cout << "tmit = " << &(*tmit) << endl;
            // Only check different points
            if ((*tmit)(0,0) != (*mit)(0,0) && (*tmit)(1,0) != (*mit)(1,0))
            {
                // This could be a new lowNeighbor
                if ((*tmit)(0,0) >= (*mit)(0,0) && (*tmit)(1,0) <= (*mit)(1,0))
                {
                    if (lowFound == 0)
                    {
                        // Now low has been found yet, so this is the best low neighbor
                        lowNeighbor(0,0) = (*tmit)(0,0);
                        lowNeighbor(1,0) = (*tmit)(1,0);

                        tempM = lowNeighbor - (*mit);
                        lowDistance = tempM.twoNormSquared();

                        lowFound = 1;
                    }
                    else
                    {
                        tempM = (*tmit) - (*mit);
                        if (tempM.twoNormSquared() < lowDistance)
                        {
                            lowNeighbor(0,0) = (*tmit)(0,0);
                            lowNeighbor(1,0) = (*tmit)(1,0);
                            lowDistance = tempM.twoNormSquared();
                        }
                    }
                }

                // This could be a new highNeighbor
                if ((*tmit)(0,0) <= (*mit)(0,0) && (*tmit)(1,0) >= (*mit)(1,0))
                {
                    if (highFound == 0)
                    {
                        // Now high has been found yet, so this is the best high neighbor
                        highNeighbor(0,0) = (*tmit)(0,0);
                        highNeighbor(1,0) = (*tmit)(1,0);

                        tempM = highNeighbor - (*mit);
                        highDistance = tempM.twoNormSquared();

                        highFound = 1;
                    }
                    else
                    {
                        tempM = (*tmit) - (*mit);
                        if (tempM.twoNormSquared() < highDistance)
                        {
                            highNeighbor(0,0) = (*tmit)(0,0);
                            highNeighbor(1,0) = (*tmit)(1,0);
                            highDistance = tempM.twoNormSquared();
                        }
                    }
                }
            }
        }
        Matrix tempM(2,1);
        // If we found a lower neighbor
        if (lowFound == 1)
        {
            double slope = ( (*mit)(1,0) - lowNeighbor(1,0) ) / ( (*mit)(0,0) - lowNeighbor(0,0) );

            for (unsigned j=lowNeighbor(1,0);j<(*mit)(1,0);j++)
            {
                unsigned tempUnsigned = ceil(lowNeighbor(0,0) + ((double)(j-lowNeighbor(1,0)))/slope);
                //cout << "(*mit)(0,0) = " << (*mit)(0,0) << "    tempUnsigned = " << tempUnsigned << endl;
                for (unsigned i=tempUnsigned;i<lowNeighbor(0,0);i++)
                //for (unsigned i=(*mit)(0,0);i<lowNeighbor(0,0);i++)
                {
                    tempM(0,0) = i;
                    tempM(1,0) = j;
                    triangularRegions.insert(tempM);
                }
            }
        }

        // If we found a higher neighbor
        if (highFound == 1)
        {
            double slope = ( (*mit)(1,0) - highNeighbor(1,0) ) / ( (*mit)(0,0) - highNeighbor(0,0) );

            for (unsigned j=(*mit)(1,0);j<highNeighbor(1,0);j++)
            {
                unsigned tempUnsigned = ceil((*mit)(0,0) + ((double)(j-(*mit)(1,0)))/slope);
                //cout << "highNeighbor(0,0) = " << highNeighbor(0,0) << "    tempUnsigned = " << tempUnsigned << endl;
                for (unsigned i=tempUnsigned;i<(*mit)(0,0);i++)
                //for (unsigned i=highNeighbor(0,0);i<(*mit)(0,0);i++)
                {
                    tempM(0,0) = i;
                    tempM(1,0) = j;
                    triangularRegions.insert(tempM);
                }
            }
        }
    }

    //cout << "BoundaryTrianglesTwoDim done." << endl;
    return triangularRegions;
}

set <unsigned> ProjBalMatroidOpt::randomBasis ()
{
    return M->randomBasis();
}

set <unsigned> ProjBalMatroidOpt::MetropolisBoltzmannUpdateFunction(set <unsigned> pivotBasis, double T)
{
    return MetropolisBoltzmannUpdateFunction(pivotBasis, T, BalanceFunction);
}

//set <unsigned> ProjBalMatroidOpt::MetropolisBoltzmannUpdateFunction(set <unsigned> pivotBasis, double T, double (*tempFct) (Matrix))
set <unsigned> ProjBalMatroidOpt::MetropolisBoltzmannUpdateFunction(set <unsigned> pivotBasis, double T, Functional *F)
{
    unsigned numNeighborsPivotBasis = 0;
    unsigned numNeighborsPotentialBasis = 0;
    set <unsigned> potentialBasis, tempBasis;

    // Count how many neighbors there are of pivotBasis.
    // We could save all neighbors, but we wish to minimize the memory footprint
    M->initializePivot(pivotBasis);
    while(M->nextPivot(potentialBasis))
    {
        numNeighborsPivotBasis++;
    }
    copy(pivotBasis.begin(), pivotBasis.end(), ostream_iterator<const int>(cout, " "));
    cout << endl;

    // Now pick a neighbor uniformly random 
    M->initializePivot(pivotBasis);
    //cout << "prob shooting for: " << (1/(double)numNeighborsPivotBasis) << endl;
    while(M->nextPivot(potentialBasis))
    {
        //cout << rand() << endl;
        //cout << RAND_MAX << endl;
        //cout << (double)rand() / (double)RAND_MAX << endl;
        if ((double)rand() / (double)RAND_MAX < (1/(double)numNeighborsPivotBasis))
        {
            cout << "   Breaking." << endl;
            break; //Break from while loop
        }
    }
    copy(potentialBasis.begin(), potentialBasis.end(), ostream_iterator<const int>(cout, " "));
    cout << endl;
    // Now potentialBasis holds the randomly chosen neighbor of pivotBasis
    
    // Count how many neighbors there are of potentialBasis.
    M->initializePivot(potentialBasis);
    while(M->nextPivot(tempBasis))
    {
        numNeighborsPotentialBasis++;
    }
    cout << "Number of Neighbors: " << numNeighborsPivotBasis << endl;
    cout << "Number of Neighbors Potential Basis: " << numNeighborsPotentialBasis << endl;

    // Now return potentialBasis with probability min(pi_j d_i / pi_i d_j, 1)
    // or return pivotBasis with probability 1 - min(pi_j d_i / pi_i d_j, 1)
    // d_i and d_j are the number of neighbors of pivotBasis and potential Basis.
    // pi_* is the probability distribution. In this case the Boltzmann distribution.

    double t1 = -(F->Evaluate((Weight.subColumns(potentialBasis)).rowSum()))/T;
    double t2 = -(F->Evaluate((Weight.subColumns(pivotBasis)).rowSum()))/T;
    cout << "t1 = " << t1 << "  exp(t1) = " << exp(t1) << endl;
    cout << "t2 = " << t2 << "  exp(t2) = " << exp(t2) << endl;
    cout << "exp(t1-t2) = " << exp(t1-t2) << endl;
    cout << "numNeighborsPivotBasis/numNeighborsPotentialBasis = " << (double)numNeighborsPivotBasis/numNeighborsPotentialBasis << endl;

    double prob1 = exp(-(F->Evaluate((Weight.subColumns(potentialBasis)).rowSum()))/T)*numNeighborsPivotBasis/(exp(-(F->Evaluate((Weight.subColumns(pivotBasis)).rowSum()))/T)*numNeighborsPotentialBasis); 
    cout << "prob1: " << prob1 << endl;

    if (prob1 > 1 || isnan(prob1)) // prob1 = min(prob1,1);
    {
        prob1 = 1;
    }
    cout << "New prob1: " << prob1 << endl;

    if ((double)rand() / (double)RAND_MAX < prob1)
    {
        return potentialBasis;
    }
    else
    {
        return pivotBasis;
    }
}
list <set <unsigned> > ProjBalMatroidOpt::SimulatedAnnealing(set <unsigned> firstBasis, list < double > temperatures, list < unsigned > times, list <set <unsigned> > &minBases)
{
    return SimulatedAnnealing(firstBasis,temperatures,times,minBases,BalanceFunction);
}

//list <set <unsigned> > ProjBalMatroidOpt::SimulatedAnnealing(set <unsigned> firstBasis, list < double > temperatures, list < unsigned > times, list <set <unsigned> > &minBases, double (*tempFct) (Matrix) )
list <set <unsigned> > ProjBalMatroidOpt::SimulatedAnnealing(set <unsigned> firstBasis, list < double > temperatures, list < unsigned > times, list <set <unsigned> > &minBases, Functional *F )
{
    if (temperatures.size () != times.size ())
    {
        cerr << "ProjBalMatroidOpt::SimulatedAnnealing called with temperatures and times of different size.\n";
        exit (0);
    }

    list <set <unsigned> > pivotBases;
    list <double>::const_iterator tempit = temperatures.begin();
    list <unsigned>::const_iterator timeit = times.begin();
    set <unsigned> currentBasis, currentMinBasis;
    double currentMin;

    currentBasis = firstBasis;

    currentMinBasis = currentBasis;
    currentMin = F->Evaluate((Weight.subColumns(currentMinBasis)).rowSum());
    minBases.push_back(currentMinBasis);

    // Go through all pairs of temperatures and times and pivot
    while (tempit != temperatures.end())
    {
        for (unsigned i=0;i<(*timeit);i++)
        {
            currentBasis = MetropolisBoltzmannUpdateFunction(currentBasis, *tempit, F);
            pivotBases.push_back(currentBasis);
            if (F->Evaluate((Weight.subColumns(currentBasis)).rowSum()) < currentMin)
            {
                currentMinBasis = currentBasis;
                currentMin = F->Evaluate((Weight.subColumns(currentMinBasis)).rowSum());
                minBases.push_back(currentMinBasis);
            }
        }
        
        tempit++;
        timeit++;
    }

    return pivotBases;
}


MinVarianceBalClustering::MinVarianceBalClustering ()
{

}

// This constructor reads in the data points. There is an implied uniform matroid structure.
// If n points are read in, there is a rank n/2 uniform matroid on n elements. 
MinVarianceBalClustering::MinVarianceBalClustering (std::istream &in)
{
    // Read in a weighting matrix which is the cluster points. 
    // Each col is a data point
    // Exit with error if number of points is not even.
    in >> Weight;
    if ((Weight.cols % 2) != 0)
    {
        cerr << "MinVarianceBalClustering Constructor: Number of data points not even." << endl;
        exit(0);
    }

    M = new UniformMatroid(Weight.cols/2,Weight.cols);
    // For our Variance functional we need the sum of all points
    //MinVarianceMatrix = Weight.rowSum();
    
    BalanceFunction = new FunctionalMinVariance(Weight.rowSum());
}

void MinVarianceBalClustering::printMathProg(std::ostream &o)
{
    o << "Minimal Variance Balanced Clustering problem." << endl;
    o << "Data points (columns)." << endl;
    o << Weight;
    o << *M;
    //string S("W");
    //Matrix M = Weight.transpose();
    //M.writeMatlab(cout,S);
}

void MinVarianceBalClustering::writeClustersMatlab(set <unsigned> S, std::ostream &o,string label)
{
    Matrix M1 = (Weight.subColumns(S)).transpose();
    M1.writeMatlab(o,label + "1");

    Matrix M2 = (Weight.subColumnsDiff(S)).transpose();
    M2.writeMatlab(o,label + "2");
}

TwoDAngleOrderedList::TwoDAngleOrderedList(set <unsigned> initialSet, ProjBalMatroidOpt *PBMO)
{
    refPBMO = PBMO;
    Matrix M = refPBMO->projectSet(initialSet);
    x = M(0,0);
    y = M(1,0);
    //cout << "Initial" << endl;
    //cout << M;
}

TwoDAngleOrderedList::~TwoDAngleOrderedList()
{
    refPBMO = 0;
}
void TwoDAngleOrderedList::insert( set <unsigned> S)
{
    Matrix M = refPBMO->projectSet(S);
    double tempx = M(0,0);
    double tempy = M(1,0);
    if (tempx == x && tempy == y)
    {
        return;
    }
    TwoDSetAngle tempTwoDSA;
    tempTwoDSA.S = S;
    tempTwoDSA.angle = atan2(tempy-y,tempx-x);
    //if (tempTwoDSA.angle < 0)
    //{
    //    tempTwoDSA.angle += M_PI*2;
    //}
    orderedRays.insert(tempTwoDSA); 
}

double TwoDAngleOrderedList::largestAngle (set <unsigned> &S1, set <unsigned> &S2)
{
    if (orderedRays.size() < 2)
    {
        //cout << "Size less than 2" << endl;
        return 0;
    }
    double tempAngle;
    if (orderedRays.size() == 2)
    {
        //cout << "Sized 2" << endl;
        set < TwoDSetAngle, ltTwoDSetAngle>::const_iterator torit = orderedRays.begin();
        S1 = (*torit).S;
        tempAngle = -((*torit).angle);
        torit++;
        S2 = (*torit).S;
        tempAngle += (*torit).angle;
        
        return tempAngle;
    }

    set < TwoDSetAngle, ltTwoDSetAngle>::const_iterator orit = orderedRays.begin();
    set < TwoDSetAngle, ltTwoDSetAngle>::const_iterator norit = orderedRays.begin();
    norit++;
    double curLargest = (*norit).angle - (*orit).angle;
    S1 = (*orit).S;
    S2 = (*norit).S;
    
    while (orit != orderedRays.end() && norit != orderedRays.end())
    {
        tempAngle = (*norit).angle - (*orit).angle; 
        if (tempAngle > curLargest)
        {
            curLargest = tempAngle;
            S1 = (*orit).S;
            S2 = (*norit).S;
            //cout << "NEW LARGEST" << endl;
            //cout << refPBMO->projectSet(S1);
            //cout << "----" << endl;
            //cout << refPBMO->projectSet(S2);
        }

        orit++;
        norit++;
    }

    // Now compare first and last;

    orit = orderedRays.begin();
    norit = orderedRays.end();
    norit--;
    tempAngle = 2*M_PI + (*orit).angle - (*norit).angle; 
    if (tempAngle > curLargest)
    {
        curLargest = tempAngle;
        S1 = (*orit).S;
        S2 = (*norit).S;
    }
    

    return curLargest;
}

void TwoDAngleOrderedList::print ()
{
    set < TwoDSetAngle, ltTwoDSetAngle>::const_iterator orit = orderedRays.begin();
    set < unsigned > S;
    int opl;

    while (orit != orderedRays.end())
    {
        S = (*orit).S;
        //copy(S.begin(), S.end(), ostream_iterator<const int>(cout, " "));
        //cout << endl;
        cout << "       Angle: " << (*orit).angle << endl; 
        opl = Matrix::printPadLength;
        Matrix::printPadLength=8;
        cout << refPBMO->projectSet(S);
        Matrix::printPadLength=opl;
        orit++;
    }
}

