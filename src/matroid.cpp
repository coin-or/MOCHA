// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev: 315 $ $Date: 2009-09-22 21:12:42 -0400 (Tue, 22 Sep 2009) $
#include "matroid.h"
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <iterator>

#define ZERO_PIVOT_THRESHOLD 0.0000001

// LAPACK declaration
extern "C" void dgesv_(int *, int *, double *, int *, int [], double *, int *, int *);

Matroid::Matroid ()
{
    matroidRank = -1;
    numElements=-1; 
}

// Do nothing. Should not be called
Matroid::Matroid (istream &in)
{
    matroidRank = -1;
    numElements=-1; 
}

Matroid::~Matroid ()
{
}

set <Matrix, ltcolvec> Matroid::calcAllBasesProj(Matrix &Weight)
{
    set <Matrix, ltcolvec> dummyReturn;
    return dummyReturn;
}

set <unsigned> Matroid::GreedyAlgorithmMax(Matrix Weight)
{
    if (Weight.rows != getNumElements() || Weight.cols != 1)
    {
        cerr << "Matroid::GreedyAlgorithmMax Weight not the appropriate dimension." << endl;
        exit (0);
    }

    set <unsigned> maxBasis;  // We will grow this one element at a time.
    set <unsigned> validWeights;
    set <unsigned>::iterator suit;
    double maxWeight;
    unsigned maxWeightIndex;

    for (unsigned i=0;i<getNumElements();i++)
    {
        validWeights.insert(i);
    }

    while(maxBasis.size() != rank())
    {
        suit = validWeights.begin();
        maxWeight = Weight(*suit,0);
        maxWeightIndex = *suit;
        suit++;
        for (;suit != validWeights.end();suit++)
        {
            if (maxWeight < Weight(*suit,0))
            {
                    maxWeight = Weight(*suit,0);
                    maxWeightIndex = *suit;
            }
        }
        maxBasis.insert(maxWeightIndex);
        if (maxBasis.size() != setRank(maxBasis))
        { // Can't add maxWeightIndex then. Remove it
            maxBasis.erase(maxBasis.find(maxWeightIndex));
        }
        validWeights.erase(validWeights.find(maxWeightIndex));

    }
    return maxBasis;
}

float Matroid::random_weight_logistic(float x){
 
  float log_rand;

  log_rand=log(x)-log(1-x);
  return(log_rand);
}

// h(t) for the logistic distribution
float Matroid::H_fn_logistic(float t){
  float Hoft;
  float interval_range;
  float left_point;
  float right_point;
  float mid_point;
  float d1, value1, d2, value2, max_value;

  left_point=0;
  right_point=1-SOLVER_ACCURACY;
  interval_range=2;     
  
  while(interval_range>SOLVER_ACCURACY){
    mid_point=(left_point+right_point)/2;
    d1=(left_point+mid_point)/2;
    d2=(right_point+mid_point)/2;
    value1=t*d1-log(PI*d1/(sin(PI*d1)));
    value2=t*d2-log(PI*d2/(sin(PI*d2)));
         if(value1>value2){
	   right_point=mid_point;
	   max_value=value1;
	 }
	 else{
	   left_point=mid_point;
	   max_value=value2;
	 }
	 interval_range=right_point-left_point;
  }
  Hoft= max_value;
  return(Hoft);
}

float Matroid::upper_logistic(float avg_GAMMA, int k){

  float upper_ans;
  
  upper_ans=avg_GAMMA;
  return(upper_ans);
}

float Matroid::lower_logistic(float avg_GAMMA, int k){

  float lower_ans;
 
  lower_ans=(k-1)*H_fn_logistic((float)avg_GAMMA/(k-1));
  return(lower_ans);
}

// Modification of rand(): since problems occur when
// rand() returns 0 or MAXIMAL_RAND
int Matroid::modified_rand(){

  int temp;

  temp=rand();
  if(temp==0){
    temp=temp+1;
  }
  else if(temp==MAXIMAL_RAND){
    temp=MAXIMAL_RAND-1;
  }
  return(temp);
}



//void Matroid::EstimateBases()
//{
//	int m, num_rows, num_cols, i, j, k, max_weight;
//	float upper_bound, lower_bound, GAMMA, avg_GAMMA;
//	string filename;
//	ifstream inFile;
//	set <unsigned> max_weight_basis;
//
//
//	cout << "This program estimates the number of bases of a matrix by assigning\n";
//	cout << "random weights to each column of the matrix and finding the maximal\n";
//	cout << "weight basis. This is done m times and the resulting average estimates\n";
//	cout << "the function GAMMA. Then upper and lower bounds for the number of bases\n";
//	cout << "are calculated.\n\n";
//
//	cout << "Enter the number of times m to sample: \n";
//	cin >> m;
//	cout << "\n";
//	while(m<=0){
//		cout << "This number must be greater than zero.\n";
//		cout << "Enter the number of times m to sample: \n";
//		cin >> m;
//		cout << "\n";
//	}
//
//	cout << "Enter the number of rows in the matrix.\n";
//	cin >> num_rows;
//	cout << "\n";
//
//	cout << "Enter the number of columns in the matrix.\n";
//	cout << "This must be more than the number of rows.\n";
//	cin >> num_cols;
//
//	int matrix[num_rows][num_cols];
//	Matrix weights = Matrix(1,num_cols);
//
//	cout << "Enter the name of the input file:\n";
//	cin >> filename;
//
//	inFile.open(filename.c_str(), ios::in);
//
//	while(!inFile){
//		cout << "File open error: '" << filename << "' ";
//		cout <<"Try again.\n";
//		inFile.close();
//		cout << "Enter the name of the input file:\n";
//		cin >> filename;
//		inFile.open(filename.c_str(), ios::in);
//	}
//
//	// read in the matrix from inFile
//	for(i=0; i<num_rows; i++){
//		for(j=0; j<num_cols; j++){
//			inFile >> matrix[i][j];
//		}
//	}
//
//	GAMMA=0;
//	srand((unsigned)time(NULL));
//	for(k=0;k<m;k++){
//		max_weight_basis.clear();
//  
//		for(i = 0; i < num_cols; i++){
//		weights(0,i)=0;
//		}
//		for(i = 0; i < num_cols; i++){
//			weights(0,i)=random_weight_logistic((float)modified_rand()/MAXIMAL_RAND);
//		}
//  
//		max_weight_basis = GreedyAlgorithmMax(weights);
//  
//		set<unsigned>::iterator it;
//		max_weight = 0;
//		for(it = max_weight_basis.begin(); it != max_weight_basis.end(); it++)
//		{
//			max_weight += *it;
//		}
//  
//		GAMMA=GAMMA+max_weight;
//  
//	}
//
//	//compute upper and lower estimates
//	avg_GAMMA = (float)GAMMA/m;
//	upper_bound = upper_logistic(avg_GAMMA, num_rows);
//	lower_bound = lower_logistic(avg_GAMMA, num_rows);
//
//	//Output results
//	cout << "\n";
//	cout << "\n";
//	cout << "Bounds for Log(X) are:\n";
//	cout << "Upper bound:\n";
//	printf("%f\n",upper_bound);
//	cout << "Lower bound:\n";
//	printf("%f\n",lower_bound);
//	cout << "GAMMA(X):\n";
//	printf("%f\n",avg_GAMMA);
//
//
//}

set <Matrix, ltcolvec> GraphicalMatroid::calcAllBasesProj(Matrix &Weight)
{
    set <Matrix, ltcolvec> projPoints;

    list <set <unsigned> > edgePartition = graphRep.NagIbar(); 

    list <vector <unsigned> > edgePartitionList;

    list <set <unsigned> >::iterator lsit = edgePartition.begin();

    Graph tempG = graphRep;

    lsit = edgePartition.begin();
    int i = 0;
    while (lsit != edgePartition.end() && (*lsit).size() > 0)
    {
        set <unsigned> someEdges = *lsit;
        cout << "Edge set " << i << endl;
        i++;

        copy(someEdges.begin(), someEdges.end(), ostream_iterator<const int>(cout, " "));
        cout << endl; 

        edgePartitionList.clear();
        unsigned n1,n2;
        vector <unsigned> someVec;
        // Create list of edges given by nodes due to indexing problem
        set <unsigned>::const_iterator sit = someEdges.begin();
        for (;sit != someEdges.end();sit++)
        {
            graphRep.edgeToNode(*sit,n1,n2);
            someVec.clear();
            someVec.push_back(n1);
            someVec.push_back(n2);
            edgePartitionList.push_front(someVec);
        }

        if (tempG.isSpanningForest(edgePartitionList) == 1)
        {
            cout << "Yes it is a spanning forest." << endl;
        }
        else 
        {
            cout << "No it is not a spanning forest." << endl;
        }
        //cout << tempG;

        tempG = tempG.subGraphDiff(edgePartitionList);

        lsit++;
        cout << endl;
    }

    cout << "Edge index: " << endl;
    cout << "calculating edgeIndex." << endl;
    Matrix edgeIndex = graphRep.calcEdgeIndex(edgePartition);
    cout << "Printing edge index." << endl;
    cout << edgeIndex;
    cout << "done Printing edge index." << endl;
    lsit = edgePartition.begin();

    set <unsigned> temps1, temps2, temps3;
    //cout << &edgeIndex << endl;
    graphRep.findChildrenSpanningTreeCount = 0;
    graphRep.findChildrenBFSLevel = 0;
    time_t startTime = time(0);
    findChildren(graphRep,*lsit, temps1, temps2, temps3, edgeIndex, Weight, projPoints, 100000,0);
    time_t endTime = time(0);

    return projPoints;
}


VectorMatroid::VectorMatroid()
{
    thisMatroidType = VECTOR_MATROID;
}

VectorMatroid::~VectorMatroid() 
{

}

GraphicalMatroid::GraphicalMatroid() 
{
    thisMatroidType = GRAPH_MATROID;

}

GraphicalMatroid::~GraphicalMatroid() 
{

}

std::ostream& operator<< (std::ostream& o, Matroid &someMatroid)
{ 
    someMatroid.printMatroid(o); 
    return o;
}

std::istream& operator>> (std::istream& in, Matroid &someMatroid)
{
    someMatroid.getMatroid(in);
    return in;
}


void VectorMatroid::printMatroid (std::ostream &o)
{
    o << "VectorMatroid     Rank: " << rank () << endl;
    //rank ();
    Matrix::printPadLength=4;
    o << "  Matrix:\n";
    o << matrixRep;
    //o << "  GE\n";
    //o << matrixRep.GE();
    Matrix::printPadLength=0;
}

void VectorMatroid::getMatroid(std::istream &in)
{
    in >> matrixRep; 
    matroidRank = matrixRep.rank ();
    numElements = matrixRep.getCols ();
}

VectorMatroid::VectorMatroid(std::istream &in)
{
    thisMatroidType = VECTOR_MATROID;
    getMatroid(in);
    //in >> matrixRep;
    //matroidRank = matrixRep.rank ();
    //numElements = matrixRep.getCols (); 
}

VectorMatroid::VectorMatroid(Matrix M)
{
    thisMatroidType = VECTOR_MATROID;
    matrixRep = M;
    matroidRank = matrixRep.rank ();
    numElements = matrixRep.getCols (); 
}

int VectorMatroid::rank ()
{
    if (matroidRank == -1)
    {
        matroidRank = matrixRep.rank();
    }
    return matroidRank;
}

int VectorMatroid::isBasis (set <unsigned> S)
{
    if (numElements == 0)
    {
        return 0;
    }
    if (S.size() > numElements)
    {
        return 0;
    }
    
    //set <int>::iterator Siter;
    //Siter = S.begin();

    Matrix M = matrixRep.subColumns(S);

    if ((M.rank()) == rank())
    {
        return 1;
    } 

    return 0;
}

int VectorMatroid::setRank (set <unsigned> S) //Returns the rank of the set
{
    if (numElements == 0)
    {
        return 0;
    }
    if (S.size() > numElements)
    {
        return 0;
    }
    
    Matrix M = matrixRep.subColumns(S);

    return M.rank();
}

set <unsigned> VectorMatroid::randomBasis () //Returns a random basis
{
    if (numElements == 0)
    {
        set <unsigned> S;
        return S;
    }

    set <unsigned> S;
    Matrix M;
    
    if (rank() == 0)
    {
        cerr << "VectorMatroid::randomBasis called with rank()==0" << endl;
        exit (0);
    }

    while ((M.rank()) != (rank()))
    {
        S.clear();
        while ((int)S.size() < rank())
        {
            S.insert(rand() % numElements);
        }
        if (S.size() == 0)
        {
            cerr << "S.size() = " << S.size() << endl;
            exit(1);
        }

        //copy(S.begin(), S.end(), ostream_iterator<const int>(cout, " "));

        M = matrixRep.subColumns(S); 
        if (S.size() == 0)
        {
            cerr << "After subColumns called, S.size() = " << S.size() << endl;
            exit(1);
        }
        //cout << "Submatrix" << endl;
        //cout << M;
        //cout << "GE" << endl;
        //Matrix::printPadLength=4;
        //cout << M.GE();
        //Matrix::printPadLength=0;
        //cout << "Rank " << M.rank() << "\n\n";
    }
    if (S.size() == 0)
    {
        cerr << "S.size() = " << S.size() << endl;
        cerr << "VectorMatroid::randomBasis S.size()==0" << endl;
        cerr << "rank() = " << rank() << endl;
        cerr << "matrixRep.rank() = " << matrixRep.rank() << endl;
        cerr << "matrixRep" << endl;
        cerr << matrixRep;
        cerr << "M.rank() = " << M.rank() << endl;
        cerr << "M" << endl;
        cerr << "M.rows = " << M.rows << endl;
        cerr << "M.cols = " << M.cols << endl;
        cerr << M ;
        exit (0);
    }

    return S;
}

int VectorMatroid::calcBases ()
{
    // First need rank
    if (matroidRank == -1)
    {
        rank ();
    }
    
    // Now go through all matroidRank sized subsets of 
    // Columns of matrixRep and see if they are rank matroidRank
    return -1;
}

void GraphicalMatroid::printMatroid (std::ostream &o)
{
    o << graphRep;
}

// This will initialize all data for enumerating
// Neighbors of initBasis. This is so specialized 
// techniques can be used for graphical matroids for example
// This uses no specialized technique. Simply checks all (i,j) and
// see if initbasis - i + j is a basis
void VectorMatroid::initializePivotGeneric(set <unsigned> initBasis)
{
    if (isBasis(initBasis) != 1)
    {
        cerr << "VectorMatroid::initializePivot called with non basis." << endl;

        cerr << "size of initBasis = " << initBasis.size() << endl;
        copy(initBasis.begin(), initBasis.end(), ostream_iterator<const int>(cout, " "));
        exit (0);
    }
    set <unsigned>::iterator si;
    int inSet;

    pivotBasis = initBasis;
    remSet.clear();
    addSet.clear();

    remSet = initBasis;
    // This is simple set difference [n] - initBasis
    for (unsigned i=0;i<(unsigned)getNumElements();i++)
    {
        inSet = 0;
        si = remSet.begin();
        remSetCount= 0;
        while (remSetCount < remSet.size())
        {
            if (*si == i)
            {
                inSet = 1;
            } 
            si++;
            remSetCount++;
        }
        if (inSet == 0)
        {
            addSet.insert(i);
        }
    }

    remSetCount = 0;
    addSetCount = 0;
    remEl = remSet.begin();
    addEl = addSet.begin();
}
// This will set &pivot to the next pivot from initBasis
// Returns 1 if &pivot is set, i.e. there is a next pivot
// Returns 0 else. A pivot is only valid if it is a basis
int VectorMatroid::nextPivotGeneric(set <unsigned> &pivot)
{
    //if (remSetCount >= remSet.size())

    int pivotIsBasis = 0;
    
    while (pivotIsBasis == 0)
    {
        if (remEl == remSet.end())
        {
            return 0;
        }
        pivotIsBasis = 0;

        si = pivotBasis.find(*remEl);
        if (si == pivotBasis.end())
        {
            cerr << "nextPivot: Can not find element to remove. Problem with pivotBasis. Probably not initialized." << endl;
            exit (0);
        }
        pivotBasis.erase(si);
        pivotBasis.insert(*addEl);   

        if (isBasis(pivotBasis) == 1)
        {
            pivot = pivotBasis;
            pivotIsBasis = 1;
        }

        si = pivotBasis.find(*addEl);
        pivotBasis.erase(si);
        pivotBasis.insert(*remEl);

        addEl++;
        addSetCount++;
        if (addSetCount < addSet.size())
        {
            // only return if we found a basis
            if (pivotIsBasis == 1)
            {
                return 1;
            }            
        }
        else
        {
            addSetCount = 0;
            addEl = addSet.begin ();
            remEl++;
            remSetCount++; 
            // only return if we found a basis
            if (pivotIsBasis == 1)
            {
                return 1;
            }            
        }
    }
    return 0;
}
void VectorMatroid::initializePivot(set <unsigned> initBasis)
{
    return initializePivotGeneric(initBasis);
    //return initializePivotLAPACK(initBasis);
}

void VectorMatroid::initializePivotLAPACK(set <unsigned> initBasis)
{
    cout << "initBasisLAPACK: " << endl;
    //copy(initBasis.begin(), initBasis.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    if (isBasis(initBasis) != 1)
    {
        cerr << "VectorMatroid::initializePivot called with non basis." << endl;
        cerr << "initBasis: ";
        copy(initBasis.begin(), initBasis.end(), ostream_iterator<const int>(cerr, " "));
        cerr << endl;
        cerr << "isBasis(initBasis) = " << isBasis(initBasis) << endl;
        cerr << matrixRep.subColumns(initBasis);
        exit (0);
    }
    //Graph   tempG = graphRep.subGraph(initBasis);

    //predMatrix = tempG.FloydWarshall();
    //cout << "Pred matrix" << endl;
    //cout << predMatrix;

    set <unsigned>::iterator si;
    int inSet;

    pivotBasis = initBasis;
    currentCycle.clear();
    addSet.clear();

    // This is simple set difference [n] - initBasis
    for (unsigned i=0;i<(unsigned)getNumElements();i++)
    {
        inSet = 0;
        si = initBasis.begin();
        remSetCount= 0;
        while (remSetCount < initBasis.size())
        {
            if (*si == i)
            {
                inSet = 1;
            } 
            si++;
            remSetCount++;
        }
        if (inSet == 0)
        {
            addSet.insert(i);
        }
    }
    //cout << "numElements = " << getNumElements() << endl;
    //cout << "addSet(" << addSet.size() << ")   ";
    //copy(addSet.begin(), addSet.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;


    remSetCount = 0;
    addSetCount = 0;
    addEl = addSet.begin();

    //unsigned n1, n2;
    //graphRep.edgeToNode(*addEl, n1, n2);

    //currentCycle = getCycle (n1, n2);
    //cout << "currentCycle(" << currentCycle.size() << ")   ";
    //copy(currentCycle.begin(), currentCycle.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
   
    // This is where we need to set currentCycle to all the x_i in Bx=A_j that are non-zero
    // where A_j is the element from the add set
    
    currentInitBasis = matrixRep.subColumns(initBasis);
    double AT[currentInitBasis.rows*currentInitBasis.cols];
    double b[currentInitBasis.rows];

    // Translate into column major order for FORTRAN LAPACK function
    for(unsigned i=0;i<currentInitBasis.cols;i++)
    {
        for(unsigned j=0;j<currentInitBasis.rows;j++)
        {
            AT[j+currentInitBasis.rows*i] = currentInitBasis(j,i);
        }
    }
    set <unsigned> oneCol;
    oneCol.insert(*addEl);
    Matrix tempM = matrixRep.subColumns(oneCol);
    for(unsigned j=0;j<currentInitBasis.rows;j++)
    {
        b[j] = tempM(j,0);
    }

    int INFO, LDA, LDB, N, NRHS;
    INFO = -1;
    int IPIV[currentInitBasis.rows];
    LDA = currentInitBasis.rows;
    LDB = currentInitBasis.rows;
    N = currentInitBasis.rows;
    NRHS = 1;

    dgesv_(&N,&NRHS,AT,&LDA,IPIV,b,&LDB,&INFO);

    if (INFO != 0)
    {
        cerr << "VectorMatroid::initializePivot  INFO returned as " << INFO << " from dgesv_" << endl;
        exit (0);
    }

    double XMAX = -1; // This will be the max over fabs(x_i)
    for(unsigned j=0;j<currentInitBasis.rows;j++)
    {
        if (fabs(b[j]) > XMAX)
        {
            XMAX = fabs(b[j]);
        }
    }
    if (XMAX == -1)
    {
        cerr << "XMAX == -1" << endl;
        exit(0);
    }
    XMAX = max((double)1,XMAX);
    //cout << "XMAX = " << XMAX << endl;

    set <unsigned>::const_iterator suci= initBasis.begin();
    // Now go through and identify all entries of x_i where |x_i| > 10^-16*XMAX
    for(unsigned i=0;i<currentInitBasis.rows;i++)
    {
        //if (fabs(b[i]) > XMAX*0.000000000000001 )
        if (fabs(b[i]) > ZERO_PIVOT_THRESHOLD )
        //if (fabs(b[i]) > XMAX*0.000000000001 )
        {
            // Need to add the element of the ground set [n] corresponding to column i to currentCycle
            currentCycle.insert(*suci); 
            //cout << "   Added " << *suci << endl;
        }
        suci++;
    }
    //cout << "done" << endl;
    //copy(currentCycle.begin(), currentCycle.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    //cout << "currentCycle.size() = " << currentCycle.size() << endl;

    

    remEl = currentCycle.begin();
    //cout << "Initialize done." << endl;
}

int VectorMatroid::nextPivot (set <unsigned> &pivot)
{
    return nextPivotGeneric(pivot);
    //return nextPivotLAPACK(pivot);
}

// This will set &pivot to the next pivot from initBasis
// Returns 1 if &pivot is set, i.e. there is a next pivot
// Returns 0 else. A pivot is only valid if it is a basis
int VectorMatroid::nextPivotLAPACK (set <unsigned> &pivot)
{
    cout << "nextPivot" << endl;
    //if (addSetCount >= addSet.size())
    if (addEl == addSet.end())
    {
        cout << "addEl == addSet, returning 0" << endl;
        return 0;
    }

    // Because we know all paths for the pivot basis we can for sure find
    // the next pivot basis.

    // Each newly added edge will create a cycle! Thus, we need only
    // add that edge and go through and remove each other edge in the cycle
    // one by one.

    // We will allow currentCycle to be empty

    si = pivotBasis.find(*remEl);
    //if (si == pivotBasis.end())
    //{
    //    cerr << "nextPivot(VectorMatroid): Can not find element to remove. Problem with pivotBasis. Probably not initialized." << endl;
    //    cerr << "*remEl = " << *remEl << endl;
    //    exit (0);
    //}
    int pivotFound = 0;
    if (si != pivotBasis.end() && currentCycle.size() != 0)
    {
        pivotBasis.erase(si);
        pivotBasis.insert(*addEl);   

        //cout << "   Current cycle ";
        //copy(currentCycle.begin(), currentCycle.end(), ostream_iterator<const int>(cout, " "));
        //cout << endl;
        //cout << "   Add element " << *addEl << endl;
        //cout << "   Rem element " << *remEl << endl;

        pivot = pivotBasis;
        // This is to check our work. We can remove this once we have 
        // verified the validity of the specialize approach
        if (isBasis(pivot) == 0)
        {
            cerr << "VectorMatroid::nextPivot pivot is not a basis." << endl;
            cerr << "*remEl = " << *remEl << endl;
            cerr << "*addEl = " << *addEl << endl;
            cerr << "pivotBasis ";
            cerr << "matrixRep.subColumns(pivotBasis)" << endl;
            cerr << "rank = " << matrixRep.subColumns(pivotBasis).rank() << endl;
            cerr << matrixRep.subColumns(pivotBasis);
            copy(pivotBasis.begin(), pivotBasis.end(), ostream_iterator<const int>(cerr, " "));
            cerr << endl;
            exit(0);
        }
        cout << "       Pivot found" << endl;
        pivotFound = 1;

        si = pivotBasis.find(*addEl);
        pivotBasis.erase(si);
        pivotBasis.insert(*remEl);

        remEl++;
    }
    else
    {
        cout << "       Pivot NOT found!" << endl;
    }

    if (remEl == currentCycle.end() || currentCycle.size() == 0)
    {
        currentCycle.clear();
        cout << "   remEl == currentCycle.end()" << endl;
        addEl++;
        if (addEl == addSet.end())
        {
            if (pivotFound == 1)
            {
                cout << "       return 1" << endl;
                return 1; // If we got to here, than return 1
            }
            else
            {
                cout << "       return 0" << endl;
                return 0;
            }
        }
        double AT[currentInitBasis.rows*currentInitBasis.cols];
        double b[currentInitBasis.rows];

        // Translate into column major order for FORTRAN LAPACK function
        for(unsigned i=0;i<currentInitBasis.cols;i++)
        {
            for(unsigned j=0;j<currentInitBasis.rows;j++)
            {
                AT[j+currentInitBasis.rows*i] = currentInitBasis(j,i);
            }
        }
        set <unsigned> oneCol;
        oneCol.insert(*addEl);
        Matrix tempM = matrixRep.subColumns(oneCol);
        for(unsigned j=0;j<currentInitBasis.rows;j++)
        {
            b[j] = tempM(j,0);
        }
        cout << "A" << endl;
        cout << currentInitBasis;
        cout << "b" << endl;
        cout << tempM << endl;


        int INFO, LDA, LDB, N, NRHS;
        INFO = -1;
        int IPIV[currentInitBasis.rows];
        LDA = currentInitBasis.rows;
        LDB = currentInitBasis.rows;
        N = currentInitBasis.rows;
        NRHS = 1;

        dgesv_(&N,&NRHS,AT,&LDA,IPIV,b,&LDB,&INFO);
        //cout << "x" << endl;
        //set <unsigned>::const_iterator suci2= pivotBasis.begin();
        //for(unsigned j=0;j<currentInitBasis.rows;j++)
        //{
        //    cout << b[j] << "   " << *suci2 << endl;
        //    suci2++;
        //}

        if (INFO != 0)
        {
            cerr << "VectorMatroid::initializePivot  INFO returned as " << INFO << " from dgesv_" << endl;
            exit (0);
        }

        double XMAX = -1; // This will be the max over fabs(x_i)
        for(unsigned j=0;j<currentInitBasis.rows;j++)
        {
            if (fabs(b[j]) > XMAX)
            {
                XMAX = fabs(b[j]);
            }
        }
        if (XMAX == -1)
        {
            cerr << "XMAX == -1" << endl;
            exit(0);
        }
        XMAX = max((double)1,XMAX);


        cout << "XMAX = " << XMAX << endl;
        cout << "x" << endl;
        cout.precision(20);
        set <unsigned>::const_iterator suci= pivotBasis.begin();
        // Now go through and identify all entries of x_i where |x_i| > 10^-16*XMAX
        for(unsigned i=0;i<currentInitBasis.rows;i++)
        {
            cout << b[i];
            //if (fabs(b[i]) > XMAX*0.000000000000001 )
            if (fabs(b[i]) > ZERO_PIVOT_THRESHOLD)
            //if (fabs(b[i]) > XMAX*0.000000000001 )
            {
                // Need to add the element of the ground set [n] corresponding to column i to currentCycle
                currentCycle.insert(*suci); 
                cout << "   " << *suci;
            }
            cout << endl;
            suci++;
        }
        remEl = currentCycle.begin();

    }

    if (pivotFound == 1)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void GraphicalMatroid::getMatroid(std::istream &in)
{
    in >> graphRep; 
    matroidRank = graphRep.rank ();
    numElements = graphRep.getNumEdges ();
}

GraphicalMatroid::GraphicalMatroid(std::istream &in)
{
    thisMatroidType = GRAPH_MATROID;
    getMatroid(in);
}

GraphicalMatroid::GraphicalMatroid(Graph G)
{
    thisMatroidType = GRAPH_MATROID;
    graphRep = G;
    matroidRank = graphRep.rank ();
    numElements = graphRep.getNumEdges ();
}

int GraphicalMatroid::rank ()
{
    if (matroidRank == -1)
    {
        matroidRank = graphRep.rank();
    }
    return matroidRank;
}

int GraphicalMatroid::isBasis (set <unsigned> someElements)
{
    if (numElements == 0)
    {
        return 0;
    }
    if (someElements.size() > numElements)
    {
        return 0;
    }
    return graphRep.isSpanningForest(someElements);
}

set <unsigned> GraphicalMatroid::randomBasis ()
{
    return graphRep.randSpanningForest ();
}

void GraphicalMatroid::initializePivot(set <unsigned> initBasis)
{
    //cout << "initBasis: ";
    //copy(initBasis.begin(), initBasis.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    if (isBasis(initBasis) != 1)
    {
        cerr << "GraphicialMatroid::initializePivot called with non basis." << endl;
        exit (0);
    }
    Graph   tempG = graphRep.subGraph(initBasis);

    predMatrix = tempG.FloydWarshall();
    //cout << "Pred matrix" << endl;
    //cout << predMatrix;

    set <unsigned>::iterator si;
    int inSet;

    pivotBasis = initBasis;
    currentCycle.clear();
    addSet.clear();

    // This is simple set difference [n] - initBasis
    for (unsigned i=0;i<(unsigned)getNumElements();i++)
    {
        inSet = 0;
        si = initBasis.begin();
        remSetCount= 0;
        while (remSetCount < initBasis.size())
        {
            if (*si == i)
            {
                inSet = 1;
            } 
            si++;
            remSetCount++;
        }
        if (inSet == 0)
        {
            addSet.insert(i);
        }
    }
    //cout << "numElements = " << getNumElements() << endl;
    //cout << "addSet(" << addSet.size() << ")   ";
    //copy(addSet.begin(), addSet.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;


    remSetCount = 0;
    addSetCount = 0;
    addEl = addSet.begin();

    unsigned n1, n2;
    graphRep.edgeToNode(*addEl, n1, n2);

    currentCycle = getCycle (n1, n2);
    //cout << "currentCycle(" << currentCycle.size() << ")   ";
    //copy(currentCycle.begin(), currentCycle.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;

    remEl = currentCycle.begin();
    //cout << "Initialize done." << endl;
}

// This will set &pivot to the next pivot from initBasis
// Returns 1 if &pivot is set, i.e. there is a next pivot
// Returns 0 else. A pivot is only valid if it is a basis
int GraphicalMatroid::nextPivot (set <unsigned> &pivot)
{
    //if (addSetCount >= addSet.size())
    if (addEl == addSet.end())
    {
        return 0;
    }

    // Because we know all paths for the pivot basis we can for sure find
    // the next pivot basis.

    // Each newly added edge will create a cycle! Thus, we need only
    // add that edge and go through and remove each other edge in the cycle
    // one by one.

    si = pivotBasis.find(*remEl);
    if (si == pivotBasis.end())
    {
        cerr << "nextPivot(GraphicialMatroid): Can not find element to remove. Problem with pivotBasis. Probably not initialized." << endl;
        exit (0);
    }
    pivotBasis.erase(si);
    pivotBasis.insert(*addEl);   

    //cout << "   Current cycle ";
    //copy(currentCycle.begin(), currentCycle.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    //cout << "   Add element " << *addEl << endl;
    //cout << "   Rem element " << *remEl << endl;

    pivot = pivotBasis;

    si = pivotBasis.find(*addEl);
    pivotBasis.erase(si);
    pivotBasis.insert(*remEl);

    remEl++;

    if (remEl == currentCycle.end())
    {
        addEl++;
        if (addEl == addSet.end())
        {
            return 1; // If we got to here, than return 1
        }
        unsigned n1, n2;
        graphRep.edgeToNode(*addEl, n1, n2);

        currentCycle = getCycle (n1, n2);
        remEl = currentCycle.begin();
    }

    return 1;
}

int GraphicalMatroid::setRank (set <unsigned> S) //Returns the rank of the set
{
    if (numElements == 0)
    {
        return 0;
    }
    if (S.size() > numElements)
    {
        return 0;
    }
    Graph tempG = graphRep.subGraph(S);

    return tempG.rank();

    return 0;
}


// This functions returns data with respect to the pivot basis
// This return the cycle formed by adding the edge (n1,n2)
// Since pivot basis is a spanning forest, there will always
// be a cycle if (n1,n2) is NOT in pivotBasis.
set <unsigned> GraphicalMatroid::getCycle (unsigned n1, unsigned n2)
{
    set <unsigned> newCycle;

    unsigned tn1 = n1;
    unsigned tn2 = n2;

    while ((int)predMatrix(tn1,tn2) != -1)
    {
    
        newCycle.insert(graphRep.edgeNumber((unsigned)predMatrix(tn1,tn2),tn2));
        tn2 = (unsigned)predMatrix(tn1,tn2);
    }

    return newCycle;
}

UniformMatroid::UniformMatroid ()
{
    thisMatroidType = UNIFORM_MATROID;
}

UniformMatroid::UniformMatroid (int rank, int elements)
{
    thisMatroidType = UNIFORM_MATROID;
    if (elements >= 0)
    {
        numElements = elements;
    }
    else
    {
        cerr << "UniformMatroid constructor called with elements < 0" << endl;
        exit (0);
    }

    if (rank >= 0 && rank <= elements)
    {
        matroidRank = rank;
    }
    else
    {
        cerr << "UniformMatroid constructor called with rank < 0 or rank > elements." << endl;
        exit (0);
    }
}

int UniformMatroid::isBasis (set <unsigned> S)
{
    if (S.size() <= matroidRank)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int UniformMatroid::setRank (set <unsigned> S)
{
    if (S.size() <= matroidRank)
    {
        return S.size();
    }
    else
    {
        return matroidRank;
    }
}

set <unsigned> UniformMatroid::randomBasis ()
{
    set <unsigned> S;

    while (S.size() < matroidRank)
    {
        S.insert(rand() % numElements);
    }
    return S;
}

void UniformMatroid::initializePivot(set <unsigned> initBasis)
{
    if (isBasis(initBasis) != 1)
    {
        cerr << "UniformMatroid::initializePivot called with non basis." << endl;
        exit (0);
    }
    set <unsigned>::iterator si;
    int inSet;

    pivotBasis = initBasis;
    remSet.clear();
    addSet.clear();

    remSet = initBasis;
    // This is simple set difference [n] - initBasis
    for (unsigned i=0;i<(unsigned)getNumElements();i++)
    {
        inSet = 0;
        si = remSet.begin();
        remSetCount= 0;
        while (remSetCount < remSet.size())
        {
            if (*si == i)
            {
                inSet = 1;
            } 
            si++;
            remSetCount++;
        }
        if (inSet == 0)
        {
            addSet.insert(i);
        }
    }

    remSetCount = 0;
    addSetCount = 0;
    remEl = remSet.begin();
    addEl = addSet.begin();
}

int UniformMatroid::nextPivot(set <unsigned> &pivot)
{

    int pivotIsBasis = 0;
    
    while (pivotIsBasis == 0)
    {
        if (remEl == remSet.end())
        {
            return 0;
        }
        pivotIsBasis = 0;

        si = pivotBasis.find(*remEl);
        if (si == pivotBasis.end())
        {
            cerr << "nextPivot: Can not find element to remove. Problem with pivotBasis. Probably not initialized." << endl;
            exit (0);
        }
        pivotBasis.erase(si);
        pivotBasis.insert(*addEl);   

        if (isBasis(pivotBasis) == 1)
        {
            pivot = pivotBasis;
            pivotIsBasis = 1;
        }

        si = pivotBasis.find(*addEl);
        pivotBasis.erase(si);
        pivotBasis.insert(*remEl);

        addEl++;
        addSetCount++;
        if (addSetCount < addSet.size())
        {
            // only return if we found a basis
            if (pivotIsBasis == 1)
            {
                return 1;
            }            
        }
        else
        {
            addSetCount = 0;
            addEl = addSet.begin ();
            remEl++;
            remSetCount++; 
            // only return if we found a basis
            if (pivotIsBasis == 1)
            {
                return 1;
            }            
        }
    }
    return 0;
}

void UniformMatroid::printMatroid (std::ostream &o) 
{
    o << "Uniform Matroid" << endl;
    o << "Rank: " << matroidRank << "   Elements: " << numElements << endl;
}

void UniformMatroid::getMatroid (std::istream &in)
{

}
