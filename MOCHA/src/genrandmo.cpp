// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev: 315 $ $Date: 2009-09-22 21:12:42 -0400 (Tue, 22 Sep 2009) $

/// \file genrandmo.cpp

#include "matroid.h"
#include "mathprog.h"
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <set>
#include <list>
#include <algorithm>
#include <iterator>
#include <string>
#include <cmath>

using namespace::std;

Matrix  *monomialMatrix;
int curCol = 0;
void nestedForRec(Matrix &M,int low, int high, int curRow)
{
    for (int i=low;i<=high;i++)
    {
        M(curRow,0)=i;
        if (curRow < (int)M.rows - 1)
        {
            nestedForRec(M,low,high,curRow+1);
        }
        else if (curRow == (int)M.rows-1)
        {
            for (int j=0;j<(int)M.rows;j++)
            {
                (*monomialMatrix)(j,curCol) = M(j,0);
            }
            curCol++;
        }
    }
}

int main ()
{
    int matroidType = -1;
    int tempInt;
    int tempInt2;
    int tempInt3;
    int tempInt4;
    char outFileName[80];
     

    cout << "Output file name? " << endl;
    cin >> outFileName;
    
    string S;
    ofstream of;
    of.open(outFileName);

    
    while (matroidType == -1)
    {
        cout << "VECTOR(1) GRAPH (2) VECTOR_ALL_TUPLES(3)?" << endl;
        cin >> tempInt;
        if (tempInt == 1 || tempInt == 2 || tempInt == 3)
        {
            matroidType = tempInt;
        }
    }

    int rows;
    int cols;
    if (matroidType == 1)
    {
        of << "VECTOR" << endl;

        cout << "Rows: ";
        cin >> rows;
        cout << "Cols: ";
        cin >> cols;
        of << rows << endl << cols << endl;
        cout << "Matrix Random Integer Low: ";
        cin >> tempInt;
        cout << "Matrix Random Integer High: ";
        cin >> tempInt2;
        Matrix M(rows,cols,tempInt,tempInt2+1);
        of << M; 
        cout << "Random weighting matrix? " ;
        cin >> S;
        if(S == "Y" || S == "y")
        {
            cout << "Number Rows Weighting Matrix: ";
            cin >> tempInt;
            cout << "Weighting Matrix Random Integer Low: ";
            cin >> tempInt2;
            cout << "Weighting Matrix Random Integer High: ";
            cin >> tempInt3;
            Matrix W(tempInt,cols,tempInt2,tempInt3+1);
            of << tempInt << endl << cols << endl;
            of << W;
        }
    }
    else if (matroidType == 2)
    {
        of << "GRAPH" << endl << "ADJACENCY" << endl;
        cout << "Number of nodes: ";
        cin >> rows;
        cout << "Number of edges: ";
        cin >> cols;
        Matrix  M(rows);
        for (int i=0;i<rows;i++)
        {
            M(i,i) = 1;
        }
        int numRand = 0;

        int row2, col2;
        srand((unsigned int)time(0));
        cout << "Generating random edges" << endl;
        while (numRand < cols)
        {
            row2 = rand () % rows; 
            col2 = rand () % (rows - row2);
            col2 += row2;
            //cout << "   row " << row2 << "  col " << col2 << endl;
            if (M(row2,col2) == 0)
            {
                M(row2,col2) = 1;
                M(col2,row2) = 1;
                numRand++;
            }
        }
        of << rows << endl << rows << endl;
        of << M;
        cout << "Random weighting matrix? " ;
        cin >> S;
        if(S == "Y" || S == "y")
        {
            cout << "Number Rows Weighting Matrix: ";
            cin >> tempInt;
            cout << "Weighting Matrix Random Integer Low: ";
            cin >> tempInt2;
            cout << "Weighting Matrix Random Integer High: ";
            cin >> tempInt3;
            Matrix W(tempInt,cols,tempInt2,tempInt3+1);
            of << tempInt << endl << cols << endl;
            of << W;
        }

    }
    else if (matroidType == 3)
    {
        int mlow,mhigh, wcols;
        cout << "Low range: " ;
        cin >> tempInt;
        cout << "High range: ";
        cin >> tempInt2;
        cout << "Dimension: ";
        cin >> tempInt3;
        cout << "Number of rows of A: ";
        cin >> tempInt4;
        wcols = (unsigned)pow((double)(tempInt2+1-tempInt),tempInt3);
        cout << "Creating matrix with of dimensions " << tempInt4 << " x " << wcols << endl;

        cout << "Matrix Random Integer Low: ";
        cin >> mlow;
        cout << "Matrix Random Integer High: ";
        cin >> mhigh;
        Matrix M(tempInt4,wcols,mlow,mhigh);
        of << "VECTOR" << endl;
        of << tempInt4 << endl << wcols << endl;
        of << M; 
        monomialMatrix = new Matrix(tempInt3,wcols);
        Matrix tempM(tempInt3,1);
        nestedForRec(tempM,tempInt,tempInt2,0);
        
        of << monomialMatrix->rows << endl << monomialMatrix->cols << endl;
        of << (*monomialMatrix);

    }

}
    
