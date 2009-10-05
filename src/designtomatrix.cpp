// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// This program will read in two matrices; a n x k exponent vector matrix B, 
// a m x k design point matrix P
//
// It then computes A_ij := \prod_{h=1}^k P_{i,h}^B_{j,h} 
// It reads in the matrices from standard input and prints the matrix to 
// standard output
#include"matrix.h"
#include<cmath>


int main (int argc, char *argv[])
{
    Matrix B, P;
    double tempPow;


    cin >> B;
    cin >> P;

    Matrix A(P.rows,B.rows);

    for (int i=0;i<A.rows;i++)
    {
        for (int j=0;j<A.cols;j++)
        {
            A(i,j) = 1;
            for (int l=0;l<B.cols;l++)
            {
                tempPow = pow(P(i,l),B(j,l));
                // I don't think pow is returning 0^0 as 1! It should!
                //if (P(i,l) == 0 && B(j,l) == 0)
                //{
                //    tempPow = 1;
                //}
                //if (P(i,l) == 0)
                //{
                //    tempPow = 1;
                //}

                A(i,j) *= tempPow;
            }
        }
    }

    cout << "VECTOR" << endl;
    cout << A.rows << endl << A.cols << endl;
    cout << A;
    cout << B.cols << endl << B.rows << endl;
    cout << B.transpose();

}
