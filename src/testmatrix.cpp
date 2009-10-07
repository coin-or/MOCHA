// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev: 315 $ $Date: 2009-09-22 21:12:42 -0400 (Tue, 22 Sep 2009) $

/// \file testmatrix.cpp

#include <iostream>
#include "matrix.h"
#include <cstdlib>
#include <ctime>
#include <set>

using namespace::std;

int main ()
{
    if ( 0 == 1)
    {
        Matrix M(2,100,0,9);
        set<unsigned> someCols;

        someCols.insert(0);
        someCols.insert(3);
        someCols.insert(4);


        //M(1,1) = 5;
        //M(7,4) = 2;
        //M(2,3) = 14;
        //M(2,4) = 9;
        //M(2,5) = 3;
        //M(5,3) = 1;
        //M(0,9) = 4;
        //M(7,4) = 19;
        //M(8,8) = 73;


        cout << M;

        cout << "Copying M to A\n";

        Matrix A(M);

        cout << A;

        Matrix B;
        
        cout << "Created empty matrix B\n";
        cout << B;

        B = A;

        cout << B;

        //cout << "Trying GE\n";

        //Matrix *C;

        //C = A.GE();
        //A.GE();

        //cout << *C;
        //cout << "Rank of C: " << C->rank() << endl;
        cout << "A is unchanged\n";
        cout << A;
        cout << "Rank of A: " << A.rank() << endl;

        Matrix D;

        D = A.subColumns(someCols);

        cout << "Sub columns 0 3 4 of A" << endl;
        cout << D;

        Matrix E;

        E = D.GE();

        cout << "Rank of D=" << D.rank() << endl;
        cout << "D.GE()" << endl;
        cout << E;
    }


    // TEST 2

    if (0==1)
    {
        Matrix   A1(500,1000,0,50);

        cout << "Computing GE/rank of A1" << endl;

        cout << A1.rank() << endl;
    }

    if (0==1)
    {
        Matrix  F(3,5,0,2);
        Matrix  G(5,7,0,2);
    
        cout << F << G;

        Matrix  H;
        H = F * G;
        cout << H; 

    }

    if (0==1)
    {
        Matrix A(3,10,0,10);
        set<unsigned> someCols;
        
        cout << A;        
        cout << A.GE();
        srand(time(0));  
        while (someCols.size() < 3)
        {
            someCols.insert(rand() % 10);
        }
        
        Matrix B;
        B = A.subColumns(someCols);
        cout << B;
        cout << B.GE();
        cout << "Rank B:" << B.rank() << endl;
    }

    if (0==1)
    {
        Matrix A;

        cin >> A;

        cout << A;

        cout << "Rank A:" << A.rank() << endl;
        cout << "GE of A:" << endl << A.GE();

    }

    if (0==0)
    {
        Matrix A1(3,7,0,10);

        cout << A1;


#ifdef HAVE_LIBGMPXX
        cout << A1.GE_gmp() << endl;
#endif
    }
}
