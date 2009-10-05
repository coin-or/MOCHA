// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev: 315 $ $Date: 2009-09-22 21:12:42 -0400 (Tue, 22 Sep 2009) $
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>

using namespace std;

int det (int **M,int row,int col)
{
    if (row != col) {
        fprintf(stderr,"row not equal to column in det\n");
        exit (0);
    }
    
    
}

int main ()
{
    int Arows, Acols, Arank;

    cin >> Arows;
    cin >> Acols;
    cin >> Arank;

    // Create A matrix
    //cout << "Matrix has " << Arows << " rows and " << Acols << " cols and " << Arank << " rank\n"; 
    int A[Arows][Acols];

    for (int i=0;i<Arows;i++){
        for (int j=0;j<Arows;j++){
            cin >> A[i][j];
        }
    }

    //cout << "Matrix:\n";
    //for (int i=0;i<Arows;i++){
    //    for (int j=0;j<Arows;j++){
    //        cout << A[i][j] << " ";
    //    }
    //    cout << "\n";
    //}

    //Find a random starting matroid basis       

    int B[Arows][Arank];

}
