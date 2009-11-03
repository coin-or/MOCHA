// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev$ $Date$
//


/// \file genmatrix.cpp

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
int main (int argc, char *argv[])
{

    if (argc != 5)
    {
        cout << "Usage: " << argv[0] << " <rows> <cols> <low> <high>" << endl;
        exit (0);
    }

    unsigned rows,cols,low,high;

    sscanf(argv[1],"%u",&rows);
    sscanf(argv[2],"%u",&cols);
    sscanf(argv[3],"%u",&low);
    sscanf(argv[4],"%u",&high);

    cout << rows << endl << cols << endl;
    Matrix M(rows,cols,low,high+1);
    cout << M;
}
