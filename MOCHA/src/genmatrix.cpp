// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev: 315 $ $Date: 2009-09-22 21:12:42 -0400 (Tue, 22 Sep 2009) $
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

    sscanf(argv[1],"%d",&rows);
    sscanf(argv[2],"%d",&cols);
    sscanf(argv[3],"%d",&low);
    sscanf(argv[4],"%d",&high);

    cout << rows << endl << cols << endl;
    Matrix M(rows,cols,low,high+1);
    cout << M;
}
