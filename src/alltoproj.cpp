// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev: 315 $ $Date: 2009-09-22 21:12:42 -0400 (Tue, 22 Sep 2009) $
// This file is meant to convert the output of findChildren printing all spanning
// trees and project them by some given weighting.
// Expects input to be sets given by unsigned ints on each line.
// First argument should be the matrix that is our weighting.
// The second argument should be the number of elements in each set

#include "matrix.h"
#include <iostream>
#include <fstream>
#include <set>
#include <list>
#include <algorithm>
#include <iterator>
#include <string>
#include <ctime>


int main (int argc, char* argv[])
{
    if (argc != 3)
    {
        cout << "Usage: " << argv[0] << " <matrix file> <set cardinality> " << endl;
        exit(0);
    }

    ifstream inputfile(argv[1]);
    if (inputfile.fail() == 1)
    {
        cout << "Can not open file" << endl;
        exit (0);
    }


    Matrix Weight;
    inputfile >> Weight;
    unsigned setSize = atoi(argv[2]);
    unsigned tempui;
    set <unsigned> currSet;

    while (cin.eof() == 0)
    {
        currSet.clear();
        for (unsigned i=0;i<setSize;i++)
        {
            cin >> tempui;
            currSet.insert(tempui);
            if (cin.eof() == 1)
            {
                exit (0);
            }
        }
        //copy(currSet.begin(), currSet.end(), ostream_iterator<const int>(cout, " "));
        //cout << endl;
        cout << ((Weight.subColumns(currSet)).rowSum()).transpose();
    }

}

