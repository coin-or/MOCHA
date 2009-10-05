// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev: 315 $ $Date: 2009-09-22 21:12:42 -0400 (Tue, 22 Sep 2009) $
// This program reads in from stdin points, and finds the pareto optimum using 
// a straightforward search

#include "matrix.h"
#include "mathprog.h"
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
    if (argc != 2)
    {
        cout << "Usage: " << argv[0] << " <dim> " << endl;
        exit(0);
    }
    unsigned dim;
    dim=atoi(argv[1]);

    ProjBalMatroidOpt MyOpt;
    Matrix  tempM(dim,1);
    int tempInt;
    set <Matrix, ltcolvec> allPoints;

    while(!cin.eof())
    {
        for (unsigned i=0;i<dim;i++)
        {
            cin >> tempInt;
            if (cin.eof())
            {
               break;
            }
            tempM(i,0) = tempInt;
        }
        allPoints.insert(tempM);
    }
    //cout << "Size allPoints.insert() = " << allPoints.size() << endl;
    set <Matrix, ltcolvec> paretoOpt = MyOpt.ParetoOptimum(allPoints);
    set <Matrix, ltcolvec>::iterator smit = paretoOpt.begin();

    for (;smit!=paretoOpt.end();smit++)
    {
        for (unsigned j=0;j<(*smit).rows;j++)
        {
            cout << (*smit)(j,0);
            if (j < (*smit).rows - 1)
            {
                cout << " ";
            }
        }
        cout << endl;
    }
}
