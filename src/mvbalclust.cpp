// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

/// \file mvbalclust.cpp

#include "matroid.h"
#include "mathprog.h"
#include <iostream>
#include <fstream>
#include <set>
#include <list>
#include <algorithm>
#include <iterator>
#include <string>
#include <ctime>

using namespace::std;

int main (int argc, char *argv[])
{
    srand(time(0));

    ifstream inputFile;
    if (argc >= 2)
    {
        inputFile.open(argv[1]);
    }
    else {
        cout << "Filename of cluster points (points as columns): ";
        char  fileName[80];
        cin >> fileName;
        inputFile.open(fileName);
    }

    if (inputFile.fail() == 1)
    {
        cout << "Can not open file " << endl;
        exit (0);
    }
    unsigned pivotLimit;

    MinVarianceBalClustering MyCluster(inputFile);
    cout << MyCluster;

    int solveMethod = 0; // 0 - localsearch, 1 - tabusearch, 2 - boundary, 3 - firstcomefirstserve
    if (argc >= 3)
    {
        solveMethod = atoi(argv[2]);
    }
    else {
        cout << "Enter solve method. 0 (localsearch) 1 (tabusearch) 2 (boundary) 3 (first come first serve)" << endl;
        cin >> solveMethod;
    }
    time_t startTime;
    time_t endTime;

    if (solveMethod == 0){ // localsearch
        startTime = time(0);
        cout << "Performing Local Search." << endl;
        list <set <unsigned> > ls = MyCluster.LocalSearchRandomStart();
        cout << "Done." << endl;
        endTime = time(0);
        cout << "Total seconds: " << endTime - startTime << endl;


        cout << "Size of ls: " << ls.size () << endl;
       
        list <set <unsigned> >::const_iterator lsit = ls.end();

        lsit--;

        string CLLS("CLLS");

        MyCluster.writeClustersMatlab(*lsit,cout, CLLS);
        cout << "plot(CLLS1(:,1),CLLS1(:,2),'ro');" << endl;
    }
    else if (solveMethod == 1){ // tabusearch
        if (argc >= 4){
            pivotLimit = atoi(argv[3]);
        }
        else {
            cout << "Tabu pivot limit: ";
            cin >> pivotLimit;
        }
        list <set <unsigned> > ts;
        
        cout << "Performing Tabu Search with tabu limit " << pivotLimit << "." << endl;
        startTime = time(0);
        MyCluster.TabuSearchHeuristicRandomStart(pivotLimit,ts);
        cout << "Done." << endl;
        endTime = time(0);
        cout << "Total seconds: " << endTime - startTime << endl;
        cout << "Size of ts: " << ts.size () << endl;
        list <set <unsigned> >::const_iterator tsit = ts.end();

        tsit--;

        string CLTS("CLTS");

        MyCluster.writeClustersMatlab(*tsit,cout, CLTS);
        cout << "hold on;" << endl;
        cout << "plot(CLTS1(:,1),CLTS1(:,2),'ro');" << endl;
    }
    else if(solveMethod == 2){ // boundary
        string CLTS("CLTS");
        set <Matrix, ltcolvec> PBP;

        cout << "Performing Boundary Calculation." << endl;
        startTime = time(0);
        list <set <unsigned> > BP = MyCluster.Boundary(PBP);
        cout << "Done." << endl;
        endTime = time(0);
        cout << "Total seconds: " << endTime - startTime << endl;

        string CLCH("CLCH"); 
        MyCluster.writeClustersMatlab(MyCluster.FindMin(BP),cout, CLTS);
        cout << "hold on;" << endl;
        cout << "plot(CLCH1(:,1),CLCH1(:,2),'ro');" << endl;
    }
    if (solveMethod == 3){ // first come first serve
        startTime = time(0);
        cout << "Performing First Come First Serve." << endl;
        list <set <unsigned> > ls = MyCluster.FirstComeFirstServeRandomStart();
        cout << "Done." << endl;
        endTime = time(0);
        cout << "Total seconds: " << endTime - startTime << endl;


        cout << "Size of ls: " << ls.size () << endl;
       
        list <set <unsigned> >::const_iterator lsit = ls.end();

        lsit--;

        string CLLS("CLFCFS");

        MyCluster.writeClustersMatlab(*lsit,cout, CLLS);
        cout << "plot(CLFCFS1(:,1),CLFCFS1(:,2),'ro');" << endl;
    }
    else {
        cout << "solveMethod not correct. Must be 0,1,2 or 3." << endl;
        exit(1);
    }
}
