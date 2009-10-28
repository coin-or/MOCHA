// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev: 315 $ $Date: 2009-09-22 21:12:42 -0400 (Tue, 22 Sep 2009) $

/// \file nagibatest.cpp

#include "graph.h"
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
    ifstream inputfile;
    char fileName[80];
    int printMod = 10000;
    int printTrees = 0;

    if (argc >= 2)
    {
        inputfile.open(argv[1]);
    }
    else {

        cout << "Input file: ";
        cin >> fileName;

        inputfile.open(fileName);
    }

    if (inputfile.fail() == 1)
    {
        cout << "Can not open file" << endl;
        exit (0);
    }

    if (argc >= 3)
    {
        sscanf(argv[2],"%d",&printMod);
    }
    if (argc >= 4)
    {
        sscanf(argv[3],"%d",&printTrees);
        if (printTrees > 1)
        {
            printTrees = 0;
        }
    }

    Graph G(inputfile);
    cout << G;

    list <set <unsigned> > edgePartition = G.NagIbar(); 

    list <vector <unsigned> > edgePartitionList;

    list <set <unsigned> >::iterator lsit = edgePartition.begin();

    Graph tempG = G;

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
            G.edgeToNode(*sit,n1,n2);
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
    Matrix edgeIndex = G.calcEdgeIndex(edgePartition);
    cout << "Printing edge index." << endl;
    cout << "done Printing edge index." << endl;
    lsit = edgePartition.begin();

    set <unsigned> temps1, temps2, temps3;
    //cout << &edgeIndex << endl;
    G.findChildrenSpanningTreeCount = 0;
    G.findChildrenBFSLevel = 0;
    time_t startTime = time(0);
    findChildren(G,*lsit, temps1, temps2, temps3, edgeIndex,printMod,printTrees);
    time_t endTime = time(0);
    cout << "Total seconds: " << endTime - startTime << endl;
    cout << "Number of spanning trees: " << G.findChildrenSpanningTreeCount << endl;
}
