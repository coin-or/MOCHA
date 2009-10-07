// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev: 315 $ $Date: 2009-09-22 21:12:42 -0400 (Tue, 22 Sep 2009) $

/// \file graphtest.cpp

#include "graph.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <set>

using namespace::std;

int main ()
{
    Graph   G;

    //G.readGraph(cin);
    cin >> G;
    cout << G;
    cout << "Connected components" << endl;
    cout << G.connComp();
    set <unsigned> SF = G.randSpanningForest();
    copy(SF.begin(), SF.end(), ostream_iterator<const int>(cout, " "));
    cout << endl;
    
    int cyclefree = 0;
    G.randSpanningForest(cyclefree);

    cout << "Cycle free: " << cyclefree << endl;

    set <unsigned> S = G.randSpanningForest ();
    cout << "Random spanning forest: ";
    copy(S.begin(), S.end(), ostream_iterator<const int>(cout, " "));
    cout << endl;
    
    cout << "Is spanning tree: " << G.isSpanningForest(S) << endl;

    cout << "Initializing pivot" << endl;


    cout << "Floyd output" << endl;
    cout << G.FloydWarshall ();

    G.printVertexEdgeMatrix ();
}
