// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev: 315 $ $Date: 2009-09-22 21:12:42 -0400 (Tue, 22 Sep 2009) $
#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include "matrix.h"
#include <set>
#include <map>
#include <vector>
#include <list>


//set <unsigned> setUnionUnsigned (set <unsigned> &, set <unsigned> &);
//set <unsigned> setDifferenceUnsigned (set <unsigned> &, set <unsigned> &);


// The primary data structure of the graph will be its adjacency matrix
// Will not change underlying graph data. Using copying and functions
// for such things.

class VertexBucket
{
public:
    list <unsigned> vertices;
    int r;
    VertexBucket(int rvalue);
protected:
};

class Graph 
{
public:
    Graph ();
    ~Graph ();
    Graph (const Graph &G); 
    Graph (Matrix &M); // Creates a graph based off an adjacency matrix M
    Graph   &operator = (const Graph &G);
    Graph (std::istream &in);
     
    friend std::ostream &operator<< (std::ostream &o, const Graph &G);
    friend std::istream &operator>> (std::istream &in, Graph &G);

    // This is an implementation of Matsui's reverse-search spanning tree enumeration algorithm.
    // It takes as input a spanning tree  T and H(T) (see Matsui), and an index of the edges as 
    // a Matrix. We use a modified version for better space complexity.  
    // If edge e has incident vertices i and j, then (i,j) is the index of e.
    // The ordering of the edge index is assumed to be given by Nagamachi and Ibaraki (Allison's part).
    // It is a BFS procedure, so it will recurse on itself.
    // Currently it will print out all spanning trees to standard output.
    // G is the underlying graph, initTree is the first tree given by order edgeIndex.
    // The current tree we care about is initTree\g union f.
    // initTree, f, g are all edges indexed according to G.
    // printMod specifies the mod value to print tree count. 0 none
    friend void findChildren(Graph &G, set <unsigned> &initTree, set <unsigned> &deltaf, set <unsigned> &deltag, set <unsigned> deltaH, Matrix &edgeIndex, unsigned printMod, unsigned printTrees);

    friend void findChildren(Graph &G, set <unsigned> &initTree, set <unsigned> &deltaf, set <unsigned> &deltag, set <unsigned> deltaH, Matrix &edgeIndex, Matrix &Weight, set <Matrix, ltcolvec> &projTrees, unsigned printMod, unsigned printTrees);


    // This will return the ``bottom'' edge of this graph as given by edgeIndex
    // Returns edge with lowest index. Returns this graphs edge index, NOT index provided.
    unsigned MatsuiBottom(Matrix &edgeIndex);
    unsigned MatsuiBottom(Matrix &edgeIndex, set <unsigned> &someEdges);

    // This will return the ``top'' edge of this graph as given by edgeIndex
    // Returns edge with highest index. Returns this graphs edge index, NOT index provided.
    unsigned MatsuiTop(Matrix &edgeIndex);
    unsigned MatsuiTop(Matrix &edgeIndex, set <unsigned> &someEdges);

    // returns the edge index of someEdge given by the matrix edgeIndex
    unsigned getEdgeIndex(unsigned someEdge, Matrix &edgeIndex);

    // Returns a subgraph with same number of nodes but only
    // Those edges in someEdges
    Graph   subGraph(set <unsigned> &someEdges);
    Graph   subGraph(list <vector <unsigned> > &someEdges);

    // Returns a subgraph with same number of nodes but only
    // Those edges NOT in someEdges
    Graph   subGraphDiff(set <unsigned> &someEdges);
    Graph   subGraphDiff(list <vector <unsigned> > &someEdges);

    // This returns a matrix of same size of adjMatrix where
    // (i,j) = 1 if there exists a path between i and j
    // Computes the transitive closure
    // Saves an internal copy of transitive closure
    Matrix  connComp ();

    // Returns 1 if the two nodes are connected on this graph.
    // Uses transitive closure and computes it if first call
    int nodesConnected(int, int);

    void    readGraph (std::istream &in);

    // Returns the edges of a random spanning forest
    // Uses BFS
    // This also computes this graphs number of connected components
    set <unsigned> randSpanningForest();

    // Same as above but it will set cycleFree to 1 if
    // the underlying graph is cycle free. 0 else. 
    set <unsigned> randSpanningForest(int &cycleFree);

    // Returns 1 if someEdges is a spanning forest. 0 Else.
    int isSpanningForest (list <vector <unsigned> > someEdges);

    // Returns 1 if someEdges is a spanning forest. 0 Else.
    int isSpanningForest (set <unsigned> someEdges);

    // Returns 1 if someEdges is cycle free. 0 Else.
    int isCycleFree (set <unsigned> someEdges);

    // Returns the shortest path (as edge subset)
    // Uses FloydWarshall. Computes predMatrix needed.
    set <unsigned> shortestPath (unsigned n1, unsigned n2);

    // Returns the shortest path (as edge subset)
    // Uses FloydWarshall. Computes predMatrix needed.
    list <vector <unsigned> > shortestPathList (unsigned n1, unsigned n2);

    // Converts a list of edges given by end nodes to a set of their edge numbers
    set <unsigned> listEdgesToSet (list <vector <unsigned> > &);

    int getNumConnComponents ();
    int rank();
    int getNumEdges ();
    void deleteDoubleArray(int **, int);

    //  This performs Floyd-Warshall's algorithm to find all pairs shortest
    //paths. It returns the predecessor matrix 
    Matrix FloydWarshall ();

    // Nagamochi Ibariki Algorithm. To be implemented by Allison
    list <set <unsigned> > NagIbar();

    // This takes the output of NagIbar and returns a matrix which indexes
    // the edges given by the input. The index of each set is orbitally chosen
    Matrix calcEdgeIndex(list <set <unsigned> > &);

    // This compares if a <= b + c where we assume a,b,c are non-negative and -1 means infinity
    // Returns 1 if a <= b + c, 0 else. If a and (b or c) is infinity return -1
    int     leftLessThanEqualRight (double a, double b, double c);

    void edgeToNode(unsigned someEdge, unsigned &, unsigned &);

    void static edgeToNodeEdgeIndex(unsigned someEdge, unsigned &, unsigned &, Matrix &edgeIndex);

    // Returns edge number between nodes (a,b)
    // Returns -1 if (a,b) is not an edge
    unsigned    edgeNumber(int a, int b);

    // Print out the vertex edge matrix. Prints with arbitrary directions to edges
    void printVertexEdgeMatrix();

    unsigned findChildrenSpanningTreeCount;
    unsigned findChildrenBFSLevel;
protected:
    // Returns a random adjacent node to startNode.
    // If no node exists, returns startNode
    //unsigned randAdjNode(unsigned startNode);
    void    initGraph ();
    // The edges are ordered from upper left to lower right as read on the adjMatrix
    Matrix  adjMatrix;
    Matrix  nodesToEdgeNumber; // This simplifies retrieving edge numbers based on nodes
    map <unsigned, vector <unsigned> >  edges;
    int     numEdges;
    int     numNodes;
    int     numConnComponents;
    Matrix  transClosure;
    int     transClosureComputed;
    int     predMatrixComputed;
    Matrix  predMatrix;
};
#endif
