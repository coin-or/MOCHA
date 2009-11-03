// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev$ $Date$


/// \file graph.cpp
#include "graph.h"
#include <iostream>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <deque>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <iterator>

VertexBucket::VertexBucket(int rvalue) : r(rvalue)
{
  
}

void Graph::initGraph ()
{
    numEdges = 0;
    numNodes = 0;
    numConnComponents = 0;
    transClosureComputed = 0;
    predMatrixComputed = 0;
}

Graph::Graph ()
{
}

Graph::~Graph ()
{
}

Graph::Graph (const Graph &G)
{
    adjMatrix = G.adjMatrix;
    nodesToEdgeNumber = G.nodesToEdgeNumber;
    edges = G.edges;
    numEdges = G.numEdges;
    numNodes = G.numNodes;
    numConnComponents = G.numConnComponents;
}

Graph::Graph (std::istream &in)
{
    readGraph(in);
}

Graph &Graph::operator = (const Graph &G)
{
    adjMatrix = G.adjMatrix;
    nodesToEdgeNumber = G.nodesToEdgeNumber;
    edges = G.edges;
    numEdges = G.numEdges;
    numNodes = G.numNodes;
    numConnComponents = G.numConnComponents;
    predMatrixComputed = G.predMatrixComputed;
    transClosureComputed = G.transClosureComputed;
    return *this; 
}

std::ostream &operator<< (std::ostream &o, const Graph &G)
{
    o << "GRAPH     Nodes: " << G.numNodes << "   Edges: " << G.numEdges << endl;
    o << "Number of connected components: " << G.numConnComponents << endl;
    o << "Adjacency matrix:" << endl << G.adjMatrix;
    map <unsigned, vector <unsigned> >::const_iterator EI = (G.edges).begin();
    //vector <unsigned> V;
    //o << "Edges:" << endl;
    //for (int i=0;i<G.numEdges;i++)
    //{
    //    V = (*EI).second;
    //    o << "[" << i << "]  " << V[0] << " " << V[1] << endl;
    //    EI++; 
    //}
    o << "nodesToEdgeNumber:" << endl;
    o << G.nodesToEdgeNumber;
    //o << "predMatrixComputed = " << G.predMatrixComputed << endl;
    //o << "nodesToEdgeNumber:" << endl;
    //o << G.nodesToEdgeNumber;
    return o;
}

std::istream &operator >> (std::istream &in, Graph &G)
{
    G.readGraph(in);
    return in;
}

void Graph::readGraph (std::istream &in)
{
    string S;
    unsigned count = 0;
    vector <unsigned> oneEdge;

    in >> S;
    if (S == "ADJACENCY")
    {
        in >> adjMatrix;
        numNodes = adjMatrix.getRows(); 
        numEdges = 0;
        Matrix  nten(numNodes);
        for (int i=0;i<numNodes;i++)
        {
            for (int j=0;j<numNodes;j++)
            {
                nten(i,j)=-1;
            }
        }
        for (int i=0;i<numNodes;i++)
        {
            for (int j=i+1;j<numNodes;j++)
            {
                if (((adjMatrix(i,j)) != 0))
                {
                    nten(i,j) = numEdges; 
                    oneEdge.clear();
                    oneEdge.push_back(i);  
                    oneEdge.push_back(j);  
                    edges[numEdges] = oneEdge; 
                    numEdges++;
                }
                else {
                    //oneEdge.clear();
                    //edges[count] = oneEdge; 
                }

                count++;
            }
        }
        nodesToEdgeNumber = nten;
        randSpanningForest (); //Run this to calculate the number of connected components
    }
    else if (S == "INCIDENCE")
    {
        // I dont handle this yet
    }
    else 
    {
        return;
    }
}

Matrix  Graph::connComp ()
{
    Matrix  M = adjMatrix;

    for(int i = 0; i < numNodes; i++)
    {
        for(int j = 0; j < numNodes; j++)
        {
            if((M(i,j)) == 1)
            {
                for(int k = 0; k < numNodes; k++)
                {
                    if((M(i,k)) == 1)
                    {
                        M(j,k) = 1;
                    }
                }
            }
        }
    } 

    if (transClosureComputed == 0)
    {
        transClosure = M;
        transClosureComputed = 1;
    }
    return M;
}

set <unsigned> Graph::randSpanningForest ()
{
    int dontCare;
    return randSpanningForest(dontCare);
}

set <unsigned> Graph::randSpanningForest (int &cycleFree)
{
    set <unsigned> SF;
    unsigned startNode;
    unsigned randIndex;
    unsigned nodesVisited=0;
    int u;
    deque <unsigned> nodeQueue;
    deque <unsigned>::iterator nqi;

    // Pick random starting node
    startNode = rand() % numNodes;
    //cout << "startNode = " << startNode << endl;

    // Mark each element white=0
    vector <unsigned> color(numNodes,0);
    vector <int> pred(numNodes,-1);

    // Mark the startNode as grey=1
    color[startNode] = 1;

    nodeQueue.push_back(startNode);
    numConnComponents=0;
    numConnComponents++;

    cycleFree = 1; // Initialize. If we see a vertex with non-white color, not cycle free
    while (nodeQueue.size() != 0)
    {
        nqi = nodeQueue.begin();
        randIndex = (rand() % (unsigned)nodeQueue.size());
        u = *(nqi + randIndex);
        nodeQueue.erase((nqi + randIndex)); 
        //cout << "   u = " << u << endl;
        // Iteratre through all adjacent nodes of u
        for (int i=0;i<numNodes;i++)
        {
            if (adjMatrix(u,i) == 1 && u != i)
            {
                if (color[i] != 0 && (pred[u] != i || pred[u] == -1))
                {
                    //cout << "       Detected a cycle! u=" << u << "  i=" << i << endl; 
                    //cout << "       color[i]=" << color[i] << endl;
                    //cout << "       pred[i] = " << pred[i] << endl;
                    //cout << "       pred[u] = " << pred[u] << endl;
                    // We have detected a cycle!
                    cycleFree = 0;
                }
                if (color[i] == 0)
                {
                    //cout << "       marking " << i << " as grey(1)." << endl;
                    pred[i] = u;
                    color[i] = 1;
                    nodeQueue.push_back(i);
                    // Add edges (u,i)
                    SF.insert(edgeNumber(u,i));
                }
            }
        }
        color[u]=2;
        nodesVisited++;

        if (nodeQueue.size () == 0 && (int)nodesVisited < numNodes)
        {
            // We havn't visited all nodes
            for (int i=0;i<numNodes;i++)
            {
                if (color[i] == 0)
                {
                    nodeQueue.push_back(i);
                    i=numNodes;    
                } 
            }
            numConnComponents++;
        }
    }

    return SF;
}

unsigned    Graph::edgeNumber(int a, int b)
{
    //cout << "(a,b) = (" << a << "," << b << ")" << endl;
    unsigned    edgeNum = 0;
    int x,y;
    if (a == b)
    {
        cerr << "Graph::edgeNumber called with a==b." << endl;
        exit (0);
    }
    if ( a < b )
    {
        x = a;
        y = b;
    }
    else
    {
        x = b;
        y = a;
    }
    
    //cout << "nodesToEdgeNumber:" << endl;
    //cout << nodesToEdgeNumber;
    return (unsigned)nodesToEdgeNumber(x,y);

    //// Lazy calculation
    //for (int i=0;i<x;i++)
    //{
    //    edgeNum += numNodes - (1 + i);
    //}
    //edgeNum += y - (x + 1);

    return edgeNum;
}

int Graph::isCycleFree (set <unsigned> someEdges)
{
    Graph   G;

    G = subGraph(someEdges);
    //cout << "Edges of subgraph: ";
    //copy(someEdges.begin(), someEdges.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    //cout << G;
    

    int cycleFree;

    G.randSpanningForest(cycleFree);

    if (cycleFree == 1)
    {
        return 1;
    }
    
    return 0; // Error handling
}

int Graph::isSpanningForest (list <vector <unsigned> > someEdges)
{
    Graph   G;

    G = subGraph(someEdges);
    //cout << "Edges of subgraph: ";
    //copy(someEdges.begin(), someEdges.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    //cout << G;
    

    int cycleFree;

    G.randSpanningForest(cycleFree);

    if (cycleFree == 1 && (int)someEdges.size () == rank())
    {
        return 1;
    }
    
    return 0; // Error handling
}

int Graph::isSpanningForest (set <unsigned> someEdges)
{
    Graph   G;

    G = subGraph(someEdges);
    //cout << "Edges of subgraph: ";
    //copy(someEdges.begin(), someEdges.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    //cout << G;
    

    int cycleFree;

    G.randSpanningForest(cycleFree);

    if (cycleFree == 1 && (int)someEdges.size () == rank())
    {
        return 1;
    }
    
    return 0; // Error handling
}

Graph Graph::subGraph(list <vector <unsigned> > &someEdges)
{
    Matrix newAdjMatrix(numNodes);
    for (int i=0;i<numNodes;i++)
    {
        newAdjMatrix(i,i) = 1;
    }

    list <vector <unsigned> >::const_iterator lvi;
    vector <unsigned> V;
    lvi = someEdges.begin();
    while (lvi != someEdges.end())
    {
        V = *lvi;
        if (V.size() != 0)
        {
            //cout << "(*sei): " << (*sei) << ". V[0]=" << V[0] << "  V[1]=" << V[1] << endl;
            newAdjMatrix(V[0],V[1])=1;
            newAdjMatrix(V[1],V[0])=1;
        }

        lvi++;
    }
    //cout << "   Subgraph" << endl;
    //cout << newAdjMatrix;
    Graph   newGraph(newAdjMatrix); 
    return newGraph;
}

Graph Graph::subGraph(set <unsigned> &someEdges)
{
    Matrix newAdjMatrix(numNodes);
    for (int i=0;i<numNodes;i++)
    {
        newAdjMatrix(i,i) = 1;
    }

    set <unsigned>::const_iterator sei;
    vector <unsigned> V;
    sei = someEdges.begin();
    while (sei != someEdges.end())
    {
        V = edges[(*sei)];
        if (V.size() != 0)
        {
            //cout << "(*sei): " << (*sei) << ". V[0]=" << V[0] << "  V[1]=" << V[1] << endl;
            newAdjMatrix(V[0],V[1])=1;
            newAdjMatrix(V[1],V[0])=1;
        }

        sei++;
    }
    //cout << "   Subgraph" << endl;
    //cout << newAdjMatrix;
    Graph   newGraph(newAdjMatrix); 
    return newGraph;
}

Graph Graph::subGraphDiff(list <vector <unsigned> >&someEdges)
{
    //cout << "subGraphDiff(list <vector <unsigned> > & called" << endl;
    Matrix newAdjMatrix = adjMatrix;

    list <vector <unsigned> >::const_iterator lvi;
    vector <unsigned> V;
    lvi = someEdges.begin();
    while (lvi != someEdges.end())
    {
        V = *lvi;
        if (V.size() == 2)
        {
            //cout << "(*sei): " << (*sei) << ". V[0]=" << V[0] << "  V[1]=" << V[1] << endl;
            newAdjMatrix(V[0],V[1])=0;
            newAdjMatrix(V[1],V[0])=0;
        }

        lvi++;
    }
    //cout << "   Subgraph" << endl;
    //cout << newAdjMatrix;
    Graph   newGraph(newAdjMatrix); 
    //cout << "subGraphDiff(list <vector <unsigned> > & done" << endl;
    return newGraph;
}

Graph Graph::subGraphDiff(set <unsigned> &someEdges)
{
    set <unsigned> setDiff;
    set <unsigned>::const_iterator sit = someEdges.begin();

    for (unsigned i=0;i<(unsigned)numEdges;i++)
    {
        if (i == *sit && sit != someEdges.end())
        {
            sit++;
        }
        else
        {
            setDiff.insert(i);
        }
    }

    //cout << "   someEdges: ";
    //copy(someEdges.begin(), someEdges.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    //cout << "   setDiff  : ";
    //copy(setDiff.begin(), setDiff.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;

    return subGraph(setDiff);
}

Graph::Graph (Matrix &M)
{
    initGraph();
    unsigned count = 0;
    vector <unsigned> oneEdge;

    adjMatrix = M;
    numNodes = adjMatrix.getRows(); 
    numEdges = 0;
    Matrix  nten(numNodes);
    for (int i=0;i<numNodes;i++)
    {
        for (int j=0;j<numNodes;j++)
        {
            nten(i,j)=-1;
        }
    }
    for (int i=0;i<numNodes;i++)
    {
        for (int j=i+1;j<numNodes;j++)
        {
            if (i != j && ((adjMatrix(i,j)) != 0))
            {
                oneEdge.clear();
                oneEdge.push_back(i);  
                oneEdge.push_back(j);  
                edges[count] = oneEdge; 
                numEdges++;
                nten(i,j) = count; 
                count++;
            }
        }
    }
    nodesToEdgeNumber = nten;
    randSpanningForest (); //Run this to calculate the number of connected components
}

int Graph::getNumConnComponents ()
{
    return numConnComponents;
}

int Graph::rank()
{
    return numNodes - numConnComponents;
}

int Graph::getNumEdges ()
{
    return numEdges;
}

/*  This performs Floyd-Warshall's algorithm to find all pairs shortests
    paths. It returns the predecessor matrix */
Matrix Graph::FloydWarshall (){
    Matrix  returnMatrix;
    Matrix  *predMatrix = new Matrix(numNodes);
    Matrix  *distMatrix = new Matrix(numNodes);// We assume all distances are 1 
    Matrix  *tempDistMatrix = new Matrix(numNodes);
    Matrix  *tempPredMatrix = new Matrix(numNodes);

    Matrix *tempMatrix;

    int     i,j,k;

    //printf("FloydWarshall called\n");

    for (i=0;i<numNodes;i++){
        for (j=0;j<numNodes;j++) {
            if (i == j) {
                (*distMatrix)(i,j) = 0;
            }
            else {
                if (adjMatrix(i,j) == 1) {
                    (*distMatrix)(i,j) = adjMatrix(i,j);
                }
                else {
                    // -1 means infinity
                    (*distMatrix)(i,j) = -1;
                }
            }

            if (i != j && adjMatrix(i,j) == 1) {
                (*predMatrix)(i,j) = i+1;
            }
            else {
                (*predMatrix)(i,j) = 0;
            }
            (*tempPredMatrix)(i,j) = 0;
            (*tempDistMatrix)(i,j) = 0;
        }
    }

    // tempPredMatrix will hold the kth iteration. predMatrix will hold the k-1 iteration 
    // tempDistMatrix will hold the kth iteration. distMatrix will hold the k-1 iteration 
    for(k=0;k<numNodes;k++){
        for (i=0;i<numNodes;i++) {
            for (j=0;j<numNodes;j++) {
                if (leftLessThanEqualRight((*distMatrix)(i,j),(*distMatrix)(i,k),(*distMatrix)(k,j))) {
                    (*tempPredMatrix)(i,j) = (*predMatrix)(i,j);
                    (*tempDistMatrix)(i,j) = (*distMatrix)(i,j);
                }
                else {
                    (*tempPredMatrix)(i,j) = (*predMatrix)(k,j);
                    (*tempDistMatrix)(i,j) = (*distMatrix)(i,k) + (*distMatrix)(k,j);
                }
            }
        }
        tempMatrix = predMatrix;
        predMatrix = tempPredMatrix;
        tempPredMatrix = tempMatrix;

        tempMatrix = distMatrix;
        distMatrix = tempDistMatrix;
        tempDistMatrix = tempMatrix;
    }

    delete tempPredMatrix;
    delete tempDistMatrix;
    delete distMatrix;


    // Need to correct for proper indexing
    for (int i=0;i<predMatrix->getRows();i++)
    {
        for (int j=0;j<predMatrix->getRows();j++)
        {
            (*predMatrix)(i,j) = (*predMatrix)(i,j) - 1;
        } 
    }
    returnMatrix = *predMatrix;
    delete predMatrix;
    
    return returnMatrix;
}

void    Graph::deleteDoubleArray(int **someArray, int numRows){
    int i;

    for(i=0;i<numRows;i++){
        delete[] someArray[i];
    }
    delete[] someArray;
}
                
// This compares if a <= b + c where we assume a,b,c are non-negative and -1 means infinity
// Returns 1 if a <= b + c, 0 else. If a and (b or c) is infinity return -1
int     Graph::leftLessThanEqualRight (double a, double b, double c) {
    //Handle the cases
    if ( (a == -1) && ( (b == -1) || (c == -1) ) ) {
        return -1;
    }

    if ( (a == -1) && ( (b != -1) && (c != -1) ) ) {
        return 0;
    }

    if ( (a != -1) && ( (b == -1) || (c == -1) ) ) {
        return 1;
    }

    if ( (a != -1) && (b != -1) && (c != -1) ) {
        if ( a <= b + c) {
            return 1;
        }
        else{
            return 0;
        }
    }
    return -2;
}

void Graph::edgeToNode(unsigned someEdge, unsigned &n1, unsigned &n2)
{
    vector <unsigned> V;

    V = edges[someEdge]; 
    
    if (V.empty() != 1)
    {
        if (V[0] < V[1])
        {
            n1 = V[0];
            n2 = V[1];
        } else {
            n2 = V[0];
            n1 = V[1];
        }
    }
}

void Graph::edgeToNodeEdgeIndex(unsigned someEdge, unsigned &n1, unsigned &n2, Matrix &edgeIndex)
{
    for (unsigned i=0;i< edgeIndex.rows;i++)
    {
        for (unsigned j=0;j< edgeIndex.cols;j++)
        {
            if (edgeIndex(i,j) == someEdge)
            {
                n1 = i;
                n2 = j;
                return;
            }
        }
    }
    cerr << "Should not get here in edgeToNodeEdgeIndex. " << endl;
    cerr << "someEdge = " << someEdge << endl;
    cerr << "edgeIndex" << endl;
    cerr << edgeIndex;
    exit (0);
}

void Graph::printVertexEdgeMatrix()
{
    vector <unsigned>   V;
    Matrix  VEM(numNodes,numEdges); // By default this is filled with zeros
    for (int i=0;i<numEdges;i++)
    {
        V = edges[i]; 

        VEM(V[0],i) = 1;
        VEM(V[1],i) = -1;

    }

    cout << "Vertex Edge matrix." << endl;
    cout << VEM;
}

int Graph::nodesConnected(int a, int b)
{
    // Compute the transitive closure if it hasn't been
    if (transClosureComputed == 0)
    {
        connComp();
    }
    //cout << &transClosure << endl;
    //cout << transClosure;
    if (transClosure(a,b) == 1)
    {
        return 1;
    }
    else 
    {
        return 0;
    }
}

list <vector <unsigned> > Graph::shortestPathList (unsigned n1, unsigned n2)
{
    //cout << "predMatrixComputed = " << predMatrixComputed << endl;
    if (predMatrixComputed == 0)
    {
        //cout << "Computing FloydWarshall" << endl;
        predMatrix = FloydWarshall();
        predMatrixComputed = 1;
    }
    //else 
    //{
    //    cout << "Not Computing FloydWarshall" << endl;
    //    cout << predMatrix;
    //    cout << "That was it." << endl;
    //}

    list <vector <unsigned> > newCycle;
    vector <unsigned> V;

    unsigned tn1 = n1;
    unsigned tn2 = n2;
    //cout << "predMatrix:" << endl;
    //cout << predMatrix;

    while ((int)predMatrix(tn1,tn2) != -1)
    {
        V.push_back((unsigned)predMatrix(tn1,tn2));
        V.push_back(tn2);
        newCycle.push_back(V);
        V.clear();
        tn2 = (unsigned)predMatrix(tn1,tn2);
    }

    return newCycle;
}

// Converts a list of edges given by end nodes to a set of their edge numbers
set <unsigned> Graph::listEdgesToSet (list <vector <unsigned> > &someEdges)
{
    set <unsigned> newEdges;

    list <vector <unsigned> >::iterator lvit = someEdges.begin();

    vector <unsigned> V;

    for (;lvit != someEdges.end();lvit++)
    {
        V = *lvit;
        newEdges.insert(edgeNumber(V[0],V[1]));
    }
    return newEdges;
}

set <unsigned> Graph::shortestPath (unsigned n1, unsigned n2)
{
    //cout << "predMatrixComputed = " << predMatrixComputed << endl;
    if (predMatrixComputed == 0)
    {
        //cout << "Computing FloydWarshall" << endl;
        predMatrix = FloydWarshall();
        predMatrixComputed = 1;
    }
    //else 
    //{
    //    cout << "Not Computing FloydWarshall" << endl;
    //    cout << predMatrix;
    //    cout << "That was it." << endl;
    //}

    set <unsigned> newCycle;

    unsigned tn1 = n1;
    unsigned tn2 = n2;
    //cout << "predMatrix:" << endl;
    //cout << predMatrix;

    while ((int)predMatrix(tn1,tn2) != -1)
    {
        newCycle.insert(edgeNumber((unsigned)predMatrix(tn1,tn2),tn2));
        tn2 = (unsigned)predMatrix(tn1,tn2);
    }

    return newCycle;
}

list <set <unsigned> > Graph::NagIbar ()
{
    
    //declare list of edge partitions to be returned
    set <unsigned> emptySet;
    list <set <unsigned> > edgePartitions(numEdges, emptySet);

    //matrix to know whether edges are scanned or not
    int **Unscanned = new int*[numNodes];
    for (int i=0;i<numNodes;i++){
        Unscanned[i] = new int[numNodes];
    }
    
    //initialize the matrix to say none of the edges are scanned
    int i, j;
    for(i = 0; i < numNodes; i++)
      for(j = 0; j < numNodes; j++)
      {	
	if(j == i)
	  Unscanned[i][j] = 0;
	else
	  Unscanned[i][j] = (int)adjMatrix(i,j);
      }

    //make buckets and put all vertex numbers in 0th bucket
    VertexBucket allVertices(0);
    for(i = 0; i < numNodes; i++)
      allVertices.vertices.push_front(i);
    list <VertexBucket> Buckets(1, allVertices);
   
    //array to hold r-values for each vertex
    int *rvalues = new int[numNodes];
    
    //initialize each r-value to 0
    for(i = 0; i < numNodes; i++)
      rvalues[i] = 0;
    
    int v, r, node;
    list <set <unsigned> >::iterator lsit;
    lsit = edgePartitions.begin();
    list <VertexBucket>::iterator lsit2;
    lsit2 = Buckets.begin();
    list <VertexBucket>::iterator lsit3;
    lsit3 = Buckets.begin();
    set <unsigned> itrSet;
    
    //as long as we still have a vertex in a bucket
    while(!Buckets.empty())
    {
      lsit2 = Buckets.begin();
      v = (*lsit2).vertices.front();
      (*lsit2).vertices.pop_front();
      if((*lsit2).vertices.empty())
	Buckets.erase(lsit2);
      for(node = 0; node < numNodes; node++)
      {
        if(Unscanned[v][node] == 1)
        {
          r = rvalues[node] + 1;
          lsit = edgePartitions.begin();
	  advance(lsit, r-1);
          (*lsit).insert(edgeNumber(v,node));
          if(rvalues[v] == rvalues[node])
            rvalues[v] += 1;
          rvalues[node] += 1;
          //move vertex node to new bucket
          lsit3 = Buckets.begin();
          for( ; (*lsit3).r != rvalues[node] && lsit3 != Buckets.end() ; lsit3++);
          if((*lsit3).r == rvalues[node] && lsit3 != Buckets.end())
          {
            (*lsit3).vertices.push_front(node); // This was causing seg faults on g1.ni
            ++lsit3;
            (*lsit3).vertices.remove(node);
            if((*lsit3).vertices.empty())
              Buckets.erase(lsit3);
          }
          else
          {
            lsit3 = Buckets.begin();
            for( ; (*lsit3).r > rvalues[node]; lsit3++);
            VertexBucket newVertexBucket(rvalues[node]); 
            newVertexBucket.vertices.push_front(node);
            Buckets.insert(lsit3, newVertexBucket);         
	    (*lsit3).vertices.remove(node);
	    if((*lsit3).vertices.empty())
	      Buckets.erase(lsit3); 
	  }
          Unscanned[v][node] = 0;
	  Unscanned[node][v] = 0;
        }
      }
    }
    for (int i=0;i<numNodes;i++){
        delete Unscanned[i];
    }
    delete Unscanned;
    delete rvalues;
    return edgePartitions;
}

// This takes the output of NagIbar and returns a matrix which indexes
// the edges given by the input. The index of each set is orbitally chosen
Matrix Graph::calcEdgeIndex(list <set <unsigned> > &edgePartition)
{
    Matrix edgeIndex(numNodes,numNodes); //Creates a numNodes^2 matrix filled with 0's
    list <set <unsigned> >::const_iterator lsit = edgePartition.begin();
    unsigned labelIndex = 0;
    unsigned n1,n2;
    for (unsigned i=0;i<(unsigned)numNodes;i++)
    {
        for (unsigned j=0;j<(unsigned)numNodes;j++)
        {
            edgeIndex(i,j) = -1;
        }
    }

    lsit = edgePartition.begin();
    set <unsigned> thisSet;
    for(;lsit != edgePartition.end();lsit++)
    {
        thisSet = *lsit;
        set <unsigned>::const_iterator  sit = thisSet.begin(); 
        for(;sit != thisSet.end();sit++)
        {
            // Find each edge and label appropriately
            edgeToNode(*sit,n1,n2);
            edgeIndex(n1,n2) = labelIndex;
            labelIndex++;
        }
    }
    return edgeIndex;
}

unsigned Graph::MatsuiBottom(Matrix &edgeIndex, set <unsigned> &someEdges)
{
    unsigned n1, n2;
    vector <unsigned> someEdge;
    unsigned currentMax;
    set <unsigned>::const_iterator sit = someEdges.begin();
    someEdge = edges[*sit];
    currentMax = (unsigned)edgeIndex(someEdge[0],someEdge[1]);
    n1 = someEdge[0];
    n2 = someEdge[1];

    while (sit != someEdges.end())
    {
        someEdge = edges[*sit];
        if (edgeIndex(someEdge[0],someEdge[1]) > currentMax)
        {
            currentMax = (unsigned)edgeIndex(someEdge[0],someEdge[1]);
            n1 = someEdge[0];
            n2 = someEdge[1];
        }

        sit++;
    }

    return edgeNumber(n1,n2);
}

unsigned Graph::MatsuiBottom(Matrix &edgeIndex)
{
    unsigned n1, n2;
    vector <unsigned> someEdge;
    unsigned currentMax;
    map <unsigned, vector <unsigned> >::const_iterator mvit = edges.begin();
    someEdge = (*mvit).second;
    currentMax = (unsigned)edgeIndex(someEdge[0],someEdge[1]);
    n1 = someEdge[0];
    n2 = someEdge[1];

    while (mvit != edges.end())
    {
        someEdge = (*mvit).second;
        if (edgeIndex(someEdge[0],someEdge[1]) > currentMax)
        {
            currentMax = (unsigned)edgeIndex(someEdge[0],someEdge[1]);
            n1 = someEdge[0];
            n2 = someEdge[1];
        }

        mvit++;
    }

    return edgeNumber(n1,n2);
}

unsigned Graph::MatsuiTop(Matrix &edgeIndex, set <unsigned> &someEdges)
{
    unsigned n1, n2;
    vector <unsigned> someEdge;
    unsigned currentMin;
    set <unsigned>::const_iterator sit = someEdges.begin();
    someEdge = edges[*sit];
    currentMin = (unsigned)edgeIndex(someEdge[0],someEdge[1]);
    n1 = someEdge[0];
    n2 = someEdge[1];

    while (sit != someEdges.end())
    {
        someEdge = edges[*sit];
        if (edgeIndex(someEdge[0],someEdge[1]) < currentMin)
        {
            currentMin = (unsigned)edgeIndex(someEdge[0],someEdge[1]);
            n1 = someEdge[0];
            n2 = someEdge[1];
        }

        sit++;
    }

    return edgeNumber(n1,n2);
}

unsigned Graph::MatsuiTop(Matrix &edgeIndex)
{
    unsigned n1, n2;
    vector <unsigned> someEdge;
    unsigned currentMin;
    map <unsigned, vector <unsigned> >::const_iterator mvit = edges.begin();
    someEdge = (*mvit).second;
    currentMin = (unsigned)edgeIndex(someEdge[0],someEdge[1]);
    n1 = someEdge[0];
    n2 = someEdge[1];

    while (mvit != edges.end())
    {
        someEdge = (*mvit).second;
        if (edgeIndex(someEdge[0],someEdge[1]) < currentMin)
        {
            currentMin = (unsigned)edgeIndex(someEdge[0],someEdge[1]);
            n1 = someEdge[0];
            n2 = someEdge[1];
        }

        mvit++;
    }

    return edgeNumber(n1,n2);
}

unsigned Graph::getEdgeIndex(unsigned someEdge, Matrix &edgeIndex)
{
    unsigned n1,n2;
   
    edgeToNode(someEdge,n1,n2);
    return (unsigned)edgeIndex(n1,n2);
}

void findChildren(Graph &G, set <unsigned> &initTree, set <unsigned> &deltaf, set <unsigned> &deltag, set <unsigned> deltaH, Matrix &edgeIndex, Matrix &Weight, set <Matrix, ltcolvec> &projTrees, unsigned printMod, unsigned printTrees)
{
    G.findChildrenSpanningTreeCount++;
    G.findChildrenBFSLevel++;
    set <unsigned> thisT;
    set <unsigned> tempUnsignedSet;
    set <unsigned> thisTremH_TSet;
    set <unsigned> H;

    if (printMod != 0)
    {
        if ((G.findChildrenSpanningTreeCount % printMod) == 0)
        {
            cout << G.findChildrenSpanningTreeCount << endl;
        }
    }

    //cout << "*************** findChildren ***************" << endl;
    //cout << "findChildrenBFSLevel = " << G.findChildrenBFSLevel << endl;
    //cout << "deltaf: ";
    //copy(deltaf.begin(), deltaf.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    //cout << "deltag: ";
    //copy(deltag.begin(), deltag.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    //cout << "deltaH: ";
    //copy(deltaH.begin(), deltaH.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    set_difference(initTree.begin(), initTree.end(), deltag.begin(), deltag.end(), insert_iterator< set <unsigned> >(tempUnsignedSet, tempUnsignedSet.begin()));
    
    set_union(tempUnsignedSet.begin(), tempUnsignedSet.end(), deltaf.begin(), deltaf.end(), insert_iterator< set <unsigned> >(thisT, thisT.begin()));

    tempUnsignedSet.clear(); // Clear for better memory
    if (printTrees == 1)
    {
        cout << "Tree: ";
        copy(thisT.begin(), thisT.end(), ostream_iterator<const int>(cout, " "));
        cout << endl;
    }

    projTrees.insert((Weight.subColumns(thisT)).rowSum());
    //cout << "thisT: ";
    //copy(thisT.begin(), thisT.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    //cout << "G.isSpanningForest(thisT) = " << G.isSpanningForest(thisT) << endl;

    // Output this spanning tree
    //cout << G.subGraph(thisT);

    // j, k and jprime are in terms of edgeIndex. Thus we need to be
    // carefull when finding e_j
    unsigned k = G.getEdgeIndex(G.MatsuiBottom(edgeIndex, thisT),edgeIndex); 
    //unsigned k = G.MatsuiBottom(edgeIndex, thisT); 
    unsigned jprime = k;
    // Add function to properly ''increment'' according to edgeIndex
    unsigned j = k + 1;
    unsigned n1,n2;

    //cout << "k = " << k << endl;
    //cout << "j = " << j << endl;
    //cout << "jprime = " << jprime << endl;
    //cout << "jprime + 2*G.numNodes - 3 = " << jprime + 2*G.numNodes - 3 << endl;
    //cout << "G.numEdges = " << G.numEdges << endl << endl;

    set_difference(initTree.begin(), initTree.end(), deltaH.begin(), deltaH.end(), insert_iterator< set <unsigned> >(H, H.begin()));
    //cout << "H: ";
    //copy(H.begin(), H.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    set_difference(thisT.begin(), thisT.end(), H.begin(), H.end(), insert_iterator< set <unsigned> >(thisTremH_TSet, thisTremH_TSet.begin()));

    // while (j <= min(jprime + 2n - 3,m) )
    while (j <= jprime + 2*G.numNodes - 3 && j < (unsigned)G.numEdges)
    {
        //cout << "   j = " << j << endl;

        Graph *thisTremH_T = new Graph;
        *thisTremH_T = G.subGraph(thisTremH_TSet);
        //cout << "T - H" << endl;
        //cout << *thisTremH_T << endl;
        
        // Find first edge j such that edge j is not connected on subgraph TremH_T
        //cout << "Calling edgeToNodeEdgeIndex 1" << endl;
        Graph::edgeToNodeEdgeIndex(j,n1,n2,edgeIndex); 
        //cout << "T-H: n1 = " << n1 << endl;
        //cout << "T-H: n2 = " << n2 << endl;
        while ((j <= jprime + 2*G.numNodes - 3 && j < (unsigned)G.numEdges) && thisTremH_T->nodesConnected(n1,n2) == 1)
        {
            // Add function to properly ''increment'' according to edgeIndex
            //cout << "j++" << endl;
            j++;
            if (j <= jprime + 2*G.numNodes - 3 && j < (unsigned)G.numEdges)
            {
                //cout << "Calling edgeToNodeEdgeIndex 2" << endl;
                Graph::edgeToNodeEdgeIndex(j,n1,n2,edgeIndex); 
            }
        }
        //cout << "Done with j++" << endl;
        delete thisTremH_T;

        if (j <= jprime + 2*G.numNodes - 3 && j < (unsigned)G.numEdges)
        {
            //cout << "   j <= jprime + 2*G.numNodes - 3 || j <= G.numEdges" << endl;
            //cout << "   j = " << j << endl;
            // f should be indexed in terms of G
            //unsigned f = j;
            //cout << "Calling edgeToNodeEdgeIndex 1" << endl;
            Graph::edgeToNodeEdgeIndex(j,n1,n2,edgeIndex);
            unsigned f = G.edgeNumber(n1,n2); 
            jprime = j;
            // Add function to properly ''increment'' according to edgeIndex
            j++;
            unsigned g = G.numNodes;

            //Find Cycle(thisT,f) \cap H 
            Graph *thisTGraph = new Graph;
            *thisTGraph = G.subGraph(thisT);
            //cout << "*thisTGraph" << endl;
            //cout << *thisTGraph;
            G.edgeToNode(f,n1,n2);
            //cout << "n1 = " << n1 << endl;
            //cout << "n2 = " << n2 << endl;
            // This cycle is not correct. This returns edges with 
            // index relative to the local graph, not the edge index
            // relative to the original graph G.
            //set <unsigned> cycle = thisTGraph->shortestPath(n1,n2);
            list <vector <unsigned> > cycleList = thisTGraph->shortestPathList(n1,n2);
            set <unsigned > cycle = G.listEdgesToSet(cycleList);
            //cout << "   cycle: ";
            //copy(cycle.begin(), cycle.end(), ostream_iterator<const int>(cout, " "));
            //cout << endl;

            delete thisTGraph; // Save memory

            // Now need cycle \cap H
            // H is simply initT - deltag
            // Already computed above. Called H
            //set_intersection(initTree.begin(), initTree.end(), cycle.begin(), cycle.end(), insert_iterator< set <unsigned> >(tempUnsignedSet, tempUnsignedSet.begin()));
            
            set <unsigned> D;
            set_intersection(H.begin(), H.end(), cycle.begin(), cycle.end(), insert_iterator< set <unsigned> >(D, D.begin()));
            tempUnsignedSet.clear();

            set <unsigned>::iterator sit;
            //cout << "D: ";
            //copy(D.begin(), D.end(), ostream_iterator<const int>(cout, " "));
            //cout << endl;

            // Now go through D by decreasing order given by the edgeIndex structure
            // These are the elements we add to deltaH, but need to remove when this
            // loop is done.
            set <unsigned> remH; 
            remH.clear ();
            while (D.empty() != 1)
            {
                //Find largest D according to edgeIndex
                sit = D.begin();
                unsigned currentMax,currentEdge;
                G.edgeToNode(*sit,n1,n2);
                currentMax = (unsigned)edgeIndex(n1,n2);
                currentEdge = *sit;
                sit++;
                for (;sit != D.end();sit++)
                {
                    G.edgeToNode(*sit,n1,n2);
                    if (edgeIndex(n1,n2) > currentMax)
                    {
                        currentMax = (unsigned)edgeIndex(n1,n2);
                        currentEdge = *sit;
                    }
                }
                g = currentEdge;        
                D.erase(g);
            
                deltag.insert(g);
                deltaf.insert(f);
                deltaH.insert(g);
                remH.insert(g);
                //int blahblah;
                //cin >> blahblah;
                //cout << "Calling findChildren with g = " << g << "  f = " << f << endl;
                findChildren(G,initTree,deltaf,deltag,deltaH,edgeIndex,Weight, projTrees, printMod,printTrees);
                deltag.erase(g);
                deltaf.erase(f);
            }
            set_difference(deltaH.begin(), deltaH.end(), remH.begin(), remH.end(), insert_iterator< set <unsigned> >(tempUnsignedSet, tempUnsignedSet.begin()));
            deltaH = tempUnsignedSet;

        }
        //else 
        //{
        //    cout << "   j > jprime + 2*G.numNodes - 3 || j <= G.numEdges" << endl;
        //    cout << "   j = " << j << endl;
        //}
    }
    G.findChildrenBFSLevel--;
}

void findChildren(Graph &G, set <unsigned> &initTree, set <unsigned> &deltaf, set <unsigned> &deltag, set <unsigned> deltaH, Matrix &edgeIndex, unsigned printMod, unsigned printTrees)
{
    G.findChildrenSpanningTreeCount++;
    G.findChildrenBFSLevel++;
    set <unsigned> thisT;
    set <unsigned> tempUnsignedSet;
    set <unsigned> thisTremH_TSet;
    set <unsigned> H;

    if (printMod != 0)
    {
        if ((G.findChildrenSpanningTreeCount % printMod) == 0)
        {
            cout << G.findChildrenSpanningTreeCount << endl;
        }
    }

    //cout << "*************** findChildren ***************" << endl;
    //cout << "findChildrenBFSLevel = " << G.findChildrenBFSLevel << endl;
    //cout << "deltaf: ";
    //copy(deltaf.begin(), deltaf.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    //cout << "deltag: ";
    //copy(deltag.begin(), deltag.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    //cout << "deltaH: ";
    //copy(deltaH.begin(), deltaH.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    set_difference(initTree.begin(), initTree.end(), deltag.begin(), deltag.end(), insert_iterator< set <unsigned> >(tempUnsignedSet, tempUnsignedSet.begin()));
    
    set_union(tempUnsignedSet.begin(), tempUnsignedSet.end(), deltaf.begin(), deltaf.end(), insert_iterator< set <unsigned> >(thisT, thisT.begin()));

    tempUnsignedSet.clear(); // Clear for better memory
    if (printTrees == 1)
    {
        cout << "Tree: ";
        copy(thisT.begin(), thisT.end(), ostream_iterator<const int>(cout, " "));
        cout << endl;
    }

    //cout << "thisT: ";
    //copy(thisT.begin(), thisT.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    //cout << "G.isSpanningForest(thisT) = " << G.isSpanningForest(thisT) << endl;

    // Output this spanning tree
    //cout << G.subGraph(thisT);

    // j, k and jprime are in terms of edgeIndex. Thus we need to be
    // carefull when finding e_j
    unsigned k = G.getEdgeIndex(G.MatsuiBottom(edgeIndex, thisT),edgeIndex); 
    //unsigned k = G.MatsuiBottom(edgeIndex, thisT); 
    unsigned jprime = k;
    // Add function to properly ''increment'' according to edgeIndex
    unsigned j = k + 1;
    unsigned n1,n2;

    //cout << "k = " << k << endl;
    //cout << "j = " << j << endl;
    //cout << "jprime = " << jprime << endl;
    //cout << "jprime + 2*G.numNodes - 3 = " << jprime + 2*G.numNodes - 3 << endl;
    //cout << "G.numEdges = " << G.numEdges << endl << endl;

    set_difference(initTree.begin(), initTree.end(), deltaH.begin(), deltaH.end(), insert_iterator< set <unsigned> >(H, H.begin()));
    //cout << "H: ";
    //copy(H.begin(), H.end(), ostream_iterator<const int>(cout, " "));
    //cout << endl;
    set_difference(thisT.begin(), thisT.end(), H.begin(), H.end(), insert_iterator< set <unsigned> >(thisTremH_TSet, thisTremH_TSet.begin()));

    // while (j <= min(jprime + 2n - 3,m) )
    while (j <= jprime + 2*G.numNodes - 3 && j < (unsigned)G.numEdges)
    {
        //cout << "   j = " << j << endl;

        Graph *thisTremH_T = new Graph;
        *thisTremH_T = G.subGraph(thisTremH_TSet);
        //cout << "T - H" << endl;
        //cout << *thisTremH_T << endl;
        
        // Find first edge j such that edge j is not connected on subgraph TremH_T
        //cout << "Calling edgeToNodeEdgeIndex 1" << endl;
        Graph::edgeToNodeEdgeIndex(j,n1,n2,edgeIndex); 
        //cout << "T-H: n1 = " << n1 << endl;
        //cout << "T-H: n2 = " << n2 << endl;
        while ((j <= jprime + 2*G.numNodes - 3 && j < (unsigned)G.numEdges) && thisTremH_T->nodesConnected(n1,n2) == 1)
        {
            // Add function to properly ''increment'' according to edgeIndex
            //cout << "j++" << endl;
            j++;
            if (j <= jprime + 2*G.numNodes - 3 && j < (unsigned)G.numEdges)
            {
                //cout << "Calling edgeToNodeEdgeIndex 2" << endl;
                Graph::edgeToNodeEdgeIndex(j,n1,n2,edgeIndex); 
            }
        }
        //cout << "Done with j++" << endl;
        delete thisTremH_T;

        if (j <= jprime + 2*G.numNodes - 3 && j < (unsigned)G.numEdges)
        {
            //cout << "   j <= jprime + 2*G.numNodes - 3 || j <= G.numEdges" << endl;
            //cout << "   j = " << j << endl;
            // f should be indexed in terms of G
            //unsigned f = j;
            //cout << "Calling edgeToNodeEdgeIndex 1" << endl;
            Graph::edgeToNodeEdgeIndex(j,n1,n2,edgeIndex);
            unsigned f = G.edgeNumber(n1,n2); 
            jprime = j;
            // Add function to properly ''increment'' according to edgeIndex
            j++;
            unsigned g = G.numNodes;

            //Find Cycle(thisT,f) \cap H 
            Graph *thisTGraph = new Graph;
            *thisTGraph = G.subGraph(thisT);
            //cout << "*thisTGraph" << endl;
            //cout << *thisTGraph;
            G.edgeToNode(f,n1,n2);
            //cout << "n1 = " << n1 << endl;
            //cout << "n2 = " << n2 << endl;
            // This cycle is not correct. This returns edges with 
            // index relative to the local graph, not the edge index
            // relative to the original graph G.
            //set <unsigned> cycle = thisTGraph->shortestPath(n1,n2);
            list <vector <unsigned> > cycleList = thisTGraph->shortestPathList(n1,n2);
            set <unsigned > cycle = G.listEdgesToSet(cycleList);
            //cout << "   cycle: ";
            //copy(cycle.begin(), cycle.end(), ostream_iterator<const int>(cout, " "));
            //cout << endl;

            delete thisTGraph; // Save memory

            // Now need cycle \cap H
            // H is simply initT - deltag
            // Already computed above. Called H
            //set_intersection(initTree.begin(), initTree.end(), cycle.begin(), cycle.end(), insert_iterator< set <unsigned> >(tempUnsignedSet, tempUnsignedSet.begin()));
            
            set <unsigned> D;
            set_intersection(H.begin(), H.end(), cycle.begin(), cycle.end(), insert_iterator< set <unsigned> >(D, D.begin()));
            tempUnsignedSet.clear();

            set <unsigned>::iterator sit;
            //cout << "D: ";
            //copy(D.begin(), D.end(), ostream_iterator<const int>(cout, " "));
            //cout << endl;

            // Now go through D by decreasing order given by the edgeIndex structure
            // These are the elements we add to deltaH, but need to remove when this
            // loop is done.
            set <unsigned> remH; 
            remH.clear ();
            while (D.empty() != 1)
            {
                //Find largest D according to edgeIndex
                sit = D.begin();
                unsigned currentMax,currentEdge;
                G.edgeToNode(*sit,n1,n2);
                currentMax = (unsigned)edgeIndex(n1,n2);
                currentEdge = *sit;
                sit++;
                for (;sit != D.end();sit++)
                {
                    G.edgeToNode(*sit,n1,n2);
                    if (edgeIndex(n1,n2) > currentMax)
                    {
                        currentMax = (unsigned)edgeIndex(n1,n2);
                        currentEdge = *sit;
                    }
                }
                g = currentEdge;        
                D.erase(g);
            
                deltag.insert(g);
                deltaf.insert(f);
                deltaH.insert(g);
                remH.insert(g);
                //int blahblah;
                //cin >> blahblah;
                //cout << "Calling findChildren with g = " << g << "  f = " << f << endl;
                findChildren(G,initTree,deltaf,deltag,deltaH,edgeIndex,printMod,printTrees);
                deltag.erase(g);
                deltaf.erase(f);
            }
            set_difference(deltaH.begin(), deltaH.end(), remH.begin(), remH.end(), insert_iterator< set <unsigned> >(tempUnsignedSet, tempUnsignedSet.begin()));
            deltaH = tempUnsignedSet;

        }
        //else 
        //{
        //    cout << "   j > jprime + 2*G.numNodes - 3 || j <= G.numEdges" << endl;
        //    cout << "   j = " << j << endl;
        //}
    }
    G.findChildrenBFSLevel--;
}

//set <unsigned> setUnionUnsigned (set <unsigned> &, set <unsigned> &)
//{
//
//}
//
//set <unsigned> setDifferenceUnsigned (set <unsigned> &, set <unsigned> &)
//{
//
//}
