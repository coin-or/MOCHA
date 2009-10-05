// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

#include "matroid.h"
#include "mathprog.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <set>
#include <list>
#include <algorithm>
#include <iterator>
#include <string>
#include <ctime>
#include <sstream>



using namespace::std;

Matrix MyFctMatrix;

double MyFct (Matrix M)
{
    double returnValue=0;
    
    for (int i=0;i<MyFctMatrix.rows && i < M.rows;i++)
    {
        returnValue += (MyFctMatrix(i,0)-M(i,0))*(MyFctMatrix(i,0)-M(i,0));
    }
    return returnValue; 
}   

int main (int argc, char *argv[])
{

   ifstream inputfile;
    char fileName[80];

    srand(time(0));
    if (argc >= 2)
    {
        inputfile.open(argv[1]);
    }
    else {
        cout << "Input file: ";
        cin >> fileName;
        inputfile.open(fileName);
    }

    ofstream outputfile;
    if (argc <= 2)
    {
        cout << "Enter output filename: ";
        cin >> fileName;
        cout << "Writing to file: " << fileName << endl;
        outputfile.open(fileName);
    }
    else if (argc >= 3)
    {
        cout << "Writing to file: " << argv[2] << endl;
        outputfile.open(argv[2]);
    }

    cout << "Selecting balance function MyFct" << endl;
    ProjBalMatroidOpt MyOpt(inputfile,&MyFct);
    cout << MyOpt;

    int LocalSearches;
    if (argc >= 4)
    {
        LocalSearches = atoi(argv[3]);
        cout << "Number of Local searches: ";
        cout << LocalSearches << endl;
        
    }
    else 
    {
        cout << "Number of Local searches: ";
        cin >> LocalSearches;
    }
    cout << LocalSearches << " local searches." << endl;
    if (argc >= 5)
    {
        ifstream MyFctMatrixFile;
        cout << "Reading file " << argv[4] << " for MyFctMatrix." << endl;
        MyFctMatrixFile.open(argv[4]);
        MyFctMatrixFile >> MyFctMatrix;
        cout << "MyFctMatrix" << endl;
        cout << MyFctMatrix << endl;
    }
    else
    {
        MyFctMatrix = Matrix(2,1);
        MyFctMatrix(0,0) = 0;
        MyFctMatrix(1,0) = 0;
    }


    list <set <unsigned> > LocalSearchPivots;


    int tenPercent = (int) (LocalSearches / 10);
    if (tenPercent == 0)
    {
        tenPercent = 1;
    }

    time_t startTime = time(0);
    time_t innerLoopStartTime = time(0);
    for (int i=0;i<LocalSearches;i++)
    {
        LocalSearchPivots.clear();


        LocalSearchPivots = MyOpt.LocalSearchRandomStart();

        std::stringstream out;
        out << i;

        string LSPLabel = "LSPIVOTS";
        LSPLabel += out.str();
        MyOpt.ProjBalMatroidOpt::writePivotsMatlab(LocalSearchPivots,outputfile,LSPLabel);
        outputfile << "plot(" << LSPLabel <<  "(:,1)," << LSPLabel <<  "(:,2),'kx-');" << endl;
        outputfile << "plot(" << LSPLabel <<  "(1,1)," << LSPLabel <<  "(1,2),'g^');" << endl;
        outputfile << "plot(" << LSPLabel <<  "(end,1)," << LSPLabel <<  "(end,2),'r*');" << endl;

        outputfile << "hold on;" << endl;
    
        if ((i % tenPercent) == 0 && i != 0)
        {
            time_t innerLoopEndTime = time(0);
            cout << i/tenPercent << "0% ";
            cout << innerLoopEndTime - innerLoopStartTime << " seconds." << endl;
            innerLoopStartTime = time(0);
        }
    }
    time_t endTime = time(0);
    outputfile << "% TOTAL SECONDS " << endTime - startTime << endl;
    outputfile << LocalSearches << " searches. " << endl;
    cout << "100%" << endl;

}

