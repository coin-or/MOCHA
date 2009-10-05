// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE

// $Rev: 315 $ $Date: 2009-09-22 21:12:42 -0400 (Tue, 22 Sep 2009) $
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

using namespace::std;


double MyFct (Matrix M)
{
    double x,y;

    x = M(0,0); 
    y = M(1,0);

    return x*x + y*y; 
}   

double MyFct2 (Matrix M)
{
    double x,y;

    x = M(0,0); 
    y = M(1,0);

    return (x-35)*(x-35) + (y-85)*(y-85); 
}   

double MyFct3 (Matrix M)
{
    double x,y;

    x = M(0,0); 
    y = M(1,0);

    return (x-50)*(x-50) + (y-60)*(y-60); 
}   

int main (int argc, char *argv[])
{
    ifstream inputfile;
    list <set <unsigned> > MyOptPivots;
    list <set <unsigned> > TabuSearchPivots;
    list <set <unsigned> > TabuSearchBases;
    list <set <unsigned> > ConvexComparePivots;
    char fileName[80];
    set <unsigned> someRandomBasis;
    set < Matrix, ltcolvec> someProjPoints;
    set < Matrix, ltcolvec> someProjPointsTS;
    set < Matrix, ltcolvec> AutoBOXsomeProjPoints;
    set < Matrix, ltcolvec> AutoBOXsomeProjPointsTS;
    set < Matrix, ltcolvec> BFSProjPoints;
    set < Matrix, ltcolvec> POBFSProjPoints;
    set < Matrix, ltcolvec> ParetoOpt;
    set < Matrix, ltcolvec> CHParetoOpt;
    set < Matrix, ltcolvec> TriangularRegions;
    set < Matrix, ltcolvec> TriangularRegionsPrint;
    set < Matrix, ltcolvec> newPointsAutoBoxLS;
    set < Matrix, ltcolvec> newPointsAutoBoxTS;
    set < Matrix, ltcolvec> newPointsBoxLS;
    set < Matrix, ltcolvec> newPointsBoxTS;
    set < Matrix, ltcolvec> allBasesProj;
    set <Matrix, ltcolvec> CH;
    string openedFileName, S;
    int TabuSearch = 0;
    unsigned TabuSearchIterationLimit = 0;
    int BFSEnumerate = 0;
    int SimAnn = 0;
    int numAnnealings = 0;
    int paretoBoundary = 0;
    int paretoBoundaryTriangle = 0;
    int pivotTestType; // 1 is Local Search, 2 is Tabu Search 
    int numPivotTests = 0;
    int paretoBFS = 0;
    int numPivots;
    int pivotLimit;
    int paretoBoundaryTriangleBFS = 0;
    int calcAll = 0;
    list <double> temperatures;
    list <unsigned> times;
    list <set <unsigned> > SimAnnPivots;
    list <set <unsigned> > SimAnnMinBases;
    int AutoBOXTest = 0;
    int AutoBOXTestSearches = 0;
    int AutoBOXTestTS = 0;
    int AutoBOXTestSearchesTS = 0;
    int AutoBOXTestTSLimit = 0;

    int BOXTest = 0;
    int BOXTestSearches = 0;
    int BOXTestTS = 0;
    int BOXTestSearchesTS = 0;
    int BOXTestTSLimit = 0;
    int ConvexHull = 0;
    int MATLABOutput = 0;
    int numSearches;    
    int BFSDepth;
    int retryLimit;
    int BoundaryRetryLimit;
    int POnumSearches;    
    int POBFSDepth;
    int POretryLimit;
    int POBoundaryRetryLimit;


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


    if (inputfile.fail() == 1)
    {
        cout << "Can not open file " << fileName << endl;
        exit (0);
    }
    openedFileName = fileName;

    cout << "Selecting balance function MyFct3" << endl;
    ProjBalMatroidOpt MyOpt(inputfile,&MyFct3);

    //MyOpt.setBalFct(&MyFct);
    
    cout << "*************** Matroid Optimization Program ***************" << endl;
    cout << MyOpt;
    
    cout << "Different Fiber BFS Enumerations? (y/n) ";
    cin >> S;
    if (S == "Y" || S == "y" || S == "yes" || S == "YES")
    {
        cout << "Number of searches? ";
        cin >> numSearches;
        cout << "BFS depth? ";
        cin >> BFSDepth;
        cout << "Interior Random retry limit? ";
        cin >> retryLimit;
        cout << "Boundary random retry limit? ";
        cin >> BoundaryRetryLimit;
        BFSEnumerate = 1;
    }
    cout << "Run Auto Box Pivot Heuristic (Local Search)? (y/n) ";
    cin >> S;
    if (S == "Y" || S == "y" || S == "yes" || S == "YES")
    { 
        cout << "Number of attempts per point: ";
        cin >> BOXTestSearches;
        
        AutoBOXTest = 1;
    }
    cout << "Run Auto Box Pivot Heuristic (Tabu Search)? (y/n) ";
    cin >> S;
    if (S == "Y" || S == "y" || S == "yes" || S == "YES")
    { 
        cout << "Number of attempts per point: ";
        cin >> BOXTestSearchesTS;
        cout << "Tabu Search Limit: ";
        cin >> BOXTestTSLimit;
        
        AutoBOXTestTS = 1;
    }


    cout << "Run Box Pivot Heuristic (Local Search)? (y/n) ";
    cin >> S;
    Matrix lowerCornerLS(MyOpt.projDim(),1);
    Matrix upperCornerLS(MyOpt.projDim(),1);
    if (S == "Y" || S == "y" || S == "yes" || S == "YES")
    { 
        cout << "Number of attempts per point: ";
        cin >> BOXTestSearches;
        
        cout << "Enter lowerCorner Coordinates: " << endl;
        int tempInt;
        for (int i=0;i<(int)MyOpt.projDim();i++)
        {
            cout << " Row[" << i << "] = ";
            cin >> tempInt; 
            lowerCornerLS(i,0) = tempInt;
        }
        cout << "Enter upperCorner Coordinates: " << endl;
        for (int i=0;i<(int)MyOpt.projDim();i++)
        {
            cout << " Row[" << i << "] = ";
            cin >> tempInt; 
            upperCornerLS(i,0) = tempInt;
        }
        BOXTest = 1;
    }
    cout << "Run Box Pivot Heuristic (Tabu Search)? (y/n) ";
    cin >> S;
    Matrix lowerCornerTS(MyOpt.projDim(),1);
    Matrix upperCornerTS(MyOpt.projDim(),1);
    if (S == "Y" || S == "y" || S == "yes" || S == "YES")
    { 
        cout << "Number of attempts per point: ";
        cin >> BOXTestSearchesTS;
        cout << "Tabu Search Limit: ";
        cin >> BOXTestTSLimit;
        
        cout << "Enter lowerCorner Coordinates: " << endl;
        int tempInt;
        for (int i=0;i<(int)MyOpt.projDim();i++)
        {
            cout << " Row[" << i << "] = ";
            cin >> tempInt; 
            lowerCornerTS(i,0) = tempInt;
        }
        cout << "Enter upperCorner Coordinates: " << endl;
        for (int i=0;i<(int)MyOpt.projDim();i++)
        {
            cout << " Row[" << i << "] = ";
            cin >> tempInt; 
            upperCornerTS(i,0) = tempInt;
        }
        BOXTestTS = 1;
    }
    cout << "Run Boundary Calculation (dim=2 only)? (y/n) ";
    cin >> S;
    if (S == "Y" || S == "y" || S == "yes" || S == "YES")
    { 
        ConvexHull = 1;
    }

    cout << "Run Pareto Optimum test? " ;
    cin >> S;
    if (S == "Y" || S == "y" || S == "yes" || S == "YES")
    {
        int choiceCount = 0;
        cout << "Run Boundary only Pareto search? ";
        cin >> S;
        if (S == "Y" || S == "y" || S == "yes" || S == "YES")
        {
            paretoBoundary = 1;
            choiceCount++;
        }
        if (choiceCount == 0)
        {
            cout << "Run Boundary Triangular Region Pareto search? ";
            cin >> S;
            if (S == "Y" || S == "y" || S == "yes" || S == "YES" && choiceCount == 0)
            {
                cout << "Number of tests for each point? ";
                cin >> numPivots;
                cout << "Local Search(1) or Tabu Search(2). ";
                cin >> S;
                if (S == "1")  
                {
                    pivotTestType = 1;
                }
                else if (S == "2")
                {
                    pivotTestType = 2;
                    cout << "Tabu Search limit? ";
                    cin >> pivotLimit;
                }
                else 
                {
                    cerr << "Must enter 1 or 2." << endl;
                    exit(0);
                }
                paretoBoundaryTriangle = 1;
                choiceCount++;
            }
        }
        if (choiceCount == 0)
        {
            cout << "Run BFS Pareto search? ";
            cin >> S;
            if (S == "Y" || S == "y" || S == "yes" || S == "YES" && choiceCount == 0)
            {
                cout << "Number of searches? ";
                cin >> POnumSearches;
                cout << "BFS depth? ";
                cin >> POBFSDepth;
                cout << "Boundary random retry limit? ";
                cin >> POBoundaryRetryLimit;
                cout << "Interior random retry limit? ";
                cin >> POretryLimit;
                paretoBFS = 1;
                choiceCount++;
            }
        }
        if (choiceCount == 0)
        {
            cout << "Run BFS and Boundary Triangular Region Pareto search? ";
            cin >> S;
            if (S == "Y" || S == "y" || S == "yes" || S == "YES" && choiceCount == 0)
            {
                cout << "Number of searches? ";
                cin >> POnumSearches;
                cout << "BFS depth? ";
                cin >> POBFSDepth;
                cout << "Boundary random retry limit? ";
                cin >> POBoundaryRetryLimit;
                cout << "Random retry limit? ";
                cin >> POretryLimit;

                cout << "Number of tests for each point? ";
                cin >> numPivots;
                cout << "Local Search(1) or Tabu Search(2). ";
                cin >> S;
                if (S == "1")  
                {
                    pivotTestType = 1;
                }
                else if (S == "2")
                {
                    pivotTestType = 2;
                    cout << "Tabu Search limit? ";
                    cin >> pivotLimit;
                }
                else 
                {
                    cerr << "Must enter 1 or 2." << endl;
                    exit(0);
                }

                paretoBoundaryTriangleBFS = 1;
                choiceCount++;
            }
        }
        if (choiceCount > 1)
        {
            cerr << "Can not choose more than one pareto option." << endl;
            exit(0);
        }
    }

    //cout << "Run Tabu Search? (y/n) ";
    //cin >> S;
    //if (S == "Y" || S == "y" || S == "yes" || S == "YES")
    //{
    //    cout << "TabuSearchIterationLimit? ";
    //    cin >> TabuSearchIterationLimit;
    //    
    //    TabuSearch = 1;
    //}

    //cout << "Run Simulated Annealing? (y/n) ";
    //cin >> S;
    //if (S == "Y" || S == "y" || S == "yes" || S == "YES")
    //{
    //    double tempTemperature;
    //    unsigned tempTime;
    //    cout << "Number of annealings?" << endl;
    //    cin >> numAnnealings;
    //    for (int i=0;i< numAnnealings;i++)
    //    {
    //        cout << "Iteration " << i + 1 << endl;
    //        cout << "Temperature? ";
    //        cin >> tempTemperature;
    //        cout << "Time? ";
    //        cin >> tempTime;
    //        temperatures.push_back(tempTemperature);
    //        times.push_back(tempTime);
    //    }

    //    SimAnn = 1;
    //}
    
    cout << "Calculate all projected bases using brute-force enumeration. (Graphical only)? (y/n) ";
    cin >> S;
    if (S == "Y" || S == "y" || S == "yes" || S == "YES")
    {
        calcAll = 1;
    }

    if (argc <= 2)
    {
        cout << "Output matlab file? (y/n) ";
        cin >> S;
        if (S == "Y" || S == "y" || S == "yes" || S == "YES")
        { 
            MATLABOutput = 1;
            cout << "Output file: ";
            cin >> fileName;
            S = fileName;
            if (S == openedFileName) {
                cout << "Can not save over " << S << endl;
                exit (0);
            }
        }
    }
    else if (argc >= 3)
    {
        MATLABOutput = 1;
    }

    cout << "*#*#*#*#*#*#*#*#*#*#*#*# PERFORMING CALCULATIONS *#*#*#*#*#*#*#*#*#*#*#*#" << endl;
    
    time_t startTime = time(0);
    if (BFSEnumerate == 1)
    { 
        cout << "Camputing BFSRandomStarts" << endl;
        MyOpt.MultiBFSRandomStarts(numSearches,BFSDepth,BoundaryRetryLimit,1,retryLimit,BFSProjPoints); 
        //cout << "Computing Pareto Optimum" << endl;
        //ParetoOpt = MyOpt.ParetoOptimum(BFSProjPoints); 
    }
    if (AutoBOXTest == 1)
    {
        cout << "Running AutoBoundsPivotTest (Local Search)." << endl;
        AutoBOXsomeProjPoints = BFSProjPoints;    
        MyOpt.AutoBoundsPivotTestLocalSearch(AutoBOXsomeProjPoints,BOXTestSearches,newPointsAutoBoxLS);
    }
    if (AutoBOXTestTS == 1)
    {
        cout << "Running AutoBoundsPivotTest (Tabu Search)." << endl;
        AutoBOXsomeProjPointsTS = BFSProjPoints;    
        MyOpt.AutoBoundsPivotTestTabuSearch(AutoBOXsomeProjPoints,BOXTestSearchesTS,newPointsAutoBoxTS,BOXTestTSLimit);
    }


    if (BOXTest == 1)
    {
        cout << "Running BoxPivotTest (Local Search)." << endl;
        someProjPoints = BFSProjPoints;    
        MyOpt.BoxPivotTestLocalSearch(lowerCornerLS,upperCornerLS,someProjPoints,BOXTestSearches,newPointsBoxLS);
    }
    if (BOXTestTS == 1)
    {
        cout << "Running BoxPivotTest (Tabu Search)." << endl;
        someProjPointsTS = BFSProjPoints;    
        MyOpt.BoxPivotTestTabuSearch(lowerCornerTS,upperCornerTS,someProjPoints,BOXTestSearchesTS,newPointsBoxTS,BOXTestTSLimit);
    }
    if (ConvexHull == 1)
    {
        cout << "Calculating Boundary." << endl;
        MyOpt.Boundary(CH);
        //CH = MyOpt.ParetoOptimum(CH);
        //CH = MyOpt.BoundaryTrianglesTwoDim (CH);

        //cout << "Convex Hull" << endl;
        //set <Matrix, ltcolvec>::const_iterator tmit = CH.begin();
        //while (tmit != CH.end())
        //{
        //    cout << (*tmit);
        //    cout << "----" << endl;
        //    tmit++;
        //}
    }
    //if (TabuSearch == 1)
    //{
    //    cout << "Running Tabu Search." << endl;
    //    someRandomBasis = MyOpt.randomBasis ();
    //    TabuSearchPivots = MyOpt.TabuSearchHeuristic(someRandomBasis,TabuSearchIterationLimit,TabuSearchBases);
    //    ConvexComparePivots = MyOpt.LocalSearch(someRandomBasis);
    //}
    if (SimAnn == 1)
    {
        cout << "Running Simulated Annealing." << endl;
        someRandomBasis = MyOpt.randomBasis ();
        SimAnnPivots = MyOpt.SimulatedAnnealing(someRandomBasis,temperatures,times, SimAnnMinBases);
    }
    if (paretoBoundary == 1)
    {
        cout << "Finding boundary and computing pareto optimum." << endl;
        MyOpt.Boundary(ParetoOpt);
        ParetoOpt = MyOpt.ParetoOptimum(ParetoOpt);
    }
    if (paretoBoundaryTriangle == 1)
    {
        cout << "Finding boundary pareto points." << endl;
        MyOpt.Boundary(CHParetoOpt);
        CHParetoOpt = MyOpt.ParetoOptimum(CHParetoOpt);
        cout << "Finding Pareto Optimum." << endl;

        cout << "Computing triangular regions." << endl;
        TriangularRegions = MyOpt.BoundaryTrianglesTwoDim (CHParetoOpt);
        TriangularRegionsPrint = TriangularRegions;

        if (pivotTestType == 1) // Local Search
        {
            TriangularRegions = MyOpt.PivotTestLocalSearch(TriangularRegions,numPivots);
        }
        if (pivotTestType == 2) // Tabu Search
        {
            TriangularRegions = MyOpt.PivotTestTabuSearch(TriangularRegions,numPivots,pivotLimit);
        }

        // Union the boundary pareto points (ParetoOpt) with the resulting tests on triangular Regions

        set <Matrix, ltcolvec>::const_iterator mit = CHParetoOpt.begin();
        for (;mit != CHParetoOpt.end();mit++)
        {
            TriangularRegions.insert(*mit);
        }

        ParetoOpt = MyOpt.ParetoOptimum(TriangularRegions);
    }
    if (paretoBFS == 1)
    {
        cout << "Camputing BFSRandomStarts" << endl;
        MyOpt.MultiBFSRandomStarts(POnumSearches,POBFSDepth,POBoundaryRetryLimit,1,POretryLimit,POBFSProjPoints); 
        ParetoOpt = MyOpt.ParetoOptimum(POBFSProjPoints);

    }
    if (paretoBoundaryTriangleBFS == 1)
    {
        cout << "Finding boundary points." << endl;
        MyOpt.Boundary(CHParetoOpt);

        cout << "Camputing BFSRandomStarts" << endl;
        MyOpt.MultiBFSRandomStarts(POnumSearches,POBFSDepth,POBoundaryRetryLimit,1,POretryLimit,CHParetoOpt); 

        cout << "Finding Pareto Optimum." << endl;
        CHParetoOpt = MyOpt.ParetoOptimum(CHParetoOpt);

        cout << "Computing triangular regions." << endl;
        TriangularRegions = MyOpt.BoundaryTrianglesTwoDim (CHParetoOpt);

        if (pivotTestType == 1) // Local Search
        {
            TriangularRegions = MyOpt.PivotTestLocalSearch(TriangularRegions,numPivots);
        }
        if (pivotTestType == 2) // Tabu Search
        {
            TriangularRegions = MyOpt.PivotTestTabuSearch(TriangularRegions,numPivots,pivotLimit);
        }

        // Union the boundary pareto points (ParetoOpt) with the resulting tests on triangular Regions

        set <Matrix, ltcolvec>::const_iterator mit = CHParetoOpt.begin();
        for (;mit != CHParetoOpt.end();mit++)
        {
            TriangularRegions.insert(*mit);
        }

        ParetoOpt = MyOpt.ParetoOptimum(TriangularRegions);
    }
    if (calcAll == 1)
    {
        allBasesProj = MyOpt.calcAllProjBases();
        cout << "Number of projected bases: " << allBasesProj.size() << endl;
    } 




    time_t endTime = time(0);
    cout << "Total seconds: " << endTime - startTime << endl;
    if (MATLABOutput == 1)
    {
        ofstream outputfile;
        if (argc <= 2)
        {
            cout << "Writing to file: " << fileName << endl;
            outputfile.open(fileName);
        }
        else if (argc >= 3)
        {
            cout << "Writing to file: " << argv[2] << endl;
            outputfile.open(argv[2]);
        }
        outputfile << "%Total seconds: " << endTime - startTime << endl;

        outputfile << "hold on;" << endl;
        if (AutoBOXTest == 1)
        {
            string BOXLabel = "AUTOBOXLS";
            ProjBalMatroidOpt::writeBFSListMatlab(newPointsAutoBoxLS,outputfile,BOXLabel);
            outputfile << "plot(" << BOXLabel <<  "(:,1)," << BOXLabel <<  "(:,2),'r.');" << endl;
        }
        if (AutoBOXTestTS == 1)
        {
            string BOXLabel = "AUTOBOXTS";
            ProjBalMatroidOpt::writeBFSListMatlab(newPointsAutoBoxTS,outputfile,BOXLabel);
            outputfile << "plot(" << BOXLabel <<  "(:,1)," << BOXLabel <<  "(:,2),'r.');" << endl;
        }

        if (BOXTest == 1)
        {
            string BOXLabel = "BOXLS";
            ProjBalMatroidOpt::writeBFSListMatlab(newPointsBoxLS,outputfile,BOXLabel);
            outputfile << "plot(" << BOXLabel <<  "(:,1)," << BOXLabel <<  "(:,2),'r.');" << endl;
        }
        if (BOXTestTS == 1)
        {
            string BOXLabel = "BOXTS";
            ProjBalMatroidOpt::writeBFSListMatlab(newPointsBoxTS,outputfile,BOXLabel);
            outputfile << "plot(" << BOXLabel <<  "(:,1)," << BOXLabel <<  "(:,2),'r.');" << endl;
        }
        if (BFSEnumerate == 1)
        {
            string BFSLabel = "BFS";
            ProjBalMatroidOpt::writeBFSListMatlab(BFSProjPoints,outputfile,BFSLabel);
            outputfile << "plot(" << BFSLabel <<  "(:,1)," << BFSLabel <<  "(:,2),'k.');" << endl;
            //string POLabel = "PO";
            //ProjBalMatroidOpt::writeBFSListMatlab(ParetoOpt,outputfile,POLabel);
            //outputfile << "plot(" << POLabel <<  "(:,1)," << POLabel <<  "(:,2),'bo');" << endl;
        }
        if (ConvexHull == 1)
        {
            string CHLabel = "CH";
            ProjBalMatroidOpt::writeBFSListMatlab(CH,outputfile,CHLabel);
            outputfile << "plot(" << CHLabel <<  "(:,1)," << CHLabel <<  "(:,2),'gx');" << endl;

            //CH = MyOpt.ParetoOptimum(CH);
            //// Pareto Optimum code
            //CH = MyOpt.BoundaryTrianglesTwoDim (CH);
            //string POLabel = "PO";

            //ProjBalMatroidOpt::writeBFSListMatlab(CH,outputfile,POLabel);
            //outputfile << "hold on;" << endl;
            //outputfile << "plot(" << POLabel <<  "(:,1)," << POLabel <<  "(:,2),'r^');" << endl;

            //set <Matrix, ltcolvec> RPO = MyOpt.PivotTestLocalSearch(CH,1);
            //string RPOLabel = "RPO";
            //ProjBalMatroidOpt::writeBFSListMatlab(RPO,outputfile,RPOLabel);
            //outputfile << "hold on;" << endl;
            //outputfile << "plot(" << RPOLabel <<  "(:,1)," << RPOLabel <<  "(:,2),'bo');" << endl;

            //set <Matrix, ltcolvec> TPO = MyOpt.ParetoOptimum(CH);
            //string TPOLabel = "TPO";
            //ProjBalMatroidOpt::writeBFSListMatlab(TPO,outputfile,TPOLabel);
            //outputfile << "hold on;" << endl;
            //outputfile << "plot(" << TPOLabel <<  "(:,1)," << TPOLabel <<  "(:,2),'b+');" << endl;
        }
        //if (TabuSearch == 1)
        //{
        //    string TSPLabel = "TABUPIVOTS";
        //    MyOpt.ProjBalMatroidOpt::writePivotsMatlab(TabuSearchPivots,outputfile,TSPLabel);
        //    outputfile << "plot(" << TSPLabel <<  "(:,1)," << TSPLabel <<  "(:,2),'kx-');" << endl;
        //    string TSBLabel = "TABUBASES";
        //    MyOpt.ProjBalMatroidOpt::writePivotsMatlab(TabuSearchBases,outputfile,TSBLabel);
        //    outputfile << "plot(" << TSBLabel <<  "(:,1)," << TSBLabel <<  "(:,2),'bx-');" << endl;
        //    string CMPLabel = "CONVCOMP";
        //    MyOpt.ProjBalMatroidOpt::writePivotsMatlab(ConvexComparePivots,outputfile,CMPLabel);
        //    outputfile << "plot(" << CMPLabel <<  "(:,1)," << CMPLabel <<  "(:,2),'rx-');" << endl;
        //}
        if (SimAnn == 1)
        {
            string SALabel = "SIMANN";
            MyOpt.ProjBalMatroidOpt::writePivotsMatlab(SimAnnPivots,outputfile,SALabel);
            outputfile << "plot(" << SALabel <<  "(:,1)," << SALabel <<  "(:,2),'r.-');" << endl;
        }
        if (paretoBoundary == 1)
        {
            string POLabel = "PO";
            ProjBalMatroidOpt::writeBFSListMatlab(ParetoOpt,outputfile,POLabel);
            outputfile << "hold on;" << endl;
            outputfile << "plot(" << POLabel <<  "(:,1)," << POLabel <<  "(:,2),'gx');" << endl;
        }
        if (paretoBoundaryTriangle == 1)
        {
            string POLabel = "PO";
            ProjBalMatroidOpt::writeBFSListMatlab(ParetoOpt,outputfile,POLabel);
            outputfile << "hold on;" << endl;
            outputfile << "plot(" << POLabel <<  "(:,1)," << POLabel <<  "(:,2),'gx');" << endl;
            string TRLabel = "TO";
            ProjBalMatroidOpt::writeBFSListMatlab(TriangularRegionsPrint,outputfile,TRLabel);
            outputfile << "hold on;" << endl;
            outputfile << "plot(" << TRLabel <<  "(:,1)," << TRLabel <<  "(:,2),'b^');" << endl;
        }
        if (paretoBFS == 1)
        {
            string POLabel = "PO";
            ProjBalMatroidOpt::writeBFSListMatlab(ParetoOpt,outputfile,POLabel);
            outputfile << "hold on;" << endl;
            outputfile << "plot(" << POLabel <<  "(:,1)," << POLabel <<  "(:,2),'gx');" << endl;
        }
        if (paretoBoundaryTriangleBFS == 1)
        {
            string POLabel = "PO";
            ProjBalMatroidOpt::writeBFSListMatlab(ParetoOpt,outputfile,POLabel);
            outputfile << "hold on;" << endl;
            outputfile << "plot(" << POLabel <<  "(:,1)," << POLabel <<  "(:,2),'gx');" << endl;
        }
        if (calcAll == 1)
        {
            string APBLabel = "APB";
            ProjBalMatroidOpt::writeBFSListMatlab(allBasesProj, outputfile,APBLabel);
            outputfile << "hold on;" << endl;
            outputfile << "plot(" << APBLabel <<  "(:,1)," << APBLabel <<  "(:,2),'g*');" << endl;
        } 
    }
    else {
        //ProjBalMatroidOpt::writeBFSListMatlab(BFSProjPoints,cout);
        //if (BOXTest == 1)
        //{
        //    ProjBalMatroidOpt::writeBFSListMatlab(someProjPoints,cout);
        //}
        //
        //cout << "hold on;" << endl;
        //if (BOXTest == 1)
        //{
        //    cout << "plot(BFSLVL2(:,1),BFSLVL2(:,2),'r.');" << endl;
        //}
        //cout << "plot(BFSLVL1(:,1),BFSLVL1(:,2),'k.');" << endl;
    }
}
