// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE


/// \file estimatebases.cpp

#include "matroid.h"
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>


int main (int argc, char *argv[])
{
	int m, num_rows, i, k, max_weight;
	float upper_bound, lower_bound, GAMMA, avg_GAMMA;
	string filename;
    string matroidType;
	ifstream inFile;
	set <unsigned> max_weight_basis;

    cout << "argc = " << argc << endl;


	cout << "This program estimates the number of bases of a matrix by assigning\n";
	cout << "random weights to each column of the matrix and finding the maximal\n";
	cout << "weight basis. This is done m times and the resulting average estimates\n";
	cout << "the function GAMMA. Then upper and lower bounds for the number of bases\n";
	cout << "are calculated.\n\n";

    if (argc > 1)
    {
        sscanf(argv[1],"%d",&m);
        cout << "m = " << m << endl;
    }
    else
    {
        cout << "Enter the number of times m to sample: \n";
        cin >> m;
        cout << "\n";
    }
    if (m <= 0)
    {
        cout << "m <= 0     exiting" << endl;
        exit(0);
    }

	//while(m<=0){
	//	cout << "This number must be greater than zero.\n";
	//	cout << "Enter the number of times m to sample: \n";
	//	cin >> m;
	//	cout << "\n";
	//}

	//cout << "Enter the number of rows in the matrix.\n";
	//cin >> num_rows;
	//cout << "\n";

	//cout << "Enter the number of columns in the matrix.\n";
	//cout << "This must be more than the number of rows.\n";
	//cin >> num_cols;

	//int matrix[num_rows][num_cols];

    Matroid *M = 0;

    if (argc > 2)
    {
        filename = argv[2];
    }
    else 
    {
        cout << "Enter the name of the input file:\n";
        cin >> filename;
    }

	inFile.open(filename.c_str(), ios::in);

	while(!inFile){
		cout << "File open error: '" << filename << "' ";
		cout <<"Try again.\n";
		inFile.close();
		cout << "Enter the name of the input file:\n";
		cin >> filename;
		inFile.open(filename.c_str(), ios::in);
	}
    if (argc > 3)
    {
        //sscanf(argv[3],"%d",&matroidType);
        matroidType = argv[3];
    }
    else
    {
        cout << "Matroid type (1) Vector, (2) Graphical: ";
        cin >> matroidType;
    }
    if (matroidType == "1")
    {
        M = new VectorMatroid(inFile);
    }
    else if (matroidType == "2")
    {
        M = new GraphicalMatroid(inFile);
    }
    cout << *M;
    num_rows = M->getNumElements();
	Matrix weights(num_rows,1);


	// read in the matrix from inFile
	//for(i=0; i<num_rows; i++){
	//	for(j=0; j<num_cols; j++){
	//		inFile >> matrix[i][j];
	//	}
	//}
    

	GAMMA=0;
	srand((unsigned)time(NULL));

    int ten_percent = (int)floor(m*0.10);
    ten_percent = (int)max((double)1,(double)ten_percent);
    cout << "ten_percent = " << ten_percent << endl;
	for(k=0;k<m;k++){
        if (k % ten_percent == 0 && k != 0)
        {
            cout << "   " << (k/ten_percent)*10 << "%" << endl;
        }
		max_weight_basis.clear();
  
        //cout << "Zeroing weights." << endl;
		for(i = 0; i < num_rows; i++){
		    weights(i,0)=0;
		}
        //cout << "calling random_weight_logistic with modified_rand." << endl;
		for(i = 0; i < num_rows; i++){
            // I'm guessing you want the input to logistic to be (0,1)
			//weights(i,0)=M->random_weight_logistic((float)M->modified_rand()/MAXIMAL_RAND);
			weights(i,0)=M->random_weight_logistic((float) ((rand() % 99999) + 1) / 100000);
		}
  
        //cout << "Calling GreedyAlgorithmMax." << endl;
		max_weight_basis = M->GreedyAlgorithmMax(weights);
  
		set<unsigned>::iterator it;
		max_weight = 0;
        //cout << "Adding to max_weight." << endl;
		for(it = max_weight_basis.begin(); it != max_weight_basis.end(); it++)
		{
			max_weight += (int)weights(*it,0);
		}
  
		GAMMA=GAMMA+(float)max_weight;
	}
    cout << "  100%" << endl;

	//compute upper and lower estimates
	avg_GAMMA = (float)GAMMA/(float)m;
	upper_bound = M->upper_logistic(avg_GAMMA, num_rows);
	lower_bound = M->lower_logistic(avg_GAMMA, num_rows);

	//Output results
	cout << "\n";
	cout << "\n";
	cout << "Bounds for Log(X) are:\n";
	cout << "Upper bound:\n";
	printf("%f\n",upper_bound);
	cout << "Lower bound:\n";
	printf("%f\n",lower_bound);
	cout << "GAMMA(X):\n";
	printf("%f\n",avg_GAMMA);

}
