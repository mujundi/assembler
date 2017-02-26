/************************************
* ACO.cpp
* Ant Colony Optimization algorithm
*
* Definition file for the ACO class
*************************************/

#include "ACO.h"

using namespace std;


ACO::ACO(int nAnts, int nNodes, double a, double b,
			double q, double r, double t, 
			int start, vector<vector<int>> X, int k): 
							numOfAnts(nAnts),
							numOfNodes(nNodes),
							startingNode(start),
							alpha(a), beta(b), Q(q), rho(r), tau(t),
							graph(X), trial(k)
{
	srand(time(NULL));
	long int seed = 5 + rand()%46;
	randoms = new Randoms(seed);
}

ACO::~ACO() 
{ 
	delete randoms; 
}


void ACO::init() 
{
	pheromones.resize(numOfNodes, vector<double>(numOfNodes, 0.0));
	deltaPheromones.resize(numOfNodes, vector<double>(numOfNodes, 0.0));
	probabilites.resize(numOfNodes, vector<double>(2, -1.0));
	paths.resize(numOfAnts, vector<int>(numOfNodes, -1));
	bestPath.resize(numOfNodes, -1);
	bestLength = (double) INT_MIN;

	for (int i=0; i<numOfNodes; i++)
	{
		for (int j=0; j<numOfNodes; j++)
		{
			if ( !(graph[i][j] == INT_MIN) )
			{
				pheromones[i][j] = randoms -> Uniform() * tau;
			}

		}
	}
}


void ACO::printResults() 
{
	cout << " BEST ROUTE:" << endl;

	for (int i=0; i<numOfNodes; i++) 
	{
		cout << bestPath[i] << " ";
	}

	cout << endl << "Length: " << bestLength << endl;
}


bool ACO::visited(int antk, int node) 
{
	for (int i=0; i<numOfNodes; i++) 
	{
		if (paths[antk][i] == -1) { break; }
		if (paths[antk][i] == node) { return true; }
	}

	return false;
}


bool ACO::exists(int nodei, int nodej)
{
	return (!(graph[nodei][nodej] == INT_MIN));
}


double ACO::PHI(int nodei, int nodej, int antk) 
{
	double eta_ij = (double) pow(graph[nodei][nodej], beta);
	double tau_ij = (double) pow(pheromones[nodei][nodej], alpha);

	double sum = 0.0;

	for (int c=0; c<numOfNodes; c++) 
	{
		if (exists(nodei, c)) 
		{
			if (!visited(antk, c)) 
			{
				double eta = (double) pow(graph[nodei][c], beta);
				double tau = (double) pow(pheromones[nodei][c], alpha);
				sum += eta * tau;
			}	
		}
	}
	return (eta_ij * tau_ij) / sum;
}


double ACO::length(int antk) 
{
	double sum = 0.0;

	for (int j=0; j<numOfNodes-1; j++) 
	{
		sum += graph[paths[antk][j]][paths[antk][j+1]];	
	}

	return sum;
}


int ACO::valid(int antk, int iteration) 
{
	for(int i=0; i<numOfNodes-1; i++) 
	{
		int nodei = paths[antk][i];
		int nodej = paths[antk][i+1];

		if (nodei < 0 || nodej < 0) { return -1; }
		if (!exists(nodei, nodej)) { return -2; }

		for (int j=0; j<i-1; j++) 
		{
			if (paths[antk][i] == paths[antk][j]) { return -3; }	
		}
	}
	
	if (!exists (startingNode, paths[antk][numOfNodes-1])) 
	{
		return -4;
	}
	
	return 0;
}


int ACO::node() 
{
	double xi = randoms -> Uniform();
	int i = 0;
	double sum = probabilites[i][0];

	while (sum < xi)
	{
		i++;
		sum += probabilites[i][0];
	}

	return (int) probabilites[i][1];
}


void ACO::route(int antk) 
{
	paths[antk][0] = startingNode;

	for (int i=0; i<numOfNodes-1; i++) 
	{
		int nodei = paths[antk][i];
		int count = 0;

		for (int c=0; c<numOfNodes; c++) 
		{
			if (nodei == c) 
			{
				continue;	
			}

			if (exists(nodei, c)) 
			{
				if (!visited(antk, c)) 
				{

					probabilites[count][0] = PHI(nodei, c, antk);
					probabilites[count][1] = (double) c;
					count++;
				}
			}
		}
		
		// deadlock
		if (0 == count) { return; }
		
		paths[antk][i+1] = node();
	}	
}


void ACO::updatePheromones() 
{
	for (int k=0; k<numOfAnts; k++) 
	{
		double rlength = length(k);
		for (int r=0; r<numOfNodes-1; r++) 
		{
			int nodei = paths[k][r];
			int nodej = paths[k][r+1];

			deltaPheromones[nodei][nodej] += rlength / Q;
			deltaPheromones[nodej][nodei] += rlength / Q;
		}
	}

	for (int i=0; i<numOfNodes; i++) 
	{
		for (int j=0; j<numOfNodes; j++) 
		{
			pheromones[i][j] = (1 - rho) * pheromones[i][j] + deltaPheromones[i][j];
			deltaPheromones[i][j] = 0.0;
		}	
	}
}


void ACO::optimize (int N) 
{
	string fileName;
	if (trial == 1) fileName = "ConvergenceData1.csv";
	else if (trial == 2) fileName = "ConvergenceData2.csv";
	else if (trial == 3) fileName = "ConvergenceData3.csv";

	ofstream csvfile(fileName);

	for (int iterations=1; iterations<=N; iterations++) 
	{
		cout << flush;
		cout << "ITERATION " << iterations << " HAS STARTED!" << endl << endl;

		for (int k=0; k<numOfAnts; k++) 
		{
			cout << " : ant " << k+1 << " has been released!" << endl;
			
			while (0 != valid(k, iterations)) 
			{
				for (int i=0; i<numOfNodes; i++) 
				{
					paths[k][i] = -1;	
				}

				route(k);
			}
			
			for (int i=0; i<numOfNodes; i++) 
			{
				cout << paths[k][i] << " ";	
			}
			
			cout << endl << "  :: route done" << endl;

			double rlength = length(k);

			if (rlength > bestLength) 
			{
				bestLength = rlength;

				for (int i=0; i<numOfNodes; i++) 
				{
					bestPath[i] = paths[k][i];
				}
			}
			cout << " ::: ant " << k+1 << " has ended!" << endl;				
		}	

		cout << endl << "updating PHEROMONES . . .";
		updatePheromones ();
		cout << " done!" << endl << endl;
		
		for (int i=0; i<numOfAnts; i++) 
		{
			for (int j=0; j<numOfNodes; j++) 
			{
				paths[i][j] = -1;
			}
		}

		writeToCSV(csvfile, iterations, N);
		cout << endl << "ITERATION " << iterations << " HAS ENDED!" << endl << endl;
	}

	csvfile.close();
}


vector<int> ACO::getBestPath()
{
	return bestPath;
}


void ACO::writeToCSV(ofstream& csvfile, int currentIteration, int iterations)
{

	if (currentIteration == 1)
	{
		csvfile << "Number of Ants,Iterations,Alpha,Beta,Q,Rho,Tau" << endl;
		csvfile << numOfAnts << "," << iterations << ",";
		csvfile << alpha << "," << beta << "," << Q << "," << rho << "," << tau << endl;
		csvfile << endl;

		csvfile << "Iteration,Best Overlap" << endl;
		csvfile << currentIteration << "," << bestLength << endl;
	}

	else
	{
		csvfile << currentIteration << "," << bestLength << endl;
	}

}