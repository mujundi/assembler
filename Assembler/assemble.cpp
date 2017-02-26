/****************************
* SHOTGUN SEQUENCING (ACO).
*
* @author Musa Jundi
* @version 0.0 2/7/2017
*
* C++11
*
*****************************/


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <omp.h>
#include "DNA.h"
#include "ACO.h"

using namespace std;

vector<vector<int>> createGraph(vector<vector<int>>);
void deleteHamNode(vector<vector<int>>& adjMatrix);
void printGraph(vector<vector<int>> x, ofstream& outfile);
void simplePrint(vector<vector<int>> x);

int main(void)
{

	#pragma omp parallel for schedule(dynamic, 1) ordered
	for (int trial=1; trial<4; trial++)
	{
		DNA dna;

		dna.findOverlaps(); 
		vector<vector<int>> overlaps = dna.getOverlaps();
		vector<vector<int>> adjMatrix = createGraph(overlaps);

		int nAnts = 100;
		int nNodes = sqrt(overlaps.size()) + 1;
		double alpha = 0.8;
		double beta = 1.0;
		double Q = 4700;
		double rho = 0.15;
		double tau = 3;
		int startingNode = 0;
		int iterations = 1000;


		ACO ants(nAnts,nNodes,alpha,beta,Q,rho,tau,startingNode,adjMatrix,trial);
		ants.init();

		auto start = chrono::high_resolution_clock::now();
		ants.optimize(iterations);
		auto finish = chrono::high_resolution_clock::now();
		chrono::duration<double> elapsed = finish - start;


		ants.printResults();

		vector<int> path = ants.getBestPath();
		dna.setSequenceOrder(path);
		string mergedSequence = dna.combineFragments(adjMatrix);

		#pragma omp ordered
		{
			ofstream outfile("mergedDNA.txt", ios::app);
			if (outfile.is_open()) 
			{
				outfile << "TRIAL " << trial << endl;
				outfile << "-------" << endl;
				outfile << "Elapsed time: " << elapsed.count() << " s" << endl;
				outfile << "Sequence length: " << mergedSequence.length() << endl << endl;
				outfile << endl << mergedSequence << endl;
				outfile << endl << endl;   
			}

			outfile.close();
		}
	}

	return 0;
}



vector<vector<int>> createGraph(vector<vector<int>> x)
{
	vector<vector<int>> overlaps = x;
	int N = sqrt(overlaps.size());

	vector<vector<int>> adjMatrix;

	adjMatrix.resize(N, vector<int>(N, INT_MIN));

	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			adjMatrix[i][j] = overlaps[i*N+j][2];
		}
	}

	vector<int> a;
	a.resize(N, 1);
	adjMatrix.insert(adjMatrix.begin(), a);

	for (int i=0; i<N+1; i++)
	{
		adjMatrix[i].insert(adjMatrix[i].begin(), 1);
	}

	return adjMatrix;
}



void deleteHamNode(vector<vector<int>>& adjMatrix)
{
	int N = adjMatrix.size();
	adjMatrix.erase(adjMatrix.begin());
	for (int i=0; i<N; i++)
		adjMatrix[i].erase(adjMatrix[i].begin());
}



void simplePrint(vector<vector<int>> x)
{
	int N = x.size();
	int M = x[0].size();
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<M; j++)
		{
			if (x[i][j]==INT_MIN) cout << "x ";
			else cout << x[i][j] << " ";
		}
		cout << endl;
	}
}



void printGraph(vector<vector<int>> x, ofstream& outfile)
{
	int N = x.size();
	outfile << " GRAPH: " << endl;
	outfile << "  | ";

	for( int i=0; i<N; i++) {
		outfile << i << " ";
	}

	outfile << endl << "- | ";

	for (int i=0; i<N; i++) {
		outfile << "- ";
	}

	outfile << endl;
	int count = 0;

	for (int i=0; i<N; i++) {
		outfile << i << " | ";
		for (int j=0; j<N; j++) {
			if(i == j) {
				outfile << "x ";	
			}
			else if (x[i][j] == INT_MIN)
			{
				outfile << "x ";
			}
			else {
				outfile << x[i][j] << " ";	
			}
			if ( !(x[i][j] == INT_MIN) ) {
				count++;	
			}
		}
		outfile << endl;
	}

	outfile << endl;
	outfile << "Number of connections: " << count << endl << endl << endl;
}
