/************************************
* ACO.h
* Ant Colony Optimization algorithm
*
* Header file for the ACO class
*************************************/

#ifndef ACO_H
#define ACO_H
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>
#include <limits>
#include <climits>
#include "randoms.cpp"


class ACO
{
private:
	int numOfAnts, numOfNodes, startingNode;
	int trial;
	double alpha, beta, Q, rho, tau;
	double bestLength;

	std::vector<int> bestPath;
	std::vector<std::vector<int>> graph;
	std::vector<std::vector<int>> paths;
	std::vector<std::vector<double>> pheromones, deltaPheromones, probabilites;
	Randoms* randoms;

	void updatePheromones();
	int valid(int antk, int iteration);
	void route(int antk);
	int node();
	double length(int antk);
	double PHI(int nodei, int nodej, int antk);
	bool visited(int antk, int c);
	bool exists(int nodei, int nodej);
	void writeToCSV(std::ofstream& csvfile, int currentIteration, int iterations);

public:
	ACO(int nAnts, int nNodes, double alpha, double beta,
			double Q, double rho, double tau, int startingNode, std::vector<std::vector<int>> graph, int trial);
	~ACO();
	void init();
	void printPheromones();
	void printResults();
	void optimize(int iterations);
	std::vector<int> getBestPath();
};