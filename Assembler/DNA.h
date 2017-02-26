/************************************
* DNA.h
*
* Header file for the DNA class
*************************************/

#ifndef DNA_H
#define DNA_H
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <climits>


class DNA
{
private:
	std::string sequence;
	std::vector<int> sequenceOrder;
	std::vector<std::string> fragments;
	std::vector<std::vector<int>> overlaps;
	int N;

	void readfile();
	int LCS(std::string a, std::string b);

public:
	DNA();
	std::string combineFragments(std::vector<std::vector<int>> X);
	void findOverlaps();
	void setSequenceOrder(std::vector<int> path);
	std::vector<std::vector<int>> getOverlaps();
};

