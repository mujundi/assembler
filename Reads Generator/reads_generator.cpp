/************************************
* reads_generator.cpp
*
* Generates fragments from genome FASTA file
*
*************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <vector>
#include <functional>
#include <algorithm>
#include <chrono>


using namespace std;


int main(void)
{
	auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
	
	mt19937 rng(seed);
	uniform_int_distribution<int> uni(5,95);
	auto random_int = bind(uni, rng);


	string sequence;
	ifstream infile("ecoli.fasta");

	cout << "Reading genome data from 'ecoli.fasta'..." << endl;

	if (infile.is_open())
	{
		getline(infile,sequence);	
		getline(infile,sequence);

		if (!sequence.empty() && *sequence.rbegin() == '\r') 
		{
    		sequence.erase( sequence.length()-1, 1);
		}		

		long int index = 0;
		bool anyBasesLeft = true;
		vector<string> reads;

		cout << "Generating simulated reads..." << endl;

		while (anyBasesLeft)
		{
			if (index >= sequence.size() - 100)
			{
				index = sequence.size() - 100;
				anyBasesLeft = false;
			}

			reads.push_back(sequence.substr(index,100));
			index += random_int();
		}

		random_shuffle(reads.begin(), reads.end());

		ofstream readsfile("ecoli_reads.fasta");

		if (readsfile.is_open())
		{
			for (long int i=0; i<reads.size(); i++)
			{
				if (i == reads.size()-1)
					readsfile << reads[i];
				else
					readsfile << reads[i] << endl;
			}
		}

		cout << "Output stored in 'ecoli_reads.fasta'." << endl;

		readsfile.close();
	}


	infile.close();



	return 0;
}