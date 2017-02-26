/************************************
* DNA.cpp
*
* Definition file for the DNA class
*************************************/

#include "DNA.h"

using namespace std;


DNA::DNA()
{
	#pragma omp critical
	readfile();
	
	N = fragments.size();
	int numRows = fragments.size() * fragments.size();
	overlaps.resize(numRows, vector<int>(3, INT_MIN));
}


string DNA::combineFragments(vector<vector<int>> X)
{
	int a;
	int b;
	int c;
	string a_seq;
	string b_seq;

	for(int i=0; i<N; i++)
	{
		if (i==0) 
		{
			a_seq = "";
			b = sequenceOrder[i];
		}
		else
		{
			a = sequenceOrder[i-1];
			b = sequenceOrder[i];

			c = X[a][b];

			a_seq.erase(a_seq.length() - c);
		}

		b_seq = fragments[b-1];

		sequence = a_seq + b_seq;
		a_seq = sequence;

	}

	return sequence;	

}


int DNA::LCS(string a, string b)
{
	int m = a.length() - 1;
	int n = b.length() - 1;
	int o = (n < m) ? n : m;

	int counter = 0;
	int maxLength = INT_MIN;


	for (int i=0; i <= o; i++)
	{
		for (int j=0; j <= i; j++)
		{
			if (a[m-i+j] == b[j]) { counter += 1; }
			else
			{
				counter = 0;
				break;
			}
		}

		if (counter > maxLength) { maxLength = counter; }

		counter = 0;
	}

	return maxLength;

}


void DNA::findOverlaps()
{
	int x;
	int lcs;
	int max = fragments.size();

	for (int i=0; i < max; i++)
	{
		for (int j=0; j < max; j++)
		{
			if (i==j) overlaps[i*max+j] = {i, j, INT_MIN};
			else
			{
				lcs = LCS(fragments[i],fragments[j]);
				x = (lcs > 0) ? lcs : INT_MIN;

				overlaps[i*max+j] = {i, j, x};
			}
		}
	}
}


void DNA::readfile()
{
	string line;

	ifstream infile("ecoli_reads.fasta");

	if (infile.is_open())
	{
		while (getline(infile,line))
		{
			if (!line.empty() && *line.rbegin() == '\r') 
			{
    			line.erase( line.length()-1, 1);
			}			
			
			fragments.push_back(line);
		}
	}

	infile.close();
}


void DNA::setSequenceOrder(vector<int> path)
{
	sequenceOrder.assign(path.begin()+1, path.end());
}

vector<vector<int>> DNA::getOverlaps()
{
	return overlaps;
}