#ifndef _GA_H_
#define _GA_H_

#include <iostream>
#include <vector>
//#include <ctime>
#include <cstdlib>
#include <numeric>
//#include <functional>

using namespace std;

// Functions for GA
int optimizeDecisivnessGA (vector < vector <int> > const& data,
	vector < vector <int> > const& legal, int maxGenes);
int countLociPresent (vector < vector <int> > const& data);
void mutate (vector < vector <int> > & data, vector < vector <int> > const& legal,
	bool const& add, bool const& remove);
void xover (vector < vector <int> > & recipient, vector < vector <int> > const& donor,
	vector < vector < int > > legal, int maxGenes);
double calcFit (vector < vector <int> > const& individual);
int stepAdd (vector < vector < int > > & data, vector < vector < int > > const& legal,
	int maxGenes);
int stepReplace (vector < vector < int > > & data, vector < vector < int > > const& legal,
	int maxGenes);
void printData (vector < vector <int> > const& data);
vector < vector <int> > getLegalMatrix (vector < vector <int> > const& data);

#endif /* _GA_H_ */