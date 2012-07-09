#ifndef _GA_H_
#define _GA_H_

using namespace std;

// Functions for GA
vector < vector <int> > GAHandler (int const& numAddGA, vector < vector <int> > const& data,
	double & GADecisiveness, bool const& referenceTaxonPresent, int const& numProcs);
int optimizeDecisivnessGA (vector < vector <int> > & data,
	vector < vector <int> > const& legal, int const& numAdd, bool const& referenceTaxonPresent,
	int const& numProcs);
int countLociPresent (vector < vector <int> > const& data);
void mutate (vector < vector <int> > & data, vector < vector <int> > const& legal,
	bool const& add, bool const& remove);
void xover (vector < vector <int> > & recipient, vector < vector <int> > const& donor,
	vector < vector < int > > legal, int maxGenes);
//double calcFit (vector < vector <int> > const& individual);
double calcFit (bool const& referenceTaxonPresent, int & numTrees,
	vector < vector <int> > const& data, int const& numProcs);
int stepAdd (vector < vector < int > > & data, vector < vector < int > > const& legal,
	int maxGenes);
int stepReplace (vector < vector < int > > & data, vector < vector < int > > const& legal,
	int const& maxGenes, bool const& referenceTaxonPresent, int const& numProcs);
void printData (vector < vector <int> > const& data);
vector < vector <int> > getLegalMatrix (vector < vector <int> > const& data);
void printGADataToFile (vector < vector <int> > const& data, vector <string> const& taxonNames,
	string const& nexusFileName, int & numAddGA);

#endif /* _GA_H_ */