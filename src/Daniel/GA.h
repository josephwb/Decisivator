#ifndef _GA_H_
#define _GA_H_

using namespace std;

// Functions for GA
int optimizeDecisivenessGA (vector < vector <int> > const& data,
	vector < vector <int> > const& legal, int maxGenes);
int countLociPresent (vector < vector <int> > const& data);
void mutate ( vector < vector <int> > & data, vector < vector <int> > const& legal,
	bool add, bool remove );
void xover ( vector < vector <int> > & recipient, vector < vector <int> > const& donor);
double calcFit( vector < vector <int> > const& individual);



#endif /* _GA_H_ */
