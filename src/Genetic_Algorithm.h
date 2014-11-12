#ifndef _GENETIC_ALGORITHM_H_
#define _GENETIC_ALGORITHM_H_

// Functions for GA
int optimizeDecisivenessGA (vector < vector <int> > const& data,
    vector < vector <int> > const& legal, int maxNumGenes);
int countLociPresent (vector < vector <int> > const& data);
void mutate ( vector < vector <int> > & data, vector < vector <int> > const& legal,
    bool add, bool remove );
void xOver ( vector < vector <int> > & recipient, vector < vector <int> > const& donor);
double calculateFitness( vector < vector <int> > const& individual);



#endif /* _GENETIC_ALGORITHM_H_ */
