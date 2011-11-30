#ifndef _MATRIX_SCAN_H_
#define _MATRIX_SCAN_H_

bool searchForReferenceTaxon (vector < vector <int> > & data, vector <int> & referenceTaxa,
	vector <string> & taxonNames);
void searchForAllTriplets (vector < vector <int> > const& data, vector < vector <int> > & triplets,
	vector < vector <int> > & tripletLocations, vector < vector <int> > & missingTriplets);
void searchForAllQuartets (vector < vector <int> > const& data, vector < vector <int> > & missingQuartets);
void searchForQuartetsWithReference (vector < vector <int> > const& data, vector < vector <int> > & missingQuartets,
	vector <int> const& referenceTaxa, vector <string> const& taxonNames);
void searchForQuartetsWithReferenceBIG (vector < vector <int> > const& data, vector <int> & missingQuartetsByTaxa,
	vector <int> const& referenceTaxa, vector <string> const& taxonNames);
void whichTaxaProblematic (vector < vector <int> > const& missingGroups, vector <string> const& taxonNames,
	string const& grouping, vector <int> const& referenceTaxa);
void whichTaxaProblematicBIG (vector <int> const& missingQuartetsByTaxa, vector <string> const& taxonNames,
	string const& grouping, vector <int> const& referenceTaxa);
void getCoverage (vector < vector <int> > & data, double & taxonCoverage);
double calculatePartialDecisiveness (bool const& referenceTaxonPresent, int & numTrees,
	vector < vector <int> > const& data, bool const& findAll);
int searchEdgePartitionsMinimum (vector < vector <int> > const& data, vector <int> const& left,
	vector <int> const& right, vector <int> const& sib, vector <int> const& upper,
	bool const& referenceTaxonPresent, int & partitionMatched);
int searchEdgePartitionsAll (vector < vector <int> > const& data, vector <int> const& left,
	vector <int> const& right, vector <int> const& sib, vector <int> const& upper,
	bool const& referenceTaxonPresent);
bool testCompleteDecisivness (vector < vector <int> > const& data, bool const& referenceTaxonPresent,
	vector <int> const& referenceTaxa, vector <string> const& taxonNames, vector < vector <int> > & missingQuartets,
	vector < vector <int> > & triplets, vector < vector <int> > & tripletLocations, vector < vector <int> > & missingTriplets);
vector < vector <double> > determineDecisivenessUserTree (vector < vector <int> > const& data,
	vector < vector < vector <bool> > > & userTrees, vector < vector <int> > const& treeTaxonOrdering,
	vector <string> const& taxonNames, vector <double> const& locusWeights, vector <double> const& taxonWeights);
void printBipartitionTable (vector < vector <bool> > const& tree, vector <double> const& decisiveness,
	vector <int> const& numSatisfied, vector <int> const& numPossible, int const& numTaxa,
	int const& numTrees, int const& treeNumber);
int searchEdgeStoredClades (vector <int> const& newClade,
	vector < vector <int> > const& storedClades);
double calculatePartialDecisivenessNEW (bool const& referenceTaxonPresent, int & numTrees,
	vector < vector <int> > const& data, bool const& findAll);
bool checkCladeSubset (vector <int> const& newClade, vector <int> const& storedClade);
bool searchEdgeOnePartition (vector < vector <int> > const& data, vector <int> const& clade,
	int const& partitionToMatch);
string getPartitionMatches (vector < vector <int> > const& data, vector <int> const& left,
	vector <int> const& right, vector <int> const& sib, vector <int> const& upper,
	int const& partitionToMatch, bool & leftValid, bool & rightValid, bool & sibValid,
	bool & upperValid);
	
#endif /* _MATRIX_SCAN_H_ */