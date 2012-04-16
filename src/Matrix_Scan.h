#ifndef _MATRIX_SCAN_H_
#define _MATRIX_SCAN_H_

//extern bool debuggering;

void checkForMissingTaxa (vector < vector <int> > const& data, vector <string> const& taxonNames);
bool searchForReferenceTaxon (vector < vector <int> > & data, vector <int> & referenceTaxa,
	vector <string> & taxonNames);
void searchForAllTriplets (vector < vector <int> > const& data, vector < vector <int> > & triplets,
	vector < vector <int> > & tripletLocations, vector < vector <int> > & missingTriplets);
void searchForAllQuartets (vector < vector <int> > const& data, vector < vector <int> > & missingQuartets);
void searchForQuartetsWithReference (vector < vector <int> > const& data, vector < vector <int> > & missingQuartets,
	vector <int> const& referenceTaxa, vector <string> const& taxonNames);
//void searchForQuartetsWithReferenceBIG (vector < vector <int> > const& data, vector <int> & missingQuartetsByTaxa,
//	vector <int> const& referenceTaxa, vector <string> const& taxonNames);
void whichTaxaProblematic (vector < vector <int> > const& missingGroups, vector <string> const& taxonNames,
	string const& grouping, vector <int> const& referenceTaxa);
//void whichTaxaProblematicBIG (vector <int> const& missingQuartetsByTaxa, vector <string> const& taxonNames,
//	string const& grouping, vector <int> const& referenceTaxa);
void getCoverage (vector < vector <int> > const& data, double & taxonCoverage);
double calculatePartialDecisiveness (bool const& referenceTaxonPresent, int & numTrees,
	vector < vector <int> > const& data, bool const& findAll, int const& numProcs);
int searchEdgePartitions (vector < vector <int> > const& data, vector <int> const& left,
	vector <int> const& right, vector <int> const& sib, vector <int> const& upper, bool const& findAll,
	bool const& referenceTaxonPresent);
bool testCompleteDecisivness (vector < vector <int> > const& data, bool const& referenceTaxonPresent,
	vector <int> const& referenceTaxa, vector <string> const& taxonNames,
	vector < vector <int> > & missingQuartets, vector < vector <int> > & triplets,
	vector < vector <int> > & tripletLocations, vector < vector <int> > & missingTriplets);
vector < vector <double> > determineDecisivenessUserTree (vector < vector <int> > const& data,
	vector < vector < vector <bool> > > & userTrees, vector < vector <int> > const& treeTaxonOrdering,
	vector <string> const& taxonNames, vector <double> const& locusWeights, vector <double> const& taxonWeights);
void printBipartitionTable (vector < vector <bool> > const& tree, vector <double> const& decisiveness,
	vector <int> const& numSatisfied, vector <int> const& numPossible, int const& numTaxa,
	int const& numTrees, int const& treeNumber);
double calculatePartialDecisivenessSinglePartition (bool const& referenceTaxonPresent, int & numTrees,
	vector < vector <int> > const& data, bool const& findAll, int const& partitionID,
	vector <string> const& locusNames);

#endif /* _MATRIX_SCAN_H_ */