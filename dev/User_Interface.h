#ifndef _USER_INTERFACE_H_
#define _USER_INTERFACE_H_

void printProgramInfo();
void printHelp ();
void printMatrix (vector < vector <int> > const& data, vector <string> const& taxonNames, 
	vector <double> const& locusWeights, vector <double> const& taxonWeights);
void printProgamOptions (bool & addGenes, bool & merge, bool & exclude, bool & deleteGenes, bool & revert,
	bool & quit, bool & print, bool & reweightLoci, bool & reweightTaxa, bool & partialTreewise,
	bool & partialBranchwise, bool & summarize, bool & testCompleteDeciveness,
	bool & writeCurrentMatrix, bool & testUserTree);
void printSummaryInformation (vector <string> const& locusNames, vector <string> const& taxonNames,
	vector < vector <int> > const& data, double const& taxonCoverage, vector <int> const& referenceTaxa, bool const& matrixDecisive,
	double const& treewiseDecisiveness, double const& branchwiseDecisiveness, bool const& completeDecisivenessDetermined, string const& nexusFileName,
	int const& numRandomTrees, int const& numUserTrees);
void processCommandLineArguments(int argc, char *argv[], string & matrixFileName,
	string & nexusFileName, string & locusWeightFileName, string & taxonWeightFileName,
	string & treeFileName, int & burnin, int & thinning);
void writeMatrix (vector <string> const& taxonNames, int const& numChar, 
	vector < vector <string> > const& taxaAlignment, vector < vector <int> > const& includedLocusRanges,
	vector <string> const& locusNames);
void printMatrixToFile (vector < vector <int> > const& data, vector <string> const& taxonNames, 
	vector <double> const& locusWeights, vector <double> const& taxonWeights);
void writeNexus (int const& numTaxa, int const& numChar, vector <string> const& taxonNames,
	string & fileName, vector < vector <string> > const& taxaAlignment, vector < vector <int> > const& includedLocusRanges,
	vector <string> const& locusNames);
void writePhylip (int const& numTaxa, int const& numChar, vector <string> const& taxonNames,
	string & fileName, vector < vector <string> > const& taxaAlignment);

#endif /* _USER_INTERFACE_H_ */