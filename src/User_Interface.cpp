#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

#include "User_Interface.h"
#include "General.h"
#include "Parse_Nexus.h"
#include "Matrix_Scan.h"

extern bool debugging; // print out extra junk to screen
extern double version;
extern string month;
extern int year;

void printProgramInfo() {
	cout << endl << 
	"************************************************" << endl <<
	"            Decisivator version " << version <<      endl <<
	"                A PuRGe Product"                  << endl <<
	"              University of Idaho"                << endl <<
	"       Department of Biological Sciences"         << endl <<
	"        Complaints: josephwb@umich.edu"           << endl <<
	"                  " << month <<", " << year <<       endl << 
	"************************************************" << endl << endl;
}

void printProgamOptions (bool & addGenes, bool & merge, bool & exclude, bool & deleteGenes, bool & revert,
	bool & quit, bool & print, bool & reweightLoci, bool & reweightTaxa, bool & partialTreewise,
	bool & partialBranchwise, bool & summarize, bool & testCompleteDeciveness,
	bool & writeCurrentMatrix, bool & testUserTree, bool & partialIndividualPartition,
	bool & printRefTaxa, bool & useGA)
{
	bool validChoice = false;
	char userChoice;
	
	while (!validChoice) {
		cout << endl
		<< "Programs options:" << endl
		<< " [A]dd virtual taxon-character(s) to matrix manually" << endl
		<< " [M]erge taxa (i.e create a chimeric taxon)" << endl
		<< " [E]xclude taxa" << endl
		<< " [D]elete partition(s)" << endl
		<< " [P]rint current matrix to the screen" << endl
		<< " [C]alculate partial decisiveness using random trees" << endl
		<< " [T]est of complete decisiveness (can take a while for large taxon samples)" << endl
		<< " [I]nvestigate decisiveness on a provided user-tree" << endl
		<< " [O]ptimally add virtual taxon-character(s) to matrix using a genetic algorithm" << endl
		<< " [S]ummarize current status" << endl
		<< " [L]ist current reference taxa (i.e. have data for all partitions)" << endl
		<< " [W]rite current matrix" << endl
		<< " [R]evert to original matrix" << endl
		<< " [Q]uit" << endl
		<< endl
		<< "Enter desired option: ";
		
		cin >> userChoice;
		cin.ignore(200, '\n');
		
		if (checkCharValue(userChoice,'a')) {
			addGenes = true;
			validChoice = true;
			continue;
		} else if (checkCharValue(userChoice,'m')) {
			merge = true;
			validChoice = true;
			continue;
		} else if (checkCharValue(userChoice,'e')) {
			exclude = true;
			validChoice = true;
			continue;
		} else if (checkCharValue(userChoice,'d')) {
			deleteGenes = true;
			validChoice = true;
			continue;
		} else if (checkCharValue(userChoice,'r')) {
			revert = true;
			validChoice = true;
			continue;
		} else if (checkCharValue(userChoice,'p')) {
			print = true;
			validChoice = true;
			continue;
		} else if (checkCharValue(userChoice,'c')) {
			while (!validChoice) {
				cout << endl
				<< "Calculate partial decisiveness:" << endl
				<< " [T]ree-wise" << endl
				<< " [B]ranch-wise" << endl
				<< " [I]ndividual partition" << endl;
				cin >> userChoice;
				cin.ignore(200, '\n');
				
				if (checkCharValue(userChoice,'t'))
				{
					partialTreewise = true;
					validChoice = true;
				} else if (checkCharValue(userChoice,'b')) {
					partialBranchwise = true;
					validChoice = true;
				} else if (checkCharValue(userChoice,'i')) {
					partialIndividualPartition = true;
					validChoice = true;
				} else {
					cout << endl << "Invalid input option (" << userChoice << "). Try again." << endl << endl;
				}
			}
			continue;
		} else if (checkCharValue(userChoice,'s')) {
			summarize = true;
			validChoice = true;
			continue;
		} else if (checkCharValue(userChoice,'t')) {
			testCompleteDeciveness = true;
			validChoice = true;
			continue;
		} else if (checkCharValue(userChoice,'l')) {
			printRefTaxa = true;
			validChoice = true;
			continue;
		} else if (checkCharValue(userChoice,'w')) {
			writeCurrentMatrix = true;
			validChoice = true;
			continue;
		} else if (checkCharValue(userChoice,'i')) {
			testUserTree = true;
			validChoice = true;
			continue;
		} else if (checkCharValue(userChoice,'q')) {
			quit = true;
			validChoice = true;
			continue;
		} else if (checkCharValue(userChoice,'o')) {
			useGA = true;
			validChoice = true;
			continue;
		} else {
			cout << "Invalid input option (" << userChoice << "). Try again." << endl << endl;
		}
	}
}

// use getopt here instead
void processCommandLineArguments(int argc, char *argv[], string & matrixFileName,
	string & nexusFileName, string & locusWeightFileName, string & taxonWeightFileName,
	string & treeFileName, int & burnin, int & thinning, int & numProcs)
{
	if (argc == 1) { // Assume that it will be Nexus
		nexusFileName = getFileName();
	} else {
		for (int i = 1; i < argc; i++) {
			string temp = argv[i];
			
			if (temp == "-h" || temp == "-help") {
				printHelp();
				exit(0);  
			} else if (temp == "-m") {
				i++;
				matrixFileName = argv[i];
				checkValidFile(matrixFileName);
				continue;
			} else if (temp == "-d") {
				i++;
				nexusFileName = argv[i];
				checkValidFile(nexusFileName);
				continue;
			} else if (temp == "-l") {
				i++;
				locusWeightFileName = argv[i];
				checkValidFile(locusWeightFileName);
				continue;
			} else if (temp == "-w") {
				i++;
				taxonWeightFileName = argv[i];
				checkValidFile(taxonWeightFileName);
				continue;
			} else if (temp == "-t") {
				i++;
				treeFileName = argv[i];
				checkValidFile(treeFileName);
				continue;
			} else if (temp == "-b") {
				i++;
				burnin = convertStringtoInt(argv[i]);
				continue;
			} else if (temp == "-n") {
				i++;
				thinning = convertStringtoInt(argv[i]);
				continue;
			} else if (temp == "-np") {
				i++;
				numProcs = convertStringtoInt(argv[i]);
				continue;
			} else if (temp == "-debug") {
				debugging = true;
				cout << endl << "Printing additional information to screen for debugging purposes." << endl << endl;
				continue;
			} else {
				cout
				<< "Unknown command-line argument '" << argv[i] << "' encountered." << endl
				<< endl;
				printHelp();
				exit(0);
			}
			cout << endl;
		}
	}
}

void printHelp ()
{
	cout << endl
	<< "Program description: Calculates phylogenetic 'decisiveness' sensu Sanderson and Steel." << endl
	<< endl
	<< "To compile, type the following in a unix prompt:" << endl
	<< endl
	<< "make" << endl
	<< endl
	<< "Usage:" << endl
	<< endl
	<< "./Decisivator [-d data_file] [-m taxon-gene_matrix] [-t tree_file] [-b burnin] [-n thinning]" << endl
	<< "   [-w taxon_weights] [-l locus_weights] [-np num_proc]" << endl
	<< endl
	<< "where:" << endl
	<< endl
	<< "'data_file' is a simple 'vanilla' Nexus file containing sequences and defined CHARSETs." << endl
	<< "  - PLEASE NOTE! Only simple CHARSETs are currently supported." << endl
	<< "    - e.g. contiguous (X-Y) or interval (e.g. codon: X-Y\\3) data are fine." << endl
	<< "    - CHARSET referencing is NOT allowed at present (but will be!)" << endl
	<< endl
	<< "'taxon-gene_matrix' is a table listing taxa (rows) and genes (columns). This is a legacy format." << endl
	<< "  '1' indicates a cell (taxon-partition) has data, while '0' indicates it does not." << endl
	<< "  First row should give locus names. First column should give taxon names." << endl
	<< endl
	<< "'tree_file' contains user tree(s) in Nexus format to evaluate decisiveness upon. Tree(s) must be" << endl
	<< "fully bifurcating." << endl
	<< endl
	<< "'burnin' is the number of trees to ignore. Only makes sense with a distribution of trees." << endl
	<< endl
	<< "'thinning' is the interval between sampling trees (i.e. where every nth tree sample will be retained)." << endl
	<< "  Only makes sense with a distribution of trees." << endl
	<< endl
	<< "'taxon_weights' is a two-column (taxon, weight; with headers) file listing weights for taxa" << endl
	<< "  based on some arbitrary accessibility criterion." << endl
	<< endl
	<< "'locus_weights' is the analogous two-column (locus, weight; with headers) file for locus weights." << endl
	<< endl
	<< "'num_proc' is the number of processors to use (default is all available). Only with OPENMP version." << endl
	<< endl
	<< "NOTE: The taxon and locus weight files need not be complete. All weights are 1.0 by default." << endl
	<< "  Enter only weights which should be changed from the default. Or do this shit within the program itself." << endl
	<< endl;
}

void printSummaryInformation (vector <string> const& locusNames, vector <string> const& taxonNames,
	vector < vector <int> > const& data, double const& taxonCoverage, vector <int> const& referenceTaxa,
	bool const& matrixDecisive, double const& treewiseDecisiveness, double const& branchwiseDecisiveness,
	bool const& completeDecisivenessDetermined, string const& nexusFileName, int const& numRandomTrees,
	int const& numUserTrees, int const& numProcs)
{
	cout << endl << endl
	<< "*******************************" << endl
	<< "*** SUMMARY OF CURRENT DATA ***" << endl
	<< "*******************************" << endl;
	
	if (numProcs == 1) {
		cout << endl << "1 processor available for analysis." << endl;
	} else {
		cout << endl << numProcs << " processors available for analysis." << endl;
	}
	
	cout << "Input file: '" << nexusFileName << "'." << endl;
	cout << "A total of " << locusNames.size() << " partitions read for " << taxonNames.size() << " taxa." << endl;
	
	if (!referenceTaxa.empty()) {
		if (referenceTaxa.size() == 1) {
			cout << "Taxon '" << taxonNames[referenceTaxa[0]] << "' serves as a reference taxon (i.e. has data for all partitions)." << endl;
		} else {
			cout << referenceTaxa.size() << " reference taxa found (i.e. have data for all partitions)." << endl;
		}
	} else {
		cout << endl << "No reference taxa found (i.e. have data for all partitions)." << endl;
	}
	if (debugging) {
		cout << "A total of " << locusNames.size() << " partitions read:" << endl << endl;
		printVectorAsList(locusNames);
		cout << endl << "A total of " << taxonNames.size() << " taxa read:" << endl << endl;
		printVectorAsList(taxonNames);
		cout << endl;
	}
	
	if (numUserTrees != 0) {
		if (numUserTrees == 1) {
			cout << numUserTrees << " user tree in memory." << endl;
		} else {
			cout << numUserTrees << " user trees in memory." << endl;
		}
	} else {
		cout << "No user trees are in memory." << endl;
	}
	cout << "Matrix coverage is currently at: " << taxonCoverage << endl;
	if (!completeDecisivenessDetermined) {
		cout << "Complete decisiveness has not yet been determined for this taxon-character matrix." << endl;
	} else {
		if (matrixDecisive) {
			cout << "Taxon-character matrix is decisive for all possible trees!" << endl;
		} else {
			cout << "Taxon-character matrix is NOT currently decisive for all possible trees." << endl;
		}
	}
	if (treewiseDecisiveness != 0) {
		cout << "Partial tree-wise decisiveness for taxon-character matrix is currently: " << treewiseDecisiveness
		<< " (determined from " << numRandomTrees << " random trees)." << endl;
	} else {
		cout << "Partial tree-wise decisiveness has not yet been determined for this taxon-character matrix." << endl;
	}
	if (branchwiseDecisiveness != 0) {
		cout << "Partial branch-wise decisiveness for taxon-character matrix is currently: " << branchwiseDecisiveness
		<< " (determined from " << numRandomTrees << " random trees)." << endl;
	} else {
		cout << "Partial branch-wise decisiveness has not yet been determined for this taxon-character matrix." << endl;
	}
	checkForMissingTaxa (data, taxonNames);
}

void printReferenceTaxa (vector <int> const& referenceTaxa, vector <string> const& taxonNames)
{
	if (referenceTaxa.size() != 0) {
		cout << endl << referenceTaxa.size() << " reference taxa found (i.e. have data for all partitions):" << endl;
		for (int i = 0; i < (int)referenceTaxa.size(); i++) {
			cout << " " << i + 1 << ". " << taxonNames[referenceTaxa[i]] << endl;
		}
	} else {
		cout << endl << "No reference taxa found (i.e. have data for all partitions)." << endl;
	}
	cout << endl;
}


// print our proportion of data present
void printMatrix (vector < vector <int> > const& data, vector <string> const& taxonNames, 
	vector <double> const& locusWeights, vector <double> const& taxonWeights)
{
	cout.precision(3);
	cout.setf(ios::fixed,ios::floatfield);
	
	int numTaxa = (int)data.size();
	int numLoci = (int)data[0].size();
	vector <double> proportionTaxonDataPresent;
	vector <double> proportionPartitionDataPresent;
	
	string maxString;
	int longestName = 0;
	
// get proportion of data present for taxa
	for (int i = 0; i < numTaxa; i++) {
		double temp = 0.0;
		for (int j = 0; j < numLoci; j++) {
			if (data[i][j]) {
				temp ++;
			}
		}
		temp = temp / (double)numLoci;
		proportionTaxonDataPresent.push_back(temp);
	}
// now for loci
	for (int i = 0; i < numLoci; i++) {
		double temp = 0.0;
		for (int j = 0; j < numTaxa; j++) {
			if (data[j][i]) {
				temp ++;
			}
		}
		temp = temp / (double)numTaxa;
		proportionPartitionDataPresent.push_back(temp);
	}
	
	for (int i = 0; i < numTaxa; i++) {
		string currentString = taxonNames[i];
		if (currentString.size() > maxString.size()) {
			maxString = currentString;
		}
	}
// Determine length of longest name. um, seems like I can just use size here...
	for (string::const_iterator iterCharacters = maxString.begin(); iterCharacters < maxString.end(); ++iterCharacters) {
		longestName++;
	}
	cout << endl << "Current matrix:" << endl;
	for (int i = 0; i < (longestName + numLoci +2); i++) {
		cout << " ";
	}
	cout << "   Prop.";
	cout << "   Weight" << endl;
	for (int i = 0; i < numTaxa; i++) {
		cout << " " ;
// Print out leading spaces
		string tempName = taxonNames[i];
		if (tempName.size() < maxString.size()) {
			string::size_type tempDiffSize;
			tempDiffSize = maxString.size() - tempName.size();
			for (string::size_type iterSpaces = 0; iterSpaces < tempDiffSize; iterSpaces++) {
				cout << " ";
			}
		}
		cout << tempName << " ";
		for (int j = 0; j < numLoci; j++) {
			cout << data[i][j];
		}
		cout << "   " << proportionTaxonDataPresent[i];
		cout << "   " << taxonWeights[i] << endl;
	}
	cout << endl;
	
	cout << "Partition coverage and weights:" << endl << endl;
	printVectorAsList (proportionPartitionDataPresent, locusWeights, "Part.", "Prop.", "Weight");
	cout << endl;
}




// USE THIS!!!
void printMatrixToFile (vector < vector <int> > const& data, vector <string> const& taxonNames, 
	vector <double> const& locusWeights, vector <double> const& taxonWeights)
{
	ofstream log;
	log.open("Decisivator.log", ios::app);
	
	int numTaxa = (int)data.size();
	int numLoci = (int)data[0].size();
	
	string maxString;
	int longestName = 0;
	
	for (int i = 0; i < numTaxa; i++) {
		string currentString = taxonNames[i];
		if (currentString.size() > maxString.size()) {
			maxString = currentString;
		}
	}
// Determine length of longest name
	for (string::const_iterator iterCharacters = maxString.begin(); iterCharacters < maxString.end(); ++iterCharacters) {
		longestName++;
	}
	log << endl << "Currently active matrix:" << endl;
	for (int i = 0; i < (longestName + numLoci +3); i++) {
		log << " ";
	}
	log << "Weight" << endl;
	for (int i = 0; i < numTaxa; i++) {
		log << " " ;
// Print out leading spaces
		string tempName = taxonNames[i];
		if (tempName.size() < maxString.size()) {
			string::size_type tempDiffSize;
			tempDiffSize = maxString.size() - tempName.size();
			for (string::size_type iterSpaces = 0; iterSpaces < tempDiffSize; iterSpaces++) {
				log << " ";
			}
		}
		log << tempName << " ";
		for (int j = 0; j < numLoci; j++) {
			log << data[i][j];
		}
		log << "   " << taxonWeights[i] << endl;
	}
	log << endl;
	for (int i = 0; i < longestName - 6; i++) {
		log << " ";
	}
	log << "Weight: ";
	for (int i = 0; i < numLoci; i++) {
		log << locusWeights[i];
	}
	log << endl;
}





void writeMatrix (vector <string> const& taxonNames, int const& numChar, 
	vector < vector <string> > const& taxaAlignment, vector < vector <int> > const& includedLocusRanges,
	vector <string> const& locusNames)
{
	int numTaxa = (int)taxonNames.size();
	bool done = false;
	string fileName;
	
	while (!done) {
		char userChoice;
		cout << endl << endl << "WRITE CURRENT MATRIX" << endl << endl
		<< "Write:" << endl
		<< " [N]exus format" << endl
		<< " [P]hylip format" << endl
		<< " [B]oth" << endl
		<< " [R]eturn to main menu" << endl
		<< endl
		<< "Enter desired option: ";
		
		cin >> userChoice;
		cin.ignore(200, '\n');
		
		if (checkCharValue(userChoice,'n') || checkCharValue(userChoice,'b')) {
			cout << "Enter name for Nexus-formatted data: ";
			cin >> fileName;
			cin.ignore(200, '\n');
			checkValidOutputFile(fileName);
			writeNexus(numTaxa, numChar, taxonNames, fileName, taxaAlignment, includedLocusRanges, locusNames);
			done = true;
		}
		if (checkCharValue(userChoice,'p') || checkCharValue(userChoice,'b')) {
			cout << "Enter name for PHYLIP-formatted data: ";
			cin >> fileName;
			cin.ignore(200, '\n');
			checkValidOutputFile(fileName);
			writePhylip(numTaxa, numChar, taxonNames, fileName, taxaAlignment);
			done = true;
		} else if (checkCharValue(userChoice,'r')) {
			done = true;
		} else if (!checkCharValue(userChoice,'p') && !checkCharValue(userChoice,'b') && !checkCharValue(userChoice,'n') && !checkCharValue(userChoice,'r')) {
			cout << "Invalid input option (" << userChoice << "). Try again." << endl;
		}
	}
}

void writeNexus (int const& numTaxa, int const& numChar, vector <string> const& taxonNames,
	string & fileName, vector < vector <string> > const& taxaAlignment, vector < vector <int> > const& includedLocusRanges,
	vector <string> const& locusNames)
{
	ofstream outFile;
	outFile.open(fileName.c_str());
	int numLoci = (int)locusNames.size();
	string maxString;
	// int longestName = 0;
	
	for (int i = 0; i < numTaxa; i++) {
		string currentString = taxonNames[i];
		if (currentString.size() > maxString.size()) {
			maxString = currentString;
		}
	}
	// longestName = (int)maxString.size();
	
	outFile << "#NEXUS" << endl << endl
	<< "[ *** File generated by Decisivator v0.3 ***]" << endl
	<< "[ *** via PuRGe University of Idaho 2011 ***]"<< endl << endl
	<< "Begin data;" << endl
	<< "	Dimensions ntax=" << numTaxa << " nchar=" << numChar << ";" << endl
	<< "	Format datatype=dna missing=? interleave=yes;" << endl
	<< "	Matrix" << endl << endl;
	
	for (int i = 0; i < numTaxa; i++) {
		bool match = false;
		int j = 0;
		while (!match) {
 			if (taxaAlignment[j][0] == taxonNames[i]) {
 				outFile << taxonNames[i];
// Print out trailing spaces
				string tempName = taxonNames[i];
				if (tempName.size() < maxString.size()) {
					string::size_type tempDiffSize;
					tempDiffSize = maxString.size() - tempName.size();
					for (string::size_type iterSpaces = 0; iterSpaces < tempDiffSize; iterSpaces++) {
						outFile << " ";
					}
				}
				outFile << "   " << taxaAlignment[j][1] << endl;
 				match = true;
 			}
			j++;
		}
	}
	outFile << ";" << endl
	<< "End;" << endl;
	
// write CHARSET information
	outFile << endl << "BEGIN ASSUMPTIONS;" << endl << endl;
	for (int partIter = 0; partIter < numLoci; partIter++) {
		if (includedLocusRanges[partIter][0] == 0) { // simple range e.g. '1-566'
			outFile << "CHARSET " << locusNames[partIter] << " = " << includedLocusRanges[partIter][1]
				<< "-" << includedLocusRanges[partIter][2] << ";" << endl;
		} else if (includedLocusRanges[partIter][0] == 1) { // interval range e.g. '1-566\3'
			outFile << "CHARSET " << locusNames[partIter] << " = " << includedLocusRanges[partIter][1]
				<< "-" << includedLocusRanges[partIter][2] << "\\" << includedLocusRanges[partIter][2]
				<< ";" << endl;
		}
	}
	outFile << endl << "End;";
	
	outFile.close();
}

void writePhylip (int const& numTaxa, int const& numChar, vector <string> const& taxonNames,
	string & fileName, vector < vector <string> > const& taxaAlignment)
{
	ofstream outFile;
	outFile.open(fileName.c_str());
	string maxString;
	// int longestName = 0;
	
	for (int i = 0; i < numTaxa; i++) {
		string currentString = taxonNames[i];
		if (currentString.size() > maxString.size()) {
			maxString = currentString;
		}
	}
	// longestName = (int)maxString.size();
	
	outFile << numTaxa << " " << numChar << endl;
	
	for (int i = 0; i < numTaxa; i++) {
		bool match = false;
		int j = 0;
		while (!match) {
 			if (taxaAlignment[j][0] == taxonNames[i]) {
 				outFile << taxonNames[i];
 				string tempName = taxonNames[i];
				if (tempName.size() < maxString.size()) {
					string::size_type tempDiffSize;
					tempDiffSize = maxString.size() - tempName.size();
					for (string::size_type iterSpaces = 0; iterSpaces < tempDiffSize; iterSpaces++) {
						outFile << " ";
					}
				}
 				outFile << "   " << taxaAlignment[j][1] << endl;
 				match = true;
 			}
			j++;
		}
	}
	outFile.close();
	cout  << endl << "File '" << fileName << "' successfully created." << endl;
}

int selectPartition (vector < vector <int> > const& data, vector <string> const& locusNames)
{
	int partitionID;
	int numLoci = data[0].size();
	bool done = false;
	while (!done) {
		cout << "Partitions available:" << endl << endl;
		printVectorAsList(locusNames);
		partitionID = checkValidIntInput("Select partition to query: ");
		if (partitionID < 1 || partitionID > numLoci) {
			cout << "Invalid choice. Must be between 1 and " << numLoci << ". Try again." << endl;
		} else {
			partitionID--;
			done = true;
		}
	}
	return(partitionID);
}
