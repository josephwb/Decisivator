#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <numeric>

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
#endif

using namespace std;

#include "General.h"
#include "Matrix_Scan.h"
#include "Trees_Edges.h"
#include "User_Interface.h"

extern bool debuggering;

void checkForMissingTaxa (vector < vector <int> > const& data, vector <string> const& taxonNames)
{
	vector <string> missingData;
	for (int i = 0; i < (int)data.size(); i++) // numTaxa
	{
		bool match = false;
		for (int j = 0; j < (int)data[0].size(); j++) // numPartitions
		{
			if (data[i][j])
			{
				match = true;
			}
		}
		if (!match)
		{
			missingData.push_back(taxonNames[i]);
		}
	}
	if (!missingData.empty())
	{
		if (missingData.size() > 1)
		{
			cout << endl << " *** WARNING! MULTIPLE TAXA WITH NO DATA!! ***" << endl << endl;
			printVectorAsList(missingData);
		}
		else
		{
			cout << endl << "*** WARNING! TAXON WITH NO DATA!! ***" << endl << endl;
			cout << "Taxon: " << missingData[0] << "." << endl << endl;
		}
	}
}

bool searchForReferenceTaxon (vector < vector <int> > & data, vector <int> & referenceTaxa,
	vector <string> & taxonNames)
{
	bool referencePresent = false;
	referenceTaxa.clear();
	int numTaxa = (int)data.size();
	int numLoci = (int)data[0].size();
	int numReference = 0;
	
	cout << endl << "SEARCHING FOR REFERENCE TAXA..." << endl << endl;
	
	vector <int> completeTaxon (numLoci,1);
	for (int i = 0; i < numTaxa; i++)
	{
		if (data[i] == completeTaxon)
		{
			if (!referencePresent)
			{
				if (debuggering) {cout << endl << "Reference taxon (i.e. sequenced for all genes)." << endl;}
			}
			referencePresent = true;
			referenceTaxa.push_back(i);
			if (debuggering) {cout << " " << numReference + 1 << ". " << taxonNames[i] << "." << endl;}
			cout << " " << numReference + 1 << ". " << taxonNames[i] << "." << endl;
			numReference++;
		}
	}
	if (numReference > 0)
	{
		cout << numReference << " reference taxa observed." << endl;
		
// take first reference taxon found, put at the bottom; speeds up satisfying edges
		if (referenceTaxa[numReference - 1] != (int)taxonNames.size() - 1) // i.e. last taxon is NOT a reference taxon
		{
			int ref = referenceTaxa[0];
			
			if (debuggering) {cout << "Moving reference taxon '" << taxonNames[ref] << "' to bottom of matrix. Don't be alarmed!" << endl;}
			
// need to fix reference taxon indexing; moving taxon i to end changes index for all taxa i+1 -> N
			if (debuggering) {cout << "Raw:" << endl;
			printVectorAsList(referenceTaxa);}
			
			referenceTaxa.push_back(numTaxa);
			referenceTaxa.erase(referenceTaxa.begin()+0);
			
			for (int i = 0; i < (int)referenceTaxa.size(); i++)
			{
				referenceTaxa[i] = referenceTaxa[i] - 1;
			}
			
			if (debuggering) {cout << "Fixed:" << endl;
			printVectorAsList(referenceTaxa);}
			
			data.push_back(data[ref]);
			data.erase(data.begin()+ref);
			
			taxonNames.push_back(taxonNames[ref]);
			taxonNames.erase(taxonNames.begin()+ref);
		}
	}
	else
	{
		cout << "No reference taxa observed." << endl;
	}
	return referencePresent;
}

// We DO want locations for triplets when we are looping over missing quartets
void searchForAllTriplets (vector < vector <int> > const& data, vector < vector <int> > & triplets,
	vector < vector <int> > & tripletLocations, vector < vector <int> > & missingTriplets)
{
// Clear vectors, as dimensions may have changed
	missingTriplets.clear();
	tripletLocations.clear();
	triplets.clear();
	
	int numTaxa = (int)data.size();
	int numLoci = (int)data[1].size();
	unsigned long tripletCounter = 0;
	unsigned long numPresent = 0;
	unsigned long numPossibleTriplets = choose(numTaxa, 3);
	
	cout << endl << "Searching for presence of all possible taxon triplets..." << endl;
	for (int i = 0; i < numTaxa - 2; i++)
	{
		vector <int> tempI;
		tempI.push_back(i);
		for (int j = i + 1; j < numTaxa - 1; j++)
		{
			vector <int> tempJ = tempI;
			tempJ.push_back(j);
			for (int k = j + 1; k < numTaxa; k++)
			{
				vector <int> tempK = tempJ;
				tempK.push_back(k);
				vector <int> locations;
				
				bool tripletFound = false;
				for (int l = 0; l < numLoci; l++)	// Loop over loci
				{
					if(data[i][l] == 1 && data[j][l] == 1 && data[k][l] == 1)
					{
						tripletFound = true;
						locations.push_back(l);
					}
				}
				triplets.push_back(tempK);
				tripletLocations.push_back(locations);
				if (!tripletFound)
				{
					missingTriplets.push_back(tempK);
				}
				else
				{
					numPresent++;
				}
				tempK.clear();
				locations.clear();
				tripletCounter++;
                if ((tripletCounter & CHUNKSIZE) == CHUNKSIZE) {
                    printProgress("Triplet", tripletCounter, numPossibleTriplets);
                }
			}
		}
	}
	cout << endl << "Counted " << tripletCounter << " total triplets, " << numPresent << " of which were observed." << endl;
	if (tripletCounter != numPresent)
	{
		cout << "Matrix is NOT decisive for all possible trees." << endl;
	}
	else
	{
		cout << " Woo-hoo! All possible taxon triplets observed. Matrix is (probably) decisive for all possible trees!" << endl;
	}
}

// May not be so useful; if not all triplets are present, neither will all quaduplets
// Therefore, only call if all triplets ARE present (a necessary but not sufficient condition)
void searchForAllQuartets (vector < vector <int> > const& data, vector < vector <int> > & missingQuartets)
{
// Clear vector, as dimensions may have changed
	missingQuartets.clear();
	
	int numTaxa = (int)data.size();
	int numLoci = (int)data[1].size();
	unsigned long quartetCounter = 0;
	unsigned long numPresent = 0;
	unsigned long numPossibleQuartets = choose(numTaxa, 4);
	
	cout << endl << "Searching for presence of all possible taxon quartets..." << endl;
	for (int i = 0; i < numTaxa - 3; i++)
	{
		vector <int> tempI;
		tempI.push_back(i);
		for (int j = i + 1; j < numTaxa - 2; j++)
		{
			vector <int> tempJ = tempI;
			tempJ.push_back(j);
			for (int k = j + 1; k < numTaxa - 1; k++)
			{
				vector <int> tempK = tempJ;
				tempK.push_back(k);
				for (int l = k + 1; l < numTaxa; l++)
				{
					vector <int> tempL = tempK;
					tempL.push_back(l);
					
					bool quartetFound = false;
					for (int m = 0; m < numLoci; m++)	// Loop over loci
					{
						if(data[i][m] == 1 && data[j][m] == 1 && data[k][m] == 1 && data[l][m] == 1)
						{
							quartetFound = true;
							continue;
						}
					}
					if (!quartetFound)
					{
						missingQuartets.push_back(tempL);
					}
					else
					{
						numPresent++;
					}
					tempL.clear();
					quartetCounter++;
                    if ((quartetCounter & CHUNKSIZE) == CHUNKSIZE) {
                        printProgress("Quartet", quartetCounter, numPossibleQuartets);
                    }
				}
			}
		}
	}
	
	cout << endl << "Counted " << quartetCounter << " total quartets, " << numPresent << " of which were observed." << endl;
	
	if (quartetCounter != numPresent)
	{
		cout << "Matrix is NOT decisive for all possible trees." << endl;
	}
	else
	{
		cout << " Woo-hoo! All possible taxon quartets observed. Matrix IS decisive for all possible trees!" << endl;
	}
}


// *** make this shit multi-threaded ***
void searchForQuartetsWithReference (vector < vector <int> > const& data, vector < vector <int> > & missingQuartets,
	vector <int> const& referenceTaxa, vector <string> const& taxonNames)
{
// Clear vector, as dimensions may have changed
	missingQuartets.clear();
	
	int numTaxa = (int)data.size();
	int numLoci = (int)data[1].size();
	int numReferenceTaxa = (int)referenceTaxa.size();
	unsigned long quartetCounter = 0;
	unsigned long numPresent = 0;
	int numTaxaToConsider = 0;
	unsigned long numPossibleQuartets = choose(numTaxa - numReferenceTaxa, 3);
	
	vector <int> availableTaxa; // Don't waste time considering reference taxa == always present
	for (int i = 0; i < numTaxa; i++)
	{
		bool match = false;
		for (int j = 0; j < numReferenceTaxa; j++)
		{
			if (i == referenceTaxa[j])
			{
				match = true;
			}
		}
		if (!match)
		{
			availableTaxa.push_back(i);		
		}
	}
	numTaxaToConsider = (int)availableTaxa.size();
	
	if (availableTaxa.size() == 0)
	{
		cout << "Matrix is complete, and so decisive for all possible trees. Lucky you." << endl;
	}
	else
	{
		if (numReferenceTaxa > 0)
		{
			cout << endl << "Searching for presence of all possible taxon quartets containing (arbitrary) reference taxon '" << taxonNames[referenceTaxa[0]] << "'..." << endl;	
		}
		else
		{
			cout << endl << "Searching for presence of all possible taxon quartets containing reference taxon '" << taxonNames[referenceTaxa[0]] << "'..." << endl;	
		}
	}
// Essentially just looking for triples with N = numTaxa-1 (a little better if several reference taxa exist), since reference taxon is a given
	vector <int> reference;
	reference.push_back(referenceTaxa[0]);
	
	if (debuggering)
	{
		cout << "numTaxaToConsider = " << numTaxaToConsider << ":" << endl;
		for (int i = 0; i < numTaxaToConsider; i++)
		{
			cout << "  " << i + 1 << ". " << taxonNames[availableTaxa[i]] << endl;
		}
	}
	
// If there are <= 3 taxa missing data, the search-loop is not entered, and things fuck up
// Add in arbitrary reference taxa
	if (numTaxaToConsider < 3)
	{
		int diff = 3 - numTaxaToConsider;
		for (int i = 0; i < diff; i++)
		{
			numTaxaToConsider++;
			availableTaxa.push_back(referenceTaxa[i+1]);
		}
	}
	
	for (int i = 0; i < numTaxaToConsider - 2; i++)
	{
		vector <int> tempI = reference;
		tempI.push_back(availableTaxa[i]);
		for (int j = i + 1; j < numTaxaToConsider - 1; j++)
		{
			vector <int> tempJ = tempI;
			tempJ.push_back(j);
			for (int k = j + 1; k < numTaxaToConsider; k++)
			{
				vector <int> tempK = tempJ;
				tempK.push_back(k);
				vector <int> locations;
				
				bool quartetFoundFound = false;
				for (int l = 0; l < numLoci; l++)	// Loop over loci
				{
//					if(data[i][l] == 1 && data[j][l] == 1 && data[k][l] == 1)
					if(data[availableTaxa[i]][l] == 1 && data[availableTaxa[j]][l] == 1 && data[availableTaxa[k]][l] == 1)
					{
						quartetFoundFound = true;
						continue; // got it; move on
					}
				}
				if (!quartetFoundFound)
				{
					missingQuartets.push_back(tempK);
				}
				else
				{
					numPresent++;
				}
				tempK.clear();
				quartetCounter++;
                if ((quartetCounter & CHUNKSIZE) == CHUNKSIZE) {
                    printProgress("Quartet", quartetCounter, numPossibleQuartets);
                }
			}
		}
	}
	
	if (quartetCounter == numPresent && quartetCounter == 0) // should be obsolete
	{
		cout << endl << "Zero quartets considered. Happens when < 4 'available' taxa present." << endl;
	}
	
	if (quartetCounter != numPresent)
	{
		cout << "Matrix is NOT decisive for all possible trees." << endl;
	}
	else if (quartetCounter == numPresent && numPresent > 0)
	{
		cout << " Woo-hoo! All possible taxon quartets observed. Matrix IS decisive for all possible trees!" << endl;
	}
}

// Summarize which taxa are missing from triplets/quartets
void whichTaxaProblematic (vector < vector <int> > const& missingGroups, vector <string> const& taxonNames,
	string const& grouping, vector <int> const& referenceTaxa)
{
	vector <int> missingTaxa;
	int numTaxa = (int)taxonNames.size();
	int numGroups = (int)missingGroups.size();
	int numReferenceTaxa = (int)referenceTaxa.size();
	
	for (int i = 0; i < numTaxa; i++)
	{
		int numMissing = 0;
		bool reference = false;
		for (int j = 0; j < numReferenceTaxa; j++)	// Don't count reference taxa in missing taxa
		{
			if (i == referenceTaxa[j])
			{
				reference = true;
			}
		}
		if (!reference)
		{
			for (int j = 0; j < numGroups; j++)
			{
				for (int k = 0; k < int(missingGroups[j].size()); k++)
				{
					if (i == missingGroups[j][k])
					{
						numMissing++;
					}
				}
			}
		}
		missingTaxa.push_back(numMissing);
	}
	
	cout << endl << "Missing " << grouping << " summary:" << endl << endl;
	for (int i = 0; i < numTaxa; i++)
	{
		if (missingTaxa[i] > 0)
		{
			cout << " Taxon '" << taxonNames[i] << "' missing from (" << missingTaxa[i] << ") taxon " << grouping << "." << endl;
		}
	}
}

// No need to parallelize here
void getCoverage (vector < vector <int> > const& data, double & taxonCoverage)
{
	int taxonGeneCount = 0;
	taxonCoverage = 0;
	taxonGeneCount = 0;
	int numTaxa = (int)data.size();
	int numLoci = (int)data[1].size();
	
	for (int i = 0; i < numTaxa; i++)
	{
		for (int j = 0; j < numLoci; j++)
		{
			if(data[i][j] == 1)
			{
				taxonGeneCount++;
			}
		}
	}
	taxonCoverage = float(taxonGeneCount)/(float(numTaxa) * float(numLoci));
}




// Eventually pass in more stuff e.g. keep track of missing clades, etc.
// multithreading implemented via openmp
double calculatePartialDecisiveness (bool const& referenceTaxonPresent, int & numTrees,
	vector < vector <int> > const& data, bool const& findAll, int const& numProcs)
{
	double partialDecisiveness = 0.0; // keep track with running mean to avoid possible overflow
	double numInternalEdges = data.size() - 3;
	int numTaxa = data.size();
	
	double foo = 0.0;
	
// should this be fixed or flexible?
	if (numTrees == 0)
	{
		numTrees = checkValidIntInput("Specify number of trees to simulate for calculating partial decisiveness: ");
	}
	
	cout << endl << "Simulating " << numTrees << " trees and testing for satisfaction of internal edges..." << endl;
	
// Running mean: ((((count - 1)/count) * runningMean) + (currentResult/count)); // prevent overflow
	
	#pragma omp parallel
	{
		vector < vector <bool> > tree;
		vector < vector <int> > sibNodes;
		vector <int> left;
		vector <int> right;
		vector <int> sib;
		vector <int> upper;
		double privateTreeCount = 0.0; // needed for multithreading
		//double privateDecisiveness = 0.0;
		
		int thread_id = omp_get_thread_num();
		
//		double temp = 0.0;
		double sum = 0.0;
		
		#pragma omp for
		for (int j = 0; j < numTrees; j++)
		{
			if (numProcs == 1) {printProgress("Tree", j + 1, numTrees);}
			privateTreeCount++;
			double currentDecisiveness = 0.0;
			double numEdgesSatisfied = 0.0;
			tree = fastBinaryTree(numTaxa, sibNodes, referenceTaxonPresent);
			if (debuggering) {printTree(tree);}
			
// For debugging:
//		printTree(tree);
//		sort(tree.begin()+numTaxa,tree.end());
					
			if (debuggering) {cout << endl << "Edges:" << endl;}
			
			for (int i = 0; i < numInternalEdges; ++i) // Walk through all internal edges
			{
				getEdges(i, tree, sibNodes, referenceTaxonPresent, left, right, sib, upper);
				
// *** Scan taxon-gene matrix here ***
				if (searchEdgePartitions(data, left, right, sib, upper, findAll, referenceTaxonPresent))
				{
					numEdgesSatisfied++;
				}
// empty vectors for next edge
				left.clear();
				right.clear();
				sib.clear();
				upper.clear();
			}
// empty vectors for next tree
			tree.clear();
			sibNodes.clear();
			currentDecisiveness = numEdgesSatisfied/numInternalEdges;
//			privateDecisiveness = ((((privateTreeCount - 1)/privateTreeCount) * privateDecisiveness)
//				+ (currentDecisiveness/privateTreeCount));
			
//			temp = privateDecisiveness;
			sum += currentDecisiveness;
		}
		
		#pragma omp critical
		{
//			cout << "Thread " << thread_id << " analyzed (" << privateTreeCount << ") trees. Yay!"
//				<< " Partial decisiveness for this thread is: " << temp << endl;
			cout << "Partial decisiveness for thread " << thread_id << " is: " << sum / privateTreeCount << endl;
//			partialDecisiveness += privateDecisiveness;
		//	foo += sum / privateTreeCount;
			foo += sum;
		}
	}
	cout << endl << "Done." << endl << endl;
	
//	partialDecisiveness = partialDecisiveness / numProcs;
	
//	cout << "Counted " << treeCount << " trees." << endl;
	
//	cout << "Partial decisiveness is currently: " << partialDecisiveness << endl;
//	cout << "Foo decisiveness is currently: " << foo / numProcs << endl;
	cout << "Foo decisiveness is currently: " << foo / numTrees << endl;
	
	return(partialDecisiveness);
}

double calculatePartialDecisivenessSinglePartition (bool const& referenceTaxonPresent, int & numTrees,
	vector < vector <int> > const& data, bool const& findAll, int const& partitionID,
	vector <string> const& locusNames)
{
	double partialDecisiveness = 0.0; // keep track with running mean to avoid possible overflow
	vector < vector <bool> > tree;
	vector < vector <int> > sibNodes;
	vector <int> left;
	vector <int> right;
	vector <int> sib;
	vector <int> upper;
	double treeCount = 0;
	double numInternalEdges = data.size() - 3;
	int numTaxa = data.size();
	vector < vector <int> > partitionToQuery;
	
	for (int i = 0; i < numTaxa; i++) // extract partition of interest
	{
		vector <int> temp;
		temp.push_back(data[i][partitionID]);
		partitionToQuery.push_back(temp);
	}
	
// should this be fixed or flexible?
	if (numTrees == 0)
	{
		numTrees = checkValidIntInput("Specify number of trees to simulate for calculating partial decisiveness: ");
	}
	
	cout << endl << "Simulating " << numTrees << " trees and testing for satisfaction of internal edges..." << endl;
	
	for (int j = 0; j < numTrees; j++)
	{
		printProgress("Tree", j + 1, numTrees);
		treeCount++;
		double currentDecisiveness = 0.0;
		double numEdgesSatisfied = 0;
		tree = fastBinaryTree(numTaxa, sibNodes, referenceTaxonPresent);
		if (debuggering) {printTree(tree);}
		
		if (debuggering) {cout << endl << "Edges:" << endl;}
		for (int i = 0; i < numInternalEdges; ++i) // Walk through all internal edges
		{
			getEdges(i, tree, sibNodes, referenceTaxonPresent, left, right, sib, upper);
			
// *** Scan taxon-gene matrix here ***
			if (searchEdgePartitions(partitionToQuery, left, right, sib, upper, findAll, referenceTaxonPresent))
			{
				numEdgesSatisfied++;
			}
			
// empty vectors for next edge
			left.clear();
			right.clear();
			sib.clear();
			upper.clear();
		}
// empty vectors for next tree
		tree.clear();
		sibNodes.clear();
		currentDecisiveness = numEdgesSatisfied/numInternalEdges;
		partialDecisiveness = ((((treeCount - 1)/treeCount) * partialDecisiveness) + (currentDecisiveness/treeCount));
	}
	cout << endl << "Done." << endl << endl;
	
	cout << "Partial decisiveness for partition '" << locusNames[partitionID] << "' is: " << partialDecisiveness << endl;
	
	return(partialDecisiveness);
}

// findAll below switches between boolean and count
int searchEdgePartitions (vector < vector <int> > const& data, vector <int> const& left,
	vector <int> const& right, vector <int> const& sib, vector <int> const& upper, bool const& findAll,
	bool const& referenceTaxonPresent)
{
	int numLoci = (int)data[0].size();
	int numSatisfied = 0;
	int numUpper = (int)upper.size();
	
	if (!findAll && referenceTaxonPresent) // only consider a single reference taxon for upper; if it fails, all will
	{
		numUpper = 1;
	}
	
// Need to loop over vectors left, right, sib, other
//	cout << endl << "Searching for presence of all possible taxon quartets..." << endl;
	for (int l = 0; l < numUpper; l++)
	{
		for (int i = 0; i < int(left.size()); i++)
		{
			for (int j = 0; j < int(right.size()); j++)
			{
				for (int k = 0; k < int(sib.size()); k++)
				{
					for (int m = 0; m < numLoci; m++)	// Loop over loci
					{
						if (data[left[i]][m] == 1 && data[right[j]][m] == 1 && data[sib[k]][m] == 1 && data[upper[l]][m] == 1)
						{
							if (!findAll)
							{
								numSatisfied++;
								if (debuggering) {cout << "Quartet (" << left[i] << "," << right[j] << "," << sib[k] << "," << upper[l] << ") found at locus " << m << "!" << endl;}
								return(numSatisfied); // exit: one is all that is needed
							}
							else
							{
								numSatisfied++;
								m = numLoci; // exit this quartet, go on to next
							}
						}
					}
				}
			}
		}
	}
	return(numSatisfied);
}

bool testCompleteDecisivness (vector < vector <int> > const& data, bool const& referenceTaxonPresent,
	vector <int> const& referenceTaxa, vector <string> const& taxonNames, vector < vector <int> > & missingQuartets,
	vector < vector <int> > & triplets, vector < vector <int> > & tripletLocations, vector < vector <int> > & missingTriplets)
{
	bool matrixDecisive = false;
	
	if (referenceTaxonPresent)
	{
		searchForQuartetsWithReference(data, missingQuartets, referenceTaxa, taxonNames);
		if (missingQuartets.size() > 0)
		{
			whichTaxaProblematic(missingQuartets, taxonNames, "quartets", referenceTaxa);
		}
		else
		{
			matrixDecisive = true;
		}
	}
	else
	{
		searchForAllTriplets(data, triplets, tripletLocations, missingTriplets); // test triplets first; far fewer
		if (missingTriplets.size() > 0)
		{
			whichTaxaProblematic(missingTriplets, taxonNames, "triplets", referenceTaxa);
		}
		else
		{
			searchForAllQuartets(data, missingQuartets);
			if (missingQuartets.size() > 0)
			{
				whichTaxaProblematic(missingQuartets, taxonNames, "quartets", referenceTaxa);
			}
		}
	}
	
	return (matrixDecisive);
}

vector < vector <double> > determineDecisivenessUserTree (vector < vector <int> > const& data,
	vector < vector < vector <bool> > > & userTrees, vector < vector <int> > const& treeTaxonOrdering,
	vector <string> const& taxonNames, vector <double> const& locusWeights, vector <double> const& taxonWeights)
{
	vector < vector <double> > result;
	vector < vector <bool> > rawTree;
	vector < vector <bool> > formattedTree;
	vector < vector <int> > sibNodes;
	vector <bool> currentClade;
	vector <int> left;
	vector <int> right;
	vector <int> sib;
	vector <int> upper;
	double treeCount = 0;
	int numTaxa = (int)userTrees[0][0].size();
	int numTrees = (int)userTrees.size();
	bool findAll = true;
	
	printMatrixToFile (data, taxonNames, locusWeights, taxonWeights);
	
	double numInternalEdges = userTrees[0][0].size() - 3; // *should* be the same for all trees in a file...
	
	if (debuggering) {cout << "Here I am! numTrees = " << numTrees << ". numInternalEdges = " << numInternalEdges
		<< ". numTaxa = " << numTaxa << "." << endl;}
	
	for (int j = 0; j < numTrees; j++)
	{
		treeCount++;
		vector <double> decisivenessCurrentTree;
		vector <int> taxonOrdering;
		vector <int> numSatisfied;
		vector <int> numPossible;
		
/* ** Need to map taxon ordering to alignment via treeTaxonOrdering ***
options:
1. rearrange tree to reflect that ordering (easy but ugly)
2. map (pretty but harder); this will be fixed when put into OO format
*/
		rawTree = userTrees[j];
		taxonOrdering = treeTaxonOrdering[j];
		if (debuggering)
		{
//			cout << "Raw tree:" << endl;
//			printTree(rawTree);
			
			cout << endl << "Translation:" << endl;
			printVectorAsList(taxonOrdering);
		}
		for (int k = 0; k < (int)rawTree.size(); k++)
		{
			vector <bool> clade (numTaxa, false);
			for (int l = 0; l < numTaxa; l++)
			{
				if (rawTree[k][l])
				{
					clade[taxonOrdering[l]] = true;
				}
			}
			formattedTree.push_back(clade);
			clade.clear();
		}
		
//			cout << "Formatted tree:" << endl;
//			printTree(formattedTree);
		
		sibNodes = getSibNodes(formattedTree); // get sibling node relationships; expected downstream
		
		if (debuggering) {cout << endl << "Edges:" << endl;}
		for (int i = 0; i < numInternalEdges; ++i) // this is the important bit; *** NEED TO INCORPORATE TAXON CODING ***
		{
			currentClade = formattedTree[numTaxa + i];
			double currentDecisiveness = 0.0;
			int numEdgesSatisfied = 0;
			
			getEdges(i, formattedTree, sibNodes, 0, left, right, sib, upper);
			int numPossibleQuartets = (int)left.size() * (int)right.size() * (int)sib.size() * (int)upper.size();
			
// *** Scan taxon-gene matrix here ***
			 numEdgesSatisfied = searchEdgePartitions(data, left, right, sib, upper, findAll, 0);
			 
			if (debuggering)
			{
				cout << "Considering " << numPossibleQuartets << " ways to satisfy current node " << i + numTaxa << ":" << endl;
				printClade(currentClade);
				cout << endl;
				cout << numEdgesSatisfied << " instances satisfied of a possible " << numPossibleQuartets
			 	<< " for node " << i + numTaxa << " (= " << ((double)numEdgesSatisfied/(double)numPossibleQuartets)*100 << "%)." << endl;
			}
			currentDecisiveness = (double)numEdgesSatisfied / (double)numPossibleQuartets;
			decisivenessCurrentTree.push_back(currentDecisiveness);
			numSatisfied.push_back(numEdgesSatisfied);
			numPossible.push_back(numPossibleQuartets);
			 
// empty vectors for next edge
			left.clear();
			right.clear();
			sib.clear();
			upper.clear();
		}
		
		printBipartitionTable(formattedTree, decisivenessCurrentTree, numSatisfied,
			numPossible, numTaxa, numTrees, j);
		result.push_back(decisivenessCurrentTree);
		
		if (numTrees > 1)
		{
			cout << "Average decisiveness for tree " << j + 1 << " is: " << sum(decisivenessCurrentTree)/(double)decisivenessCurrentTree.size() << endl;
		}
		else
		{
			cout << "Average decisiveness for this tree is: " << sum(decisivenessCurrentTree)/(double)decisivenessCurrentTree.size() << endl;
		}
// empty vectors for next tree
		rawTree.clear();
		formattedTree.clear();
		sibNodes.clear();
		decisivenessCurrentTree.clear();
		numSatisfied.clear();
		numPossible.clear();
	//	partialDecisiveness = ((((treeCount - 1)/treeCount) * partialDecisiveness) + (currentDecisiveness/treeCount));
	}
	cout << "Done." << endl << endl;
	
	return(result);
}

// PAUP-style
void printBipartitionTable (vector < vector <bool> > const& tree, vector <double> const& decisiveness,
	vector <int> const& numSatisfied, vector <int> const& numPossible, int const& numTaxa,
	int const& numTrees, int const& treeNumber)
{
	ofstream log;
	log.open("Decisivator.log",ios::app);
	
	if (numTrees > 1)
	{
		cout << endl << "Bipartition Decisiveness Scores for tree #" << treeNumber + 1 << ":" << endl << endl;
		log << endl << "Bipartition Decisiveness Scores for tree #" << treeNumber + 1 << ":" << endl << endl;
	}
	else
	{
		cout << endl << endl << "Bipartition Decisiveness Scores:" << endl << endl;
		log << endl << endl << "Bipartition Decisiveness Scores:" << endl << endl;
	}
	
	int counter = 0;
	int divisor = 0;
	if (numTaxa > 999)
	{
		counter = 1;
		divisor = 1000;
		for (int i = 1; i <= numTaxa; i++)
		{
			if (i % divisor == 0)
			{
				cout << counter;
				log << counter;
				counter ++;
			}
			else
			{
				cout << " ";
				log << " ";
			}
		}
		cout << endl;
		log << endl;
	}
	if (numTaxa > 99)
	{
		counter = 1;
		divisor = 100;
		for (int i = 1; i <= numTaxa; i++)
		{
			if (i % divisor == 0)
			{
				cout << counter;
				log << counter;
				counter ++;
			}
			else
			{
				cout << " ";
				log << " ";
			}
		}
		cout << endl;
		log << endl;
	}
	if (numTaxa > 9)
	{
		counter = 1;
		divisor = 10;
		for (int i = 1; i <= numTaxa; i++)
		{
			if (i % divisor == 0)
			{
				cout << counter;
				log << counter;
				counter ++;
			}
			else
			{
				cout << " ";
				log << " ";
			}
		}
		cout << endl;
		log << endl;
	}
	counter = 1;
	for (int i = 1; i <= numTaxa; i++)
	{
		if (counter == 10)
		{
			counter = 0;
			cout << counter;
			counter ++;
		}
		else
		{
			cout << counter;
			log << counter;
			counter ++;
		}
	}
	cout << "	Freq.	Poss.	 %" << endl;
	log << "	Freq.	Poss.	 %" << endl;
	
	for (int i = 1; i <= numTaxa; i++)
	{
		cout << "-";
		log << "-";
	}
	cout << "---------------------------" << endl;
	log << "---------------------------" << endl;
	
// Last 'node' defines entire tree, penultimate node defines root (all other nodes already taken care of)
	counter = 0;
	for (int i = numTaxa; i < (int)tree.size() - 2; i++)
	{
		vector <bool> currentClade = tree[i];
		for (int j = 0; j < numTaxa; j++)
		{
			if (currentClade[j])
			{
				cout << "*";
				log << "*";
			}
			else
			{
				cout << ".";
				log << ".";
			}
		}
		cout << "	" << numSatisfied[counter] << "	" << numPossible[counter] << "	" << decisiveness[counter]*100 << endl;
		log << "	" << numSatisfied[counter] << "	" << numPossible[counter] << "	" << decisiveness[counter]*100 << endl;
		counter++;
	}
	cout << endl;
	log << endl;
}


// *** DEPRECATED CODE *** //

// Just uses custom int type that is larger than standard; deprecated
void searchForQuartetsWithReferenceBIG (vector < vector <int> > const& data, vector <int> & missingQuartetsByTaxa,
	vector <int> const& referenceTaxa, vector <string> const& taxonNames)
{
	int numTaxa = (int)data.size();
	int numLoci = (int)data[1].size();
	int numReferenceTaxa = (int)referenceTaxa.size();
	vector <int> temp (numTaxa, 0);
	missingQuartetsByTaxa = temp;
	
	unsigned long quartetCounter = 0;
	unsigned long numPresent = 0;
	int numTaxaToConsider = 0;
	unsigned long numPossibleQuartets = choose(numTaxa - numReferenceTaxa, 3);
	
	vector <int> availableTaxa; // Don't waste time considering reference taxa == always present
	for (int i = 0; i < numTaxa; i++)
	{
		bool match = false;
		for (int j = 0; j < numReferenceTaxa; j++)
		{
			if (i == referenceTaxa[j])
			{
				match = true;
			}
		}
		if (!match)
		{
			availableTaxa.push_back(i);		
		}
	}
	numTaxaToConsider = (int)availableTaxa.size();
	
	if (availableTaxa.size() == 0)
	{
		cout << "Matrix is complete, and so decisive for all possible trees. Lucky you." << endl;
	}
	else
	{
		if (numReferenceTaxa > 0)
		{
			cout << endl << "Searching for presence of all possible taxon quartets containing (arbitrary) reference taxon '" << taxonNames[referenceTaxa[0]] << "'..." << endl;	
		}
		else
		{
			cout << endl << "Searching for presence of all possible taxon quartets containing reference taxon '" << taxonNames[referenceTaxa[0]] << "'..." << endl;	
		}
	}
	
	if (debuggering)
	{
		cout << "numTaxaToConsider = " << numTaxaToConsider << ":" << endl;
		for (int i = 0; i < numTaxaToConsider; i++)
		{
			cout << "  " << i + 1 << ". " << taxonNames[availableTaxa[i]] << endl;
		}
	}
	
// If there are <= 3 taxa missing data, the search-loop is not entered, and things fuck up
// Add in arbitrary reference taxa
	if (numTaxaToConsider < 3)
	{
		int diff = 3 - numTaxaToConsider;
		for (int i = 0; i < diff; i++)
		{
			numTaxaToConsider++;
			availableTaxa.push_back(referenceTaxa[i+1]);
		}
	}
	
	for (int i = 0; i < numTaxaToConsider - 2; i++)
	{
		for (int j = i + 1; j < numTaxaToConsider - 1; j++)
		{
			for (int k = j + 1; k < numTaxaToConsider; k++)
			{
				bool quartetFoundFound = false;
				for (int l = 0; l < numLoci; l++)	// Loop over loci
				{
					if(data[availableTaxa[i]][l] == 1 && data[availableTaxa[j]][l] == 1 && data[availableTaxa[k]][l] == 1)
					{
						quartetFoundFound = true;
						continue; // got it; move on
					}
				}
				if (!quartetFoundFound)
				{
					missingQuartetsByTaxa[i]++;
					missingQuartetsByTaxa[j]++;
					missingQuartetsByTaxa[k]++;
				}
				else
				{
					numPresent++;
				}
				quartetCounter++;
                if ((quartetCounter & CHUNKSIZE) == CHUNKSIZE) {
                    printProgress("Quartet", quartetCounter, numPossibleQuartets);
                }
			}
		}
	}
	
	if (quartetCounter == numPresent && quartetCounter == 0) // should be obsolete
	{
		cout << endl << "Zero quartets considered. Happens when < 4 'available' taxa present." << endl;
	}
	if (quartetCounter != numPresent)
	{
		cout << "Matrix is NOT decisive for all possible trees." << endl;
	}
	else if (quartetCounter == numPresent && numPresent > 0)
	{
		cout << " Woo-hoo! All possible taxon quartets observed. Matrix IS decisive for all possible trees!" << endl;
	}
}

void whichTaxaProblematicBIG (vector <int> const& missingQuartetsByTaxa, vector <string> const& taxonNames,
	string const& grouping, vector <int> const& referenceTaxa)
{
	vector <int> missingTaxa;
	int numTaxa = (int)taxonNames.size();
	int currentMax = 0;
	
	for (int i = 0; i < numTaxa; i++)
	{
		if (missingQuartetsByTaxa[i] > currentMax)
		{
			missingTaxa.clear();
			currentMax = missingQuartetsByTaxa[i];
			missingTaxa.push_back(i);
		}
		else if (missingQuartetsByTaxa[i] == currentMax)
		{
			missingTaxa.push_back(i);
		}
	}
	
	if (currentMax > 0)
	{
		cout << endl << "Missing " << grouping << " summary:" << endl << endl;
		for (int i = 0; i < int(missingTaxa.size()); i++)
		{
			cout << " Taxon '" << taxonNames[missingTaxa[i]] << "' missing from (" << currentMax << ") taxon " << grouping << "." << endl;
		}
	}
}