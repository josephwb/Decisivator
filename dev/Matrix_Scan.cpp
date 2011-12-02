#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <numeric>

using namespace std;

#include "Matrix_Scan.h"
#include "Trees_Edges.h"
#include "General.h"
#include "User_Interface.h"

extern bool DEBUG;

bool searchForReferenceTaxon (vector < vector <int> > & data, vector <int> & referenceTaxa,
	vector <string> & taxonNames)
{
	bool referencePresent = false;
	referenceTaxa.clear();
	int numTaxa = data.size();
	int numLoci = data[0].size();
	int numReference = 0;
	
	cout << endl << "SEARCHING FOR REFERENCE TAXA..." << endl << endl;
	
	vector <int> completeTaxon (numLoci,1);
	for (int i = 0; i < numTaxa; i++)
	{
		if (data[i] == completeTaxon)
		{
			if (!referencePresent)
			{
				if (DEBUG) {cout << endl << "Reference taxon (i.e. sequenced for all genes)." << endl;}
			}
			referencePresent = true;
			referenceTaxa.push_back(i);
			if (DEBUG) {cout << " " << numReference + 1 << ". " << taxonNames[i] << "." << endl;}
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
			
			if (DEBUG) {cout << "Moving reference taxon '" << taxonNames[ref] << "' to bottom of matrix. Don't be alarmed!" << endl;}
			
// need to fix reference taxon indexing; moving taxon i to end changes index for all taxa i+1 -> N
			if (DEBUG) {cout << "Raw:" << endl;
			printVectorAsList(referenceTaxa);}
			
			referenceTaxa.push_back(numTaxa);
			referenceTaxa.erase(referenceTaxa.begin()+0);
			
			for (int i = 0; i < (int)referenceTaxa.size(); i++)
			{
				referenceTaxa[i] = referenceTaxa[i] - 1;
			}
			
			if (DEBUG) {cout << "Fixed:" << endl;
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
	
	int numTaxa = data.size();
	int numLoci = data[1].size();
	BigInt tripletCounter = 0;
	BigInt numPresent = 0;
	BigInt numPossibleTriplets = choose(numTaxa, 3);
	
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
				printProgress("Triplet", tripletCounter, numPossibleTriplets);
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
	
	int numTaxa = data.size();
	int numLoci = data[1].size();
	BigInt quartetCounter = 0;
	BigInt numPresent = 0;
	BigInt numPossibleQuartets = choose(numTaxa, 4);
	
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
					printProgress("Quartet", quartetCounter, numPossibleQuartets);
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

// Just uses custom int type that is larger than standard
void searchForQuartetsWithReferenceBIG (vector < vector <int> > const& data, vector <int> & missingQuartetsByTaxa,
	vector <int> const& referenceTaxa, vector <string> const& taxonNames)
{
	int numTaxa = data.size();
	int numLoci = data[1].size();
	int numReferenceTaxa = referenceTaxa.size();
	vector <int> temp (numTaxa, 0);
	missingQuartetsByTaxa = temp;
	
	BigInt quartetCounter = 0;
	BigInt numPresent = 0;
	int numTaxaToConsider = 0;
	BigInt numPossibleQuartets = choose(numTaxa - numReferenceTaxa, 3);
	
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
	numTaxaToConsider = availableTaxa.size();
	
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
	
	if (DEBUG)
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
				printProgress("Quartet", quartetCounter, numPossibleQuartets);
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

void searchForQuartetsWithReference (vector < vector <int> > const& data, vector < vector <int> > & missingQuartets,
	vector <int> const& referenceTaxa, vector <string> const& taxonNames)
{
// Clear vector, as dimensions may have changed
	missingQuartets.clear();
	
	int numTaxa = data.size();
	int numLoci = data[1].size();
	int numReferenceTaxa = referenceTaxa.size();
	BigInt quartetCounter = 0;
	BigInt numPresent = 0;
	int numTaxaToConsider = 0;
	BigInt numPossibleQuartets = choose(numTaxa - numReferenceTaxa, 3);
	
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
	numTaxaToConsider = availableTaxa.size();
	
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
	
	if (DEBUG)
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
				printProgress("Quartet", quartetCounter, numPossibleQuartets);
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
	int numTaxa = taxonNames.size();
	int numGroups = missingGroups.size();
	int numReferenceTaxa = referenceTaxa.size();
	
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

void whichTaxaProblematicBIG (vector <int> const& missingQuartetsByTaxa, vector <string> const& taxonNames,
	string const& grouping, vector <int> const& referenceTaxa)
{
	vector <int> missingTaxa;
	int numTaxa = taxonNames.size();
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

void getCoverage (vector < vector <int> > & data, double & taxonCoverage)
{
	int taxonGeneCount = 0;
	taxonCoverage = 0;
	taxonGeneCount = 0;
	int numTaxa = data.size();
	int numLoci = data[1].size();
	
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

int searchEdgePartitionsAll (vector < vector <int> > const& data, vector <int> const& left,
	vector <int> const& right, vector <int> const& sib, vector <int> const& upper,
	bool const& referenceTaxonPresent)
{
	int numLoci = data[0].size();
	int numSatisfied = 0;
	int numUpper = upper.size();
	
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
						if(data[left[i]][m] == 1 && data[right[j]][m] == 1 && data[sib[k]][m] == 1 && data[upper[l]][m] == 1)
						{
							numSatisfied++;
						}
					}
				}
			}
		}
	}
	return numSatisfied;
}

int searchEdgePartitionsMinimum (vector < vector <int> > const& data, vector <int> const& left,
	vector <int> const& right, vector <int> const& sib, vector <int> const& upper,
	bool const& referenceTaxonPresent, int & partitionMatched)
{
	int numLoci = data[0].size();
	int numSatisfied = 0;
	int numUpper = upper.size();
	
	if (referenceTaxonPresent) // only consider a single reference taxon for upper; if it fails, all will
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
						if(data[left[i]][m] == 1 && data[right[j]][m] == 1 && data[sib[k]][m] == 1 && data[upper[l]][m] == 1)
						{
							numSatisfied++;
							partitionMatched = m; // may want to know which partition is involved
							if (DEBUG) {cout << "Quartet (" << left[i] << "," << right[j] << "," << sib[k] << "," << upper[l] << ") found at locus " << m << "!" << endl;}
							return(numSatisfied); // exit: one is all that is needed
						}
					}
				}
			}
		}
	}
	return numSatisfied;
}

// searching clade for a cache of a given partition
bool searchEdgeOnePartition (vector < vector <int> > const& data, vector <int> const& clade,
	int const& partitionToMatch)
{
	bool matched = true;
	for (int i = 0; i < (int)clade.size(); i++)
	{
		if (data[clade[i]][partitionToMatch] != 1) // exit, not a match
		{
			return false;
		}
	}
	return matched;
}

string getPartitionMatches (vector < vector <int> > const& data, vector <int> const& left,
	vector <int> const& right, vector <int> const& sib, vector <int> const& upper,
	int const& partitionToMatch, bool & leftValid, bool & rightValid, bool & sibValid,
	bool & upperValid)
{
	string newMatches;
	
// 3 valid
	if (leftValid && rightValid && sibValid && !upperValid) // L, R, S
	{
		upperValid = searchEdgeOnePartition (data, upper, partitionToMatch);
		if (upperValid)
		{
			if (newMatches.empty())
			{
				cout << "Here I am, ashole!" << endl;
				char foo = 'U';
				string temp (1, foo);
				newMatches = temp;
				cout << "Hmm. That seemed to have worked..." << endl;
			}
			else
			{
				newMatches += 'U';
			}
		}
	}
	else if (leftValid && rightValid && !sibValid && upperValid) // L, R, U
	{
		sibValid = searchEdgeOnePartition (data, sib, partitionToMatch);
		if (sibValid)
		{
			if (newMatches.empty())
			{
				char foo = 'S';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'S';
			}
		}
	}
	else if (leftValid && !rightValid && sibValid && upperValid) // L, S, U
	{
		rightValid = searchEdgeOnePartition (data, right, partitionToMatch);
		if (rightValid)
		{
			if (newMatches.empty())
			{
				char foo = 'R';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'R';
			}
		}
	}
	else if (!leftValid && rightValid && sibValid && upperValid) // R, S, U
	{
		leftValid = searchEdgeOnePartition (data, left, partitionToMatch);
		if (leftValid)
		{
			if (newMatches.empty())
			{
				char foo = 'L';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'L';
			}
		}
	}

// 2 valid
	else if (leftValid && rightValid && !sibValid && !upperValid) // L, R
	{
		upperValid = searchEdgeOnePartition (data, upper, partitionToMatch);
		if (upperValid)
		{
			if (newMatches.empty())
			{
				char foo = 'U';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'U';
			}
		}
		sibValid = searchEdgeOnePartition (data, sib, partitionToMatch);
		if (sibValid)
		{
			if (newMatches.empty())
			{
				char foo = 'S';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'S';
			}
		}
	}
	else if (leftValid && !rightValid && sibValid && !upperValid) // L, S
	{
		upperValid = searchEdgeOnePartition (data, upper, partitionToMatch);
		if (upperValid)
		{
			if (newMatches.empty())
			{
				char foo = 'U';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'U';
			}
		}
		rightValid = searchEdgeOnePartition (data, right, partitionToMatch);
		if (rightValid)
		{
			if (newMatches.empty())
			{
				char foo = 'R';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'R';
			}
		}
	}
	else if (leftValid && !rightValid && !sibValid && upperValid) // L, U
	{
		rightValid = searchEdgeOnePartition (data, right, partitionToMatch);
		if (rightValid)
		{
			if (newMatches.empty())
			{
				char foo = 'R';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'R';
			}
		}
		sibValid = searchEdgeOnePartition (data, sib, partitionToMatch);
		if (sibValid)
		{
			if (newMatches.empty())
			{
				char foo = 'S';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'S';
			}
		}
	}
	else if (!leftValid && rightValid && sibValid && !upperValid) // R, S
	{
		upperValid = searchEdgeOnePartition (data, upper, partitionToMatch);
		if (upperValid)
		{
			if (newMatches.empty())
			{
				char foo = 'U';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'U';
			}
		}
		leftValid = searchEdgeOnePartition (data, left, partitionToMatch);
		if (leftValid)
		{
			if (newMatches.empty())
			{
				char foo = 'L';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'L';
			}
		}
	}
	else if (!leftValid && rightValid && !sibValid && upperValid) // R, U
	{
		leftValid = searchEdgeOnePartition (data, left, partitionToMatch);
		if (leftValid)
		{
			if (newMatches.empty())
			{
				char foo = 'L';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'L';
			}
		}
		sibValid = searchEdgeOnePartition (data, sib, partitionToMatch);
		if (sibValid)
		{
			if (newMatches.empty())
			{
				char foo = 'S';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'S';
			}
		}
	}
	else if (!leftValid && !rightValid && sibValid && upperValid) // S, U
	{
		leftValid = searchEdgeOnePartition (data, left, partitionToMatch);
		if (leftValid)
		{
			if (newMatches.empty())
			{
				char foo = 'L';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'L';
			}
		}
		rightValid = searchEdgeOnePartition (data, right, partitionToMatch);
		if (rightValid)
		{
			if (newMatches.empty())
			{
				char foo = 'R';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'R';
			}
		}
	}

// 1 valid
	else if (leftValid && !rightValid && !sibValid && !upperValid) // L
	{
		upperValid = searchEdgeOnePartition (data, upper, partitionToMatch);
		if (upperValid)
		{
			if (newMatches.empty())
			{
				char foo = 'U';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'U';
			}
		}
		rightValid = searchEdgeOnePartition (data, right, partitionToMatch);
		if (rightValid)
		{
			if (newMatches.empty())
			{
				char foo = 'R';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'R';
			}
		}
		sibValid = searchEdgeOnePartition (data, sib, partitionToMatch);
		if (sibValid)
		{
			if (newMatches.empty())
			{
				char foo = 'S';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'S';
			}
		}
	}
	else if (!leftValid && rightValid && !sibValid && !upperValid) // R
	{
		upperValid = searchEdgeOnePartition (data, upper, partitionToMatch);
		if (upperValid)
		{
			if (newMatches.empty())
			{
				char foo = 'U';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'U';
			}
		}
		leftValid = searchEdgeOnePartition (data, left, partitionToMatch);
		if (leftValid)
		{
			if (newMatches.empty())
			{
				char foo = 'L';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'L';
			}
		}
		sibValid = searchEdgeOnePartition (data, sib, partitionToMatch);
		if (sibValid)
		{
			if (newMatches.empty())
			{
				char foo = 'S';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'S';
			}
		}
	}
	else if (!leftValid && !rightValid && sibValid && !upperValid) // S
	{
		upperValid = searchEdgeOnePartition (data, upper, partitionToMatch);
		if (upperValid)
		{
			if (newMatches.empty())
			{
				char foo = 'U';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'U';
			}
		}
		leftValid = searchEdgeOnePartition (data, left, partitionToMatch);
		if (leftValid)
		{
			if (newMatches.empty())
			{
				char foo = 'L';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'L';
			}
		}
		else
		{
			newMatches += 'L';
		}
		rightValid = searchEdgeOnePartition (data, right, partitionToMatch);
		if (rightValid)
		{
			if (newMatches.empty())
			{
				char foo = 'R';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'R';
			}
		}
	}
	else if (!leftValid && !rightValid && !sibValid && upperValid) // U
	{
		leftValid = searchEdgeOnePartition (data, left, partitionToMatch);
		if (leftValid)
		{
			if (newMatches.empty())
			{
				char foo = 'L';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'L';
			}
		}
		rightValid = searchEdgeOnePartition (data, right, partitionToMatch);
		if (rightValid)
		{
			if (newMatches.empty())
			{
				char foo = 'R';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'R';
			}
		}
		sibValid = searchEdgeOnePartition (data, sib, partitionToMatch);
		if (sibValid)
		{
			if (newMatches.empty())
			{
				char foo = 'S';
				string temp (1, foo);
				newMatches = temp;
			}
			else
			{
				newMatches += 'S';
			}
		}
	}
	
	else
	{
		cout << "WTF?!?" << endl;
	}
	return newMatches;
}

// search stored satisfied clades; controlling function loops over partitions
int searchEdgeStoredClades (vector <int> const& newClade,
	vector < vector <int> > const& storedClades)
{
	int numClades = storedClades.size();
	
//	cout << endl << "Searching for presence of all possible taxon quartets..." << endl;
	for (int i = 0; i < numClades; i++)
	{
		if (checkCladeSubset(newClade, storedClades[i]))
		{
			return true;
		}
	}
	return false;
}

bool checkCladeSubset (vector <int> const& newClade, vector <int> const& storedClade)
{
	bool valid = false;
	
	vector <int> difference;
	
	valid = subtractVector(newClade, storedClade, difference);
	
	if (!valid) // if stored clade has taxa not in new clade, abandon
	{
		return false;
	}
	else
	{
		if ((accumulate(newClade.begin(),newClade.end(),0) -
			accumulate(storedClade.begin(),storedClade.end(),0))
			== accumulate(difference.begin(),difference.end(),0))
		{
			return true;
		}
	}
	return false;
}


// Eventually pass in more stuff e.g. keep track of missing clades, etc.
double calculatePartialDecisivenessOLD (bool const& referenceTaxonPresent, int & numTrees,
	vector < vector <int> > const& data, bool const& findAll)
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
//	int numTaxa = data.size();
	
	
// should this be fixed or flexible?
	if (numTrees == 0)
	{
		numTrees = checkValidIntInput("Specify number of trees to simulate for calculating partial decisiveness: ");
	}
	
	cout << endl << "Simulating " << numTrees << " trees and testing for satisfaction of internal edges..." << endl;
	
// Running mean: ((((count - 1)/count) * runningMean) + (currentResult/count));
	
	for (int j = 0; j < numTrees; j++)
	{
		printProgress("Tree", j + 1, numTrees);
		treeCount++;
		double currentDecisiveness = 0.0;
		double numEdgesSatisfied = 0;
		tree = fastBinaryTree(data.size(), sibNodes, referenceTaxonPresent);
		if (DEBUG) {printTree(tree);}
		
//		printTree(tree);
//		sort(tree.begin()+numTaxa,tree.end());
//		printTree(tree);
		
//		for (int i = numTaxa; i < (int)tree.size(); i++)
//		{
//			cout << "Clade " << i << " sums to: " << accumulate(tree[i].begin(),tree[i].end(),0) << endl;;
//		}
		
		if (DEBUG) {cout << endl << "Edges:" << endl;}
		for (int i = 0; i < numInternalEdges; ++i) // Walk through all internal edges
		{
			getEdges(i, tree, sibNodes, referenceTaxonPresent, left, right, sib, upper);
			
// *** Scan taxon-gene matrix here ***
			if (findAll)
			{
				if (searchEdgePartitionsAll(data, left, right, sib, upper, referenceTaxonPresent))
				{
					numEdgesSatisfied++;
				}
			}
			else
			{
				int partitionMatched = 0;
				if (searchEdgePartitionsMinimum(data, left, right, sib, upper, referenceTaxonPresent, partitionMatched))
				{
					numEdgesSatisfied++;
				}
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
	
	cout << "Partial decisiveness is currently: " << partialDecisiveness << endl;
	
	return partialDecisiveness;
}


// this version caches locus-clades that have been previously satisfied; for use in future clade lookups
double calculatePartialDecisiveness (bool const& referenceTaxonPresent, int & numTrees,
	vector < vector <int> > const& data, bool const& findAll)
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
	int numPartitions = data[0].size();
	
// should this be fixed or flexible?
	if (numTrees == 0)
	{
		numTrees = checkValidIntInput("Specify number of trees to simulate for calculating partial decisiveness: ");
	}
	
	cout << endl << "Simulating " << numTrees << " trees and testing for satisfaction of internal edges..." << endl;
	
// Running mean: ((((count - 1)/count) * runningMean) + (currentResult/count));
	
	for (int j = 0; j < numTrees; j++)
	{
		printProgress("Tree", j + 1, numTrees);
		treeCount++;
		vector < vector <int> > temp;
		vector <int> foo;
		vector < vector < vector <int> > > cacheSatisfiedClades (numPartitions, temp);
		
		double currentDecisiveness = 0.0;
		double numEdgesSatisfied = 0;
		tree = fastBinaryTree(data.size(), sibNodes, referenceTaxonPresent);
		if (DEBUG) {printTree(tree);}
		
		if (DEBUG) {cout << endl << "Edges:" << endl;}
		for (int i = 0; i < numInternalEdges; ++i) // Walk through all internal edges
		{
			getEdges(i, tree, sibNodes, referenceTaxonPresent, left, right, sib, upper);
			
			bool rightValid = false;
			bool leftValid = false;
			bool sibValid = false;
			bool upperValid = false;
			bool match = false;
			string newMatches;
			
// search through cached edges first; could be much faster
// no use pushing back consistent clades here = redundant
			for (int k = 0; k < numPartitions; k++)
			{
				vector < vector <int> > partitionClades = cacheSatisfiedClades[k];
				if (partitionClades.empty())
				{
					cout << "partitionClades empty for edge #" << i << ", partition #" << k << endl;
				}
				
				
				if (!partitionClades.empty())
				{
					cout << "partitionClades not empty for edge #" << i << ", partition #" << k << endl;
					
					leftValid = searchEdgeStoredClades (left, partitionClades);
					rightValid = searchEdgeStoredClades (right, partitionClades);
					sibValid = searchEdgeStoredClades (sib, partitionClades);
					upperValid = searchEdgeStoredClades (upper, partitionClades);
					
					
					if (leftValid && rightValid && sibValid && upperValid)
					{
						match = true;
						k = numPartitions;
						cout << "Satisfied all four clades in one fell swoop for edge #" << i << endl;
						continue;
					}
					
					cout << "leftValid = " << leftValid << "; rightValid = " << rightValid <<
						"; sibValid = " << sibValid << "; upperValid = " << upperValid << endl;
						
					if (!leftValid && !rightValid && !sibValid && !upperValid)
					{
						if (DEBUG) {cout << "No valid clades. Gotta search'em all." << endl;}
						continue;
					}
					else
					{
						newMatches = getPartitionMatches(data, left, right, sib,
							upper, k, leftValid, rightValid, sibValid, upperValid);
						if (!newMatches.empty())
						{	
							if (DEBUG) {cout << "Caching new clades." << endl;}
							// loop over size of newMatches
							// push into cache
							for (string::iterator iter = newMatches.begin(); iter != newMatches.end(); iter++)
							{
		// *** Dying here ***
		// initialization must be wrong...
								if (DEBUG) {cout << "Working with character '" << *iter << "'." << endl;}
								if (*iter == 'L')
								{
									cacheSatisfiedClades[k].push_back(left);
								}
								else if (*iter == 'R')
								{
									cacheSatisfiedClades[k].push_back(right);
								}
								else if (*iter == 'S')
								{
									cacheSatisfiedClades[k].push_back(sib);
								}
								else if (*iter == 'U')
								{
									cacheSatisfiedClades[k].push_back(upper);
								}
							}
						}
					}
					if (leftValid && rightValid && sibValid && upperValid)
					{
						match = true;
						k = numPartitions;
						continue;
					}
				
				
				
				
				
				
				}
			}
// *** Scan taxon-gene matrix here ***
			if (match)
			{
				numEdgesSatisfied++;
			}
			else
			{
				int partitionMatched = 0;
				if (searchEdgePartitionsMinimum(data, left, right, sib, upper,
					referenceTaxonPresent, partitionMatched))
				{
					numEdgesSatisfied++;
					
					
					
					
					
					
					
					
					
					
					
// push back satisfied clades HERE
					cacheSatisfiedClades[partitionMatched].push_back(left);
					cacheSatisfiedClades[partitionMatched].push_back(right);
					cacheSatisfiedClades[partitionMatched].push_back(sib);
					cacheSatisfiedClades[partitionMatched].push_back(upper);
					
					
					
					
					
					
					
					
					
					
					
				}
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
		cacheSatisfiedClades.clear();
		currentDecisiveness = numEdgesSatisfied/numInternalEdges;
		
// running mean to avoid potential overflow
		partialDecisiveness = ((((treeCount - 1)/treeCount) * partialDecisiveness) + (currentDecisiveness/treeCount));
	}
	cout << endl << "Done." << endl << endl;
	
	cout << "Partial decisiveness is currently: " << partialDecisiveness << endl;
	
	return partialDecisiveness;
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
	
	return matrixDecisive;
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
	int numTaxa = userTrees[0][0].size();
	int numTrees = userTrees.size();
	
	printMatrixToFile (data, taxonNames, locusWeights, taxonWeights);
	
	double numInternalEdges = userTrees[0][0].size() - 3; // *should* be the same for all trees in a file...
	
	if (DEBUG) {cout << "Here I am! numTrees = " << numTrees << ". numInternalEdges = " << numInternalEdges
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
		if (DEBUG)
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
		
		if (DEBUG)
		{
//			cout << "Formatted tree:" << endl;
//			printTree(formattedTree);
		}
		
		sibNodes = getSibNodes(formattedTree); // get sibling node relationships; expected downstream
		
		if (DEBUG) {printTree(formattedTree);}
		if (DEBUG) {cout << endl << "Edges:" << endl;}
		for (int i = 0; i < numInternalEdges; ++i) // this is the important bit; *** NEED TO INCORPORATE TAXON CODING ***
		{
			currentClade = formattedTree[numTaxa + i];
			double currentDecisiveness = 0.0;
			int numEdgesSatisfied = 0;
			
			getEdges(i, formattedTree, sibNodes, 0, left, right, sib, upper);
			int numPossibleQuartets = (int)left.size() * (int)right.size() * (int)sib.size() * (int)upper.size();
			
// *** Scan taxon-gene matrix here ***
			 numEdgesSatisfied = searchEdgePartitionsAll(data, left, right, sib, upper, 0);
			 
			if (DEBUG)
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
		printBipartitionTable(formattedTree, decisivenessCurrentTree, numSatisfied, numPossible, numTaxa, numTrees, j);
		result.push_back(decisivenessCurrentTree);
		
		if (numTrees > 1)
		{
			cout << "Average decisiveness for tree " << j + 1 << " is: " << sum(decisivenessCurrentTree)/decisivenessCurrentTree.size() << endl;
		}
		else
		{
			cout << "Average decisiveness for this tree is: " << sum(decisivenessCurrentTree)/decisivenessCurrentTree.size() << endl;
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
	
	return result;
}

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