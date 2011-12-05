#include <iostream>
#include <vector>
#include <sstream>
#include <numeric>
#include <cstdlib>
#include <algorithm>

using namespace std;

#include "Manipulate_Matrix.h"
#include "General.h"
#include "Matrix_Scan.h"

extern bool debuggering;

void addTaxonGeneToMatrix (vector < vector <int> > & data, vector <string> const& taxonNames,
	vector <string> & locusNames, vector <double> & locusWeights, vector <double> & taxonWeights)
{
	bool taxaDone = false;
	bool locusDone = false;
	vector <int> userInput;
	bool validIntEntry = false;
	int numEntries = 0;
	int taxonSelection;
	vector <string> availableTaxa;
	vector <string> availableLoci;
	vector <int> availableTaxonIndex; // need to map to data from reduced vector
	vector <int> availableLocusIndex; // need to map to data from reduced vector
	
	for (int i = 0; i < int(taxonNames.size()); i ++)
	{
		if (taxonWeights[i] != 0)
		{
			availableTaxa.push_back(taxonNames[i]);
			availableTaxonIndex.push_back(i);
		}
	}
	for (int i = 0; i < int(locusNames.size()); i ++)
	{
		if (locusWeights[i] != 0)
		{
			availableLoci.push_back(locusNames[i]);
			availableLocusIndex.push_back(i);
		}
	}
	
	while (!taxaDone)
	{
		cout << endl << "Taxa available (those with non-zero weight):" << endl;
		printVectorAsList(availableTaxa);
		
		cout << endl << "Enter indicies of taxa you would like to add virtual genes (separated by spaces), or 0 to exit: ";
		numEntries = 0;	
		string tempString;
		getline(cin,tempString);
		
// Ability to read in multiple inputs
		validIntEntry = false;
		istringstream tempStream(tempString);
		while (tempStream >> taxonSelection)
		{
			validIntEntry = true;
			userInput.push_back(taxonSelection);
			numEntries++;
		}
		if (taxonSelection == 0)
		{
			taxaDone = true;
			continue;
		}
		if (!validIntEntry)
		{
			cout << endl << "*** Invalid input. Integer must be between 0 and " << taxonNames.size() << " ***" << endl;
			cin.clear();
			userInput.clear();
			continue;
		}
		if (validIntEntry)
		{
			if (cin.fail() || taxonSelection < 0 || taxonSelection > int(availableTaxa.size()))
			{
				cout << endl << "*** Invalid input. Integer must be between 0 and " << availableTaxa.size() << " ***" << endl;
				cin.clear();
				userInput.clear();
				continue;
			}
			else
			{
				sort(userInput.begin(), userInput.end());
// Order from 'largest' to 'smallest' index; otherwise need keep track of changing indices
				reverse(userInput.begin(), userInput.end());
				
				while (!locusDone)
				{
					cout << endl << "Loci available (those with non-zero weight):" << endl;
					printVectorAsList(availableLoci);
					int newGene;
					cout << endl << "Enter index of desired gene, or 0 to exit: ";
					cin >> newGene;
					
					if (cin.fail() || newGene < 0 || newGene > int(availableLoci.size()))
					{
						cout << endl << "*** Invalid input. Integer must be between 0 and " << availableLoci.size() << " ***" << endl;
						cin.clear();
						continue;
					}
					else
					{
						if (newGene == 0) // user changed mind
						{
							locusDone = true;
							taxaDone = true;
							continue;
						}
						else
						{
							for (int editIter = 0; editIter < numEntries; editIter++)
							{
// Need to fix this, using availableTaxa and availableLoci, availableTaxonIndex and availableLocusIndex
								data[availableTaxonIndex[userInput[editIter] - 1]][availableLocusIndex[newGene - 1]] = 1;
							}
							locusDone = true;
						}
					}
				}
			}
		}
		userInput.clear();
		taxaDone = true;
	}
}

// Um, I guess not useful
void deleteGenesFromMatrix (vector < vector <int> > & data, vector <string> & locusNames, vector <double> & locusWeights)
{
	bool done = false;
	vector <int> userInput;
	bool validIntEntry = false;
	int numTaxa = (int)data.size();
	int numEntries = 0;
	int geneSelection;
	
	while (!done)
	{
		cout << endl << "Loci available:" << endl;
		printVectorAsList(locusNames);
		
		cout << endl << "Enter indicies of loci you would like to delete (if multiple, separated by spaces), or 0 to exit: ";
		numEntries = 0;	
		string tempString;
		getline(cin,tempString);
		
// Ability to read in multiple inputs
		validIntEntry = false;
		istringstream tempStream(tempString);
		while (tempStream >> geneSelection)
		{
			validIntEntry = true;
			userInput.push_back(geneSelection);
			numEntries++;
		}
		if (geneSelection == 0)
		{
			done = true;
			continue;
		}
		if (!validIntEntry)
		{
			cout << endl << "*** Invalid input. Integer must be between 0 and " << locusNames.size() << " ***" << endl;
			cin.clear();
			userInput.clear();
			done = false;
			continue;
		}
		if (validIntEntry)
		{
			if (cin.fail() || geneSelection < 0 || geneSelection > int(locusNames.size()))
			{
				cout << endl << "*** Invalid input. Integer must be between 0 and " << locusNames.size() << " ***" << endl;
				cin.clear();
				userInput.clear();
				done = false;
				continue;
			}
			else
			{
				sort(userInput.begin(), userInput.end());
// Order from 'largest' to 'smallest' index; otherwise need keep track of changing indices
				reverse(userInput.begin(), userInput.end());
				
				for (int editIter = 0; editIter < numEntries; editIter++)
				{
					int currentLocusID = userInput[editIter] - 1;
					for (int i = 0; i < numTaxa; i++)
					{
						data[i].erase(data[i].begin()+currentLocusID);
					}
					locusNames.erase(locusNames.begin()+currentLocusID);
					locusWeights.erase(locusWeights.begin()+currentLocusID);
				}
				userInput.clear();
			}
		}
		done = true;
	}
}

void excludeTaxa (vector < vector <int> > & data, vector <string> & taxonNames, vector <double> & taxonWeights,
	double & revisedCoverage, vector <string> const& locusNames)
{
	bool done = false;
	vector <int> userInput;
	bool validIntEntry = false;
	int numEntries = 0;
	int numPartitions = (int)locusNames.size();
	int userSelection;
	string taxonName;
	vector <string> excludedTaxa;
	
	while (!done)
	{
		userInput.clear();
		bool validChoice = false;
		while (!validChoice)
		{
			char userChoice;
			int numTaxa = (int)taxonNames.size();
			getCoverage (data, revisedCoverage);
			cout << endl << endl << "EXCLUDE TAXA FROM MATRIX" << endl << endl;
			cout << "Matrix currently contains " << numTaxa << " taxa." << endl
			<< "Matrix coverage is currently at: " << revisedCoverage << endl << endl
			<< "Exclude:" << endl
			<< " Taxa by [I]ndex" << endl
			<< " Taxa by [N]ame" << endl
			<< " Taxa missing [E]xactly N partitions" << endl
			<< " Taxa missing N or [M]ore partitions" << endl
			<< " Taxa possessing [O]nly N partitions" << endl
			<< " Taxa possessing N or [F]ewer partitions" << endl
			<< " Taxa possessing data for only a [S]pecific partition" << endl
			<< " Taxa exhibiting minimal [P]artition overlap" << endl
			<< " [B]ack to main menu" << endl
			<< endl
			<< "Enter desired option: ";
			
			cin >> userChoice;
			cin.ignore(200, '\n');
			
			if (checkCharValue(userChoice,'i')) // by index
			{
				validChoice = true;
				cout << endl << "Taxa available:" << endl;
				printVectorAsList(taxonNames);
				cout << endl << "Enter indicies of taxa you would like to exclude (separated by spaces), or 0 to exit: ";
				numEntries = 0;	
				string tempString;
				getline(cin,tempString);
				
		// Ability to read in multiple inputs
				validIntEntry = false;
				istringstream tempStream(tempString);
				while (tempStream >> userSelection)
				{
					validIntEntry = true;
					userInput.push_back(userSelection);
					numEntries++;
				}
				if (userSelection == 0)
				{
					done = true;
					continue;
				}
				if (!validIntEntry)
				{
					cout << endl << "*** Invalid input. Integer must be between 0 and " << taxonNames.size() << " ***" << endl;
					cin.clear();
					userInput.clear();
					done = false;
					continue;
				}
				if (validIntEntry)
				{
					if (cin.fail() || userSelection < 0 || userSelection > int(taxonNames.size()))
					{
						cout << endl << "*** Invalid input. Integer must be between 0 and " << taxonNames.size() << " ***" << endl;
						cin.clear();
						userInput.clear();
						done = false;
						continue;
					}
					else
					{
						sort(userInput.begin(), userInput.end());
		// Order from 'largest' to 'smallest' index; otherwise need keep track of changing indices
						reverse(userInput.begin(), userInput.end());
						
						for (int editIter = 0; editIter < numEntries; editIter++)
						{
							int currentTaxonID = userInput[editIter] - 1;
							data.erase(data.begin()+currentTaxonID);
							taxonNames.erase(taxonNames.begin()+currentTaxonID);
							taxonWeights.erase(taxonWeights.begin()+currentTaxonID);
						}
						userInput.clear();
					}
				}
			}
			else if (checkCharValue(userChoice,'n')) // by name
			{
				cout << endl << "Taxa available:" << endl;
				printVectorAsList(taxonNames);
				cout << endl << "Enter name(s) of taxa you would like to exclude (separated by spaces), or 0 to exit: ";
				numEntries = 0;	
				string tempString;
				getline(cin,tempString);
				vector <string> tempVector;
				
		// Ability to read in multiple inputs
				istringstream tempStream(tempString);
				while (tempStream >> taxonName)
				{
					tempVector.push_back(taxonName);
					numEntries++;
				}
				if (userSelection == 0)
				{
					done = true;
					continue;
				}
				else
				{
					vector <int> matched;
					for (int editIter = 0; editIter < numEntries; editIter++)
					{
						string toMatch = tempVector[0];
						bool match = false;
						for (int taxonIter = 0; taxonIter < numTaxa; taxonIter++)
						{
							if (caseInsensitiveStringCompare(toMatch, taxonNames[taxonIter]))
							{
								match = true;
								matched.push_back(taxonIter);
								taxonIter = numTaxa;
								continue;
							}
						}
						if (!match)
						{
							cout << endl << "OOPS!!! No match found for entry '" << toMatch << "'. Typo?" << endl;
						}
					}
					if (matched.size() > 0)
					{
						sort (matched.begin(), matched.end());
						reverse(matched.begin(), matched.end());
						
						for (int editIter = 0; editIter < int(matched.size()); editIter++)
						{
							data.erase(data.begin()+editIter);
							taxonNames.erase(taxonNames.begin()+editIter);
							taxonWeights.erase(taxonWeights.begin()+editIter);
						}
					}
				}
				tempVector.clear();
			}
			else if (checkCharValue(userChoice,'e')) // exactly N missing
			{
				int partitionsMissing = checkValidIntInput ("Enter exact number of missing partitions to consider: ");
				
// Put in error-checking regarding number of actual genes
				
				excludeTaxaMissingNGenes (partitionsMissing, data, taxonNames, taxonWeights, true);
			}
			else if (checkCharValue(userChoice,'m')) // N or more missing
			{
				int partitionsMissing = checkValidIntInput ("Enter minimum number of missing partitions to consider: ");
				excludeTaxaMissingNGenes (partitionsMissing, data, taxonNames, taxonWeights, false);
			}
			else if (checkCharValue(userChoice,'o')) // possess only N
			{
				int partitionsPossessed = checkValidIntInput ("Enter exact number of possessed partitions to consider: ");
				
// Put in error-checking regarding number of actual genes
				
				excludeTaxaPossessingNGenes (partitionsPossessed, data, taxonNames, taxonWeights, true);
			}
			else if (checkCharValue(userChoice,'f')) // possess N or fewer
			{
				int partitionsPossessed = checkValidIntInput ("Enter maximum number of possessed partitions to consider: ");
				excludeTaxaPossessingNGenes (partitionsPossessed, data, taxonNames, taxonWeights, false);
			}
			else if (checkCharValue(userChoice,'s')) // only present for a specific partition
			{
				validChoice = true;
				
				cout << endl << "Partitions available:" << endl;
				printVectorAsList(locusNames);
				
				cout << endl << "Enter indicies of partitions you would like to exclude (separated by spaces), or 0 to exit: ";
				numEntries = 0;	
				string tempString;
				getline(cin,tempString);
				
		// Ability to read in multiple inputs
				validIntEntry = false;
				istringstream tempStream(tempString);
				while (tempStream >> userSelection)
				{
					validIntEntry = true;
					userInput.push_back(userSelection);
					numEntries++;
				}
				if (userSelection == 0)
				{
					done = true;
					continue;
				}
				if (!validIntEntry)
				{
					cout << endl << "*** Invalid input. Integer must be between 0 and " << numPartitions << " ***" << endl;
					cin.clear();
					userInput.clear();
					done = false;
					continue;
				}
				if (validIntEntry)
				{
					if (cin.fail() || userSelection < 0 || userSelection > numPartitions)
					{
						cout << endl << "*** Invalid input. Integer must be between 0 and " << numPartitions << " ***" << endl;
						cin.clear();
						userInput.clear();
						done = false;
						continue;
					}
					else
					{
						vector <int> excludeTaxa;
						int numTaxa = (int)taxonNames.size();
						
						for (int editIter = 0; editIter < numEntries; editIter++)
						{
							int curPart = userInput[editIter] - 1;
							for (int taxIter = 0; taxIter < numTaxa; taxIter++)
							{
								int sum = 0;
								if (data[taxIter][curPart] == 1)
								{
									for (int partIter = 0; partIter < numPartitions; partIter++)
									{
										if (data[taxIter][partIter] == 1)
										{
											sum +=1;
										}
									}
								}
								if (sum == 1) // only present for focal gene
								{
									excludeTaxa.push_back(taxIter);
									excludedTaxa.push_back(taxonNames[taxIter]);
								}
							}
						}
						if (excludeTaxa.size() != 0)
						{
							reverse(excludeTaxa.begin(), excludeTaxa.end());
							for (int editIter = 0; editIter < int(excludeTaxa.size()); editIter++)
							{
								int currentTaxon = excludeTaxa[editIter];
								data.erase(data.begin()+currentTaxon);
								taxonNames.erase(taxonNames.begin()+currentTaxon);
								taxonWeights.erase(taxonWeights.begin()+currentTaxon);
							}
						}
						else
						{
							cout << "No currently implemented taxa conform to condition" << endl;
						}
						userInput.clear();
					}
				}
				cout << endl << "Excluded taxa:" << endl;
				printVectorAsList (excludedTaxa);
				excludedTaxa.clear();
			}
			else if (checkCharValue(userChoice,'p'))
			{
				excludeTaxaMinimalOverlap (data, taxonNames, taxonWeights);
			}
			else if (checkCharValue(userChoice,'b'))
			{
				validChoice = true;
				done = true;
				continue;
			}
			else
			{
				cout << "Invalid input option (" << userChoice << "). Try again." << endl;
			}
		}
	}
}

void excludeTaxaMissingNGenes (int const& partitionsMissing, vector < vector <int> > & data,
	vector <string> & taxonNames, vector <double> & taxonWeights, bool const& missingExact)
{
	vector <string> excludedTaxa;
	vector <int> excludedIndices;
	int numTaxa = (int)data.size();
	int numPartitions = (int)data[0].size();
	
	if (missingExact)
	{
		for (int taxonIter = 0; taxonIter < numTaxa; taxonIter++)
		{
			int sum = 0;
			for (int partitionIter = 0; partitionIter < numPartitions; partitionIter++)
			{
				sum += data[taxonIter][partitionIter];
			}
			if (sum == (numPartitions - partitionsMissing))
			{
				excludedIndices.push_back(taxonIter);
				excludedTaxa.push_back(taxonNames[taxonIter]);
			}
		}
	}
	else // minimum bound
	{
		for (int taxonIter = 0; taxonIter < numTaxa; taxonIter++)
		{
			int sum = 0;
			for (int partitionIter = 0; partitionIter < numPartitions; partitionIter++)
			{
				sum += data[taxonIter][partitionIter];
			}
			if (sum <= (numPartitions - partitionsMissing))
			{
				excludedIndices.push_back(taxonIter);
				excludedTaxa.push_back(taxonNames[taxonIter]);
			}
		}
	}
	if (excludedIndices.size() != 0)
	{
		int numExcluded = (int)excludedTaxa.size();
		cout << endl << "Excluded taxa:" << endl;
		printVectorAsList (excludedTaxa);
		
// reverse-sort excludedTaxa vector
		reverse(excludedIndices.begin(),excludedIndices.end());
		for (int i = 0; i < numExcluded; i++)
		{
			data.erase(data.begin()+excludedIndices[i]);
			taxonNames.erase(taxonNames.begin()+excludedIndices[i]);
			taxonWeights.erase(taxonWeights.begin()+excludedIndices[i]);
		}
	}
	else
	{
		cout << "No currently implemented taxa conform to condition." << endl;
	}
}

void excludeTaxaPossessingNGenes (int const& partitionsPossessed, vector < vector <int> > & data,
	vector <string> & taxonNames, vector <double> & taxonWeights, bool const& possessingExact)
{
	vector <string> excludedTaxa;
	vector <int> excludedIndices;
	int numTaxa = (int)data.size();
	int numPartitions = (int)data[0].size();
	
	if (possessingExact)
	{
		for (int taxonIter = 0; taxonIter < numTaxa; taxonIter++)
		{
			int sum = 0;
			for (int partitionIter = 0; partitionIter < numPartitions; partitionIter++)
			{
				sum += data[taxonIter][partitionIter];
			}
			if (sum == partitionsPossessed)
			{
				excludedIndices.push_back(taxonIter);
				excludedTaxa.push_back(taxonNames[taxonIter]);
			}
		}
	}
	else // minimum bound
	{
		for (int taxonIter = 0; taxonIter < numTaxa; taxonIter++)
		{
			int sum = 0;
			for (int partitionIter = 0; partitionIter < numPartitions; partitionIter++)
			{
				sum += data[taxonIter][partitionIter];
			}
			if (sum <= partitionsPossessed)
			{
				excludedIndices.push_back(taxonIter);
				excludedTaxa.push_back(taxonNames[taxonIter]);
			}
		}
	}
	if (excludedIndices.size() != 0)
	{
		int numExcluded = (int)excludedTaxa.size();
		cout << endl << "Excluded taxa:" << endl;
		printVectorAsList (excludedTaxa);
		
// reverse-sort excludedTaxa vector
		reverse(excludedIndices.begin(),excludedIndices.end());
		for (int i = 0; i < numExcluded; i++)
		{
			data.erase(data.begin()+excludedIndices[i]);
			taxonNames.erase(taxonNames.begin()+excludedIndices[i]);
			taxonWeights.erase(taxonWeights.begin()+excludedIndices[i]);
		}
	}
	else
	{
		cout << "No currently implemented taxa conform to condition." << endl;
	}
}

void excludeTaxaMinimalOverlap (vector < vector <int> > & data, vector <string> & taxonNames,
	vector <double> & taxonWeights)
{
	int numTaxa = (int)data.size();
	int numPartitions = (int)data[0].size();
	vector <int> counts;
	int minCount = numTaxa;
	vector <string> excludedTaxa;
	double minOverlap = 0.0;
	bool invokeThreshold = false;
	char userChoice;
	double thresholdValue = 0.0;
	int thresholdCount = 0;
	
	cout << "Exclude [M]inimally overlapping taxa (i.e. just the worst) or implement [T]hreshold: ";
	cin >> userChoice;
	
	if (checkCharValue(userChoice,'t'))
	{
		invokeThreshold = true;
		cout << "Enter threshold value (minimum proportion of taxonomic overlap): ";
		cin >> thresholdValue;
		
		thresholdCount = (int)thresholdValue * numTaxa;
		
		if (debuggering) {cout << "thresholdCount calculated to be: " << thresholdCount << endl;}
	}
	
// Count
	for (int taxonIter = 0; taxonIter < numTaxa; taxonIter++)
	{
		int sum = 0;
		vector <int> temp (numTaxa, 0);
		
		for (int partitionIter = 0; partitionIter < numPartitions; partitionIter++)
		{
			if (data[taxonIter][partitionIter] == 1)
			{
				for (int taxIter = 0; taxIter < numTaxa; taxIter++)
				{
					if (data[taxIter][partitionIter] == 1)
					{
						temp[taxIter] = 1;
					}
				}
			}
		}
		
		sum = accumulate(temp.begin(), temp.end(), 0);
		
		if (debuggering) {cout << "Overlap count for taxon '" << taxonNames[taxonIter] << "' is: " << sum << endl;}
		
		counts.push_back(sum);
		if (sum < minCount)
		{
			minCount = sum;
		}
	}
	
	if (debuggering) {cout << "minCount = " << minCount << endl;}
	
	minOverlap = double(minCount) / double(numTaxa);
	
	if (minCount < numTaxa)
	{
		if (!invokeThreshold) // just take the lowest motherfuckers
		{
			cout << endl << "Minimum taxon overlap is: " << minOverlap << endl;
			for (int taxonIter = numTaxa-1; taxonIter >= 0; taxonIter--)
			{
				if (debuggering) {cout << "taxon '" << taxonNames[taxonIter] << "' has count: " << counts[taxonIter] << endl;}
				if (counts[taxonIter] == minCount)
				{
					excludedTaxa.push_back(taxonNames[taxonIter]);
					if (debuggering) {cout << "DELETE! taxonIter = " << taxonIter << "; taxon = " << taxonNames[taxonIter] << endl;}
					data.erase(data.begin()+taxonIter);
					taxonNames.erase(taxonNames.begin()+taxonIter);
					taxonWeights.erase(taxonWeights.begin()+taxonIter);
				}
			}
			cout << endl << "Excluded taxa:" << endl;
			printVectorAsList (excludedTaxa);
		}
		else
		{
			cout << endl << "Threshold taxon overlap is: " << thresholdValue << endl;
			for (int taxonIter = numTaxa-1; taxonIter >= 0; taxonIter--)
			{
				if (debuggering) {cout << "taxon '" << taxonNames[taxonIter] << "' has count: " << counts[taxonIter] << endl;}
				if (counts[taxonIter] <= thresholdCount)
				{
					excludedTaxa.push_back(taxonNames[taxonIter]);
					if (debuggering) {cout << "DELETE! taxonIter = " << taxonIter << "; taxon = " << taxonNames[taxonIter] << endl;}
					data.erase(data.begin()+taxonIter);
					taxonNames.erase(taxonNames.begin()+taxonIter);
					taxonWeights.erase(taxonWeights.begin()+taxonIter);
				}
			}
			if (!excludedTaxa.empty())
			{
				cout << endl << "Excluded taxa:" << endl;
				printVectorAsList (excludedTaxa);
			}
			else
			{
				cout << endl << "No taxon meets threshold criterion. Your data are better than that. Yay!" << endl;
			}
		}
	}
	else
	{
		cout << endl << "Taxon overlap is complete; no taxa to delete. Matrix is completely decisive, dude." << endl;
	}
}

// create chimeric taxa
void mergeTaxa (vector < vector <int> > & data, vector <string> & taxonNames, vector <double> & revisedTaxonWeights)
{
	bool done = false;
	vector <int> userInput;
	bool validIntEntry = false;
	int numEntries = 0;
	int taxonSelection;
	
	while (!done)
	{
		cout << endl << "Taxa available:" << endl;
		printVectorAsList(taxonNames);
		
		cout << endl << "Enter indicies of taxa you would like to merge (separated by spaces), or 0 to exit: ";
		numEntries = 0;	
		string tempString;
		getline(cin,tempString);
		
// Ability to read in multiple inputs
		validIntEntry = false;
		istringstream tempStream(tempString);
		while (tempStream >> taxonSelection)
		{
			validIntEntry = true;
			userInput.push_back(taxonSelection);
			numEntries++;
		}
		if (taxonSelection == 0)
		{
			done = true;
			continue;
		}
		if (!validIntEntry)
		{
			cout << endl << "*** Invalid input. Integer must be between 0 and " << taxonNames.size() << " ***" << endl;
			cin.clear();
			userInput.clear();
			done = false;
			continue;
		}
		if (validIntEntry)
		{
			if (cin.fail() || taxonSelection < 0 || taxonSelection > int(taxonNames.size()))
			{
				cout << endl << "*** Invalid input. Integer must be between 0 and " << taxonNames.size() << " ***" << endl;
				cin.clear();
				userInput.clear();
				done = false;
				continue;
			}
			else if (numEntries < 2)
			{
				cout << endl << "Must reference at least 2 taxa. Try again." << endl;
				cin.clear();
				userInput.clear();
				done = false;
				continue;
			}
			else
			{
				sort(userInput.begin(), userInput.end());
// Order from 'largest' to 'smallest' index; otherwise need keep track of changing indices
				reverse(userInput.begin(), userInput.end());
				
				string mergedTaxonName;
				vector < vector <int> > tempIntVector;
				vector <double> tempWeightVector;
				for (int editIter = 0; editIter < numEntries; editIter++)
				{
					int currentTaxonID = userInput[editIter] - 1;
					if (editIter == 0)
					{
						mergedTaxonName = taxonNames[currentTaxonID];
					}
					else
					{
						mergedTaxonName = mergedTaxonName + "_+_" + taxonNames[currentTaxonID];
					}
					tempIntVector.push_back(data[currentTaxonID]);
					tempWeightVector.push_back(revisedTaxonWeights[currentTaxonID]);
					data.erase(data.begin()+currentTaxonID);
					taxonNames.erase(taxonNames.begin()+currentTaxonID);
					revisedTaxonWeights.erase(revisedTaxonWeights.begin()+currentTaxonID);
				}
	// Long names can be cumbersome when printing; offer rename
				bool renameMergedTaxon = false;
				string query = "Default name for chimeric taxon is '" + mergedTaxonName + "'. " + 
					"Enter (0) to keep this name, or (1) to provide a new name: ";
				renameMergedTaxon = checkValidBoolInput(query);
				if (renameMergedTaxon)
				{
					cout << "Enter name for chimeric taxon: ";
					cin >> mergedTaxonName;
				}
				sort(tempWeightVector.begin(), tempWeightVector.end());
				reverse(tempWeightVector.begin(), tempWeightVector.end());
				
				int numLoci = (int)tempIntVector[0].size();
				vector <int> mergedData (numLoci, 0);
				for (int i = 0; i < numLoci; i++) // Merge data
				{
					for (int j = 0; j < int(tempIntVector.size()); j++)
					{
						if (tempIntVector[j][i] == 1)
						{
							mergedData[i] = 1;
							j = int(tempIntVector.size()); // Match found; exit
						}
					}
				}
				data.push_back(mergedData);
				taxonNames.push_back(mergedTaxonName);
				revisedTaxonWeights.push_back(tempWeightVector[0]); // Apply highest weight from merged taxa
				userInput.clear();
				tempWeightVector.clear();
			}
		}
		done = true;
	}
}
