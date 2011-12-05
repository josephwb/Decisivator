#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>

using namespace std;

#include "General.h"
#include "Parse_Nexus.h"
#include "Parse_Data.h"
#include "User_Interface.h"

extern bool debuggering;

/*
  NOTE: deals only with the vanilla-iest of Nexus file formats at the moment.
*/

void parseNexus (string const& nexusFileName, vector < vector <int> > & data, vector <string> & taxonNames,
	vector <string> & locusNames, int & numChar, vector < vector <string> > & taxaAlignment,
	vector < vector <int> > & includedLocusRanges)
{
//	vector < vector <int> > includedLocusRanges;
	
	bool interleavedData = false;
	string fileName;
	int numTaxa = 0;
	
	getNumTaxaChar(nexusFileName, numTaxa, numChar, interleavedData);
	cout << endl << "Parsing Nexus file of " << numTaxa << " taxa and " << numChar << " characters..." << endl;
	taxaAlignment = collectTaxaAlignment(nexusFileName, numTaxa, numChar, interleavedData, taxonNames);
	locusNames = collectCharsets(nexusFileName, locusNames, includedLocusRanges, numChar);
	printCollectedCharsets(locusNames, includedLocusRanges);
	constructMatrix (taxaAlignment, includedLocusRanges, data, locusNames);
}


// Problems parsing when spaces surround '='
void getNumTaxaChar (string fileName, int & numTaxa, int & numChar, bool & interleavedData)
{
	ifstream inputUserFile;
	bool commentLine = false;
	bool whiteSpaceOnly = false;
	bool numTaxaEncountered = false;
	bool numCharEncountered = false;
	bool equalSignEncountered = false;
	bool semicolonEncountered = false;
	bool done = false; // know to stop looking
	
	inputUserFile.open(fileName.c_str());
	string line;
	
// Looking for pattern like 'dimensions ntax=53 nchar=16620;'
// 	- can be in either order, but must be stated on same line (for now)
// 	- no spaces allowed next to equal sign (for now)

// Looking for pattern like 'Format datatype=dna [gap=-] [missing=?] {[interleave=yes] or [interleave]};'
// 	- no spaces allowed next to equal sign (for now)
	
	while (getline(inputUserFile,line) && !done)
	{
		if (debuggering) {cout << "Current line: " << line << endl;}
		int stringPosition = 0;
		commentLine = checkCommentLine(line);
		whiteSpaceOnly = checkWhiteSpaceOnly(line);
		if (line.empty() || commentLine || whiteSpaceOnly)
		{
			continue;
		}
		else
		{
			if (checkStringValue(line, "matrix", stringPosition))	// Done - won't find the information further down; really only used for 'interleave'
			{
				done = true;
				if (debuggering) {cout << "Encountered 'matrix'" << endl;}
				continue;
			}
			else if (checkStringValue(line, "dimensions", stringPosition))
			{
				if (debuggering) {cout << "Encountered 'dimensions'" << endl;}
				stringPosition = 0;
				while (!numTaxaEncountered || !numCharEncountered)
				{
					stringPosition++;
					string tempString = removeStringSuffix(extractStringElement(line,stringPosition), '=', equalSignEncountered);
					if (checkStringValue(tempString, "ntax", 0))
					{
						tempString = removeStringPrefix(extractStringElement(line,stringPosition), '=');
						tempString = removeStringSuffix(tempString, ';', numTaxaEncountered);
						numTaxa = convertStringtoInt(tempString);
						if (debuggering) {cout << "NTax = " << numTaxa << endl;}
						numTaxaEncountered = true;
					}
					if (checkStringValue(tempString, "nchar", 0))
					{
						tempString = removeStringPrefix(extractStringElement(line,stringPosition), '=');
						tempString = removeStringSuffix(tempString, ';', numCharEncountered);
						numChar = convertStringtoInt(tempString);
						if (debuggering) {cout << "NChar = " << numChar << endl;}
						numCharEncountered = true;
					}
				}
			}
			else if (checkStringValue(line, "format", stringPosition)) // check to see if interleaved
			{
				stringPosition = 0;
				while (!semicolonEncountered)
				{
					stringPosition++;
					string tempString = removeStringSuffix(extractStringElement(line,stringPosition), '=', equalSignEncountered); // where '=' is used i.e. 'interleave=yes;'
					if (checkStringValue(tempString, "interleave", 0))
					{
						tempString = removeStringPrefix(extractStringElement(line,stringPosition), '=');
						tempString = removeStringSuffix(tempString, ';', semicolonEncountered);
						if (checkStringValue(tempString, "yes", 0))
						{
							interleavedData = true;
							semicolonEncountered = true;
							cout << "Data are in interleaved format." << endl;
							continue;
						}
						else if (checkStringValue(tempString, "no", 0))
						{
							interleavedData = false;
							semicolonEncountered = true;
							cout << "Data are not in interleaved format." << endl;
							continue;
						}
					}
					tempString = removeStringSuffix(extractStringElement(line,stringPosition), ';', semicolonEncountered); // where '=' is NOT used i.e. 'interleave;'
					if (checkStringValue(tempString, "interleave", 0))											  // OR a space follows interleave declaration i.e. 'interleave [=];'
					{
						interleavedData = true;
						semicolonEncountered = true;
						cout << "Data are in interleaved format." << endl;
						continue;
					}
				}
			}
		}
	}
	inputUserFile.close();
	if (debuggering) {cout << endl;}
}

vector <string> collectCharsets (string charsetFileName, vector <string> inputCharsets,
	vector < vector <int> > & includedLocusRanges, int const& numChar)
{
	ifstream declaredCharsets;
	vector <string> collectedCharsets;
	vector <string> tempStringVector;
	vector <int> tempIntVector;
	bool commentLine = false;
	bool whiteSpaceOnly = false;
	declaredCharsets.open(charsetFileName.c_str());
	bool charsetsEncountered = false;
	string line;
	
// Read in every non-empty (or non-whitespace), non-commented-out line
	while (getline(declaredCharsets,line))
	{
		commentLine = checkCommentLine(line);
		whiteSpaceOnly = checkWhiteSpaceOnly(line);
		if (line.empty() || commentLine || whiteSpaceOnly)
		{
			continue;
		}
		else
		{
			tempStringVector.push_back(line);
		}
		line.clear();
	}
	
// Test each collected line for presence of CHARSET
// COLLECT CHARSET RANGES - START WITH SIMPLE X-Y DECLARATIONS
// will take the form of: 'CHARSET ALD = 1-566;'
	
	for (vector <string>::iterator lineIter = tempStringVector.begin(); lineIter < tempStringVector.end(); lineIter++)
	{
		string charsetName;
		string temp;
		int stringPosition = 0;
		bool semiColonEncountered = false;	// Check that charset declaration is complete
		string tempString;
		
		if (checkStringValue(*lineIter, "charset", stringPosition))
		{
			charsetsEncountered = true;
			bool equalSignEncountered = false;
			stringPosition++;
			
			charsetName = extractStringElement(*lineIter, stringPosition);
			stringPosition++;
			
			while (!equalSignEncountered)
			{
				equalSignEncountered = checkStringValue(*lineIter, "=", stringPosition);
				stringPosition++;
			}
// Parse and collect start/stop information
			while (!semiColonEncountered)
			{
				string charsetDeclaration = removeStringSuffix(extractStringElement(*lineIter,stringPosition), ';', semiColonEncountered);
				int start = 0;
				int stop = 0;
				int interval = 0;
// Need to check format i.e. '1-566;' vs '1-566\3;'
				bool separatorEncountered = false;
				tempString = removeStringSuffix(*lineIter, '\\', separatorEncountered);
				if (separatorEncountered)
				{
					if (debuggering) {cout << "Interval CHARSET encountered: '" << *lineIter << "." << endl;}
					extractIntervalRange(charsetDeclaration, start, stop, interval);
					tempIntVector.push_back(1);		// Code for interval range
					tempIntVector.push_back(start);
					tempIntVector.push_back(stop);
					tempIntVector.push_back(interval);
					continue;
				}
				else
				{
					extractSimpleCharsetRanges(charsetDeclaration, start, stop);
					tempIntVector.push_back(0);		// Code for simple range
					tempIntVector.push_back(start);
					tempIntVector.push_back(stop);
				}
			}
// Need to test for duplicates e.g. if CHARSETS are declared in both PAUP* and MrBayes blocks
			bool matchFound = false;
			for (vector <string>::const_iterator matchIter = collectedCharsets.begin(); matchIter < collectedCharsets.end(); matchIter++)
			{
				if (charsetName == *matchIter)
				{
					matchFound = true;
				}
			}
			if (!matchFound)
			{
				collectedCharsets.push_back(charsetName);
				includedLocusRanges.push_back(tempIntVector);
			}
			tempIntVector.clear();
		}
	}
// If no charsets detected, construct a single partition
	if (!charsetsEncountered)
	{
		collectedCharsets.push_back("All");
		tempIntVector.push_back(0);		// Code for simple range
		tempIntVector.push_back(1);
		tempIntVector.push_back(numChar);
		includedLocusRanges.push_back(tempIntVector);
	}
	declaredCharsets.close();
	return collectedCharsets;
}

void extractIntervalRange (string charsetDeclaration, int & start, int & stop, int & interval)
{
// This simple function assumes a very strict format: e.g. '1-566\3'
// - need to remove middle '-', separator '\', store interval
	string temp;
	vector<char> tempVector;
	int charCounter = 0;
	int dashPosition = 0;
	int separatorPosition = 0;
	
	for (string::iterator charIter = charsetDeclaration.begin(); charIter < charsetDeclaration.end(); charIter++ )
	{
		tempVector.push_back(*charIter);
		if (*charIter == '-')
		{
			dashPosition = charCounter;
		}
		if (*charIter == '\\')
		{
			separatorPosition = charCounter;
		}
		charCounter++;
	}
// Start position
	for (int charIter = 0; charIter < dashPosition; charIter++)
	{
		temp += tempVector[charIter];
	}
	start = convertStringtoInt(temp);
	temp.clear();
	
// Stop position
	for (int charIter = dashPosition + 1; charIter < charCounter; charIter++)
	{
		temp += tempVector[charIter];
	}
	stop = convertStringtoInt(temp);
	temp.clear();
	
// Interval position
	for (int charIter = separatorPosition + 1; charIter < charCounter; charIter++)
	{
		temp += tempVector[charIter];
	}
	interval = convertStringtoInt(temp);
}

void extractSimpleCharsetRanges (string charsetDeclaration, int & start, int & stop)
{
// This simple function assumes a very strict format: e.g. '1-566'
// - need to remove middle '-'
	string temp;
	vector<char> tempVector;
	int charCounter = 0;
	int dashPosition = 0;
	
	for (string::iterator charIter = charsetDeclaration.begin(); charIter < charsetDeclaration.end(); charIter++ )
	{
		tempVector.push_back(*charIter);
		if (*charIter == '-')
		{
			dashPosition = charCounter;
		}
		charCounter++;
	}
// Start position
	for (int charIter = 0; charIter < dashPosition; charIter++)
	{
		temp += tempVector[charIter];
	}
	start = convertStringtoInt(temp);
	temp.clear();
// Stop position
	for (int charIter = dashPosition + 1; charIter < charCounter; charIter++)
	{
		temp += tempVector[charIter];
	}
	stop = convertStringtoInt(temp);
}

int getIntervalPartitionLength (int const& partitionID, vector < vector <int> > const& includedLocusRanges)
{
	int partitionLength = 0;
	int start = includedLocusRanges[partitionID][1];
	int stop = includedLocusRanges[partitionID][2];
	int interval = includedLocusRanges[partitionID][3];
	for (int iter = start; iter <= stop; iter += interval)
	{
		partitionLength++;
	}
	return partitionLength;
}

void printCollectedCharsets (vector <string> const& collectedCharsets, vector < vector <int> > const& includedLocusRanges)
{
	int counter = 0;
	int sum = 0;
	cout << "READING IN CHARSETS..." << endl << endl;
	for (vector <string>::const_iterator charsetIter = collectedCharsets.begin(); charsetIter < collectedCharsets.end(); charsetIter++)
	{
		if (counter < 9)
		{
			cout << " ";
		}
		if (includedLocusRanges[counter][0] == 0)	// Code for simple ranges
		{
			int tempLength = includedLocusRanges[counter][2] - includedLocusRanges[counter][1] + 1;
			
			cout << "(" << counter + 1 << "): " << *charsetIter << ", positions " << includedLocusRanges[counter][1] << "-" << includedLocusRanges[counter][2] <<
			" (" << tempLength;
			if (tempLength == 1)
			{
				cout << " character)" << endl;
			}
			else
			{
				cout << " characters)" << endl;
			}
			sum += tempLength;
		}
		else if (includedLocusRanges[counter][0] == 1)	// Code for interval ranges
		{
			int tempLength = getIntervalPartitionLength(counter, includedLocusRanges);
			
			cout << "(" << counter + 1 << "): " << *charsetIter << ", positions " << includedLocusRanges[counter][1] << "-" << includedLocusRanges[counter][2] << ", every " << includedLocusRanges[counter][3]
			<< " sites (" << tempLength;
			if (tempLength == 1)
			{
				cout << " character)" << endl;
			}
			else
			{
				cout << " characters)" << endl;
			}
			sum += tempLength;
		}
		counter++;
	}
	cout << "Total number of partitions is: " << sum << endl << endl;
}

bool checkCommentLine (string stringToParse)
{
	bool commentLine = false;
	char firstCharacter = stringToParse[0];
	if (firstCharacter == '[')
	{
//		cout << "Dude, you've got yourself a comment there. Ignoring entire line, assuming comment does not extend across multiple lines." << endl;
		commentLine = true;
	}
	return commentLine;
}

string removeStringSuffix (string stringToParse, char suffixToRemove, bool & suffixEncountered)
{
	string temp;
	vector<char> tempVector;
	int charCounter = 0;
	int suffixStart = 0;
	suffixEncountered = false;
	
	for (string::iterator charIter = stringToParse.begin(); charIter < stringToParse.end(); charIter++ )
	{
		tempVector.push_back(*charIter);
		if (*charIter == suffixToRemove)
		{
			suffixStart = charCounter;
			suffixEncountered = true;
		}
		charCounter++;
	}
	if (suffixEncountered)
	{
		for (int charIter = 0; charIter < suffixStart; charIter++)
		{
			temp += tempVector[charIter];
		}
	}
	if (!suffixEncountered)
	{
		temp = stringToParse;
	}
	return temp;
}

string removeStringPrefix(string stringToParse, char characterToRemove)
{
	string tempString;
	vector<char> tempVector;
	int charCounter = 0;
	int characterPosition = 0;
	bool characterEncountered = false;
	
	for (string::iterator charIter = stringToParse.begin(); charIter < stringToParse.end(); charIter++ )
	{
		tempVector.push_back(*charIter);
		if (*charIter == characterToRemove)
		{
			characterPosition = charCounter;
			characterEncountered = true;
		}
		charCounter++;
	}
	if (characterEncountered)
	{
		for (int charIter = characterPosition + 1; charIter < charCounter; charIter++)
		{
			tempString += tempVector[charIter];
		}
	}
	if (!characterEncountered)
	{
		tempString = stringToParse;
	}
	return tempString;
}

vector < vector <string> > collectTaxaAlignment (string fileName, int const& numTaxa, int const& numChar,
	bool const& interleavedData, vector <string> & taxonNames)
{
// PLEASE NOTE: search strategy below uses very strict format assumptions - that which is exported by PAUP*
// 	- do not be surprised if this fucks up - it is probably a simple rearrangement of terms
	vector < vector <string> > taxaAlignment;
	ifstream inputAlignment;
	vector <string> tempStringVector;
	bool commentLine = false;
	bool whiteSpaceOnly = false;
	bool matrixEncountered = false;
	bool allCharacterRead = false;
	int numCharRead = 0;
	
	inputAlignment.open(fileName.c_str());
	string line;
	
	if (!interleavedData)
	{
		cout << endl << endl << "READING IN NON-INTERLEAVED DATA..." << endl;
// Ignore lines until 'matrix' is encountered
		while (!matrixEncountered)
		{
			getline(inputAlignment,line);
			commentLine = checkCommentLine(line);
			whiteSpaceOnly = checkWhiteSpaceOnly(line);
			
			if (line.empty() || commentLine || whiteSpaceOnly)
			{
				continue;
			}
			else if (checkStringValue(line, "matrix", 0))
			{
				matrixEncountered = true;
			}
			else
			{
				continue;
			}
		}
		numCharRead = 0;
// Read in every non-empty (or non-whitespace), non-commented-out line
		for (int taxonIter = 0; taxonIter < numTaxa; taxonIter++)
		{
			getline(inputAlignment,line);
			commentLine = checkCommentLine(line);
			whiteSpaceOnly = checkWhiteSpaceOnly(line);
			
			if (line.empty() || commentLine || whiteSpaceOnly)
			{
				taxonIter--;
				continue;
			}
			else
			{
// First string is taxon name, second is sequence
				if (debuggering) {cout << "Reading in taxon '" << extractStringElement(line, 0) << "'..." << endl;}
				tempStringVector.push_back(extractStringElement(line, 0));
				tempStringVector.push_back(extractStringElement(line, 1));
				taxonNames.push_back(extractStringElement(line, 0));
				
				taxaAlignment.push_back(tempStringVector);
				tempStringVector.clear();
				
// Count sites encountered - for error-checking; only checking first sequence (for now)
				if (taxonIter == 0)
				{
					int charCounter = (int)extractStringElement(line, 1).size();
					numCharRead += charCounter;
				}
				if (numCharRead != numChar)
				{
					cout << "numCharRead (" << numCharRead << ") does not equal numChar (" << numChar << ") declared in Nexus file. Likely a formatting issue (my bad)."
					<< endl << "Exiting." << endl;
					exit(1);
				}
			}
		}
		if (numCharRead == numChar)
		{
			allCharacterRead = true;
			if (debuggering) {cout << "numCharRead (" << numCharRead << ") == numChar (" << numChar << ") declared in Nexus file. Woo-hoo!" << endl;}
		}
	}
	else if (interleavedData)
	{
// Ignore lines until 'matrix' is encountered
		cout << endl << endl << "READING IN INTERLEAVED DATA..." << endl;
		while (!matrixEncountered)
		{
			getline(inputAlignment,line);
			commentLine = checkCommentLine(line);
			whiteSpaceOnly = checkWhiteSpaceOnly(line);
			
			if (line.empty() || commentLine || whiteSpaceOnly)
			{
				continue;
			}
			else if (checkStringValue(line, "matrix", 0))
			{
				matrixEncountered = true;
			}
			else
			{
				continue;
			}
		}
		bool firstPass = true;
		numCharRead = 0;
		while (!allCharacterRead)
		{
// Read in every non-empty (or non-whitespace), non-commented-out line
			for (int taxonIter = 0; taxonIter < numTaxa; taxonIter++)
			{
				getline(inputAlignment,line);
				commentLine = checkCommentLine(line);
				whiteSpaceOnly = checkWhiteSpaceOnly(line);
				
				if (line.empty() || commentLine || whiteSpaceOnly)
				{
					taxonIter--;
					continue;
				}
				else
				{
// First string is taxon name, second is sequence
					string taxonName = (extractStringElement(line, 0));
					string taxonSequence = (extractStringElement(line, 1));
					
					if (firstPass)
					{
						if (debuggering) {cout << "Reading in (interleaved) taxon '" << extractStringElement(line, 0) << "'..." << endl;}
						taxonNames.push_back(taxonName);
						tempStringVector.push_back(taxonName);		// Taxon name
						tempStringVector.push_back(taxonSequence);	// sequence
						taxaAlignment.push_back(tempStringVector);
					}
					if (!firstPass)
					{
						taxaAlignment[taxonIter][1] += taxonSequence;
					}
// Count sites encountered - for error-checking; only checking first sequence (for now)
					if (taxonIter == 0)
					{
						int charCounter = 0;
						for (string::size_type iterCharacters = 0; iterCharacters < taxonSequence.size(); iterCharacters++)
						{
							charCounter++;
						}
						numCharRead += charCounter;
					}
					tempStringVector.clear();
				}
			}
			firstPass = false;
			if (numCharRead == numChar)
			{
				allCharacterRead = true;
				if (debuggering) {cout << "numCharRead (" << numCharRead << ") == numChar (" << numChar << ") declared in Nexus file. Woo-hoo!" << endl;}
			}
// *** Need to print out some error here if not all characters are read ***
// At the moment just segfaults
		}
	}
	cout << endl;
	return taxaAlignment;
}

void constructMatrix (vector < vector <string> > const& taxaAlignment, vector < vector <int> > const& includedLocusRanges,
	vector < vector <int> > & data, vector <string> const& locusNames)
{
	int numMissing = 0;
	for (int i = 0; i < int(taxaAlignment.size()); i++)
	{
		vector <int> temp;
		for (int j = 0; j < int(includedLocusRanges.size()); j++)
		{
			bool sequencePresent = false;
			int begin = includedLocusRanges[j][1];
			int end = includedLocusRanges[j][2];
			int interval = 1;
			if (includedLocusRanges[j][0] == 1) // Code for interval range
			{
				interval = includedLocusRanges[j][3];
			}
			for (int k = begin - 1; k <= end - 1; k += interval)
			{
//				if (taxaAlignment[i][1][k] == 'A' || taxaAlignment[i][1][k] == 'C' || taxaAlignment[i][1][k] == 'G' || taxaAlignment[i][1][k] == 'T')
				if (validCharacterEncountered(taxaAlignment[i][1][k]))
				{
					sequencePresent = true;
					k = end;	// leave; data are present
					temp.push_back(1);
					continue;
				}
			}
			if (!sequencePresent)
			{
				temp.push_back(0);
				if (debuggering) {cout << "No sequence for taxon " << taxaAlignment[i][0] << " for locus " << locusNames[j] << "." << endl;}
				numMissing++;
			}
		}
		data.push_back(temp);
		temp.clear();
	}
	cout << "Data contains " << data.size() << " taxa, and " << data[0].size() << " partitions." << endl;
	cout << "There are " << numMissing << " missing taxon-partitions." << endl << endl;
}

// Made more general; can now handle any kind of data
bool validCharacterEncountered (char const& character)
{
	bool match = false;
	if (character != '?' && character != '-' && character != 'N')
	{
		match = true;
	}
	else
	{
		if(debuggering) {cout << "Character '" << character << "' is invalid." << endl;}
	}
	return match;
}