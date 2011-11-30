#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <cstdlib>

using namespace std;

#include "General.h"

extern bool DEBUG;

// Perform case-insenstive char match test
bool checkCharValue (char const& charInput, char const& charToMatch)
{
	if (toupper(charInput) == toupper(charToMatch))
	{
		return true;;
	}
	else
	{
		return false;
	}
}

string getFileName ()
{
	string fileName;
	bool validFileName = false;
	string tempFileName;
	cout << endl << "What is the name of the input file?" << endl;
	
	while (!validFileName)
	{
		cin >> tempFileName;
		validFileName = checkValidFile(tempFileName);
		if (cin.fail())
		{
			cin.clear(); 
			cin.ignore(200, '\n');
			cout << "Invalid input! Try again." << endl << endl;
			cout << "What is the name of the input file?" << endl;
		}
		else
		{
			fileName = tempFileName;
			cin.clear();
			cin.ignore(200, '\n');
		}
	}
	return fileName;
}

bool checkValidFile (string fileName)
{
	ifstream inFile;
	bool validFile = true;
	
	inFile.open(fileName.c_str());
	if (inFile.fail())
	{
		ofstream errorReport("Error.Decisivator.txt");
		errorReport << "Decisivator analysis failed." << endl << "Error: unable to open file '";
		errorReport << fileName << "'" << endl;
		errorReport.close();

		cerr << endl << "Decisivator analysis failed." << endl << "Error: unable to open file '";
		cerr << fileName << "'" <<  endl;
		exit(1);
	}
	else
	{
		cout << "Successfully opened file '" << fileName << "'." <<  endl;
		inFile.close();
		inFile.clear();
	}
	return validFile;
}

bool checkValidOutputFile (string & outputFileName)
{
	bool testOutBool = true;
	bool fileNameAcceptable = false;
	bool keepFileName = false;
	
// First, check if file already exists, so overwriting can be prevented
	fstream testIn;
	while (!fileNameAcceptable)
	{
		testIn.open(outputFileName.c_str());
		if (!testIn)
		{
			testIn.close();
			fileNameAcceptable = true;
		}
		else
		{
			testIn.close();
			cout << endl << "File exists!  Change name (0) or overwrite (1)? ";
			cin >> keepFileName;
			if (!keepFileName)
			{
				cout << "Enter new output file name: ";
				cin >> outputFileName;
			}
			else
			{
				cout << "Overwriting existing file '" << outputFileName << "'." << endl;
				fileNameAcceptable = true;
			}
		}
	}
	
	ofstream outFile;
	outFile.open(outputFileName.c_str());
	
	if (outFile.fail())
	{
		ofstream errorReport("Error.BEASTifier.txt");
		errorReport << "BEASTifier analysis failed." << endl << "Error: unable to open file '";
		errorReport << outputFileName << "'" << endl;
		errorReport.close();

		cerr << endl << "BEASTifier analysis failed." << endl << "Error: unable to open file '";
		cerr << outputFileName << "'" <<  endl;
		testOutBool = false;
		exit(1);
	}
	else
	{
		outFile.close();
		outFile.clear();
	}
	return testOutBool;
}

int checkValidIntInput (string queryString)
{
	int userInput = 0;
	bool validInput = false;
	cout << endl << queryString;
	validInput = false;
	while (!validInput)
	{
		cin >> userInput;
		if (cin.fail())
		{
			cin.clear(); 
			cin.ignore(200, '\n');
			cout << endl << "Invalid input! Must be an integer. Try again." << endl;
			cout << queryString;
		}
		else
		{
			validInput = true;
			cin.clear(); 			// Get rid of any extra characters accidentally entered
			cin.ignore(200, '\n');
		}
	}
	return userInput;
}

bool checkValidBoolInput (string queryString)
{
	int userInput = 0;
	bool validInput = false;
	cout << queryString;
	validInput = false;
	while (!validInput)
	{
		cin >> userInput;
		if (cin.fail())
		{
			cin.clear(); 
			cin.ignore(200, '\n');
			cout << endl << "Invalid input! Must be an integer. Try again." << endl;
			cout << queryString;
		}
		else if (userInput < 0 || userInput > 1)
		{
			cin.clear(); 
			cin.ignore(200, '\n');
			cout << endl << "Invalid input! Boolean value - must be 0 or 1. Try again." << endl;
			cout << queryString;
		}
		else
		{
			validInput = true;
			cin.clear(); 			// Get rid of any extra characters accidentally entered
			cin.ignore(200, '\n');
		}
	}
	return userInput;
}

bool checkStringValue (string stringToParse, string stringToMatch, int stringPosition)
{
// Performs case-insenstive string match test
	string testString = extractStringElement(stringToParse, stringPosition);
	if (testString.size() != stringToMatch.size())
	{
		return false;
	}
	if (testString == stringToMatch)
	{
		return true;
	}
	
	for (size_t i = 0; i < testString.size(); ++i)
	{
// 
		if (testString[i] == stringToMatch[i] || testString[i]+32 == stringToMatch[i] || testString[i] == stringToMatch[i]+32)
		{
			continue;
		}
		else
		{
			return false;
		}
	}
	return true;
}

BigInt choose (int const& n, int const& r)
{
	BigInt result = 0;
	if (r == 0) {return (1);}
	
	BigInt numerator = n;
	for (int i = n - 1; i > n - r; i --)
	{
		numerator *= i;
	}
	BigInt denominator = factorial(r);
	result = numerator / denominator;
	return result;
}

BigInt factorial (int const& num)
{
	if (num <= 1) return (1);
	return factorial(num - 1) * num; // recursive call
}

/*
BigInt Stirling2ndKind (int const& n, int const& k)
{
	BigInt S2K = 0;
	
	if (n == k) {return (1);}
	for (int i = 0; i < n; i++)
	{
		S2K += pow(-1,i) * choose(k,i) * pow(k-i,n);
	}
	S2K = S2K / factorial(k);
	return(S2K);
}
*/

// Given some existing string, extract (copy) the ith element and store as new string
string extractStringElement (string & stringToParse, int const& position)
{
	vector <string> tempVector = storeStringVector(stringToParse);
	string returnString = tempVector[position];
	
	return returnString;
}

// At the moment, assume matrix is filled with "1" and "0" only; generalize later (e.g. "?", "NA", etc.)
int convertStringtoInt (string stringToConvert)
{
	int tempInt = 0;
	istringstream tempStream(stringToConvert);
	tempStream >> tempInt;
	
	return tempInt;
}

double convertStringtoDouble (string stringToConvert)
{
	double tempDouble = 0;
	istringstream tempStream(stringToConvert);
	tempStream >> tempDouble;
	
	return tempDouble;
}

// Store elements of a string in a vector
vector <string> storeStringVector (string & stringToParse)
{
	vector <string> tempVector;
	istringstream tempStream(stringToParse);
	string tempString;
	string returnString;
	
	while (tempStream >> tempString)
	{
		tempVector.push_back(tempString);
	}
	
	return tempVector;
}

// Print purdy lists, aligning element names
void printVectorAsList (vector <string> const& vectorToPrint)
{
	string tempString = vectorToPrint[0];
	int numElements = vectorToPrint.size();
	string maxString;
	int longestName = 0;
	
	for (int i = 0; i < numElements; i++)
	{		
		string currentString = vectorToPrint[i];
		if (currentString.size() > maxString.size())
		{
			maxString = currentString;
		}
	}
// Determine length of longest name
	for (string::const_iterator iterCharacters = maxString.begin(); iterCharacters < maxString.end(); ++iterCharacters)
	{
		longestName++;
	}
	
	for (int i = 0; i < numElements; i++)
	{		
		cout << "  ";
		if (numElements >= 1000)
		{
			if (i + 1 < 10) {cout << "   ";}
			else if (i + 1 < 100) {cout << "  ";}
			else if (i + 1 < 1000) {cout << " ";}
		}
		else if (numElements >= 100)
		{
			if (i + 1 < 10) {cout << "  ";}
			else if (i + 1 < 100) {cout << " ";}
		}
		else if (numElements >= 10)
		{
			if (i + 1 < 10) {cout << " ";}
		}
		cout << i + 1 << ". ";
// Print out leading spaces
		string tempName = vectorToPrint[i];
		if (tempName.size() < maxString.size())
		{
			string::size_type tempDiffSize;
			tempDiffSize = maxString.size() - tempName.size();
			for (string::size_type iterSpaces = 0; iterSpaces < tempDiffSize; iterSpaces++)
			{
				cout << " ";
			}
		}
		cout << tempName << endl;
	}
}



void printVectorAsList (vector <int> const& vectorToPrint) // overloading for debugging
{
	int numElements = vectorToPrint.size();
	for (int i = 0; i < numElements; i++)
	{		
		cout << "  ";
		if (numElements >= 1000)
		{
			if (i + 1 < 10) {cout << "   ";}
			else if (i + 1 < 100) {cout << "  ";}
			else if (i + 1 < 1000) {cout << " ";}
		}
		else if (numElements >= 100)
		{
			if (i + 1 < 10) {cout << "  ";}
			else if (i + 1 < 100) {cout << " ";}
		}
		else if (numElements >= 10)
		{
			if (i + 1 < 10) {cout << " ";}
		}
		cout << i + 1 << ". ";
// Print out leading spaces
		cout << vectorToPrint[i] << endl;
	}
	cout << endl;
}



vector <string> collectData (string const& fileName)
{
	vector <string> inputData;
	ifstream inFile;
	string line;
	inFile.open(fileName.c_str());
	
	while (getline(inFile, line))
	{		
		if (line.empty())
		{
			continue;
		}
		else
		{
			inputData.push_back(line);
		}
		line.clear();
	}
	inFile.close();
	inFile.clear();
	
	return inputData;
}

void printProgress (string const& elementType, BigInt const& current, BigInt const& upper)
{
	cout << "  " << elementType << ": " << current << " of " << upper << "\r" << flush;
}

bool caseInsensitiveStringCompare (string const& str1, string const& str2)
{
	bool match = true;
	if (str1.size() != str2.size())
	{
		match = false;
	}
	else if (str1.size() == str2.size())
	{
		for (string::const_iterator c1 = str1.begin(), c2 = str2.begin(); c1 != str1.end(); ++c1, ++c2)
		{
//			cout << "Comparing " << *c1 << " and " << *c2 << endl;
			if (tolower(*c1) != tolower(*c2))
			{
				 match = false;
			}
		}
	}
	return match;
}

bool checkWhiteSpaceOnly (string stringToParse)
{
	bool whiteSpaceOnly = true;
	vector<string> tempVector;
	istringstream tempStream(stringToParse);
	string tempString;
	while (tempStream >> tempString)
	{
		if (tempString != "	" && tempString != " ")
		{
			whiteSpaceOnly = false;
		}
	}
	return whiteSpaceOnly;
}

float sum (vector <double> const& x)
{
	float total = 0.0;  // the sum is accumulated here
	for (int i = 0; i < (int)x.size(); i++)
	{
		total = total + x[i];  
    }
    return total;
}

bool subtractVector (vector <int> const& minuend, vector <int> const& subtrahend,
	vector <int> & difference)
{
	bool valid = true;
	
	for (int i = 0; i < (int)minuend.size(); i++)
	{
		difference.push_back(minuend[i] - subtrahend[i]);
		if ((minuend[i] - subtrahend[i]) < 0)
		{
			valid = false;
		}
	}
	
	return valid;
}