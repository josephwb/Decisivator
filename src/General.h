#ifndef _GENERAL_H_
#define _GENERAL_H_

#define CHUNKSIZE 0x0FFF

// #include "ttmath/ttmath.h"
// typedef ttmath::UInt<100> unsigned long ; // for large matrices

// General functions

bool checkWhiteSpaceOnly (string stringToParse);
bool checkCharValue (char const& charInput, char const& charToMatch);
string getFileName();
bool checkValidFile (string);
bool checkValidOutputFile (string & outputFileName);
int checkValidIntInput (string);
bool checkValidBoolInput (string queryString);
bool checkStringValue (string stringToParse, string stringToMatch, int stringPosition);
unsigned long choose (int const& n, int const& r);
unsigned long factorial (int const& num);
//unsigned long Stirling2ndKind (int const& n, int const& k);
string extractStringElement (string & stringToParse, int const& position);
int convertStringtoInt (string stringToConvert);
string convertIntToString (int & intToConvert);
double convertStringtoDouble (string stringToConvert);
vector <string> storeStringVector (string & stringToParse);
void printVectorAsList (vector <string> const& vectorToPrint);
void printVectorAsList (vector <int> const& vectorToPrint); // overloading for debugging
void printVectorAsList (vector <double> const& vectorToPrint); // overloading for debugging
void printVectorAsList (vector <double> const& vectorOneToPrint,
	vector <double> const& vectorTwoToPrint, string const& columnOneName,
	string const& columnTwoName, string const& columnThreeName); // overloading for debugging
void printVectorAsList (vector <string> const& vectorOneToPrint, vector <int> const& vectorTwoToPrint,
	vector <double> const& vectorThreeToPrint, string const& columnOneName, string const& columnTwoName,
	string const& columnThreeName, string const& columnFourName);
vector <string> collectData (string const& fileName);
void printProgress (string const& elementType, unsigned long const& current, unsigned long const& upper);
bool caseInsensitiveStringCompare (string const& str1, string const& str2);
double sum (vector <double> const& x);
#endif /* _GENERAL_H_ */