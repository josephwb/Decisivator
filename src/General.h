#ifndef _GENERAL_H_
#define _GENERAL_H_

#include "ttmath/ttmath.h"
typedef ttmath::UInt<100> BigInt; // for large matrices

// General functions

bool checkWhiteSpaceOnly (string stringToParse);
bool checkCharValue (char const& charInput, char const& charToMatch);
string getFileName();
bool checkValidFile (string);
bool checkValidOutputFile (string & outputFileName);
int checkValidIntInput (string);
bool checkValidBoolInput (string queryString);
bool checkStringValue (string stringToParse, string stringToMatch, int stringPosition);
BigInt choose (int const& n, int const& r);
BigInt factorial (int const& num);
//BigInt Stirling2ndKind (int const& n, int const& k);
string extractStringElement (string & stringToParse, int const& position);
int convertStringtoInt (string stringToConvert);
double convertStringtoDouble (string stringToConvert);
vector <string> storeStringVector (string & stringToParse);
void printVectorAsList (vector <string> const& vectorToPrint);
void printVectorAsList (vector <int> const& vectorToPrint); // overloading for debugging
vector <string> collectData (string const& fileName);
void printProgress (string const& elementType, BigInt const& current, BigInt const& upper);
bool caseInsensitiveStringCompare (string const& str1, string const& str2);
float sum (vector <double> const& x);
#endif /* _GENERAL_H_ */