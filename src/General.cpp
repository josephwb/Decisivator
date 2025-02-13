#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <numeric>

using namespace std;

#include "General.h"

extern bool debugging;


// Perform case-insenstive char match test
bool checkCharValue (char const& charInput, char const& charToMatch) {
    if (toupper(charInput) == toupper(charToMatch)) {
        return true;
    } else {
        return false;
    }
}


// split string by whitespace
vector <string> tokenize (string const& input) {
    vector <string> tokens;
    string temp;
    istringstream str(input);
    while (str >> temp) {
        tokens.push_back(temp);
    }
    return tokens;
}


string getFileName () {
    string fileName;
    bool validFileName = false;
    string tempFileName;
    cout << endl << "What is the name of the input file?" << endl;
    
    while (!validFileName) {
        cin >> tempFileName;
        validFileName = checkValidFile(tempFileName);
        if (cin.fail()) {
            cin.clear(); 
            cin.ignore(200, '\n');
            cout << "Invalid input! Try again." << endl << endl;
            cout << "What is the name of the input file?" << endl;
        } else {
            fileName = tempFileName;
            cin.clear();
            cin.ignore(200, '\n');
        }
    }
    return fileName;
}


bool checkValidFile (string fileName) {
    ifstream inFile;
    bool validFile = true;
    
    inFile.open(fileName.c_str());
    if (inFile.fail()) {
        ofstream errorReport("Error.Decisivator.txt");
        errorReport << "Decisivator analysis failed." << endl << "Error: unable to open file '";
        errorReport << fileName << "'" << endl;
        errorReport.close();

        cerr << endl << "Decisivator analysis failed." << endl << "Error: unable to open file '";
        cerr << fileName << "'" <<  endl;
        exit(1);
    } else {
        cout << "Successfully opened file '" << fileName << "'." <<  endl;
        inFile.close();
        inFile.clear();
    }
    return(validFile);
}


bool checkValidOutputFile (string & outputFileName) {
    bool testOutBool = true;
    bool fileNameAcceptable = false;
    bool keepFileName = false;
    
// First, check if file already exists, so overwriting can be prevented
    fstream testIn;
    while (!fileNameAcceptable) {
        testIn.open(outputFileName.c_str());
        if (!testIn) {
            testIn.close();
            fileNameAcceptable = true;
        } else {
            testIn.close();
            
            cout << endl << "File '" << outputFileName << "' exists!  ";
            
            keepFileName = checkValidBoolInput("Change name (0) or overwrite (1)? ");
            
            if (!keepFileName) {
                cout << "Enter new output file name: ";
                cin >> outputFileName;
            } else {
                cout << "Overwriting existing file '" << outputFileName << "'." << endl;
                fileNameAcceptable = true;
            }
        }
    }
    
    ofstream outFile;
    outFile.open(outputFileName.c_str());
    
    if (outFile.fail()) {
        ofstream errorReport("Error.BEASTifier.txt");
        errorReport << "BEASTifier analysis failed." << endl << "Error: unable to open file '";
        errorReport << outputFileName << "'" << endl;
        errorReport.close();

        cerr << endl << "BEASTifier analysis failed." << endl << "Error: unable to open file '";
        cerr << outputFileName << "'" <<  endl;
        testOutBool = false;
        exit(1);
    } else {
        outFile.close();
        outFile.clear();
    }
    return testOutBool;
}


int checkValidIntInput (string queryString) {
    int userInput = 0;
    bool validInput = false;
    cout << endl << queryString;
    validInput = false;
    while (!validInput) {
        cin >> userInput;
        if (cin.fail()) {
            cin.clear(); 
            cin.ignore(200, '\n');
            cout << endl << "Invalid input! Must be an integer. Try again." << endl;
            cout << queryString;
        } else {
            validInput = true;
            cin.clear();             // Get rid of any extra characters accidentally entered
            cin.ignore(200, '\n');
        }
    }
    return userInput;
}


bool checkValidBoolInput (string queryString) {
    int userInput = 0;
    bool validInput = false;
    cout << queryString;
    validInput = false;
    while (!validInput) {
        cin >> userInput;
        if (cin.fail()) {
            cin.clear(); 
            cin.ignore(200, '\n');
            cout << endl << "Invalid input! Must be an integer. Try again." << endl;
            cout << queryString;
        } else if (userInput < 0 || userInput > 1) {
            cin.clear(); 
            cin.ignore(200, '\n');
            cout << endl << "Invalid input! Boolean value - must be 0 or 1. Try again." << endl;
            cout << queryString;
        } else {
            validInput = true;
            cin.clear();             // Get rid of any extra characters accidentally entered
            cin.ignore(200, '\n');
        }
    }
    return userInput;
}


// stringToMatch should be all UPPERCASE
bool checkStringValue (string & inString, string const& stringToMatch) {
    bool match = false;
    transform(inString.begin(), inString.end(), inString.begin(), ::toupper);
    if (inString == stringToMatch) {
        return true;
    }
    return match;
}


bool checkStringValue (string stringToParse, string stringToMatch, int stringPosition) {
// Performs case-insenstive string match test
    string testString = extractStringElement(stringToParse, stringPosition);
    if (testString.size() != stringToMatch.size()) {
        return false;
    }
    if (testString == stringToMatch) {
        return true;
    }
    
    for (size_t i = 0; i < testString.size(); ++i) {
        if (testString[i] == stringToMatch[i] || testString[i]+32 == stringToMatch[i] || testString[i] == stringToMatch[i]+32) {
            continue;
        } else {
            return false;
        }
    }
    return true;
}


unsigned long choose (int const& n, int const& r) {
    unsigned long result = 0;
    if (r == 0) {return (1);}
    
    unsigned long numerator = n;
    for (int i = n - 1; i > n - r; i --) {
        numerator *= i;
    }
    unsigned long denominator = factorial(r);
    result = numerator / denominator;
    return (result);
}


unsigned long factorial (int const& num) {
    if (num <= 1) return (1);
    return factorial(num - 1) * num; // recursive call
}


// Given some existing string, extract (copy) the ith element and store as new string
string extractStringElement (string & stringToParse, int const& position) {
    vector <string> tempVector = storeStringVector(stringToParse);
    string returnString = tempVector[position];
    
    return returnString;
}


// remove trailing filename extension
string removeStringSuffix (string stringToParse, char suffixToRemove) {
    string temp;
    vector <char> tempVector;
    int charCounter = 0;
    int suffixStart = 0;
    bool suffixEncountered = false;
    for (string::iterator charIter = stringToParse.begin(); charIter < stringToParse.end(); charIter++) {
        tempVector.push_back(*charIter);
        if (*charIter == suffixToRemove) {
            suffixStart = charCounter;
            suffixEncountered = true;
        }
        charCounter++;
    }
    if (suffixEncountered) {
        for (int charIter = 0; charIter < suffixStart; charIter++) {
            temp += tempVector[charIter];
        }
    }
    if (!suffixEncountered) {
        temp = stringToParse;
    }
    return temp;
}


int convertStringtoInt (string stringToConvert) {
    int tempInt = 0;
    istringstream tempStream(stringToConvert);
    tempStream >> tempInt;
    
    return tempInt;
}


string convertIntToString (int & intToConvert) {
    string tempString;
    stringstream tempStream;
    tempStream << intToConvert;
    tempString = tempStream.str();
    
    return tempString;
}


double convertStringtoDouble (string stringToConvert) {
    double tempDouble = 0;
    istringstream tempStream(stringToConvert);
    tempStream >> tempDouble;
    
    return tempDouble;
}


// make sure 0.0 and 1.0 values are output with decimals
string outputDoublePrecision (double & val) {
    string tempString;
    stringstream tempStream;
    tempStream << val;
    tempString = tempStream.str();
    if (tempString.size() == 1) {
        tempString += ".0";
    }

    return tempString;
}


// Store elements of a string in a vector
vector <string> storeStringVector (string & stringToParse) {
    vector <string> tempVector;
    istringstream tempStream(stringToParse);
    string tempString;
    
    while (tempStream >> tempString) {
        tempVector.push_back(tempString);
    }
    
    return tempVector;
}


// Print purdy lists, aligning element names
void printVectorAsList (vector <string> const& vectorToPrint) {
    string tempString = vectorToPrint[0];
    int numElements = (int)vectorToPrint.size();
    string maxString;
    int longestName = 0;
    
    for (int i = 0; i < numElements; i++) {        
        string currentString = vectorToPrint[i];
        if (currentString.size() > maxString.size()) {
            maxString = currentString;
        }
    }
    longestName = (int)maxString.size();
    for (int i = 0; i < numElements; i++) {
        cout << "  ";
        if (numElements >= 1000) {
            if (i + 1 < 10) {cout << "   ";}
            else if (i + 1 < 100) {cout << "  ";}
            else if (i + 1 < 1000) {cout << " ";}
        } else if (numElements >= 100) {
            if (i + 1 < 10) {cout << "  ";}
            else if (i + 1 < 100) {cout << " ";}
        } else if (numElements >= 10) {
            if (i + 1 < 10) {cout << " ";}
        }
        cout << i + 1 << ". ";
// Print out leading spaces
        string tempName = vectorToPrint[i];
        if ((int)tempName.size() < longestName) {
            string::size_type tempDiffSize;
            tempDiffSize = longestName - tempName.size();
            for (string::size_type iterSpaces = 0; iterSpaces < tempDiffSize; iterSpaces++) {
                cout << " ";
            }
        }
        cout << tempName << endl;
    }
}


void printVectorAsList (vector <string> const& vectorOneToPrint, vector <int> const& vectorTwoToPrint,
    vector <double> const& vectorThreeToPrint, string const& columnOneName, string const& columnTwoName,
    string const& columnThreeName, string const& columnFourName)
{
    cout.precision(3);
    cout.setf(ios::fixed,ios::floatfield);
    
    string tempString = vectorOneToPrint[0];
    int numElements = (int)vectorOneToPrint.size();
    string maxString;
    int longestName = 0;
    
    for (int i = 0; i < numElements; i++) {
        string currentString = vectorOneToPrint[i];
        if (currentString.size() > maxString.size()) {
            maxString = currentString;
        }
    }
    longestName = (int)maxString.size();
    for (int i = -1; i < numElements; i++) {
        if (i >= 0) {
            cout << "  ";
            if (numElements >= 1000) {
                if (i + 1 < 10) {cout << "   ";}
                else if (i + 1 < 100) {cout << "  ";}
                else if (i + 1 < 1000) {cout << " ";}
            } else if (numElements >= 100) {
                if (i + 1 < 10) {cout << "  ";}
                else if (i + 1 < 100) {cout << " ";}
            } else if (numElements >= 10) {
                if (i + 1 < 10) {cout << " ";}
            }
            cout << i + 1 << ". ";
    // Print out leading spaces
            string tempName = vectorOneToPrint[i];
            if ((int)tempName.size() < longestName) {
                string::size_type tempDiffSize;
                tempDiffSize = longestName - tempName.size();
                for (string::size_type iterSpaces = 0; iterSpaces < tempDiffSize; iterSpaces++) {
                    cout << " ";
                }
            }
            cout << tempName << "   " << vectorTwoToPrint[i] << "   " << vectorThreeToPrint[i] << endl;
        } else {
            cout << " " << columnOneName;
            string::size_type tempDiffSize;
            tempDiffSize = longestName - columnTwoName.size();
            for (string::size_type iterSpaces = 0; iterSpaces < tempDiffSize; iterSpaces++) {
                cout << " ";
            }
            cout << columnTwoName << "   " << columnThreeName << "  " << columnFourName << endl;
        }
    }
}


void printVectorAsList (vector <int> const& vectorToPrint) { // overloading for debugging
    int numElements = (int)vectorToPrint.size();
    for (int i = 0; i < numElements; i++) {
        cout << "  ";
        if (numElements >= 1000) {
            if (i + 1 < 10) {cout << "   ";}
            else if (i + 1 < 100) {cout << "  ";}
            else if (i + 1 < 1000) {cout << " ";}
        } else if (numElements >= 100) {
            if (i + 1 < 10) {cout << "  ";}
            else if (i + 1 < 100) {cout << " ";}
        } else if (numElements >= 10) {
            if (i + 1 < 10) {cout << " ";}
        }
        cout << i + 1 << ". ";
// Print out leading spaces
        cout << vectorToPrint[i] << endl;
    }
    cout << endl;
}


void printVectorAsList (vector <double> const& vectorToPrint) { // overloading for debugging
    int numElements = (int)vectorToPrint.size();
    for (int i = 0; i < numElements; i++) {
        cout << "  ";
        if (numElements >= 1000) {
            if (i + 1 < 10) {cout << "   ";}
            else if (i + 1 < 100) {cout << "  ";}
            else if (i + 1 < 1000) {cout << " ";}
        } else if (numElements >= 100) {
            if (i + 1 < 10) {cout << "  ";}
            else if (i + 1 < 100) {cout << " ";}
        } else if (numElements >= 10) {
            if (i + 1 < 10) {cout << " ";}
        }
        cout << i + 1 << ". ";
// Print out leading spaces
        cout << vectorToPrint[i] << endl;
    }
    cout << endl;
}


void printVectorAsList (vector <double> const& vectorOneToPrint,
    vector <double> const& vectorTwoToPrint, string const& columnOneName,
    string const& columnTwoName, string const& columnThreeName) { // overloading for debugging
    int numElements = (int)vectorOneToPrint.size();
    for (int i = -1; i < numElements; i++) {
        if (i >= 0) {
            cout << "  ";
            if (numElements >= 1000) {
                if (i + 1 < 10) {cout << "   ";}
                else if (i + 1 < 100) {cout << "  ";}
                else if (i + 1 < 1000) {cout << " ";}
            } else if (numElements >= 100) {
                if (i + 1 < 10) {cout << "  ";}
                else if (i + 1 < 100) {cout << " ";}
            } else if (numElements >= 10) {
                if (i + 1 < 10) {cout << " ";}
            }
            cout << i + 1 << ".  ";
            cout << vectorOneToPrint[i] << "  " << vectorTwoToPrint[i] << endl;
        } else {
            if (numElements >= 1000) {
                if (i + 1 < 10) {cout << "  ";}
                else if (i + 1 < 100) {cout << " ";}
            } else if (numElements >= 100) {
                if (i + 1 < 10) {cout << " ";}
            }
            if (numElements > 10) {
                cout << columnOneName << "  " << columnTwoName << "  " << columnThreeName << endl;
            } else {
                cout << columnOneName << " " << columnTwoName << "  " << columnThreeName << endl;
            }
        }
    }
}


vector <string> collectData (string const& fileName) {
    vector <string> inputData;
    ifstream inFile;
    string line;
    inFile.open(fileName.c_str());
    
    while (getline(inFile, line)) {
        if (line.empty()) {
            continue;
        } else {
            inputData.push_back(line);
        }
        line.clear();
    }
    inFile.close();
    inFile.clear();
    
    return inputData;
}


void printProgress (string const& elementType, unsigned long const& current, unsigned long const& upper) {
    cout << "  " << elementType << ": " << current << " of " << upper << "\r" << flush;
}


bool caseInsensitiveStringCompare (string const& str1, string const& str2) {
    bool match = true;
    if (str1.size() != str2.size()) {
        match = false;
    } else if (str1.size() == str2.size()) {
        for (string::const_iterator c1 = str1.begin(), c2 = str2.begin(); c1 != str1.end(); ++c1, ++c2) {
//            cout << "Comparing " << *c1 << " and " << *c2 << endl;
            if (tolower(*c1) != tolower(*c2)) {
                 match = false;
            }
        }
    }
    return match;
}


bool checkWhiteSpaceOnly (string stringToParse) {
    bool whiteSpaceOnly = true;
    vector<string> tempVector;
    istringstream tempStream(stringToParse);
    string tempString;
    while (tempStream >> tempString) {
        if (tempString != "    " && tempString != " ") {
            whiteSpaceOnly = false;
        }
    }
    return whiteSpaceOnly;
}


double sum (vector <double> const& x) {
    double total = 0.0;  // the sum is accumulated here *** <- this is available natively ***
    
    total += accumulate(x.begin(), x.end(), 0.0);
    /*
    for (int i = 0; i < (int)x.size(); i++) {
        total = total + x[i];  
    }
    */
    return total;
}


// *** DEPRECATED CODE *** //

/*
unsigned long Stirling2ndKind (int const& n, int const& k)
{
    unsigned long S2K = 0;
    
    if (n == k) {return (1);}
    for (int i = 0; i < n; i++) {
        S2K += pow(-1,i) * choose(k,i) * pow(k-i,n);
    }
    S2K = S2K / factorial(k);
    return(S2K);
}
*/
