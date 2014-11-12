#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <time.h>

using namespace std;

#include "General.h"
#include "Matrix_Scan.h"
#include "Parse_Data.h"

// Store data separately from taxon/locus names
vector < vector <int> > storeData (vector <string> & inputData, int const& numLoci)
{
    vector < vector <int> > data;
    vector <string> tempStringVector;
    vector <int> tempIntVector;
    int numTaxa = (int)inputData.size() - 1;
    
    for (int i = 1; i <= numTaxa; i++) {
        tempStringVector = storeStringVector(inputData[i]);
        tempStringVector.erase(tempStringVector.begin());
        for (int j = 0; j < numLoci; j++) {
            tempIntVector.push_back(convertStringtoInt(tempStringVector[j]));
        }
        data.push_back(tempIntVector);
        tempIntVector.clear();
    }
    return data;
}

void getWeights (string const& weightFileName, vector <double> & weights, vector <string> const& names)
{
    vector <string> tempVector;
    
    int numElements = (int)names.size();
    weights = vector <double> (numElements, 1.0); // By default, apply weight of 1.0 to all
    
    if (!weightFileName.empty()) {
        tempVector = collectData(weightFileName);
        int numReweight = (int)tempVector.size();
        
        for (int i = 1; i < numReweight; i++) { // skip header line
            for (int j = 0; j < numElements; j++) {
                if (extractStringElement(tempVector[i],0) == names[j]) {
                    weights[j] = convertStringtoDouble(extractStringElement(tempVector[i],1));
                }
            }
        }
    }
}

vector <string> getTaxonNames(vector <string> & inputData)
{
    vector <string> names;
    int numTaxa = (int)inputData.size() - 1;
    
    for (int i = 1; i <= numTaxa; i++) {
        names.push_back(extractStringElement(inputData[i], 0));
    }
    
    return names;
}

// Gather up all useful nuggets; void function used so that multiple items can be returned simultaneously
void parseInputMatrix (string const& inputFileName, vector <string> & locusNames, vector <string> & taxonNames,
    vector < vector <int> > & data)
{
    vector <string> inputData;
    int numLoci = 1;
    
    inputData = collectData(inputFileName);
    locusNames = storeStringVector(inputData[0]);
    taxonNames = getTaxonNames(inputData);
    numLoci = (int)locusNames.size();
    data = storeData(inputData, numLoci);
}
