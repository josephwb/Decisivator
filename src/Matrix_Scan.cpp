#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <numeric>

#include <time.h> 

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

extern bool debugging;

void checkForMissingTaxa (vector < vector <int> > const& data, vector <string> const& taxonNames)
{
    vector <string> missingData;
    for (int i = 0; i < (int)data.size(); i++) { // numTaxa
        bool match = false;
        for (int j = 0; j < (int)data[0].size(); j++) { // numPartitions
            if (data[i][j]) {
                match = true;
            }
        }
        if (!match) {
            missingData.push_back(taxonNames[i]);
        }
    }
    if (!missingData.empty()) {
        if (missingData.size() > 1) {
            cout << endl << " *** WARNING! MULTIPLE TAXA WITH NO DATA!! ***" << endl << endl;
            printVectorAsList(missingData);
        } else {
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
    for (int i = 0; i < numTaxa; i++) {
        if (data[i] == completeTaxon) {
            if (!referencePresent) {
                if (debugging) {cout << endl << "Reference taxon (i.e. sequenced for all genes)." << endl;}
            }
            referencePresent = true;
            referenceTaxa.push_back(i);
            if (debugging) {cout << " " << numReference + 1 << ". " << taxonNames[i] << "." << endl;}
            cout << " " << numReference + 1 << ". " << taxonNames[i] << "." << endl;
            numReference++;
        }
    }
    if (numReference > 0) {
        cout << numReference << " reference taxa observed." << endl;
        
// take first reference taxon found, put at the bottom; speeds up satisfying edges
        if (referenceTaxa[numReference - 1] != (int)taxonNames.size() - 1) // i.e. last taxon is NOT a reference taxon
        {
            int ref = referenceTaxa[0];
            
            if (debugging) {cout << "Moving reference taxon '" << taxonNames[ref] << "' to bottom of matrix. Don't be alarmed!" << endl;}
            
// need to fix reference taxon indexing; moving taxon i to end changes index for all taxa i+1 -> N
            if (debugging) {cout << "Raw:" << endl;
            printVectorAsList(referenceTaxa);}
            
            referenceTaxa.push_back(numTaxa);
            referenceTaxa.erase(referenceTaxa.begin()+0);
            
            for (int i = 0; i < (int)referenceTaxa.size(); i++) {
                referenceTaxa[i] = referenceTaxa[i] - 1;
            }
            
            if (debugging) {cout << "Fixed:" << endl;
            printVectorAsList(referenceTaxa);}
            
            data.push_back(data[ref]);
            data.erase(data.begin()+ref);
            
            taxonNames.push_back(taxonNames[ref]);
            taxonNames.erase(taxonNames.begin()+ref);
        }
    } else {
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
    for (int i = 0; i < numTaxa - 2; i++) {
        vector <int> tempI;
        tempI.push_back(i);
        for (int j = i + 1; j < numTaxa - 1; j++) {
            vector <int> tempJ = tempI;
            tempJ.push_back(j);
            for (int k = j + 1; k < numTaxa; k++) {
                vector <int> tempK = tempJ;
                tempK.push_back(k);
                vector <int> locations;
                
                bool tripletFound = false;
                for (int l = 0; l < numLoci; l++) {    // Loop over loci
                    if (data[i][l] == 1 && data[j][l] == 1 && data[k][l] == 1) {
                        tripletFound = true;
                        locations.push_back(l);
                    }
                }
                triplets.push_back(tempK);
                tripletLocations.push_back(locations);
                if (!tripletFound) {
                    missingTriplets.push_back(tempK);
                } else {
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
    if (tripletCounter != numPresent) {
        cout << "Matrix is NOT decisive for all possible trees." << endl;
    } else {
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
    for (int i = 0; i < numTaxa - 3; i++) {
        vector <int> tempI;
        tempI.push_back(i);
        for (int j = i + 1; j < numTaxa - 2; j++) {
            vector <int> tempJ = tempI;
            tempJ.push_back(j);
            for (int k = j + 1; k < numTaxa - 1; k++) {
                vector <int> tempK = tempJ;
                tempK.push_back(k);
                for (int l = k + 1; l < numTaxa; l++) {
                    vector <int> tempL = tempK;
                    tempL.push_back(l);
                    
                    bool quartetFound = false;
                    for (int m = 0; m < numLoci; m++) {    // Loop over loci
                        if (data[i][m] == 1 && data[j][m] == 1 && data[k][m] == 1 && data[l][m] == 1) {
                            quartetFound = true;
                            continue;
                        }
                    }
                    if (!quartetFound) {
                        missingQuartets.push_back(tempL);
                    } else {
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
    
    if (quartetCounter != numPresent) {
        cout << "Matrix is NOT decisive for all possible trees." << endl;
    } else {
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
    for (int i = 0; i < numTaxa; i++) {
        bool match = false;
        for (int j = 0; j < numReferenceTaxa; j++) {
            if (i == referenceTaxa[j]) {
                match = true;
            }
        }
        if (!match) {
            availableTaxa.push_back(i);        
        }
    }
    numTaxaToConsider = (int)availableTaxa.size();
    
    if (availableTaxa.size() == 0) {
        cout << "Matrix is complete, and so decisive for all possible trees. Lucky you." << endl;
    } else {
        if (numReferenceTaxa > 0) {
            cout << endl << "Searching for presence of all possible taxon quartets containing (arbitrary) reference taxon '" << taxonNames[referenceTaxa[0]] << "'..." << endl;    
        } else {
            cout << endl << "Searching for presence of all possible taxon quartets containing reference taxon '" << taxonNames[referenceTaxa[0]] << "'..." << endl;    
        }
    }
// Essentially just looking for triples with N = numTaxa-1 (a little better if several reference taxa exist), since reference taxon is a given
    vector <int> reference;
    reference.push_back(referenceTaxa[0]);
    
    if (debugging) {
        cout << "numTaxaToConsider = " << numTaxaToConsider << ":" << endl;
        for (int i = 0; i < numTaxaToConsider; i++) {
            cout << "  " << i + 1 << ". " << taxonNames[availableTaxa[i]] << endl;
        }
    }
    
// If there are <= 3 taxa missing data, the search-loop is not entered, and things fuck up
// Add in arbitrary reference taxa
    if (numTaxaToConsider < 3) {
        int diff = 3 - numTaxaToConsider;
        for (int i = 0; i < diff; i++) {
            numTaxaToConsider++;
            availableTaxa.push_back(referenceTaxa[i+1]);
        }
    }
    
    for (int i = 0; i < numTaxaToConsider - 2; i++) {
        vector <int> tempI = reference;
        tempI.push_back(availableTaxa[i]);
        for (int j = i + 1; j < numTaxaToConsider - 1; j++) {
            vector <int> tempJ = tempI;
            tempJ.push_back(j);
            for (int k = j + 1; k < numTaxaToConsider; k++) {
                vector <int> tempK = tempJ;
                tempK.push_back(k);
                vector <int> locations;
                
                bool quartetFoundFound = false;
                for (int l = 0; l < numLoci; l++) {    // Loop over loci
//                    if(data[i][l] == 1 && data[j][l] == 1 && data[k][l] == 1)
                    if(data[availableTaxa[i]][l] == 1 && data[availableTaxa[j]][l] == 1 && data[availableTaxa[k]][l] == 1) {
                        quartetFoundFound = true;
                        continue; // got it; move on
                    }
                }
                if (!quartetFoundFound) {
                    missingQuartets.push_back(tempK);
                } else {
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
    
    if (quartetCounter == numPresent && quartetCounter == 0) { // should be obsolete
        cout << endl << "Zero quartets considered. Happens when < 4 'available' taxa present." << endl;
    }
    
    if (quartetCounter != numPresent) {
        cout << "Matrix is NOT decisive for all possible trees." << endl;
    } else if (quartetCounter == numPresent && numPresent > 0) {
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
    
    for (int i = 0; i < numTaxa; i++) {
        int numMissing = 0;
        bool reference = false;
        for (int j = 0; j < numReferenceTaxa; j++) {    // Don't count reference taxa in missing taxa
            if (i == referenceTaxa[j]) {
                reference = true;
            }
        }
        if (!reference) {
            for (int j = 0; j < numGroups; j++) {
                for (int k = 0; k < int(missingGroups[j].size()); k++) {
                    if (i == missingGroups[j][k]) {
                        numMissing++;
                    }
                }
            }
        }
        missingTaxa.push_back(numMissing);
    }
    
    cout << endl << "Missing " << grouping << " summary:" << endl << endl;
    for (int i = 0; i < numTaxa; i++) {
        if (missingTaxa[i] > 0) {
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
    
    for (int i = 0; i < numTaxa; i++) {
        for (int j = 0; j < numLoci; j++) {
            if(data[i][j] == 1) {
                taxonGeneCount++;
            }
        }
    }
    taxonCoverage = float(taxonGeneCount)/(float(numTaxa) * float(numLoci));
}


// Eventually pass in more stuff e.g. keep track of missing clades, etc.
// multithreading implemented via openmp
double calculatePartialDecisiveness (bool const& referenceTaxonPresent, int & numTrees,
    vector < vector <int> > const& data, bool const& findAll, int const& numProcs, bool const& verbose)
{
    double partialDecisiveness = 0.0; // keep track with running mean to avoid possible overflow
    int numInternalEdges = data.size() - 3;
    int numTaxa = data.size();
    
    double foo = 0.0;
    
// should this be fixed or flexible?
    if (numTrees == 0) {
        numTrees = checkValidIntInput("Specify number of trees to simulate for calculating partial decisiveness: ");
    }
    
    if (verbose) {
        cout << endl << "Simulating " << numTrees << " trees and testing for satisfaction of internal edges..." << endl;
    }
    
// Running mean: ((((count - 1)/count) * runningMean) + (currentResult/count)); // prevent overflow
    
    #pragma omp parallel num_threads(numProcs)
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
        
//        double temp = 0.0;
        double sum = 0.0;
        
        #pragma omp for
        for (int j = 0; j < numTrees; j++) {
            if (numProcs == 1 && verbose) {printProgress("Tree", j + 1, numTrees);}
            privateTreeCount++;
            double currentDecisiveness = 0.0;
            double numEdgesSatisfied = 0.0;
            tree = fastBinaryTree(numTaxa, sibNodes, referenceTaxonPresent);
            if (debugging) {printTree(tree);}
            
// For debugging:
//        printTree(tree);
//        sort(tree.begin()+numTaxa,tree.end());
                    
            if (debugging) {cout << endl << "Edges:" << endl;}
            
            for (int i = 0; i < numInternalEdges; ++i) { // Walk through all internal edges
                getEdges(i, tree, sibNodes, referenceTaxonPresent, left, right, sib, upper);
                
// *** Scan taxon-gene matrix here ***
        // treating searchEdgePartitions here as boolean, but returns an int
                if (searchEdgePartitions(data, left, right, sib, upper, findAll, referenceTaxonPresent)) {
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
            currentDecisiveness = numEdgesSatisfied / (double)numInternalEdges;
//            privateDecisiveness = ((((privateTreeCount - 1)/privateTreeCount) * privateDecisiveness)
//                + (currentDecisiveness/privateTreeCount));
            
//            temp = privateDecisiveness;
            sum += currentDecisiveness;
        }
        
        #pragma omp critical
        {
//            cout << "Thread " << thread_id << " analyzed (" << privateTreeCount << ") trees. Yay!"
//                << " Partial decisiveness for this thread is: " << temp << endl;
            if (verbose) {cout << "Partial decisiveness for thread " << thread_id << " is: "
                << sum / privateTreeCount << endl;}
        //    foo += sum / privateTreeCount;
            foo += sum;
        }
    }
    if (verbose) {cout << endl << "Done." << endl << endl;}
    
    partialDecisiveness = (foo / numTrees);
    if (verbose) {cout << "Decisiveness is currently: " << partialDecisiveness << endl;}
    
    return(partialDecisiveness);
}


// need multithreading here!
double calculatePartialDecisivenessSinglePartition (bool const& referenceTaxonPresent, int & numTrees,
    vector < vector <int> > const& data, bool const& findAll, int const& partitionID,
    bool const& verbose, int const& numProcs)
{
    double partialDecisiveness = 0.0; // keep track with running mean to avoid possible overflow
    int numInternalEdges = data.size() - 3;
    int numTaxa = data.size();
    vector < vector <int> > partitionToQuery;
    
    double foo = 0.0;
    
    for (int i = 0; i < numTaxa; i++) { // extract partition of interest
        vector <int> temp;
        temp.push_back(data[i][partitionID]);
        partitionToQuery.push_back(temp);
    }
    
// should this be fixed or flexible?
    if (numTrees == 0) {
        numTrees = checkValidIntInput("Specify number of trees to simulate for calculating partial decisiveness: ");
    }
    
    if (verbose) {
        cout << endl << "Simulating " << numTrees << " trees and testing for satisfaction of internal edges..." << endl;
    }
    
    #pragma omp parallel num_threads(numProcs)
    {
        vector < vector <bool> > tree;
        vector < vector <int> > sibNodes;
        vector <int> left;
        vector <int> right;
        vector <int> sib;
        vector <int> upper;
        double treeCount = 0;
        double privateTreeCount = 0.0;
        
        int thread_id = omp_get_thread_num();
        
        double sum = 0.0;
        
        #pragma omp for
        for (int j = 0; j < numTrees; j++) {
            if (numProcs == 1 && verbose) {printProgress("Tree", j + 1, numTrees);}
            treeCount++;
            privateTreeCount++;
            double currentDecisiveness = 0.0;
            double numEdgesSatisfied = 0;
            tree = fastBinaryTree(numTaxa, sibNodes, referenceTaxonPresent);
            if (debugging) {printTree(tree);}
            
            if (debugging) {cout << endl << "Edges:" << endl;}
            for (int i = 0; i < numInternalEdges; ++i) { // Walk through all internal edges
                getEdges(i, tree, sibNodes, referenceTaxonPresent, left, right, sib, upper);
                
    // *** Scan taxon-gene matrix here ***
                if (searchEdgePartitions(partitionToQuery, left, right, sib, upper, findAll, referenceTaxonPresent)) {
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
            currentDecisiveness = numEdgesSatisfied / (double)numInternalEdges;
            sum += currentDecisiveness;
            //partialDecisiveness = ((((treeCount - 1)/treeCount) * partialDecisiveness) + (currentDecisiveness/treeCount));
        }
        
        #pragma omp critical
        {
//            cout << "Thread " << thread_id << " analyzed (" << privateTreeCount << ") trees. Yay!"
//                << " Partial decisiveness for this thread is: " << temp << endl;
            if (verbose) {cout << "Partial decisiveness for thread " << thread_id << " is: "
                << sum / privateTreeCount << endl;}
        //    foo += sum / privateTreeCount;
            foo += sum;
        }
    }
    if (verbose) {cout << endl << "Done." << endl << endl;}
    
    partialDecisiveness = (foo / numTrees);
    //cout << "Partial decisiveness for partition '" << partitionID << "' is: " << partialDecisiveness << endl;
    
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
    
    if (!findAll && referenceTaxonPresent) { // only consider a single reference taxon for upper; if it fails, all will
        numUpper = 1;
    }
    
// Need to loop over vectors left, right, sib, other
//    cout << endl << "Searching for presence of all possible taxon quartets..." << endl;
    for (int l = 0; l < numUpper; l++) {
        for (int i = 0; i < (int)left.size(); i++) {
            for (int j = 0; j < (int)right.size(); j++) {
                for (int k = 0; k < (int)sib.size(); k++) {
                    for (int m = 0; m < numLoci; m++) {    // Loop over loci
                        if (data[left[i]][m] == 1 && data[right[j]][m] == 1 && data[sib[k]][m] == 1 && data[upper[l]][m] == 1) {
                            if (!findAll) {
                                numSatisfied++;
                                if (debugging) {cout << "Quartet (" << left[i] << "," << right[j] << "," << sib[k] << "," << upper[l] << ") found at locus " << m << "!" << endl;}
                                return(numSatisfied); // exit: one is all that is needed
                            } else {
                                numSatisfied++;
                                m = numLoci; // exit this quartet, go on to next
                            }
                        }
                    }
                }
            }
        }
    }
    return(numSatisfied); // hmm. returning int here, but treated as boolean elsewhere
}


bool testCompleteDecisivness (vector < vector <int> > const& data, bool const& referenceTaxonPresent,
    vector <int> const& referenceTaxa, vector <string> const& taxonNames, vector < vector <int> > & missingQuartets,
    vector < vector <int> > & triplets, vector < vector <int> > & tripletLocations, vector < vector <int> > & missingTriplets)
{
    bool matrixDecisive = false;
    
    if (referenceTaxonPresent) {
        searchForQuartetsWithReference(data, missingQuartets, referenceTaxa, taxonNames);
        if (missingQuartets.size() > 0) {
            whichTaxaProblematic(missingQuartets, taxonNames, "quartets", referenceTaxa);
        } else {
            matrixDecisive = true;
        }
    } else {
        searchForAllTriplets(data, triplets, tripletLocations, missingTriplets); // test triplets first; far fewer
        if (missingTriplets.size() > 0) {
            whichTaxaProblematic(missingTriplets, taxonNames, "triplets", referenceTaxa);
        } else {
            searchForAllQuartets(data, missingQuartets);
            if (missingQuartets.size() > 0) {
                whichTaxaProblematic(missingQuartets, taxonNames, "quartets", referenceTaxa);
            }
        }
    }
    
    return (matrixDecisive);
}


vector < vector <double> > determineDecisivenessUserTree (vector < vector <int> > const& data,
    vector < vector < vector <bool> > > & userTrees, vector < vector <int> > const& treeTaxonOrdering,
    vector <string> const& taxonNames, vector <double> const& locusWeights, vector <double> const& taxonWeights,
    int const& numProcs)
{
    vector < vector <double> > result;
    vector < vector <bool> > rawTree;
    vector < vector <bool> > formattedTree;
    vector < vector <int> > sibNodes;
    vector <bool> currentClade;
    double treeCount = 0;
    int numTaxa = (int)userTrees[0][0].size();
    int numTrees = (int)userTrees.size();
    
    time_t start = time(NULL);
    time_t stop = start;
    double seconds = 0;
    
    
    
    
    // *** FIX THIS ***
    //bool findAll = true; // *** this should be an option. might just want # genes that speak to each edge (not # quartets).
    bool findAll = true;
    bool debugging = true;
    
    
    printMatrixToFile (data, taxonNames, locusWeights, taxonWeights);
    int numInternalEdges = userTrees[0][0].size() - 3; // *should* be the same for all trees in a file...
    
    if (debugging) {cout << "Here I am! numTrees = " << numTrees << ". numInternalEdges = " << numInternalEdges
        << ". numTaxa = " << numTaxa << "." << endl;}
    
    if (numTrees > 1) {
        cout << endl << "Conducting analysis on " << numTrees << " provided user-trees." << endl;
    } else {
        cout << endl << "Conducting analysis on 1 provided user-tree." << endl;
    }
    
    for (int j = 0; j < numTrees; j++) {
        treeCount++;
        vector <int> taxonOrdering;
        vector <double> decisivenessCurrentTree (numInternalEdges, 0.0);
        
        
    // *** These two should be > int for large trees! ***
        vector <unsigned long int> numSatisfied (numInternalEdges, 0);
        vector <unsigned long int> numPossible (numInternalEdges, 0);
        
/* ** Need to map taxon ordering to alignment via treeTaxonOrdering ***
options:
1. rearrange tree to reflect that ordering (easy but ugly)
2. map (pretty but harder); this will be fixed when put into OO format
*/
        rawTree = userTrees[j];
        taxonOrdering = treeTaxonOrdering[j];
        if (debugging) {
//            cout << "Raw tree:" << endl;
//            printTree(rawTree);
            
            cout << endl << "Translation:" << endl;
            printVectorAsList(taxonOrdering);
        }
        
        for (int k = 0; k < (int)rawTree.size(); k++) {
            vector <bool> clade (numTaxa, false);
            for (int l = 0; l < numTaxa; l++) {
                if (rawTree[k][l]) {
                    clade[taxonOrdering[l]] = true;
                }
            }
            formattedTree.push_back(clade);
            clade.clear();
        }
        
//        cout << "Formatted tree:" << endl;
//        printTree(formattedTree);
        
        
        
    // *** YIKES! This takes a long time (64 of a total 74 second analysis). ***
        sibNodes = getSibNodes(formattedTree, numProcs); // get sibling node relationships; expected downstream
        
//        stop = time(NULL);
//        seconds = difftime(stop, start);
//        cout << endl << "getSibNodes took " << seconds << " seconds." << endl;
//        start = stop;
        
        if (debugging) {cout << endl << "Edges:" << endl;}
        
        
        
    // multi-thread this shit!!!
        #pragma omp parallel num_threads(numProcs)
        {
            vector <int> left;
            vector <int> right;
            vector <int> sib;
            vector <int> upper;
            #pragma omp for
            for (int i = 0; i < numInternalEdges; ++i) { // this is the important bit; *** NEED TO INCORPORATE TAXON CODING ***
                currentClade = formattedTree[numTaxa + i];
                double currentDecisiveness = 0.0;
                unsigned long int numEdgesSatisfied = 0;
            
                getEdges(i, formattedTree, sibNodes, 0, left, right, sib, upper);
                unsigned long int numPossibleQuartets = (int)left.size() * (int)right.size() * (int)sib.size() * (int)upper.size();
            
    // *** Scan taxon-gene matrix here ***
                numEdgesSatisfied = searchEdgePartitions(data, left, right, sib, upper, findAll, 0);
             
                if (debugging) {
                    cout << "Considering " << numPossibleQuartets << " ways to satisfy current node " << i + numTaxa << ":" << endl;
                    printClade(currentClade);
                    cout << endl;
                    cout << numEdgesSatisfied << " instances satisfied of a possible " << numPossibleQuartets
                    << " for node " << i + numTaxa << " (= " << ((double)numEdgesSatisfied/(double)numPossibleQuartets)*100 << "%)." << endl;
                }
                currentDecisiveness = (double)numEdgesSatisfied / (double)numPossibleQuartets;
                decisivenessCurrentTree[i] = currentDecisiveness;
                numSatisfied[i] = numEdgesSatisfied;
                numPossible[i] = numPossibleQuartets;
             
    // empty vectors for next edge
                left.clear();
                right.clear();
                sib.clear();
                upper.clear();
            }
        }
        
        printBipartitionTable(formattedTree, decisivenessCurrentTree, numSatisfied, numPossible, numTaxa, numTrees, j);
        result.push_back(decisivenessCurrentTree);
        
        if (numTrees > 1) {
            cout << "Average decisiveness for tree " << j + 1 << " is: " << sum(decisivenessCurrentTree)/(double)decisivenessCurrentTree.size() << endl;
        } else {
            cout << "Average decisiveness for this tree is: " << sum(decisivenessCurrentTree)/(double)decisivenessCurrentTree.size() << endl;
        }
// empty vectors for next tree
        rawTree.clear();
        formattedTree.clear();
        sibNodes.clear();
        decisivenessCurrentTree.clear();
        numSatisfied.clear();
        numPossible.clear();
    //    partialDecisiveness = ((((treeCount - 1)/treeCount) * partialDecisiveness) + (currentDecisiveness/treeCount));
    }
    
    stop = time(NULL);
    seconds = difftime(stop, start);
    
    if (numTrees > 1) {
        cout << endl << "Completed analysis of provided trees in " << seconds << " seconds. " << endl;
    } else {
        cout << endl << "Completed analysis of provided tree in " << seconds << " seconds. " << endl;
    }
    cout << "Bipartition decisiveness scores have been printed to file 'Decisivator.log'." << endl;
    return(result);
}




// PAUP-style
// Hmm. Probably don't want this to print to screen (except maybe when debugging).
// Should tailor the name of the file.
void printBipartitionTable (vector < vector <bool> > const& tree, vector <double> const& decisiveness,
    vector <unsigned long int> const& numSatisfied, vector <unsigned long int> const& numPossible, int const& numTaxa,
    int const& numTrees, int const& treeNumber)
{
    ofstream log;
    log.open("Decisivator.log",ios::app);
    
    if (numTrees > 1) {
        if (debugging) {cout << endl << "Bipartition Decisiveness Scores for tree #" << treeNumber + 1 << ":" << endl << endl;}
        log << endl << "Bipartition Decisiveness Scores for tree #" << treeNumber + 1 << ":" << endl << endl;
    } else {
        if (debugging) {cout << endl << endl << "Bipartition Decisiveness Scores:" << endl << endl;}
        log << endl << endl << "Bipartition Decisiveness Scores:" << endl << endl;
    }
    

    // TODO: there is at least 1 one-off error here; just involves index labelling at top of table

    int counter = 0;
    int divisor = 0;
    if (numTaxa > 999) {
        counter = 1;
        divisor = 1000;
        for (int i = 1; i <= numTaxa; i++) {
            if (i % divisor == 0) {
                if (debugging) {cout << counter;}
                if (counter == divisor) {
                    counter = 0;
                }
                log << counter;
                counter ++;
            } else {
                if (debugging) {cout << " ";}
                log << " ";
            }
        }
        if (debugging) {cout << endl;}
        log << endl;
    }
    if (numTaxa > 99) {
        counter = 1;
        divisor = 100;
        for (int i = 1; i <= numTaxa; i++) {
            if (i % divisor == 0) {
            	if (counter == divisor) {
            	    counter = 0;
            	}
                log << counter;
                counter ++;
            } else {
                if (debugging) {cout << " ";}
                log << " ";
            }
        }
        if (debugging) {cout << endl;}
        log << endl;
    }
    if (numTaxa > 9) {
        counter = 1;
        divisor = 10;
        for (int i = 1; i <= numTaxa; i++) {
            if (i % divisor == 0) {
                if (debugging) {cout << counter;}
                if (counter == divisor) {
                	counter = 0;
                }
                log << counter;
                counter ++;
            } else {
                if (debugging) {cout << " ";}
                log << " ";
            }
        }
        if (debugging) {cout << endl;}
        log << endl;
    }
    counter = 1;
    for (int i = 1; i <= numTaxa; i++) {
        if (counter == 10) {
            counter = 0;
            log << counter;
            if (debugging) {cout << counter;}
            counter ++;
        } else {
            if (debugging) {cout << counter;}
            log << counter;
            counter ++;
        }
    }
    if (debugging) {cout << "    Freq.    Poss.     %" << endl;}
    log << "    Freq.    Poss.     %" << endl;
    
    for (int i = 1; i <= numTaxa; i++) {
        if (debugging) {cout << "-";}
        log << "-";
    }
    if (debugging) {cout << "---------------------------" << endl;}
    log << "---------------------------" << endl;
    
// Last 'node' defines entire tree, penultimate node defines root (all other nodes already taken care of)
    counter = 0;
    for (int i = numTaxa; i < (int)tree.size() - 2; i++) {
        vector <bool> currentClade = tree[i];
        for (int j = 0; j < numTaxa; j++) {
            if (currentClade[j]) {
                if (debugging) {cout << "*";}
                log << "*";
            } else {
                if (debugging) {cout << ".";}
                log << ".";
            }
        }
        if (debugging) {cout << "    " << numSatisfied[counter] << "    " << numPossible[counter] << "    " << decisiveness[counter]*100 << endl;}
        log << "    " << numSatisfied[counter] << "    " << numPossible[counter] << "    " << decisiveness[counter]*100 << endl;
        counter++;
    }
    if (debugging) {cout << endl;}
    log << endl;
    log.close();
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
    for (int i = 0; i < numTaxa; i++) {
        bool match = false;
        for (int j = 0; j < numReferenceTaxa; j++) {
            if (i == referenceTaxa[j]) {
                match = true;
            }
        }
        if (!match) {
            availableTaxa.push_back(i);        
        }
    }
    numTaxaToConsider = (int)availableTaxa.size();
    
    if (availableTaxa.size() == 0) {
        cout << "Matrix is complete, and so decisive for all possible trees. Lucky you." << endl;
    } else {
        if (numReferenceTaxa > 0) {
            cout << endl << "Searching for presence of all possible taxon quartets containing (arbitrary) reference taxon '" << taxonNames[referenceTaxa[0]] << "'..." << endl;    
        } else {
            cout << endl << "Searching for presence of all possible taxon quartets containing reference taxon '" << taxonNames[referenceTaxa[0]] << "'..." << endl;    
        }
    }
    
    if (debugging) {
        cout << "numTaxaToConsider = " << numTaxaToConsider << ":" << endl;
        for (int i = 0; i < numTaxaToConsider; i++) {
            cout << "  " << i + 1 << ". " << taxonNames[availableTaxa[i]] << endl;
        }
    }
    
// If there are <= 3 taxa missing data, the search-loop is not entered, and things fuck up
// Add in arbitrary reference taxa
    if (numTaxaToConsider < 3) {
        int diff = 3 - numTaxaToConsider;
        for (int i = 0; i < diff; i++) {
            numTaxaToConsider++;
            availableTaxa.push_back(referenceTaxa[i+1]);
        }
    }
    
    for (int i = 0; i < numTaxaToConsider - 2; i++) {
        for (int j = i + 1; j < numTaxaToConsider - 1; j++) {
            for (int k = j + 1; k < numTaxaToConsider; k++) {
                bool quartetFoundFound = false;
                for (int l = 0; l < numLoci; l++) {    // Loop over loci
                    if(data[availableTaxa[i]][l] == 1 && data[availableTaxa[j]][l] == 1 && data[availableTaxa[k]][l] == 1) {
                        quartetFoundFound = true;
                        continue; // got it; move on
                    }
                }
                if (!quartetFoundFound) {
                    missingQuartetsByTaxa[i]++;
                    missingQuartetsByTaxa[j]++;
                    missingQuartetsByTaxa[k]++;
                } else {
                    numPresent++;
                }
                quartetCounter++;
                if ((quartetCounter & CHUNKSIZE) == CHUNKSIZE) {
                    printProgress("Quartet", quartetCounter, numPossibleQuartets);
                }
            }
        }
    }
    
    if (quartetCounter == numPresent && quartetCounter == 0) { // should be obsolete
        cout << endl << "Zero quartets considered. Happens when < 4 'available' taxa present." << endl;
    }
    if (quartetCounter != numPresent) {
        cout << "Matrix is NOT decisive for all possible trees." << endl;
    } else if (quartetCounter == numPresent && numPresent > 0) {
        cout << " Woo-hoo! All possible taxon quartets observed. Matrix IS decisive for all possible trees!" << endl;
    }
}

void whichTaxaProblematicBIG (vector <int> const& missingQuartetsByTaxa, vector <string> const& taxonNames,
    string const& grouping, vector <int> const& referenceTaxa)
{
    vector <int> missingTaxa;
    int numTaxa = (int)taxonNames.size();
    int currentMax = 0;
    
    for (int i = 0; i < numTaxa; i++) {
        if (missingQuartetsByTaxa[i] > currentMax) {
            missingTaxa.clear();
            currentMax = missingQuartetsByTaxa[i];
            missingTaxa.push_back(i);
        } else if (missingQuartetsByTaxa[i] == currentMax) {
            missingTaxa.push_back(i);
        }
    }
    
    if (currentMax > 0) {
        cout << endl << "Missing " << grouping << " summary:" << endl << endl;
        for (int i = 0; i < int(missingTaxa.size()); i++) {
            cout << " Taxon '" << taxonNames[missingTaxa[i]] << "' missing from (" << currentMax << ") taxon " << grouping << "." << endl;
        }
    }
}
