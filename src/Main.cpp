/*  
   Decisivator
   Copyright (c) 2014, PuRGe C/O JWB
   PuRGe is from:
   Department of Biological Sciences
   University of Idaho, Moscow, Idaho, USA
   JWB is currently at:
   Department of Ecology & Evolutionary Biology
   University of Michigan, Ann Arbor, Michigan, USA
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   Any feedback is very welcome.
   email: josephwb@umich.edu
*/


/*

A program to analyze a matrix of taxon-locus cells for the purpose of determining idealized phylogenetic
'decisiveness' sensu Michael Sanderson and Mike Steel.

Presence of a locus is denoted by '1', absence by '0'. User may pass in a 0-1 matrix or (more likely)
a Nexus file containing CHARSET declarations (in the latter case, the 0-1 matrix is determined
automatically). 

Searches for all taxon triplets (necessary) and quartets (sufficient) conditions (or, if a reference
taxon exists, that is, one sequenced for all loci, all quads that contain that reference taxon), i.e. 
determines decisiveness of matrix as a whole for ALL possible trees.

At the moment, deals only with a single matrix. User interacts, adding taxon-loci, merging taxa 
(forming chimeric sequences) and/or deleting taxa and loci that may be causing problems (although
deleting loci will never help, except to form reference taxa). In this way, one may (eventually,
and pretty inefficiently) attain a matrix that is phylogenetically decisive.

To determine partial decisiveness (branch-wise), the following will occur:

1. get a topology (edgelengths not required) - DONE
2. traverse tree, edge by edge - DONE
3. determine if two taxa on EACH side of the edge are sequenced for a particular locus
   
'Partial decisivess' is of interest for both the focal tree and 'all' trees. For the latter, it is not
possible to analyze every possible configuration. As an approximation, simulate multiple random trees
of the same taxon size. An important consideration may be to constrain parts of the generated trees
to display non-controversial edges.

Ultimate goal: how to add taxon-loci to the matrix to most optimally (i.e. most frugally or with fixed
effort) improve phylogenetic decisiveness.

The code below will, after final structure decisions (i.e. 'structurally decisive'! Ha!), be
reorganized into an object-oriented form.

Oh, and assumes unix line breaks.

*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_num_procs() 1
#endif

using namespace std;

#include "General.h"
#include "Trees_Edges.h"
#include "Manipulate_Matrix.h"
#include "Matrix_Scan.h"
#include "Parse_Data.h"
#include "User_Interface.h"
#include "Parse_Nexus.h"
#include "User_Tree.h"
#include "GAoptimize.h"


// version information
double version = 0.58;
string month = "March";
int year = 2014;

bool debugging = false; // should extra comments be printed to stdout?

int main(int argc, char *argv[]) {
    srand ((unsigned int)time(NULL) );     // initialize random seed
    string matrixFileName;
    string locusWeightFileName;
    string taxonWeightFileName;
    string nexusFileName;
    int numProcs = omp_get_num_procs();
    
// Trees
    string treeFileName;
    vector <string> rawTrees;
    vector < vector < vector <bool> > > userTrees; // can store multiple trees
    int numUserTrees = 0;
    int burnin = 0;
    int thinning = 1;
    vector < vector <int> > treeTaxonOrdering; // maps to alignment ordering
    vector < vector <double> > userTreeDecisiveness;
    
// Original matrix
    vector <string> locusNames;
    vector <string> taxonNames;
    vector < vector <int> > data;
    vector < vector <string> > taxaAlignment;
    vector < vector <int> > includedLocusRanges;
    vector <int> translationTable;
    vector <double> locusWeights; // default of 1.0; user can change
    vector <double> taxonWeights; // default of 1.0; user can change
    string dataType;
    
    vector < vector <int> > triplets;
    vector < vector <int> > tripletLocations;
    vector < vector <int> > missingTriplets;
    vector < vector <int> > missingQuartets;
    bool referenceTaxonPresent = false;
    vector <int> referenceTaxa;    // Possible that there are several
    vector <int> missingQuartetsByTaxa; // for very large matrices
    double taxonCoverage = 0.0;
    bool matrixDecisive = false;
    double branchwiseDecisiveness = 0.0;
    double treewiseDecisiveness = 0.0;
    int numChar = 0;
    double GADecisiveness = 0.0;
    
    printProgramInfo();
    
    processCommandLineArguments(argc, argv, matrixFileName, nexusFileName, locusWeightFileName,
        taxonWeightFileName, treeFileName, burnin, thinning, numProcs);
    
    cout << numProcs << " processors available for analyisis." << endl << endl;
    
// Read in data, store
    if (!nexusFileName.empty()) {
        parseNexus(nexusFileName, data, taxonNames, locusNames, numChar, taxaAlignment,
            includedLocusRanges, dataType);
    } else if (!matrixFileName.empty()) {
        parseInputMatrix(matrixFileName, locusNames, taxonNames, data);
    } else {
        cout << "Ergh. Something fucked up... Shit." << endl;
        exit(1);
    }
    if (!matrixFileName.empty()) {
        nexusFileName = matrixFileName;
    }

// *** If reference taxon IS present, put at bottom of taxon-locus matrix (as outgroup taxon is always last in tree)
//     Will allow faster satisfaction (but not rejection) of internal edges.
    referenceTaxonPresent = searchForReferenceTaxon(data, referenceTaxa, taxonNames);
    
    if (!treeFileName.empty()) {
        getUserTrees (treeFileName, rawTrees, userTrees, taxonNames, translationTable, burnin,
            thinning, treeTaxonOrdering);
        numUserTrees = (int)userTrees.size();
        if (debugging) {
            cout << "userTrees has size " << numUserTrees << "; treeTaxonOrdering has size "
                << treeTaxonOrdering.size() << endl;
            cout << "Collected user trees:" << endl;
            printVectorAsList(rawTrees);}
    }
    
// Get weights (if present)
// Probably change the format of this to a single vector < vector <double> >
    getWeights(locusWeightFileName, locusWeights, locusNames);
    getWeights(taxonWeightFileName, taxonWeights, taxonNames);
    getCoverage(data, taxonCoverage);
    
    int numPartitions = (int)data[0].size();
// Data structures for revised matrix; these will be deprecated when reorganized in object-oriented form
    vector < vector <int> > revisedData = data;
    vector < vector <int> > revisedTriplets;
    vector < vector <int> > revisedTripletLocations;
    vector < vector <int> > revisedMissingTriplets;
    vector < vector <int> > revisedMissingQuartets;
    vector <double> revisedLocusWeights = locusWeights;
    vector <double> revisedTaxonWeights = taxonWeights;
    vector <string> revisedLocusNames = locusNames;
    vector <string> revisedTaxonNames = taxonNames;
    bool revisedReferenceTaxonPresent = referenceTaxonPresent;
    vector <int> revisedReferenceTaxa = referenceTaxa;    // Possible that there are several
    double revisedCoverage = taxonCoverage;
    vector <double> partitionDecisiveness(numPartitions, 0.0);
    int numRandomTrees = 0;
    
    bool verbose = true;
    
// this is a placeholder; currently cannot prune taxa from a tree
    vector < vector < vector <bool> > > revisedUserTrees = userTrees; // i.e. if taxa are deleted
    bool completeDecisivenessDetermined = false; // UPDATE THIS!
        
    printSummaryInformation(locusNames, taxonNames, data, taxonCoverage, referenceTaxa, matrixDecisive,
        treewiseDecisiveness, branchwiseDecisiveness, completeDecisivenessDetermined, nexusFileName,
        numRandomTrees, numUserTrees, numProcs);
        
    bool doneEditing = false;
    // bool newMatrix = false;
    bool revisedMatrixDecisive = matrixDecisive;
    
// Where user interfaces
    while (!doneEditing) {
        bool addGenes = false;
        bool merge = false;
        bool exclude = false;
        bool deleteGenes = false;
        bool revert = false;
        bool quit = false;
        bool print = false;
        bool reweightLoci = false;
        bool reweightTaxa = false;
        bool partialBranchwise = false;
        bool partialTreewise = false;
        bool partialIndividualPartition = false;
        bool testCompleteDeciveness = false;
        bool summarize = false;
        bool printRefTaxa = false;
        bool writeCurrentMatrix = false;
        bool testUserTree = false;
        bool findAll = false;
        bool useGA = false;
        int numAddGA = 0;
        
        printProgamOptions (addGenes, merge, exclude, deleteGenes, revert, quit, print, reweightLoci,
            reweightTaxa, partialTreewise, partialBranchwise, summarize, testCompleteDeciveness,
            writeCurrentMatrix, testUserTree, partialIndividualPartition, printRefTaxa, useGA);
        
        if (revert) {    // User wants to start over manipulating original taxon-locus matrix
            revisedLocusNames = locusNames;
            revisedTaxonNames = taxonNames;
            revisedData = data;
            revisedTriplets = triplets;
            revisedTripletLocations = tripletLocations;
            revisedMissingTriplets = missingTriplets;
            revisedMissingQuartets = missingQuartets;
            revisedReferenceTaxonPresent = referenceTaxonPresent;
            revisedCoverage = taxonCoverage;
            revisedReferenceTaxa = referenceTaxa;
            revisedLocusWeights = locusWeights;
            revisedTaxonWeights = taxonWeights;
            revisedMatrixDecisive = matrixDecisive;
            treewiseDecisiveness = 0;
            branchwiseDecisiveness = 0;
            
            cout << endl << "Reverting to original matrix." << endl;
            
            completeDecisivenessDetermined = false; // NO LONGER APPROPRIATE
            
            printSummaryInformation(revisedLocusNames, revisedTaxonNames, revisedData,
                revisedCoverage, revisedReferenceTaxa, revisedMatrixDecisive,
                treewiseDecisiveness, branchwiseDecisiveness, completeDecisivenessDetermined,
                nexusFileName, numRandomTrees, numUserTrees, numProcs);
        } else if (merge) {
            mergeTaxa(revisedData, revisedTaxonNames, revisedTaxonWeights);
        } else if (deleteGenes) { // Um, not useful...
            deletePartitionsFromMatrix(revisedData, revisedLocusNames, revisedLocusWeights, revisedCoverage);
        } else if (exclude) { // get rid of shitty taxa to improve matrix decisiveness
            excludeTaxa(revisedData, revisedTaxonNames, revisedTaxonWeights, revisedCoverage, revisedLocusNames);
        } else if (addGenes) { // virtual genes i.e. for targetted sequencing
            addTaxonGeneToMatrix(revisedData, revisedTaxonNames, revisedLocusNames, revisedLocusWeights, revisedTaxonWeights);
        } else if (partialTreewise || partialBranchwise || partialIndividualPartition) {
            if (partialTreewise) {
                findAll = false;
                treewiseDecisiveness = calculatePartialDecisiveness(revisedReferenceTaxonPresent,
                    numRandomTrees, revisedData, findAll, numProcs, verbose);
            } else if (partialBranchwise) {
                findAll = true;
                branchwiseDecisiveness = calculatePartialDecisiveness(revisedReferenceTaxonPresent,
                    numRandomTrees, revisedData, findAll, numProcs, verbose);
            } else if (partialIndividualPartition) {
                findAll = true;
                int partitionID = selectPartition (revisedData, revisedLocusNames);
                partitionDecisiveness[partitionID] = calculatePartialDecisivenessSinglePartition(revisedReferenceTaxonPresent,
                    numRandomTrees, revisedData, findAll, partitionID, 0, numProcs); // set verbose to false
                cout << "Partial decisiveness for partition '" << locusNames[partitionID]
                    << "' is: " << partitionDecisiveness[partitionID] << endl;
            }
        } else if (testCompleteDeciveness) {
            if (completeDecisivenessDetermined) { // Nothing has changed; report recorded values
                cout << endl << "*** No changes have been made to the matrix since last test ***";
            } else { // some change to matrix has been made; recalculate complete deciveness
                revisedMatrixDecisive = testCompleteDecisivness(revisedData, revisedReferenceTaxonPresent, revisedReferenceTaxa, revisedTaxonNames, revisedMissingQuartets, revisedTriplets, revisedTripletLocations, revisedMissingTriplets);
                // newMatrix = false;
                completeDecisivenessDetermined = true;
            }
        } else if (useGA) {
            cout << endl;
            numAddGA = checkValidIntInput("Enter how many virtual taxon-character(s) to add to matrix: ");
            
            revisedData = GAHandler(numAddGA, revisedData, GADecisiveness, revisedReferenceTaxonPresent, numProcs);
            printGADataToFile(revisedData, taxonNames, nexusFileName, numAddGA);
            
        } else if (print) {
            printMatrix(revisedData, revisedTaxonNames, revisedLocusWeights, revisedTaxonWeights);
        } else if (printRefTaxa) {
            printReferenceTaxa (revisedReferenceTaxa, revisedTaxonNames);
        } else if (writeCurrentMatrix) { // output matrix in nexus or phylip format
            writeMatrix (revisedTaxonNames, numChar, taxaAlignment, includedLocusRanges, revisedLocusNames, dataType);
        } else if (testUserTree) { // determine decisiveness on passed-in user tree
// if doesn't yet exist, get filename from user
            if (treeFileName.size() == 0) {
                treeFileName = getFileName ();
                getUserTrees (treeFileName, rawTrees, userTrees, taxonNames, translationTable, burnin,
                    thinning, treeTaxonOrdering);
            }
            findAll = checkValidBoolInput("Search for minimal (0) or exhaustive (1) coverage? ");
            userTreeDecisiveness = determineDecisivenessUserTree (revisedData, userTrees,
                treeTaxonOrdering, revisedTaxonNames, revisedLocusWeights, revisedTaxonWeights, findAll, numProcs);
            
            writeAnnotatedTrees(rawTrees, translationTable, userTreeDecisiveness, revisedTaxonNames);
        } else if (quit) {
            doneEditing = true;
        }
// If matrix has changed in any way, reset
        if (merge || deleteGenes || exclude || addGenes) {
            treewiseDecisiveness = 0.0;
            branchwiseDecisiveness = 0.0;
            // newMatrix = true;
            completeDecisivenessDetermined = false;
            revisedReferenceTaxonPresent = searchForReferenceTaxon(revisedData, revisedReferenceTaxa, revisedTaxonNames);
            getCoverage(revisedData, revisedCoverage);
        }
        if (partialTreewise || partialBranchwise || testCompleteDeciveness || summarize || merge || deleteGenes || exclude || addGenes) {
            printSummaryInformation(revisedLocusNames, revisedTaxonNames, revisedData,
                revisedCoverage, revisedReferenceTaxa, revisedMatrixDecisive,
                treewiseDecisiveness, branchwiseDecisiveness, completeDecisivenessDetermined,
                nexusFileName, numRandomTrees, numUserTrees, numProcs);
        }
    }
    
    cout << endl << "Fin." << endl << endl;
    return 0;
}
