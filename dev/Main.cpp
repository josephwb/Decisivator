/*  
   Decisivator
   Copyright (c) 2011, PuRGe C/O JWB
   Department of Biological Sciences
   University of Idaho, Moscow, Idaho, USA
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

   Any feedback is very welcome, yo.
   email: josephwb@uidaho.edu
*/


/*

A program to analyze a matrix of taxon-gene cells for the purpose of determining idealized phylogenetic
'decisiveness' sensu Michael Sanderson and Mike Steel.

Presence of a gene is denoted by '1', absence by '0'. User may pass in a 0-1 matrix or (more likely)
a Nexus file containing CHARSET declarations (in the latter case, the 0-1 matrix is determined
automatically). 

Searches for all taxon triplets (necessary) and quartets (sufficient) conditions (or, if a reference
taxon exists, that is, one sequenced for all genes, all quads that contain that reference taxon), i.e. 
determines decisiveness of matrix as a whole for ALL possible trees.

At the moment, deals only with a single matrix. User interacts, adding taxon-genes, merging taxa 
(forming chimeric sequences) and/or deleting taxa and genes that may be causing problems (although
deleting genes will never help, except to form reference taxa). In this way, one may (eventually,
and pretty inefficiently) attain a matrix that is phylogenetically decisive.

To determine partial decisiveness (branch-wise), the following will occur:

1. get a topology (edgelengths not required) - DONE
2. traverse tree, edge by edge - DONE
3. determine if two taxa on EACH side of the edge are sequenced for a particular gene
   
'Partial decisivess' is of interest for both the focal tree and 'all' trees. For the latter, it is not
possible to analyze every possible configuration. As an approximation, simulate multiple random trees
of the same taxon size. An important consideration may be to constrain parts of the generated trees
to display non-controversial edges.

Ultimate goal: how to add taxon-genes to the matrix to most optimally (i.e. most frugally or with fixed
effort) improve phylogenetic decisiveness.

The code below will, after final structure decisions (i.e. 'structurally decisive'! Ha!), be
reorganized into an object-oriented form.

Oh, and assumes unix line breaks.

*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

using namespace std;

#include "General.h"
#include "Trees_Edges.h"
#include "Manipulate_Matrix.h"
#include "Matrix_Scan.h"
#include "Parse_Data.h"
#include "User_Interface.h"
#include "Parse_Nexus.h"
#include "User_Tree.h"

bool DEBUG = false; // Print out extra bits to screen

// version information
double version = 0.5;
string month = "November";
int year = 2011;

int main(int argc, char *argv[])
{
	srand ( time(NULL) );     // initialize random seed
	string matrixFileName;
	string locusWeightFileName;
	string taxonWeightFileName;
	string nexusFileName;
	
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
	
	vector < vector <int> > triplets;
	vector < vector <int> > tripletLocations;
	vector < vector <int> > missingTriplets;
	
	vector < vector <int> > missingQuartets;
	vector <double> locusWeights; // default of 1.0; user can change
	vector <double> taxonWeights; // default of 1.0; user can change
	bool referenceTaxonPresent = false;
	vector <int> referenceTaxa;	// Possible that there are several
	vector <int> missingQuartetsByTaxa; // for very large matrices
	double taxonCoverage = 0.0;
	bool matrixDecisive = false;
	double branchwiseDecisiveness = 0.0;
	double treewiseDecisiveness = 0.0;
	int numChar = 0;
	
	printProgramInfo();
	
	processCommandLineArguments(argc, argv, matrixFileName, nexusFileName, locusWeightFileName,
		taxonWeightFileName, treeFileName, burnin, thinning);

// Read in data, store
	if (!nexusFileName.empty())
	{
		parseNexus(nexusFileName, data, taxonNames, locusNames, numChar, taxaAlignment, includedLocusRanges);
	}
	else if (!matrixFileName.empty())
	{
		parseInputMatrix(matrixFileName, locusNames, taxonNames, data);
	}
	else
	{
		cout << "Ergh. Something fucked up... Shit." << endl;
		exit(1);
	}
	
// *** If reference taxon IS present, put at bottom of taxon-gene matrix (as outgroup taxon is always last in tree)
//     Will allow faster satisfaction (but not rejection) of internal edges.
	referenceTaxonPresent = searchForReferenceTaxon(data, referenceTaxa, taxonNames);
	
	if (!treeFileName.empty())
	{
		getUserTrees (treeFileName, rawTrees, userTrees, taxonNames, translationTable, burnin,
			thinning, treeTaxonOrdering);
		numUserTrees = userTrees.size();
		if (DEBUG)
		{
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
	vector <int> revisedReferenceTaxa = referenceTaxa;	// Possible that there are several
	double revisedCoverage = taxonCoverage;
	int numRandomTrees = 0;
	
	vector < vector < vector <bool> > > revisedUserTrees = userTrees; // i.e. if taxa are deleted
	bool completeDecisivenessDetermined = false; // UPDATE THIS!
		
	printSummaryInformation(locusNames, taxonNames, data, taxonCoverage, referenceTaxa, matrixDecisive,
		treewiseDecisiveness, branchwiseDecisiveness, completeDecisivenessDetermined, nexusFileName,
		numRandomTrees, numUserTrees);
		
	bool doneEditing = false;
	bool newMatrix = false;
	bool revisedMatrixDecisive = matrixDecisive;
	
// Where user interfaces
	while (!doneEditing)
	{
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
		bool testCompleteDeciveness = false;
		bool summarize = false;
		bool writeCurrentMatrix = false;
		bool testUserTree = false;
		bool findAll = false;
		
		printProgamOptions (addGenes, merge, exclude, deleteGenes, revert, quit, print, reweightLoci,
			reweightTaxa, partialTreewise, partialBranchwise, summarize, testCompleteDeciveness,
			writeCurrentMatrix, testUserTree);
		
		if (revert)	// User wants to start over manipulating original taxon-gene matrix
		{
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
			
			printSummaryInformation(locusNames, taxonNames, data, taxonCoverage, referenceTaxa, matrixDecisive,
				treewiseDecisiveness, branchwiseDecisiveness, completeDecisivenessDetermined, nexusFileName,
				numRandomTrees, numUserTrees);
		}
		else if (merge)
		{
			mergeTaxa(revisedData, revisedTaxonNames, revisedTaxonWeights);
		}
		else if (deleteGenes) // Um, not useful...
		{
			deleteGenesFromMatrix(revisedData, revisedLocusNames, revisedLocusWeights);
		}
		else if (exclude) // get rid of shitty taxa to improve matrix decisiveness
		{
			excludeTaxa(revisedData, revisedTaxonNames, revisedTaxonWeights, revisedCoverage, revisedLocusNames);
		}
		else if (addGenes) // virtual genes i.e. for targetted sequencing
		{
			addTaxonGeneToMatrix(revisedData, revisedTaxonNames, revisedLocusNames, revisedLocusWeights, revisedTaxonWeights);
		}
		else if (partialTreewise || partialBranchwise)
		{
			if (partialTreewise)
			{
				findAll = false;
				treewiseDecisiveness = calculatePartialDecisiveness(revisedReferenceTaxonPresent,
					numRandomTrees, revisedData, findAll);
			}
			else
			{
				findAll = true;
				branchwiseDecisiveness = calculatePartialDecisiveness(revisedReferenceTaxonPresent,
					numRandomTrees, revisedData, findAll);
			}
		}
		else if (testCompleteDeciveness)
		{
			if (completeDecisivenessDetermined) // Nothing has changed; report recorded values
			{
				cout << endl << "*** No changes have been made to the matrix since last test ***";
			}
			else // some change to matrix has been made; recalculate complete deciveness
			{
				revisedMatrixDecisive = testCompleteDecisivness(revisedData, revisedReferenceTaxonPresent, revisedReferenceTaxa, revisedTaxonNames, revisedMissingQuartets, revisedTriplets, revisedTripletLocations, revisedMissingTriplets);
				newMatrix = false;
				completeDecisivenessDetermined = true;
			}
		}
		else if (print)
		{
			printMatrix(revisedData, revisedTaxonNames, revisedLocusWeights, revisedTaxonWeights);
		}
		else if (writeCurrentMatrix) // output matrix in nexus or phylip format
		{
			writeMatrix (revisedTaxonNames, numChar, taxaAlignment, includedLocusRanges, locusNames);
		}
		else if (testUserTree) // determine decisiveness on passed-in user tree
		{
			
			
			
// WORKING
// if doesn't yet exist, get filename from user
			if (treeFileName.size() == 0)
			{
				treeFileName = getFileName ();
				getUserTrees (treeFileName, rawTrees, userTrees, taxonNames, translationTable, burnin,
					thinning, treeTaxonOrdering);
			}
			userTreeDecisiveness = determineDecisivenessUserTree (revisedData, userTrees,
				treeTaxonOrdering, revisedTaxonNames, revisedLocusWeights, revisedTaxonWeights);
			
			writeAnnotatedTrees (rawTrees, translationTable, userTreeDecisiveness, taxonNames);
			
			
			
		}
		else if (quit)
		{
			doneEditing = true;
		}
// If matrix has changed in any way, reset
		if (merge || deleteGenes || exclude || addGenes)
		{
			treewiseDecisiveness = 0.0;
			branchwiseDecisiveness = 0.0;
			newMatrix = true;
			completeDecisivenessDetermined = false;
			revisedReferenceTaxonPresent = searchForReferenceTaxon(revisedData, revisedReferenceTaxa, revisedTaxonNames);
			getCoverage(revisedData, revisedCoverage);
		}
		if (partialTreewise || partialBranchwise || testCompleteDeciveness || summarize || merge || deleteGenes || exclude || addGenes)
		{
			printSummaryInformation(locusNames, taxonNames, data, taxonCoverage, referenceTaxa, matrixDecisive,
				treewiseDecisiveness, branchwiseDecisiveness, completeDecisivenessDetermined, nexusFileName,
				numRandomTrees, numUserTrees);
		}
	}
	
	cout << endl << "Fin." << endl << endl;
	return 0;
}