#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <set>
#include <numeric>

using namespace std;

#include "Trees_Edges.h"
#include "General.h"

extern bool debugging; // print out extra junk to screen

// Adapted from LJH R code
vector < vector <bool> > fastBinaryTree (int const& numTaxa, vector < vector <int> > & sibNodes,
    bool const& revisedReferenceTaxonPresent)
{
    vector <int> anc;
    vector <int> dec;
    
    int simTaxa = 0;
    if (revisedReferenceTaxonPresent) { // if reference taxon present, simulate tree of size numTaxa - 1, add outgroup after
        simTaxa = numTaxa - 1;
    } else {
        simTaxa = numTaxa;
    }
    
    sibNodes.clear();
    vector < vector <bool> > tree (simTaxa, vector <bool> (simTaxa, false)); // empty tree representation
    vector <int> alive;       // index for available lineages
    vector <int> toCoalesce;  // randomly sampled lineages to coalesce
    int randomNum = 0;
    
    for (int i = 0; i < simTaxa; i++) { // identify tips == lineages to initially sample
        tree[i][i] = true;
        alive.push_back(i);
    }
    
    if (debugging) {cout << endl;}
    for (int i = simTaxa; i < (2 * simTaxa - 1); i++) { // there will be T-2 internal nodes (unrooted tree)
        vector <bool> newClade;
        if (debugging) {cout << "Sampled lineages (node " << i << "):";}
        for (int j = 0; j < 2; j++) { // select two available lineages to coalesce
            randomNum = rand() % (int)alive.size();
            
            if (debugging) {cout << "randomNum = " << randomNum << endl;}
            
            toCoalesce.push_back(alive[randomNum]);
            
            dec.push_back(alive[randomNum]);
            anc.push_back(i);
            
            if (debugging) {cout << " " << alive[randomNum];}
            alive.erase(alive.begin()+randomNum);
        }
        if (debugging) {cout << endl;}
        sibNodes.push_back(toCoalesce);
        for (int k = 0; k < simTaxa; k++) { // coalesce
            newClade.push_back(tree[toCoalesce[0]][k] + tree[toCoalesce[1]][k]); // just add the boolean values from clades toCoalesce[0] and toCoalesce[1]
        }
        alive.push_back(i); // new lineage (clade) available
        tree.push_back(newClade);
        newClade.clear();
        toCoalesce.clear();
    }
    if (revisedReferenceTaxonPresent) { // reference taxon present; add to tree
        vector <bool> root (numTaxa, true); // add reference taxon as root
        for (int i = 0; i < int(tree[1].size()); i++) {
            tree[i].push_back(false);
        }
        tree.push_back(root);
    }
    
    if (debugging) {
        printTree(tree);
        
        cout << "Sibnodes:" << endl;
        for (int i = 0; i < (int)sibNodes.size(); i++) {
            printClade(sibNodes[i]);
        }
        cout << endl;
    }
    
    if (debugging) {
        cout << "anc    dec" << endl;
        for (int i = 0; i < (int)anc.size(); i++) {
            cout << anc[i] << "    " << dec[i] << endl;
        }
    }
    
    return (tree);
}





// return two-column list of edges; should be more efficient
vector < vector <int> > convertTreeToApeFormat (vector < vector <bool> > const& tree)
{
    vector < vector <int> > treeEdges;
    
    
    
    
    
    
    
    
    
    
    
    
    return treeEdges;
}

/*
  1. identify internal edge (total of numTaxa - 3) - DONE
  2. left descendants - DONE
  3. right descendants - DONE
  4. sib descendants - DONE
  5. remaining taxa (i.e. all taxa on the other side of the ancestral edge) - DONE
     A slightly different approach is needed for 'rooted' vs. 'unrooted' trees.
    - this is called 'upper' below
        - may be ancestral edge (if present)
        - upperwise, one half of sib clade
  
  I added the option of having a reference taxon (one sequenced for all genes).
  In this case, simulate tree for n - 1 taxa, tack on outgroup (James Foster's idea).
  
  Reference taxon (if present) is always last in the tree; sort taxon-gene matrix accordingly. Did I do this yet?!?
  
  *** Hmm. This shit doesn't appear necessary at all. ***
  
*/

void printClade (vector <int> const& clade)
{
    for (int i = 0; i < int(clade.size()); i++) {
        cout << " " << clade[i];
    }
    cout << ";" << endl;
}

void printClade (vector <bool> const& clade) // overloaded for debugging
{
    for (int i = 0; i < int(clade.size()); i++) {
        cout << " " << clade[i];
    }
}

// Just for debugging
void printTree (vector < vector <bool> > const& tree)
{
    cout << endl << "Clades:" << endl;
    for (int i = 0; i < int(tree.size()); ++i) {
        cout << i << "    ";
        for (int j = 0; j < int(tree[1].size()); j++) {
            cout << tree[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

vector <int> gatherTips (vector <int> & includedTips, vector <bool> const& clade)
{
    for (int i = 0; i < int(clade.size()); i++) {
        if (clade[i]) {
            includedTips.push_back(i);
        }
    }
    return(includedTips);
}

vector <int> getRemainingTaxa (vector <int> & includedTips, int const& numTaxa)
{
    vector <int> remainingTaxa;
    sort(includedTips.begin(), includedTips.end());
    for (int i = 0; i < numTaxa; i++) {
        bool match = false;
        for (int j = 0; j < int(includedTips.size()); j++) {
            if (includedTips[j] == i) {
                match = true;
            }
        }
        if (!match) {
            remainingTaxa.push_back(i);
        }
    }
    return(remainingTaxa);
}




void getEdges (int const& edge, vector < vector <bool> > const& tree, vector < vector <int> > & sibNodes,
    bool const& revisedReferenceTaxonPresent, vector <int> & left, vector <int> & right, vector <int> & sib,
    vector <int> & upper)
{
    int numTaxa = (int)tree[0].size();
    int simTaxa = 0;
    vector <int> includedTips;
    
    if (revisedReferenceTaxonPresent) { // reference taxon will root the tree
        simTaxa = numTaxa - 1;
    } else {
        simTaxa = numTaxa;
    }
    
    int node = edge + simTaxa; // first N nodes are 'singletons'. Gah! who wrote this shit?!? Oh, yeah...
    int sibID = 0;
    if (debugging) {cout << endl << "Internal node #" << node << endl;}
    
// Right and left children:
    includedTips = gatherTips(includedTips, tree[sibNodes[edge][0]]);
    left = gatherTips(left, tree[sibNodes[edge][0]]);
    includedTips = gatherTips(includedTips, tree[sibNodes[edge][1]]);
    right = gatherTips(right, tree[sibNodes[edge][1]]);
    
//    cout << "Ancestral node: (";
    for (int j = edge + 1; j < int(sibNodes.size()); j++) {
        if (sibNodes[j][0] == node) {
//            cout << simTaxa + j << ");" << endl;    // ancestral node
            sibID = sibNodes[j][1];
            sib = gatherTips(sib, tree[sibNodes[j][1]]);
            includedTips = gatherTips(includedTips, tree[sibNodes[j][1]]);
        } else if (sibNodes[j][1] == node) {
//            cout << simTaxa + j << ");" << endl;    // ancestral node
            sibID = sibNodes[j][0];
            sib = gatherTips(sib, tree[sibNodes[j][0]]);
            includedTips = gatherTips(includedTips, tree[sibNodes[j][0]]);
        }
    }
    
    upper = getRemainingTaxa(includedTips, numTaxa);
    
    if (debugging) {
        if (upper.empty()) {
            cout << "*** UPPER IS EMPTY!!! ***" << endl;
            cout << "Left:";  printClade(left);
            cout << "Right:"; printClade(right);
            cout << "Sib:";   printClade(sib);
            cout << "Upper:"; printClade(upper);
            cout << endl << "Fixed:" << endl;
        }
    }
    
    if (upper.size() == 0) { // only occurs for 'unrooted' trees
        splitEdge(sib, upper, sibNodes, sibID - numTaxa, tree);
    }
    includedTips.clear();
    
// reverse the 'upper' vector, as reference taxon (if present) will be the last taxon
    reverse(upper.begin(), upper.end());
    
    if (debugging) {
        cout << "Left:";  printClade(left);
        cout << "Right:"; printClade(right);
        cout << "Sib:";   printClade(sib);
        cout << "Upper:"; printClade(upper);
    }
}

// Unrooted tree; split sib edge
void splitEdge (vector <int> & sib, vector <int> & upper, vector < vector <int> > const& sibNodes,
    int const& sibID, vector < vector <bool> > const& tree)
{
    sib.clear();
    upper.clear();
    sib = gatherTips(sib, tree[sibNodes[sibID][0]]);
    upper = gatherTips(upper, tree[sibNodes[sibID][1]]);
}





// currently a HUGE bottleneck. FIXED!!! Went from 60+ seconds to 0 seconds for Angiosperm tree
// need to use some clever set algorithm to speed things up.
// still want to use multithreading somehow... depends on whether the ordering of sibNodes is important (it probably is).
// or store a real tree and traverse it
/* Currently: tree is stored in sorted fashion, from smaller (tips) to larger clades. e.g.:
Clades:
0    0 0 0 0 0 0 0 0 0 1 
1    0 0 0 0 0 0 0 0 1 0 
2    0 0 0 0 0 0 0 1 0 0 
3    0 0 0 0 0 0 1 0 0 0 
4    0 0 0 0 0 1 0 0 0 0 
5    0 0 0 0 1 0 0 0 0 0 
6    0 0 0 1 0 0 0 0 0 0 
7    0 0 1 0 0 0 0 0 0 0 
8    0 1 0 0 0 0 0 0 0 0 
9    1 0 0 0 0 0 0 0 0 0 
10    0 0 0 0 0 0 0 0 1 1 
11    0 0 0 0 0 1 1 0 0 0 
12    0 0 0 0 1 1 1 0 0 0 
13    0 0 0 1 1 1 1 0 0 0 
14    0 0 0 1 1 1 1 1 0 0 
15    0 0 0 1 1 1 1 1 1 1 
16    0 1 1 0 0 0 0 0 0 0 
17    0 1 1 1 1 1 1 1 1 1 
18    1 1 1 1 1 1 1 1 1 1
So, guaranteed to find descendant clades above clade of interest.
*/
vector < vector <int> > getSibNodes (vector < vector <bool> > & tree, int const& numProcs)
{
    vector < vector <int> > sibNodes;
    int numTaxa = (int)tree[0].size();
    vector <int> alive (numTaxa, 0);        // index for available lineages
    vector <int> temp;
    
    vector <int> cladeSizes ((int)tree.size(), 0);    // used for avoiding certain comparisons
    for (int i = 0; i < (int)tree.size(); i++) {
        cladeSizes[i] = std::accumulate(tree[i].begin(),tree[i].end(),0);
    }
    
    for (int i = 0; i < numTaxa; i++) {
        alive[i] = i;
    }
    
    if (debugging) {cout << "Tree before sib search:" << endl; printTree(tree);}
    
//    cout << "Tree before sib search:" << endl; printTree(tree);
//    cout << "Worst case scenario: make " << (((int)tree.size() - numTaxa) * (int)alive.size() * (int)alive.size() * numTaxa) << " comparisons here." << endl;
    
// loop over nodes. this is fucking ugly. find a good algorithm to perform this.
    for (int i = numTaxa; i < (int)tree.size(); i++) {
        vector <bool> testClade;
//        if (debugging) {cout << "Looking for descendant of node: " << i << endl;}
        
        int currentCladeSize = cladeSizes[i];
    //    cout << "Dealing with clade #" << i << endl;
        
        for (int j = 0; j < (int)alive.size(); j++) {
            // check that clade j is 1) smaller than the current clade and 2) can possibly be a descendant of the current clade
            // checks like: tree[alive[j]] <= tree[i] are element-wise set operations, bailing on the first false
            if ((cladeSizes[alive[j]] < currentCladeSize) && (tree[alive[j]] <= tree[i])) {
    //            cout << "We've got a contender here." << endl;
                for (int k = 0; k < (int)alive.size(); k++) {
                // check that 1) clade j + k == current clade size and 2) clade k can possibly be a descendant of the current clade
                    if ((cladeSizes[alive[j]] + cladeSizes[alive[k]] == currentCladeSize) && (tree[alive[k]] <= tree[i])) {
                        testClade.reserve(tree[i].size());
                        transform(tree[alive[j]].begin(), tree[alive[j]].end(), tree[alive[k]].begin(), back_inserter(testClade), plus<bool>());
                        if (testClade == tree[i]) {
                            if (debugging) {cout << "Holy shit. It worked?!? Descendant nodes are: " << alive[j] << " & " << alive[k] << "." << endl;}
                            temp.push_back(alive[j]);
                            temp.push_back(alive[k]);
                            sibNodes.push_back(temp);
                    
                            alive.erase(alive.begin()+k);
                            alive.erase(alive.begin()+j);
                            if (i != ((int)tree.size() - 1)) {
                                alive.push_back(i);
                                if (debugging) {cout << "Adding node " << i << " to alive set." << endl;}
                            }
                            j = k = (int)alive.size(); // exit
                            temp.clear();
                        }
                        testClade.clear();
                    }
                }
            }
        }
        
        if (i == ((int)tree.size() - 1) && !alive.empty()) { // i.e. dealing with an unrooted tree
            int numNodes = (int)tree.size() - 1;
            temp.push_back(alive[0]); // take first two; 3rd in the last clade identified
            temp.push_back(alive[1]);
            
            if (debugging) {
                cout << "alive is of size: " << alive.size() << endl;
                printVectorAsList(alive);
                cout << "Adding the unrooted node:" << endl;
                printVectorAsList(temp);
            }
            
            sibNodes.push_back(temp);
            temp.clear();
            
            for (int l = 0; l < numTaxa; l++) {
                testClade.push_back(tree[alive[0]][l] + tree[alive[1]][l]);
            }
            
            tree.insert((tree.end()-1), testClade);
            temp.push_back(numNodes-1);
            temp.push_back(numNodes);
            sibNodes.push_back(temp);
            if (debugging) {cout << "Added sibNodes " << alive[0] << " & " << alive[1]
                << " and " << numNodes << " & " << numNodes-1 << "." << endl;}
            
            if (debugging) {printTree(tree);}
            i = (int)tree.size();
        }
    }
    
    if (debugging) {cout << "Got all sibNodes!" << endl; printTree(tree);}
    
    return (sibNodes);
}
