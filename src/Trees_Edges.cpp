#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>

using namespace std;

#include "Trees_Edges.h"
#include "General.h"

extern bool debuggering; // print out extra junk to screen

// Adapted from LJH R code
vector < vector <bool> > fastBinaryTree (int const& numTaxa, vector < vector <int> > & sibNodes,
	bool const& revisedReferenceTaxonPresent)
{
	vector <int> anc;
	vector <int> dec;
	
	int simTaxa = 0;
	if (revisedReferenceTaxonPresent) // if reference taxon present, simulate tree of size numTaxa - 1, add outgroup after
	{
		simTaxa = numTaxa - 1;
	}
	else
	{
		simTaxa = numTaxa;
	}
	
	sibNodes.clear();
	vector < vector <bool> > tree (simTaxa, vector <bool> (simTaxa, false)); // empty tree representation
	vector <int> alive;       // index for available lineages
	vector <int> toCoalesce;  // randomly sampled lineages to coalesce
	int randomNum = 0;
	
	for (int i = 0; i < simTaxa; i++) // identify tips == lineages to initially sample
	{
		tree[i][i] = true;
		alive.push_back(i);
	}
	
	if (debuggering) {cout << endl;}
	for (int i = simTaxa; i < (2 * simTaxa - 1); i++) // there will be T-2 internal nodes (unrooted tree)
	{
		vector <bool> newClade;
		if (debuggering) {cout << "Sampled lineages (node " << i << "):";}
		for (int j = 0; j < 2; j++) // select two available lineages to coalesce
		{
			randomNum = rand() % (int)alive.size();
			
			if (debuggering) {cout << "randomNum = " << randomNum << endl;}
			
			toCoalesce.push_back(alive[randomNum]);
			
			
			dec.push_back(alive[randomNum]);
			anc.push_back(i);
			
			
			if (debuggering) {cout << " " << alive[randomNum];}
			alive.erase(alive.begin()+randomNum);
		}
		if (debuggering) {cout << endl;}
		sibNodes.push_back(toCoalesce);
		for (int k = 0; k < simTaxa; k++) // coalesce
		{
			newClade.push_back(tree[toCoalesce[0]][k] + tree[toCoalesce[1]][k]); // just add the boolean values from clades toCoalesce[0] and toCoalesce[1]
		}
		alive.push_back(i); // new lineage (clade) available
		tree.push_back(newClade);
		newClade.clear();
		toCoalesce.clear();
	}
/*	
	if (revisedReferenceTaxonPresent) // reference taxon present; add to tree
	{
		vector <bool> root (numTaxa, true); // add reference taxon as root
		for (int i = 0; i < int(tree[1].size()); i++)
		{
			tree[i].push_back(false);
		}
		tree.push_back(root);
	}
*/
	
	if (debuggering) {printTree(tree);}
	
	printTree(tree);
	
	
	
	
	
	cout << "anc	dec" << endl;
	for (int i = 0; i < (int)anc.size(); i++)
	{
		cout << anc[i] << "	" << dec[i] << endl;
	}
	
	return (tree);
}





// return two-column list of edges
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
	for (int i = 0; i < int(clade.size()); i++)
	{
		cout << " " << clade[i];
	}
	cout << ";" << endl;
}

void printClade (vector <bool> const& clade) // overloaded for debugging
{
	for (int i = 0; i < int(clade.size()); i++)
	{
		cout << " " << clade[i];
	}
}


// Just for debugging
void printTree (vector < vector <bool> > const& tree)
{
	cout << endl << "Clades:" << endl;
	for (int i = 0; i < int(tree.size()); ++i)
	{
		cout << i << "	";
		for (int j = 0; j < int(tree[1].size()); j++)
		{
			cout << tree[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

vector <int> gatherTips (vector <int> & includedTips, vector <bool> const& clade)
{
	for (int i = 0; i < int(clade.size()); i++)
	{
		if (clade[i])
		{
			includedTips.push_back(i);
		}
	}
	return(includedTips);
}

vector <int> getRemainingTaxa (vector <int> & includedTips, int const& numTaxa)
{
	vector <int> remainingTaxa;
	sort(includedTips.begin(), includedTips.end());
	for (int i = 0; i < numTaxa; i++)
	{
		bool match = false;
		for (int j = 0; j < int(includedTips.size()); j++)
		{
			if (includedTips[j] == i)
			{
				match = true;
			}
		}
		if (!match)
		{
			remainingTaxa.push_back(i);
		}
	}
	return(remainingTaxa);
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

vector < vector <int> > getSibNodes (vector < vector <bool> > & tree)
{
	vector < vector <int> > sibNodes;
	int numTaxa = (int)tree[0].size();
	vector <int> alive;       // index for available lineages
	vector <int> temp;
	
	for (int i = 0; i < numTaxa; i++)
	{
		alive.push_back(i);
	}
	
	if (debuggering) {cout << "Tree before sib search:" << endl; printTree(tree);}
	
// loop over nodes
	for (int i = numTaxa; i < (int)tree.size(); i++)
	{
		vector <bool> testClade;
		if (debuggering) {cout << "Looking for descendant of node: " << i << endl;}
		
		for (int j = 0; j < (int)alive.size(); j++)
		{
			for (int k = 0; k < (int)alive.size(); k++)
			{
				for (int l = 0; l < numTaxa; l++)
				{
					testClade.push_back(tree[alive[j]][l] + tree[alive[k]][l]);
				}
				if (testClade == tree[i])
				{
					if (debuggering) {cout << "Holy shit. It worked?!? Descendant nodes are: "
						<< alive[j] << " & " << alive[k] << "." << endl;}
					temp.push_back(alive[j]);
					temp.push_back(alive[k]);
					sibNodes.push_back(temp);
					
					alive.erase(alive.begin()+k);
					alive.erase(alive.begin()+j);
					if (i != ((int)tree.size() - 1))
					{
						alive.push_back(i);
						if (debuggering) {cout << "Adding node " << i << " to alive set." << endl;}
					}
					j = k = (int)alive.size(); // exit
					temp.clear();
				}
				testClade.clear();
			}
		}
		
		if (i == ((int)tree.size() - 1) && !alive.empty()) // i.e. dealing with an unrooted tree
		{
			int numNodes = (int)tree.size() - 1;
			temp.push_back(alive[0]); // take first two; 3rd in the last clade identified
			temp.push_back(alive[1]);
			
			if (debuggering)
			{
				cout << "alive of of size: " << alive.size() << endl;
				cout << "alive of of size: " << alive.size() << endl;
				printVectorAsList(alive);
				cout << "Adding the unrooted node:" << endl;
				printVectorAsList(temp);
			}
			
			sibNodes.push_back(temp);
			temp.clear();
			
			for (int l = 0; l < numTaxa; l++)
			{
				testClade.push_back(tree[alive[0]][l] + tree[alive[1]][l]);
			}
			
			tree.insert((tree.end()-1), testClade);
			temp.push_back(numNodes-1);
			temp.push_back(numNodes);
			sibNodes.push_back(temp);
			if (debuggering) {cout << "Added sibNodes " << alive[0] << " & " << alive[1]
				<< " and " << numNodes << " & " << numNodes-1 << "." << endl;}
			
			if (debuggering) {printTree(tree);}
			i = (int)tree.size();
		}
	}
	
	if (debuggering) {cout << "Got all sibNodes!" << endl; printTree(tree);}
	
	return (sibNodes);
}

void getEdges (int const& edge, vector < vector <bool> > const& tree, vector < vector <int> > & sibNodes,
	bool const& revisedReferenceTaxonPresent, vector <int> & left, vector <int> & right, vector <int> & sib,
	vector <int> & upper)
{
	int numTaxa = (int)tree[0].size();
	int simTaxa = 0;
	vector <int> includedTips;
	
	if (revisedReferenceTaxonPresent)
	{
		simTaxa = numTaxa - 1;
	}
	else
	{
		simTaxa = numTaxa;
	}
		
	int node = edge + simTaxa; // first N nodes are 'singletons'. Gah! who wrote this shit?!? Oh, yeah...
	int sibID = 0;
	if (debuggering) {cout << endl << "Internal node #" << node << endl;}
	
// Right and left children:
	includedTips = gatherTips(includedTips, tree[sibNodes[edge][0]]);
	left = gatherTips(left, tree[sibNodes[edge][0]]);
	includedTips = gatherTips(includedTips, tree[sibNodes[edge][1]]);
	right = gatherTips(right, tree[sibNodes[edge][1]]);
	
//	cout << "Ancestral node: (";
	for (int j = edge + 1; j < int(sibNodes.size()); j++)
	{
		if (sibNodes[j][0] == node)
		{
//			cout << simTaxa + j << ");" << endl;	// ancestral node
			sibID = sibNodes[j][1];
			sib = gatherTips(sib, tree[sibNodes[j][1]]);
			includedTips = gatherTips(includedTips, tree[sibNodes[j][1]]);
		}
		else if (sibNodes[j][1] == node)
		{
//			cout << simTaxa + j << ");" << endl;	// ancestral node
			sibID = sibNodes[j][0];
			sib = gatherTips(sib, tree[sibNodes[j][0]]);
			includedTips = gatherTips(includedTips, tree[sibNodes[j][0]]);
		}
	}
	
	
	upper = getRemainingTaxa(includedTips, numTaxa);
	
	if (debuggering)
	{
		if (upper.empty())
		{
			cout << "*** UPPER IS EMPTY!!! ***" << endl;
			cout << "Left:";  printClade(left);
			cout << "Right:"; printClade(right);
			cout << "Sib:";   printClade(sib);
			cout << "Upper:"; printClade(upper);
			cout << endl << "Fixed:" << endl;
		}
	}
	
	if (upper.size() == 0) // only occurs for 'unrooted' trees
	{
		splitEdge(sib, upper, sibNodes, sibID - numTaxa, tree);
	}
	includedTips.clear();
	
// reverse the 'upper' vector, as reference taxon (if present) will be the last taxon
	reverse(upper.begin(), upper.end());
	
	if (debuggering)
	{
		cout << "Left:";  printClade(left);
		cout << "Right:"; printClade(right);
		cout << "Sib:";   printClade(sib);
		cout << "Upper:"; printClade(upper);
	}
}
