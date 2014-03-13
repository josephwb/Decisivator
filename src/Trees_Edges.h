#ifndef _TREES_AND_EDGES_H_
#define _TREES_AND_EDGES_H_

// Functions for generating binary trees and extracting internal edges
vector < vector <bool> > fastBinaryTree (int const& numTaxa, vector < vector <int> > & sibNodes,
	bool const& refTaxon);

vector < vector <int> > convertTreeToApeFormat (vector < vector <bool> > const& tree);

void printClade (vector <int> const& clade);
void printClade (vector <bool> const& clade); // overloaded for debugging
void printTree (vector < vector <bool> > const& tree);
void getEdges (int const& edge, vector < vector <bool> > const& tree, vector < vector <int> > & sibNodes,
	bool const& refTaxon, vector <int> & left, vector <int> & right, vector <int> & sib, vector <int> & upper);
vector <int> gatherTips (vector <int> & includedTips, vector <bool> const& clade);
vector <int> getRemainingTaxa (vector <int> & includedTips, int const& numTaxa);
void splitEdge (vector <int> & sib, vector <int> & upper, vector < vector <int> > const& sibNodes,
	int const& sibID, vector < vector <bool> > const& tree);
vector < vector <int> > getSibNodes (vector < vector <bool> > & tree, int const& numProcs);

#endif /* _TREES_AND_EDGES_H_ */
