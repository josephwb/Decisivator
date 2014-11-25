#ifndef _USER_TREES_H_
#define _USER_TREES_H_

void getUserTrees (string const& treeFileName, vector <string> & rawTrees,
    vector < vector < vector <bool> > > & userTrees, vector <string> const& taxonNames,
    vector <int> & translationTable, int const& burnin, int const& thinning,
    vector < vector <int> > & treeTaxonOrdering);
vector < vector <bool> > parseTree (string & tree, vector <string> const& taxonNames,
    vector <int> & translationTable, vector < vector <int> > & treeTaxonOrdering);
vector <int> collectTaxonCoding (vector <string> & rawInput, vector <string> const& taxonNames);
void writeAnnotatedTrees (string const& treeFileName, vector <string> const& rawTrees, vector <int> & translationTable,
    vector < vector <double> > userTreeDecisiveness, vector <string> const& taxonNames);

#endif /* _USER_TREES_H_ */
