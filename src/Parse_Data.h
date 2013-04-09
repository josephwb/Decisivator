#ifndef _PARSE_DATA_H_
#define _PARSE_DATA_H_

vector < vector <int> > storeData (vector <string> & inputData, int const& numLoci);
vector <string> getTaxonNames(vector <string> & inputData);
void getWeights (string const& weightFileName, vector <double> & weights, vector <string> const& names);
void parseInputMatrix (string const& inputFileName, vector <string> & locusNames, vector <string> & taxonNames,
	vector < vector <int> > & data);

#endif /* _PARSE_DATA_H_ */
