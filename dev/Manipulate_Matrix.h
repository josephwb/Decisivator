#ifndef _MANIPULATE_MATRIX_H_
#define _MANIPULATE_MATRIX_H_

void addTaxonGeneToMatrix (vector < vector <int> > & data, vector <string> const& taxonNames,
	vector <string> & locusNames, vector <double> & locusWeights, vector <double> & taxonWeights);
void deleteGenesFromMatrix (vector < vector <int> > & data, vector <string> & locusNames, vector <double> & locusWeights);
void excludeTaxa (vector < vector <int> > & data, vector <string> & taxonNames, vector <double> & taxonWeights,
	double & revisedCoverage, vector <string> const& locusNames);
void mergeTaxa (vector < vector <int> > & data, vector <string> & taxonNames, vector <double> & revisedTaxonWeights);
void excludeTaxaMissingNGenes (int const& partitionsMissing, vector < vector <int> > & data,
	vector <string> & taxonNames, vector <double> & taxonWeights, bool const& missingExact);
void excludeTaxaPossessingNGenes (int const& partitionsPossessed, vector < vector <int> > & data,
	vector <string> & taxonNames, vector <double> & taxonWeights, bool const& possessingExact);
void excludeTaxaMinimalOverlap (vector < vector <int> > & data, vector <string> & taxonNames,
	vector <double> & taxonWeights);

#endif /* _MANIPULATE_MATRIX_H_ */