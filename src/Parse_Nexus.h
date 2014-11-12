#ifndef _PARSE_NEXUS_H_
#define _PARSE_NEXUS_H_

void parseNexus (string const& nexusFileName, vector < vector <int> > & data, vector <string> & taxonNames,
    vector <string> & locusNames, int & numChar, vector < vector <string> > & taxaAlignment,
    vector < vector <int> > & includedLocusRanges, string & dataType);

bool checkCommentLine (string stringToParse);
string removeStringSuffix (string stringToParse, char suffixToRemove, bool & suffixEncountered);
string removeStringPrefix (string stringToParse, char characterToRemove);

void getAttributes (string fileName, int & numTaxa, int & numChar, bool & interleavedData,
    string & dataType);
vector <string> collectCharsets (string charsetFileName, vector<string> inputCharsets,
    vector < vector <int> > & includedLocusRanges, int const& numChar);
vector < vector <string> > collectTaxaAlignment (string fileName, int const& numTaxa, int const& numChar,
    bool const& interleavedData, vector <string> & taxonNames);
void printCollectedCharsets (vector <string> const& collectedCharsets, vector < vector <int> > const& includedLocusRanges);
int getIntervalPartitionLength (int const& partitionID, vector < vector <int> > const& includedLocusRanges);
void extractSimpleCharsetRanges (string charsetDeclaration, int & start, int & stop);
void extractIntervalRange(string charsetDeclaration, int & start, int & stop, int & interval);
void constructMatrix (vector < vector <string> > const& taxaAlignment, vector < vector <int> > const& includedLocusRanges,
    vector < vector <int> > & data, vector <string> const& locusNames);
bool validCharacterEncountered (char const& character);

#endif /* _PARSE_NEXUS_H_ */
