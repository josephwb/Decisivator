#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

using namespace std;

#include "User_Tree.h"
#include "General.h"
#include "Parse_Nexus.h"
#include "Trees_Edges.h"

extern bool DEBUG; // print out extra junk to screen

// *** Need to determine if translation table exists (it probably does). DONE. ***
// *** Log with vector <int> (mapping names to translation table). DONE. ***

void getUserTrees (string const& treeFileName, vector <string> & rawTrees,
	vector < vector < vector <bool> > > & userTrees, vector <string> const& taxonNames,
	vector <int> & translationTable, int const& burnin, int const& thinning,
	vector < vector <int> > & treeTaxonOrdering)
{
	ifstream inputTree;
	vector <string> rawInput;
	vector <string> treeNames;
	vector < vector <bool> > tempTree;
	
	bool commentLine = false;
	bool whiteSpaceOnly = false;
	
	int treeCounter = 0;
	
	inputTree.open(treeFileName.c_str());
	string line;
	
// Read in every non-empty (or non-whitespace), non-commented-out line
	while (getline(inputTree,line))
	{
		commentLine = checkCommentLine(line);
		whiteSpaceOnly = checkWhiteSpaceOnly(line);
		if (line.empty() || commentLine || whiteSpaceOnly)
		{
			continue;
		}
		else
		{
			rawInput.push_back(line);
		}
		line.clear();
	}
	inputTree.close();
	
	cout << endl << "READING IN USER TREE(S)..." << endl;
	
	translationTable = collectTaxonCoding (rawInput, taxonNames);
	
	if (translationTable.size() == 0)
	{
		cout << "No translation table found." << endl << endl;
	}
	else
	{
		cout << "Translation table found." << endl;
		cout << "Successfully read in coding for " << translationTable.size() << " taxa." << endl << endl;
		if (DEBUG) {printVectorAsList(translationTable);}
	}
	
//	'tree MUSCLE50genePart = (((((((Rhynochetos:29.585829,Eurypyga:29.585829):39.775412,Phaethon:69.361241):2.341594,((Syrrhaptes:67.522426,(Podiceps:43.000'
	for (vector <string>::iterator lineIter = rawInput.begin(); lineIter < rawInput.end(); lineIter++)
	{
		string name;
		int stringPosition = 0;
		commentLine = false;
		string treeDeclaration;
		
		if (checkStringValue(*lineIter, "tree", stringPosition))
		{
			stringPosition++;
			treeCounter++;
			
			if ((treeCounter > burnin) && (treeCounter % thinning == 0)) // don't bother with irrelevant trees
			{
				rawTrees.push_back(*lineIter);
				name = extractStringElement(*lineIter, stringPosition);
				cout << "Processing tree '" << name << "'." << endl;
				stringPosition++;
				
				bool equalSignEncountered = false;
				while (!equalSignEncountered)
				{
					equalSignEncountered = checkStringValue(*lineIter, "=", stringPosition);
					stringPosition++;
				}
				
				
				bool weAreCool = false;
				while (!weAreCool) // will skip rooting flavour indicator
				{
					if ((commentLine = checkCommentLine(extractStringElement(*lineIter, stringPosition))))
					{
						stringPosition++;
					}
					else
					{
						treeDeclaration = extractStringElement(*lineIter, stringPosition);
						weAreCool = true;
					}
				}
				if (DEBUG) {cout << "Tree '" << name << "' = "<< treeDeclaration << endl << endl;}
				
//				rawTrees.push_back(treeDeclaration);
				treeNames.push_back(name);
				tempTree = parseTree(treeDeclaration, taxonNames, translationTable, treeTaxonOrdering);
				userTrees.push_back(tempTree);
				
				
				if (DEBUG)
				{
					cout << "Tree '" << name << "' = "<< treeDeclaration << endl;
					printTree(tempTree);
				}
			}
			treeDeclaration.clear();
			name.clear();
			tempTree.clear();
		}
	}
	rawInput.clear();
	
	if (treeCounter != (int)userTrees.size())
	{
		cout << "Read in " << treeCounter << " user trees; retained " << userTrees.size() << " of them (discarded first "
		<< burnin << " trees, retained every " << thinning << " tree thereafter)." << endl << endl;
	}
	else
	{
		cout << "Read in " << userTrees.size() << " user trees." << endl;
	}
}

// *** pass in taxonNames and (if present) taxonCoding; map that shit
vector < vector <bool> > parseTree (string & tree, vector <string> const& taxonNames,
	vector <int> & translationTable, vector < vector <int> > & treeTaxonOrdering)
{
	vector < vector <bool> > formattedTree;
	bool complete = false;
	string curentName;
	int numTaxa = taxonNames.size();
	
// Hmm... Should I use Ape-style tree storage?
	
	vector < vector <int> > nodes;
	vector <bool> activeNode;
	int currentNode = -1; // this should be the largest open node
	string currentString;
	string comment;
	int numNodes = 0;
	int countTaxa = -1; // map to ordering in taxonNames
	vector < vector <int> > activeTaxa;
	
	vector <string> edgeLengths;
	vector <string> comments;
	vector <string> tipNames;
	vector <int> indexPosition;
	
	for (string::size_type iterator = 0; iterator != tree.size(); ++iterator)
	{
		char currentChar = tree[iterator];
		if (currentChar == '(')
		{
			currentNode = numNodes;
			numNodes++;
			
			vector <int> dummy;
			activeTaxa.push_back(dummy); // yeah, i know, ugly...
			activeNode.push_back(true);
			if (DEBUG) {cout << "Opening node #" << currentNode << endl;}
			continue;
		}
		else if (currentChar == ')') // close node
		{
			activeNode[currentNode] = false;
			if (DEBUG) {cout << "Closing node #" << currentNode << endl;}
			for (int iterActive = 0; iterActive < (int)activeNode.size(); iterActive++) // find largest open node
			{
				if (activeNode[iterActive])
				{
					currentNode = iterActive;
				}
			}
			continue;
		}
		else if (currentChar == '[') // comment. just fuggeddaboudit
		{
			complete = false;
			comment = "[";
			while (!complete)
			{
				iterator++;
				char temp = tree[iterator];
				if (temp == ']')
				{
					comment += temp;
					complete = true;
				}
				else
				{
					comment += temp;
				}
			}
//			comments.push_back(comment);
			comment.clear();
			continue;
		}
		else if (currentChar == ':') // edge length
		{
			complete = false;
			iterator++;
			char temp = tree[iterator];
			string edgeLength(1, temp);
			while (!complete)
			{
				iterator++;
				temp = tree[iterator];
				if (temp == ',' || temp == ')' || temp == '[' || temp == ';')
				{
					iterator--;
					complete = true;
				}
				else
				{
					edgeLength += temp;
				}
			}
//			cout << "Processed edge length '" << edgeLength << "'." << endl;
			edgeLengths.push_back(edgeLength);
			edgeLength.clear();
			continue;
		}
		else if (currentChar == ',')
		{
			continue;
		}
		else if (currentChar == ';')
		{
			if (DEBUG) {cout << "End of tree description." << endl;}
		}
		else if (currentChar != ',') // tip name; could be an integer if a translation table is being used.
		{
			complete = false;
			countTaxa++;
			
			for (int nodeIter = 0; nodeIter < (int)activeNode.size(); nodeIter++)
			{
				if (activeNode[nodeIter]) // is node active?
				{
					activeTaxa[nodeIter].push_back(countTaxa); // then add mother-fucking current taxon to it's list of taxa
				}
			}
			
			char temp = tree[iterator];
			string tipName(1, temp);
			while (!complete)
			{
				iterator++;
				temp = tree[iterator];
				if (temp == ',' || temp == ')' || temp == '[' || temp == ':')
				{
					complete = true;
					iterator--;
				}
				else
				{
					tipName += temp;
				}
			}
//			cout << "Processed tip name '" << tipName << "'." << endl;
			
			bool match = false;
// Check against taxonNames list first case: no translation table
			if (translationTable.size() == 0) // first case: no translation table
			{
				for (int taxaIter = 0; taxaIter < (int)taxonNames.size(); taxaIter++)
				{
					if (tipName == taxonNames[taxaIter])
					{
						indexPosition.push_back(taxaIter);
						match = true;
					}
				}
				if (!match) // check for mis-spelling of taxon names
				{
					cout << "ERROR: Tree taxon '" << tipName << "' not found in input matrix. Exiting." << endl << endl;
					exit(1);
				}
			}
// Okay, translation table is present...
			else
			{
				int taxonCode = convertStringtoInt(tipName);
				if (taxonCode > numTaxa || taxonCode < 1)
				{
					cout << "ERROR: Tree taxon index '" << taxonCode << "' not within valid range of active taxa. Exiting." << endl << endl;
					exit(1);
				}
				else // this be where fuck-ups exist
				{
					indexPosition.push_back(translationTable[taxonCode - 1]);
				}
			}
			tipNames.push_back(tipName);
			tipName.clear();
			continue;
		}
		else
		{
			cout << "Huh? What be this then? Don't recognize '" << currentChar << "'. Shit." << endl;
		}
	}
	
	treeTaxonOrdering.push_back(indexPosition);
	
	if (DEBUG)
	{
		cout << "Encountered " << numNodes << " nodes." << endl;
		cout << "Encountered " << countTaxa + 1 << " taxa.." << endl;
	}
	
	
// ok. so we have the tree parsed. now put it into binary tree format... fucker.
	for (int nodeIter = 0; nodeIter < numNodes; nodeIter++)
	{
// vector < vector <int> > activeTaxa <- lists active taxa for each node. 
		vector <bool> node;
		int numActive = activeTaxa[nodeIter].size();
		for (int taxaIter = 0; taxaIter < numTaxa; taxaIter++)
		{
			bool match = false;
			for (int activeIter = 0; activeIter < numActive; activeIter++)
			{
				if (activeTaxa[nodeIter][activeIter] == taxaIter)
				{
					match = true;
				}
			}
			node.push_back(match);
		}
		formattedTree.push_back(node);
	}
	
// annoyingly, functions downstream are expecting the singleton clades to be present as well. gah.
	for (int taxaIter = 0; taxaIter < numTaxa; taxaIter++)
	{
		vector <bool> node;
		for (int matchIter = 0; matchIter < numTaxa; matchIter++)
		{
			if (matchIter == taxaIter)
			{
				node.push_back(true);
			}
			else
			{
				node.push_back(false);
			}
		}
		formattedTree.push_back(node);
	}

// reverse the order of clades as the edge-reader function is expecting coalescent-ordering
	reverse(formattedTree.begin(), formattedTree.end());

// test!
	if (DEBUG) {printTree(formattedTree);}
	if (DEBUG) {cout << endl;}
	
	return formattedTree;
}

// extract coding from tranlation table (if present); map to alignment ordering
vector <int> collectTaxonCoding (vector <string> & rawInput, vector <string> const& taxonNames)
{
	bool translate = false;
	bool complete = false;
	bool semicolonEncountered = false;
	bool commaEncountered = false;
	vector <int> taxonCoding;
	
	vector <string> storeNames;	// for debugging
	vector <int> storeIndices;	// for debugging
	
	int numTaxa = taxonNames.size();
	int countTaxa = -1;
	
// Looking for 'translate' command; will almost certainly be there
	for (vector <string>::iterator lineIter = rawInput.begin(); lineIter < rawInput.end(); lineIter++)
	{		
		int index = 0;
		int stringPosition = 0;
		
		if (!translate)
		{
			if (checkStringValue(*lineIter, "translate", 0))
			{
				translate = true;
				vector <int> foo (numTaxa, 0);
				taxonCoding = foo; // initialize the MF at all 0s
				continue;
			}
			else // garbage
			{
				continue;
			}
		}
		else if (translate && !complete)
		{
// Format: index_value taxon_name(,/;/ )
			if (checkStringValue(*lineIter, ";", 0)) // ";" alone; done
			{
				semicolonEncountered = true;
				complete = true;
				continue;
			}
			else
			{
				countTaxa++;
				index = convertStringtoInt(extractStringElement (*lineIter, stringPosition));
				index--; // translation tables start indexing at 1
				
				storeIndices.push_back(index);
				
				string tipName = removeStringSuffix(extractStringElement(*lineIter, 1), ';', semicolonEncountered);
				if (semicolonEncountered)
				{
					complete = true;
				}
				else
				{
					tipName = removeStringSuffix(extractStringElement(*lineIter, 1), ',', commaEncountered);
				}
				storeNames.push_back(tipName);
				if (DEBUG) {cout << "Stored index #" << index << " and name '" << tipName << "'."<< endl;}
// Check against taxonNames list
				bool match = false;
				for (int taxaIter = 0; taxaIter < numTaxa; taxaIter++)
				{
					if (tipName == taxonNames[taxaIter])
					{
						taxonCoding[index] = taxaIter;
						match = true;
						if (DEBUG) {cout << "Matched tip '" << tipName << "' at index #" << index << " to taxon '"
							<< taxonNames[taxaIter] << "' which is in " << taxaIter << " place in the alignment." << endl;}
					}
				}
				if (!match)
				{
					cout << "ERROR: Tree taxon '" << tipName << "' not found in input matrix. Exiting." << endl << endl;
					exit(1);
				}
				if (DEBUG) {cout << "Encountered tree taxon '" << tipName << "' with index = " << index << endl;}
			}
		}
	}
	
	if (DEBUG) // print out mapping of translation table to alignment
	{
		if (taxonCoding.size() != 0)
		{
			for (int taxonIter = 0; taxonIter < numTaxa; taxonIter++)
			{
				if (numTaxa >= 99)
				{
					if (taxonIter < 99)
					{
						cout << " ";
					}
					if (taxonIter <= 9)
					{
						cout << " ";
					}
				}
				else if (numTaxa >= 9 && taxonIter <= 9)
				{
					cout << " ";
				}
				cout << "Taxon #" << taxonIter << " '" << taxonNames[taxonIter] << "' corresponds to the #"
				<< taxonCoding[taxonIter] << " taxon in the tree '" << storeNames[taxonCoding[taxonIter]] << "'." << endl;
			}
		}
	}
	
	
	cout << endl;
	return taxonCoding;
}

void writeAnnotatedTrees (vector <string> const& rawTrees, vector <int> & translationTable,
	vector < vector <double> > userTreeDecisiveness, vector <string> const& taxonNames)
{
// may or may not have a root-type declaration and/or probability:
//	tree rep.1 = ((((((((((((((((((((((57:0.100000
//	tree TREE1  = [&R] ((((1
//	tree STATE_0 [&lnP=-52901.05871008702] = [&R] (((1
	
	ofstream annotated_trees;
	
	string outTrees = "Decisivator.trees";
	checkValidOutputFile (outTrees);
	
	
	annotated_trees.open(outTrees.c_str());
	
	int numTrees = rawTrees.size();
	vector <double> currentDecisiveness;
	string tree;
	int annotationCount = 0;
	
	cout << "Preparing to print out " << numTrees << " trees." << endl;
	bool complete = false;
	
	annotated_trees << "#NEXUS" << endl << endl;
	annotated_trees << "Begin trees;" << endl;
	
// *** check for presence of a translation table ***
	if (translationTable.size() != 0)
	{
		cout << endl << "We've got a translation table to print out here." << endl;
/*
   translate
       1 Allosaurus_fragilis,
       2 Sinraptor,
       3 Dilong_paradoxus,
       4 Eotyrannus_lengi,
       5 Tyrannosaurus_rex,
       6 Gorgosaurus_libratus,
       7 Tanycolagreus_topwilsoni,
       8 Coelurus_fragilis,
       9 Ornitholestes_hermanni,
      10 Huaxiagnathus_orientalis,
      11 Sinosauropteryx_prima,
      12 Deinocheirus_mirificus,
      13 Harpymimus_okladnikovi,
      14 Pelecanimimus_polyodon,
      15 Shenzhousaurus_orientalis,
      16 Archaeornithomimus_asiaticus,
      87 Microraptor_zhaoianus,
      88 NGMC91_unnameddromaeosaurid,
      89 Buitreraptor_gonzalezorum;
*/
		annotated_trees << "Translate" << endl;
		
		for (int i = 0; i < (int)taxonNames.size(); i++)
		{
			annotated_trees << i + 1 << " " << taxonNames[translationTable[i]];
			if (i != ((int)taxonNames.size() - 1))
			{
				annotated_trees << "," << endl;
			}
			else
			{
				annotated_trees << ";" << endl;
			}
		}
		annotated_trees << endl;
	}
	
// may or may not have a root-type declaration and/or probability:
//	tree rep.1 = ((((((((((((((((((((((57:0.100000
//	tree TREE1  = [&R] ((((1
//	tree STATE_0 [&lnP=-52901.05871008702] = [&R] (((1
	
	for (int i = 0; i < numTrees; i++)
	{
		currentDecisiveness = userTreeDecisiveness[i];
		
		tree = rawTrees[i];
		reverse(currentDecisiveness.begin(), currentDecisiveness.end()); // this was flipped during searching. ugh. li.
		
		if (DEBUG) {cout << endl << tree << endl << endl;}
		
		bool start = false;
		int position = 0;
		while (!start)
		{
			string temp = extractStringElement(tree, position);
			if (temp[0] == '(')
			{
				start = true;
				tree = temp;
			}
			else
			{
				cout << temp << " ";
				annotated_trees << temp << " ";
				position++;
			}
		}
		
		for (int iterator = 0; iterator != (int)tree.size(); ++iterator)
		{
			char currentChar = tree[iterator];
			if (currentChar == '(')
			{
				cout << currentChar;
				annotated_trees << currentChar;
				
				continue;
			}
			else if (currentChar == ')') // close node; check if branch lengths or annotations exist
			{
				cout << currentChar;
				annotated_trees << currentChar;
				
				iterator++;
				
				currentChar = tree[iterator]; // end of tree; GET OUTTA THERE!!!
				if (currentChar == ';')
				{
					cout << currentChar << endl;
					annotated_trees << currentChar << endl;
					continue;
				}
				else if (currentChar == ',' || currentChar == ')') // no edgelengths present; simply a topology
				{
					cout << "[&decisiveness=" << currentDecisiveness[annotationCount] << "]";
					annotated_trees << "[&decisiveness=" << currentDecisiveness[annotationCount] << "]";
					annotationCount++;
					iterator--;
				}
				else if (currentChar == ':') // edge length; need to look for '[' following
				{
					cout << currentChar;
					annotated_trees << currentChar;
					complete = false;
					while (!complete)
					{
						iterator++;
						currentChar = tree[iterator];
						if (currentChar == ',' || currentChar == ')')
						{
							cout << "[&decisiveness=" << currentDecisiveness[annotationCount] << "]";
							annotated_trees << "[&decisiveness=" << currentDecisiveness[annotationCount] << "]";
							annotationCount++;
							iterator--;
							complete = true;
						}
						else if (currentChar == '[')
						{
							complete = false;
							cout << currentChar;
							annotated_trees << currentChar;
							while (!complete)
							{
								iterator++;
								currentChar = tree[iterator];
								if (currentChar == ']')
								{
									complete = true;
									cout << ",decisiveness=" << currentDecisiveness[annotationCount] << "]";
									annotated_trees << ",decisiveness=" << currentDecisiveness[annotationCount] << "]";
								}
								else
								{
									cout << currentChar;
									annotated_trees << currentChar;
								}
							}
							continue;
						}
						else
						{
							cout << currentChar; // print out edge length
							annotated_trees << currentChar;
						}
					}
				}
				else if (currentChar == '[') // should work whether edge lengths are present or not
				{
					complete = false;
					cout << currentChar;
					annotated_trees << currentChar;
					while (!complete)
					{
						iterator++;
						currentChar = tree[iterator];
						if (currentChar == ']')
						{
							complete = true;
							cout << ",decisiveness=" << currentDecisiveness[annotationCount] << "]";
							annotated_trees << ",decisiveness=" << currentDecisiveness[annotationCount] << "]";
						}
						else
						{
							cout << currentChar;
							annotated_trees << currentChar;
						}
					}
					continue;
				}
				else // no edgelengths present; simply a topology
				{	
					cout << "[&decisiveness=" << currentDecisiveness[annotationCount] << "]";
					annotated_trees << "[&decisiveness=" << currentDecisiveness[annotationCount] << "]";
					annotationCount++;
				}
			}
			else if (currentChar == '[') // annotation
			{
				complete = false;
				cout << currentChar;
				annotated_trees << currentChar;
				while (!complete)
				{
					iterator++;
					currentChar = tree[iterator];
					if (currentChar == ']')
					{
						complete = true;
						cout << ",&decisiveness=" << currentDecisiveness[annotationCount] << "]";
						annotated_trees << ",&decisiveness=" << currentDecisiveness[annotationCount] << "]";
					}
					else
					{
						cout << currentChar;
						annotated_trees << currentChar;
					}
				}
				continue;
			}
			else if (currentChar == ':') // edge length
			{
				cout << currentChar;
				annotated_trees << currentChar;
				complete = false;
				while (!complete)
				{
					iterator++;
					currentChar = tree[iterator];
					if (currentChar == ',' || currentChar == ')' || currentChar == '[')
					{
						iterator--;
						complete = true;
					}
					else
					{
						cout << currentChar;
						annotated_trees << currentChar;
					}
				}
				continue;
			}
			else if (currentChar == ',')
			{
				cout << currentChar;
				annotated_trees << currentChar;
				continue;
			}
			else if (currentChar == ';')
			{
				if (DEBUG) {cout << "End of tree description." << endl;}
				cout << currentChar << endl;
				annotated_trees << currentChar << endl;
			}
			else // tip name; could be an integer if a translation table is being used.
			{
				complete = false;
				cout << currentChar;
				annotated_trees << currentChar;
				
				currentChar = tree[iterator];
				while (!complete)
				{
					iterator++;
					currentChar = tree[iterator];
					if (currentChar == ',' || currentChar == ')' || currentChar == '[' || currentChar == ':' || currentChar == '(')
					{
						complete = true;
						iterator--;
					}
					else
					{
						cout << currentChar;
						annotated_trees << currentChar;
					}
				}
			}
		}
		tree.clear();
	}
	annotated_trees << "End;";
}