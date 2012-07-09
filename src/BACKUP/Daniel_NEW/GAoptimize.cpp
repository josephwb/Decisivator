// Daniel Beck
// 11/22/2011
// Genetic algorithm (and other optimization functions) to search for combination of genes that maximizes decisivness

// Fucked up by JWB

#include "GAoptimize.h"

// GA function accepts data matrix, legal gene matrix, and number of total genes wanted in final matrix
	// change to total number of genes added (more general)
int optimizeDecisivnessGA (vector < vector <int> > & data,
	vector < vector <int> > const& legal, int maxGenes) // get rid of maxGenes. use: int & numAdd
{
// initialize random number generator - already done in main
	// srand((unsigned)time(0));

// Constants. Some should probably be user specified: tournamentSize, populationSize, mutationProbability, maxgen
	// int numTaxa = data.size(); // not used
	// int numLoci = data[0].size(); // not used
	int numInitialLoci = countLociPresent(data);
	int tournamentSize = 5; // size of tournament group.
	const int populationSize = 100;
	double mutationProbability = 0.9; // probability of mutation, otherwise xover
	int maxgen = 100000; // maximum number of generations
	
// Store fitness
	double fitness[populationSize];
	
	// cout << "initialize " << endl;
// Initialize population
	vector < vector <int> > population[populationSize];
	int numAdd = maxGenes - numInitialLoci; 	// number of genes to add
	
// Function currently errors out when too many genes initially present.
// This should be changed to allow genes to be removed.
	if (numAdd < 0) // subtract genes
	{
		cout << "Genes will be removed" << endl;
		for (int i = 0; i < populationSize; i++)
		{ 
			population[i] = data;
			for (int j = numAdd; j < 0; j++)
			{
				mutate(population[i], legal, 0, 1);
			}
		}
	}
	else // add genes
	{
		cout << numAdd << " genes will be added" << endl;
		for (int i = 0; i < populationSize; i++)
		{ 
			population[i] = data;
			for (int j = 0; j < numAdd; j++) // Randomly add numAdd genes in legal positions
			{
				mutate(population[i], legal, 1, 0); // last 2 arguments are add, remove (boolean)
			}
		}
	}

	// cout << "calcfit" << endl;
// Calculate initial fitness. For each matrix, crunch through 1000 trees
	for (int i = 0; i < populationSize; i++)
	{
		fitness[i] = calcFit(population[i]);
	}

	bool done = 0;		// flag to end evolution loop
	int gen = 0; 		// counter for number of generations

// Check for out of bounds parameter values
	if (populationSize < 3)
	{ 
		cout << "Population size must be greater than 2" << endl;
		return(1);
	}
	if (tournamentSize < 2)
	{
		cout << "Tournament group size must be greater than 1" << endl;
		return(1);
	}

// Evolution loop
	cout << "Start evolution" << endl;
	// cout << "Progress:" << endl;
	while (!done)
	{
		gen = gen + 1;
	//	cout << ".";
// Choose tournament group
// initialize best and worst with random
		int random = (int)(rand()%(populationSize));
		int best = random;		// integer to hold position of best member of tournament group
		int worst = random;		// position of worst member of tournament group
		double bestFit = fitness[random];	// best fitness
		double worstFit = fitness[random];	// worst fitness
		
		for (int i = 1; i < tournamentSize; i++)
		{
			random = (int)(rand()%(populationSize));
			if (fitness[random] > bestFit)
			{
				bestFit = fitness[random];
				best = random;
			}
			if (fitness[random] < worstFit)
			{
				worstFit = fitness[random];
				worst = random;
			}
		}
		
// make sure best and worst are not the same
// if they are the same, choose new best randomly from population
		while (best == worst)
		{
			best = (int)(rand()%(populationSize));
		}

// Copy best individual to worst individual
		population[worst] = population[best];
		
// Mutate with some probability, otherwise xover
		double randmu = (double) rand() / RAND_MAX;
		if (randmu < mutationProbability)
		{
// mutate			
			mutate(population[worst], legal, 1, 1);
		}
		else
		{
// xover with random individual in population
			int randx = (int)(rand()%populationSize);
			while (randx == worst)
			{
				randx = (int)(rand()%populationSize);
			}
			xover(population[worst],population[randx], legal, maxGenes);
		}
		
// Calculate fitness of new individual	
		fitness[worst]=calcFit(population[worst]);
// Check to see if max number of generations has been reached
		if (gen == maxgen)
		{
			done = 1;
		}
	} 
// End of evolution loop

	double bestFit = 0;
	int bestInd = 0;
	for (int i = 0; i < populationSize; i++)
	{
		if (fitness[i] > bestFit) {
			bestFit = fitness[i];
			bestInd = i;
		}
	}
// copy best evolved matrix to original data matrix
	data = population[bestInd];
	cout << endl;
	
	return(0);
}


// Function to count number of loci present in matrix
// Actually, number of cells that are occupied
int countLociPresent (vector < vector <int> > const& data)
{
	int numTaxa = data.size();
	//int numLoci = data[0].size();

	int count = 0;
	
/*
	for (int i = 0; i < numTaxa; i++)
	{ 
		for (int j = 0; j < numLoci; j++)
		{ 
			if (data[i][j])
			{
				count += 1;
			}
		}
	}
*/
	     
	//cout << "count = " << count;
	
// faster option:
	for (int i = 0; i < numTaxa; i++)
	{
		count += accumulate(data[i].begin(), data[i].end(), 0);
	}

	return(count);
}


/* Function to mutate individual. Currently this function calculates number of legal
positions each time. This can be avoided by calculating this once and passing these
values to this function. add and remove are flags to determine if gene is added and/or removed. */
void mutate (vector < vector <int> > & data, vector < vector <int> > const& legal,
	bool const& add, bool const& remove)
{
	int numTaxa = data.size();
	int numLoci = data[0].size();

	int count0 = 0;		// number of legal absent positions
	int count1 = 0;		// number of legal present positions

// count number of legal positions
// probably a faster way to do this
	for (int i = 0; i < numTaxa; i++)
	{ 
		for (int j = 0; j < numLoci; j++)
		{ 
			if (legal[i][j])
			{
				if (data[i][j])
				{
					count1++;
				}
				else
				{
					count0++;
				}
			}
		}
	}	
	
// select random currently absent gene to add
	int radd=0;
	if (add)
	{
		if (!count0)
		{
			cout << "All legal loci full, cannot add" << endl;
			return;
		}
		else
		{
			radd = rand() % count0;
			cout << "Randomly adding " << radd << " taxon-genes." << endl;
		}
	}
// select random currently present gene to remove
	int rremove = 0;
	if (remove)
	{
		if (!count1)
		{
			cout << "All legal loci empty, cannot remove" << endl;
			return;
		}
		else
		{
			rremove = rand() % count1;
			cout << "Randomly removing " << rremove << " taxon-genes." << endl;
		}
	}
		
// modify data, adding random gene, removing random gene
	count0 = -1;
	count1 = -1;

	for (int i = 0; i < numTaxa; i++)
	{ 
		for (int j = 0; j < numLoci; j++)
		{ 

			if (legal[i][j])
			{
				if (data[i][j])
				{
					count1++;
				}
				if (!data[i][j])
				{
					count0++;
				}

// after finding correct position, change value
				if (remove)
				{
					if (count1 == rremove)
					{
						data[i][j] = 0; rremove = -1;
					}
				}
				if (add)
				{
					if (count0 == radd)
					{
						data[i][j] = 1; radd = -1;
					}
				}
// rremove and radd set to -1 after change to prevent multiple changes
			}
		}
	}
}


// Function to cross-over two individuals. Can take out one loop to make faster
void xover ( vector < vector <int> > & recipient, vector < vector <int> > const& donor,
	vector < vector <int> > legal, int maxGenes)
{
	int numTaxa = recipient.size();
	int numLoci = recipient[0].size();

	int lower = rand()%(numLoci);
	int upper = lower + rand()%(numLoci-lower);
	for (int i = 0; i < numTaxa; i++)
	{ 
		for (int j = lower; j <= upper; j++)
		{ 
			recipient[i][j] = donor[i][j];
		}
	}
// randomly add or remove genes if too many or too few
	int numAdd = maxGenes - countLociPresent(recipient);	
	if (numAdd > 0)
	{
		for (int i = 0; i < numAdd; i++)
		{
			mutate(recipient, legal, 1, 0);
		}
	}
	if (numAdd < 0)
	{
		for (int i=numAdd; i<0; i++)
		{
			mutate(recipient, legal, 0, 1);
		}
	}
}


// Function to calculate fitness of individual. 
// Decisiveness caclulation should go here. Fitness is expected to be between 0 and 1 with 1 being optimum.
double calcFit( vector < vector <int> > const& individual)
{
// double fit = calculateDecisivness(individual);
	
	int numTaxa = individual.size();
	int numLoci = individual[0].size();
	
	int sumT[numTaxa];
	int sumL[numLoci];
	for (int i = 0; i < numTaxa; i++) {sumT[i]=0;} // must be a better way to initialize this
	for (int i = 0; i < numLoci; i++) {sumL[i]=0;}
	
//	int sumT[numTaxa, 0];
//	int sumL[numLoci, 0];	
	
	for (int i = 0; i < numTaxa; i++)
	{
		for (int j = 0; j < numLoci; j++)
		{
			sumT[i] = sumT[i] + individual[i][j];
			sumL[j] = sumL[j] + individual[i][j];
		}
	}
	
	double prodT = 1;
	double prodL = 1;
	for (int i = 0; i < numTaxa; i++)
	{
		if (sumT[i] >= 1)
		{
			prodT = prodT * sumT[i];
		}
	}
	for (int i = 0; i < numLoci; i++)
	{
		if (sumL[i] >= 1)
		{
			prodL = prodL*sumL[i];
		}
	}
	double fit = 0;
	//cout << "prodL: " << prodL << " prodT: " << prodT << endl; 
	if ((prodT + prodL) > 0)
	{
		fit = prodT/(prodT + prodL);
	}
	if ((fit > 1) || (fit < 0))
	{
		cout << " fitness error " << endl;
		return(1);
	}

	return (fit);
}		

// Functions to add genes one by one to maximize decisivness. 
// 12/7/2011

// Stepwise addition function. Sequentially adds genes at best positions
int stepAdd (vector < vector < int > > & data, vector < vector < int > > const& legal, int maxGenes)
{
	int numTaxa = data.size();
	int numLoci = data[0].size();
	int numInitialLoci = countLociPresent(data);

	int bestTax;	// store best taxa
	int bestLoci;	// store best loci
	double bestDecisive=calcFit(data); // store best decisiveness

// Add until reaching the maximum number of genes
	for (int l = 0; l < (maxGenes-numInitialLoci); l++)
	{
// test all possible positions
		for (int i = 0; i < numTaxa; i++)
		{ 
			for (int j = 0; j < numLoci; j++)
			{ 
				if (legal[i][j] & !data[i][j])
				{
					vector < vector < int > > temp = data;
					temp[i][j] = 1;
					double tempDecisive = calcFit(temp);
					if (tempDecisive >= bestDecisive)
					{
						bestDecisive = tempDecisive;
						bestTax = i;
						bestLoci = j;
					}
				}
			}
		}
		data[bestTax][bestLoci] = 1;
	}
	return(0);
}

// 12/10/2011
/* Hill climbing approach similar to stepAdd. In this function genes are added randomly until
maxGenes reached. The genes are then sequentially moved from current location to best location
one by one. This process repeats untill no moves increases fitness. Can be stochastic result,
depends on starting condition (initialized randomly). This could result in an infinite loop.
TODO: set maximum number of loops */
int stepReplace (vector < vector < int > > & data, vector < vector < int > > const& legal,
	int maxGenes)
{
	int maxLoops = 100000; 	// maximum number of iterations through while loop. This would be better as a user defined parameter.
	int numTaxa = data.size();
	int numLoci = data[0].size();
	int numInitialLoci = countLociPresent(data);
	
	// fill matrix randomly 
	int numAdd=maxGenes-numInitialLoci; // number of genes to add

	if (numAdd < 0) // subtract genes
	{
		cout << "Genes will be removed" << endl;
		for (int j = numAdd; j < 0; j++)
		{
			mutate(data, legal, 0, 1);
		}

	}
	else // add genes
	{
		for (int j = 0; j < numAdd; j++)
		{
			mutate(data, legal, 1, 0);
		}
	}

	bool done = 0; 	// flag for end of loop
	int loops = 0;  // counter for number of times through while loop
	while (!done)
	{
		loops++;
		int moves = 0;
		for (int taxa = 0; taxa < numTaxa; taxa++)
		{
			for (int loci = 0; loci < numLoci; loci++)
			{
				if (data[taxa][loci] & legal[taxa][loci]) {	// if loci present and movable
					data[taxa][loci] = 0;	// remove loci

					int bestTax;	// store best taxa
					int bestLoci;	// store best loci
					double bestDecisive = calcFit(data); // store best decisiveness

					// test all possible positions
					for (int i = 0; i < numTaxa; i++)
					{ 
						for (int j = 0; j < numLoci; j++)
						{ 
							if (legal[i][j] & !data[i][j])
							{
								vector < vector < int > > temp = data;
								temp[i][j] = 1;
								double tempDecisive = calcFit(temp);
								if (tempDecisive >= bestDecisive)
								{
									bestDecisive = tempDecisive;
									bestTax = i;
									bestLoci = j;
								}
							}
						}
					}
					data[bestTax][bestLoci]=1; 	// add loci to best position
					if ((bestTax != taxa)||(bestLoci!=loci))
					{
						moves++;
					}
				}
			}
		}
// done when no moves made
		if (moves < 1) done = 1; 
// check for maximum number of iterations
		if (loops >= maxLoops) { cout << "Max loops reached" << endl; done=1;}
	}
	return(0);
}





// Temporary main to error check GA function

int main() {	

	srand((unsigned)time(0));	
	
	//srand((unsigned)time(0)); // this will be done upstream in main.
	
	const int nTax = 10;
	const int nLoci = 25;
	const int numAdd = 30;
	const int numAdd2 = 10;
	const int maxGenes = 40;

	vector < vector <int> > data;
	vector < vector <int> > legal;
	
	data.resize(nTax, vector <int> (nLoci, 0));
//	printData(data);
	legal.resize(nTax, vector <int> (nLoci, 1));
//	printData(legal);
		
// simulate some data
	for (int j = 0; j < numAdd; j++)
	{
		mutate(data, legal, 1, 0);
	}

	// print data
	cout << "Starting matrix" << " " << calcFit(data) << endl;
	printData(data);

	vector < vector < int > > data2 = data;
	vector < vector < int > > data3 = data;
	vector < vector < int > > data4 = data;
	
	legal = getLegalMatrix(data);
	
	// cout << "Legal matrix" << endl;
	// printData(legal);

// add some more to data, if wanted
	for (int j = 0; j < numAdd2; j++)
	{
		mutate(data, legal, 1, 0);
	}
	
	cout << "Original score:" << " " << calcFit(data)<< endl;
	optimizeDecisivnessGA(data, legal, maxGenes);

	cout << "GA" << " " << calcFit(data)<< endl;	

// stepwise addition
	stepAdd (data2, legal, maxGenes);

	cout << "stepwise addition" << " " << calcFit(data2) << endl;
	printData(data2);
	
// hill climb
	stepReplace (data3, legal, maxGenes);
	cout << "hill climb" << " " << calcFit(data3) << endl;
	printData(data3);

// stepwise addition + hill climb
	stepAdd (data4, legal, maxGenes);
	stepReplace (data4, legal, maxGenes);
	cout << "step + hill climb" << " " << calcFit(data4) << endl;

// print data
	printData(data4);
}

vector < vector <int> > getLegalMatrix (vector < vector <int> > const& data)
{
	vector < vector <int> > legal;
	int numTaxa = data.size();
	int numLoci = data[0].size();
	
	legal.resize(numTaxa, vector <int> (numLoci, 1));
	
//	cout << "First:" << endl;
//	printData(legal);
	for (int i = 0; i < data.size(); i++)
	{
		for (int j = 0; j < data[0].size(); j++)
		{
			legal[i][j] =! data[i][j]; // clever
		}
	}
//	cout << "New:" << endl;
//	printData(legal);
	return(legal);
}

void printData (vector < vector <int> > const& data)
{
	for (int i = 0; i < data.size(); i++)
	{
		for (int j = 0; j < data[0].size(); j++)
		{
			cout << data[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

}