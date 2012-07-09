// Daniel Beck
// 11/22/2011
// Genetic algorithm to search for combination of genes that maximizes decisiveness

// g++ -Wall -m64 -O3 GA.cpp -o GA

#include <iostream>
#include <vector>
//#include <ctime>

#include <cstdlib>


#include "GA.h"

using namespace std;


// GA function accepts data matrix, legal gene matrix, and number of total genes wanted in final matrix
int optimizeDecisivenessGA (vector < vector <int> > & data,
	vector < vector <int> > const& legal, int maxGenes)
{
// initialize random number generator
//	srand((unsigned)time(0));	// don't do this twice in a program
	
// Constants (some should probably be user specified)
//	int numTaxa = (int)data.size();
//	int numLoci = (int)data[0].size();
	int numInitialLoci = countLociPresent(data);
	int tournamentSize = 5; 	// size of tournament group. must be 2 or greater
	int populationSize = 100;
	double mutationProbability = 0.9; 		// probability of mutation, otherwise xover
	int maxgen = 1000;		// maximum number of generations
	
// Store fitness
	double fitness[populationSize];
	
	cout << "initialize " << endl;

// Initialize population
	vector < vector <int> > population[populationSize];
	int numAdd=maxGenes-numInitialLoci; 	// number of genes to add

// Function currently errors out when too many genes initially present.
// This should be changed to allow genes to be removed.
	if (numAdd < 0){ cout << "No Genes can be added." << endl; return(1);}

	for (int i = 0; i < populationSize; i++)
	{ 
		population[i] = data;
		for (int j = 0; j < numAdd; j++)
		{
			mutate(population[i], legal, 1, 0);
		}
	}
	
	cout << "calcfit" << endl;
// Calculate initial fitness
	for (int i = 0; i < populationSize; i++)
	{
		fitness[i]=calcFit(population[i]);
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
	cout << "start evolution"<< endl;
	while (!done)
	{
		gen = gen + 1;
	
// Choose tournament group
		int best;		// integer to hold position of best member of tournament group
		int worst;		// position of worst member of tournament group
		double bestFit=0; 	// best fitness
		double worstFit=1;	 // worst fitness

		for (int i=0; i<tournamentSize; i++)
		{
			int random = (int)(rand()%(populationSize));

			if (fitness[random] >= bestFit)
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
			best = (int)(rand()%(populationSize-1));
		}

// Copy best individual to worst individual
		population[worst]=population[best];

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
			xover(population[worst],population[(int)(rand()%(populationSize-1))]);
		}
		
// Calculate fitness of new individual
		fitness[worst]=calcFit(population[worst]);
		
// Check to see if max number of generations has been reached
		if (gen == maxgen) { done = 1; }
	} 	// End of evolution loop

	double bestFit = 0;
	int bestInd = 0;
	for (int i=0; i<populationSize; i++)
	{
		if (fitness[i] > bestFit)
		{
			bestFit = fitness[i];
			bestInd = i;
		}
	}
// copy best evolved matrix to original data matrix
	data = population[bestInd];

	return(0);
}

// Function to count number of loci present in matrix
int countLociPresent (vector < vector <int> > const& data)
{
	int numTaxa = data.size();
	int numLoci = data[0].size();

	int count = 0;

	for (int i=0; i<numTaxa; i++)
	{ 
		for (int j=0; j<numLoci; j++)
		{ 
			if (data[i][j]){ count = count + 1;}
		}
	}
	return(count);
}

// Function to mutate individual. Currently this function calculates number of legal positions each time.
// This can be avoided by calculating this once and passing these values to this function.
// add and remove are flags to determine if gene is added and/or removed.
void mutate (vector < vector <int> > & data, vector < vector <int> > const& legal,
	bool add, bool remove )
{
	int numTaxa = (int)data.size();
	int numLoci = (int)data[0].size();

	int count0 = 0;		// number of legal absent positions
	int count1 = 0;		// number of legal present positions

// count number of legal positions
	for (int i = 0; i < numTaxa; i++)
	{
		for (int j = 0; j < numLoci; j++)
		{
			if (legal[i][j])
			{
				if (data[i][j]) { count1 = count1 + 1;}
				else { count0 = count0 + 1;}
			}
		}
	}
	
// select random currently absent gene to add
	int radd = 0;
	if (add) { radd = rand() % count0 + 1;}

// select random currently present gene to remove
	int rremove=0;
	if (remove) { rremove =rand() % count1 + 1;}

// modify data, adding random gene, removing random gene
	count0 = 0;
	count1 = 0;
	for (int i = 0; i < numTaxa; i++)
	{ 
		for (int j = 0; j < numLoci; j++)
		{ 

			if (legal[i][j])
			{
				if (data[i][j]) { count1 = count1 + 1;}
				else { count0 = count0 + 1;}

				if (remove) {if (count1 == rremove){ data[i][j] = 0;}}
				if (add) {if (count0 == radd){ data[i][j] = 1;}}
			}
		}
	}

}

// Function to cross-over two individuals. Can take out one loop to make faster
void xover (vector < vector <int> > & recipient, vector < vector <int> > const& donor)
{
	int numTaxa = (int)recipient.size();
	int numLoci = (int)recipient[0].size();

	int lower = rand()%(numLoci);
	int upper = lower + rand()%(numLoci-lower);
	for (int i = 0; i < numTaxa; i++)
	{ 
		for (int j = lower; j <= upper; j++)
		{ 
			recipient[i][j]=donor[i][j];
		}
	}

}

// Function to calculate fitness of individual. 
// Decisiveness caclulation should go here. Fitness is expected to be between 0 and 1 with 1 being optimum. 
double calcFit (vector < vector <int> > const& individual)
{
	int numLoci = (int)individual[0].size();
	int sum = 0;
	
	for (int i = 0; i < numLoci; i++)
	{
		sum = sum + individual[0][i];
	}
	double fit = (double) sum/numLoci;	
	return (fit);
}		

// Temporary main to error check GA function
int main()
{	
	srand((unsigned)time(0));	

	const int nTax = 10;
	const int nLoci = 20;
	const int numAdd = 50;
	const int maxGenes = 100;

	vector < vector <int> > data;
	vector < vector <int> > legal;
	
	data.resize( nTax, vector<int> (nLoci, 0));
	legal.resize( nTax, vector<int> (nLoci, 1));

	for (int i = 0; i < nTax; i++)
	{
		for (int j = 0; j < nLoci; j++)
		{
			legal[i][j]=1;
			data[i][j]=0;
		}
	}
	
// simulate some data
	cout << "Mutate" << endl;
	for (int j = 0; j<numAdd; j++)
	{
		mutate(data, legal, 1, 0);
	}

// print data
	for (int i = 0; i < nTax; i++)
	{
		for (int j = 0; j < nLoci; j++)
		{
			cout << data[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

// all non-data positions are legal
	for (int i = 0; i < nTax; i++)
	{
		for (int j = 0; j < nLoci; j++)
		{
			legal[i][j]=!data[i][j];
		}
	}

	cout << "GA" << endl;	
	optimizeDecisivenessGA(data, legal, maxGenes);
	cout << "GA done" << endl;

	for (int i = 0; i < nTax; i++)
	{
		for (int j = 0; j < nLoci; j++)
		{
			cout << data[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	
	return 0;
}
