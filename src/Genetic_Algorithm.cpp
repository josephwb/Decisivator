// Daniel Beck
// 11/22/2011
// Genetic algorithm to search for combination of genes that maximizes decisiveness

// g++ -Wall -m64 -O3 GA.cpp -o GA

#include <iostream>
#include <vector>
#include <ctime>

using namespace std;

#include "Genetic_Algorithm.h"


// GA function accepts data matrix, legal gene matrix, and number of total genes wanted in final matrix
int optimizeDecisivenessGA (vector < vector <int> > & data,
	vector < vector <int> > const& legal, int maxNumGenes)
{
	int numInitialLoci = countLociPresent(data);
	
	
// Pass these two values in
// size of tournament group. must be 2 or greater.
	int tournamentSize = 5;
// make this data size specific?
	int populationSize = 100;
	
	
	
	double mutationProbability = 0.9;  // probability of mutation, otherwise xOver
	int maxNumGenerations = 1000;  // maximum number of generations
	int numAdd = maxNumGenes - numInitialLoci; // number of genes to add
	double fitness[populationSize]; // Store fitness

// Function currently errors out when too many genes initially present.
// This should be changed to allow genes to be removed.
	if (numAdd < 0)
	{
		cout << "No Genes can be added." << endl;
		return(1);
	}
	
// Initialize population
	cout << "Initializing population..." << endl;
	vector < vector <int> > population[populationSize];
	
// Set up population
	for (int i = 0; i < populationSize; i++)
	{ 
		population[i] = data; // take data as input, mutate
		for (int j = 0; j < numAdd; j++)
		{
			mutate(population[i], legal, 1, 0);
		}
	}
	
	cout << "Calculating fitness..." << endl;
// Calculate initial fitness
	for (int i = 0; i < populationSize; i++)
	{
		fitness[i] = calculateFitness(population[i]);
	}

	bool done = 0;		// flag to end evolution loop
	int genCounter = 0; // counter for number of generations

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
	cout << "Starting evolution..."<< endl;
	while (!done)
	{
		genCounter++;
	
// Choose tournament group
		int best;		// integer to hold position of best member of tournament group
		int worst;		// position of worst member of tournament group
		double bestFit = 0; 	// best fitness
		double worstFit = 1;	// worst fitness
		
// compare tournamentSize number of random individuals to current best fit
		for (int i = 0; i < tournamentSize; i++)
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
		population[worst] = population[best];

// Mutate with some probability, otherwise xOver
		double randmu = (double) rand() / RAND_MAX;
		if (randmu < mutationProbability) // mutate
		{			
			mutate(population[worst], legal, 1, 1);
		}
		else // xOver with random individual in population
		{
			xOver(population[worst],population[(int)(rand()%(populationSize-1))]);
		}
		
// Calculate fitness of new individual
		fitness[worst] = calculateFitness(population[worst]);
		
// Check to see if max number of generations has been reached
		if (genCounter == maxNumGenerations)
		{
			done = true;
		}
	} 	// End of evolution loop

	double bestFit = 0;
	int bestInd = 0;
	for (int i = 0; i < populationSize; i++)
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
	int numTaxa = (int)data.size();
	int numLoci = (int)data[0].size();

	int count = 0;

	for (int i = 0; i < numTaxa; i++)
	{ 
		for (int j = 0; j < numLoci; j++)
		{ 
			if (data[i][j])
			{
				count = count + 1;
			}
		}
	}
	return(count);
}

// Function to mutate individual. Currently this function calculates number of legal positions each time.
// This can be avoided by calculating this once and passing these values to this function.
// add and remove are flags to determine if gene is added and/or removed.
void mutate (vector < vector <int> > & data, vector < vector <int> > const& legal,
	bool add, bool remove)
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
				if (data[i][j])
				{
					count1 = count1 + 1;
				}
				else
				{
					count0 = count0 + 1;
				}
			}
		}
	}
	
// select random currently absent gene to add
	int randomAdd = 0;
	if (add)
	{
		randomAdd = rand() % count0 + 1;
	}

// select random currently present gene to remove
	int randomRemove = 0;
	if (remove)
	{
		randomRemove = rand() % count1 + 1;
	}

// modify data, adding random gene, removing random gene
	count0 = 0;
	count1 = 0;
	
	for (int i = 0; i < numTaxa; i++)
	{ 
		for (int j = 0; j < numLoci; j++)
		{ 
			if (legal[i][j])
			{
				if (data[i][j])
				{
					count1 = count1 + 1;
				}
				else
				{
					count0 = count0 + 1;
				}
				if (remove)
				{
					if (count1 == randomRemove)
					{
						data[i][j] = 0;
					}
				}
				if (add)
				{
					if (count0 == randomAdd)
					{
						data[i][j] = 1;
					}
				}
			}
		}
	}
}

// Function to cross-over two individuals. Can take out one loop to make faster
void xOver (vector < vector <int> > & recipient, vector < vector <int> > const& donor)
{
	int numTaxa = (int)recipient.size();
	int numLoci = (int)recipient[0].size();

	int lower = rand()%(numLoci);
	int upper = lower + rand()%(numLoci - lower);
	for (int i = 0; i < numTaxa; i++)
	{ 
		for (int j = lower; j <= upper; j++)
		{ 
			recipient[i][j] = donor[i][j];
		}
	}
}

// Function to calculate fitness of individual. 
// Decisiveness caclulation should go here. Fitness is expected to be between 0 and 1 with 1 being optimum. 
double calculateFitness (vector < vector <int> > const& individual)
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

/* 
Temporary main to error check GA function
int main()
{	
	srand((unsigned)time(0));	

	const int nTax = 10;
	const int nLoci = 20;
	const int numAdd = 50;
	const int maxNumGenes = 100;

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
	
simulate some data
	cout << "Mutate" << endl;
	for (int j = 0; j<numAdd; j++)
	{
		mutate(data, legal, 1, 0);
	}

print data
	for (int i = 0; i < nTax; i++)
	{
		for (int j = 0; j < nLoci; j++)
		{
			cout << data[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

all non-data positions are legal
	for (int i = 0; i < nTax; i++)
	{
		for (int j = 0; j < nLoci; j++)
		{
			legal[i][j]=!data[i][j];
		}
	}

	cout << "GA" << endl;	
	optimizeDecisivenessGA(data, legal, maxNumGenes);
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

 */
