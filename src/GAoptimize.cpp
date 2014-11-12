// Daniel Beck
// 11/22/2011
// Genetic algorithm (and other optimization functions) to search for combination of genes that maximizes decisivness

// Fucked up by JWB

#include <iostream>
#include <vector>
#include <cstdlib>
#include <numeric>
#include <fstream>



#include "GAoptimize.h"
#include "Matrix_Scan.h"
#include "General.h"

extern bool debugging;

vector < vector <int> > GAHandler (int const& numAddGA, vector < vector <int> > const& data,
    double & GADecisiveness, bool const& referenceTaxonPresent, int const& numProcs)
{
    vector < vector <int> > GAmatrix = data;
// allowed cells. later will add functionality to preclude hard-to-get cells.
    vector < vector < int > > legal = getLegalMatrix(data);
    
    optimizeDecisivnessGA (GAmatrix, legal, numAddGA, referenceTaxonPresent, numProcs);
        
    return(GAmatrix);
}









// Calculate fitness i.e. decisiveness of a particular matrix.
// Could implement different fitness metrics e.g. average over loci.
double calcFit (bool const& referenceTaxonPresent, int & numTrees,
    vector < vector <int> > const& data, int const& numProcs)
{
    double fit = 0.0;
    bool verbose = false;
    bool includeDiversity = false;
    
// 0 is for !findAll
    fit = calculatePartialDecisiveness(referenceTaxonPresent, numTrees, data, 0, numProcs, verbose);
    
    
    if (includeDiversity) { // should distribution of decisiveness across loci be taken into account?
        double weightedFit = 0.0;
        double harmonicMean = 0.0;
        double inverse_sum = 0.0;
        double weightOverall = 0.5;
        double arithmeticMean = 0.0;
        
        for (int i = 0; i < (int)data[0].size(); i++) {
            double tempFit = calculatePartialDecisivenessSinglePartition(referenceTaxonPresent,
                numTrees, data, 0, i, verbose, numProcs); // set verbose to false
            inverse_sum += (1/tempFit);
            arithmeticMean += tempFit / (double)data[0].size();
        //    cout << "inverse_sum = " << inverse_sum << "." << endl;
            
        }
        harmonicMean = (double)data[0].size() / inverse_sum;
        //cout << "Harmonic mean = " << harmonicMean << "." << endl;
        
        
        weightedFit = fit + harmonicMean; // can this be > 1?
        //weightedFit = fit * harmonicMean;
        //weightedFit = (weightOverall * fit) + ((1 - weightOverall) * harmonicMean);
        //cout << "weightedFit = " << weightedFit << "." << endl;
        
        //fit = weightedFit;
        
        fit = fit + arithmeticMean;
    }
    return(fit);
}



// double calculatePartialDecisivenessSinglePartition (bool const& referenceTaxonPresent, int & numTrees,
//     vector < vector <int> > const& data, bool const& findAll, int const& partitionID,
//     vector <string> const& locusNames, bool const& verbose, int const& numProcs)








// GA function accepts data matrix, legal gene matrix, and number of total genes wanted in final matrix
// New stopping criterion: when PD does not increase for maxGenStalled generations
int optimizeDecisivnessGA (vector < vector <int> > & data,
    vector < vector <int> > const& legal, int const& numAdd, bool const& referenceTaxonPresent,
    int const& numProcs)
{
// initialize random number generator - already done in main
    // srand((unsigned)time(0));

// Constants. Some should probably be user specified: tournamentSize, populationSize, mutationProbability, maxgen
    // int numTaxa = data.size(); // not used
    // int numLoci = data[0].size(); // not used
    int numInitialLoci = countLociPresent(data);
    int tournamentSize = 5; // size of tournament group.
    const int populationSize = 500;
    double mutationProbability = 0.9; // probability of mutation, otherwise xover
    //double stepReplaceProbability = 0.01;
    int maxgen = 100000; // maximum number of generations; ramp up to 100000 or whatever
    int numTrees = 100; // ramp up to 1000
    
    double bestSolution = 0.0;
    vector < vector <int> > bestConfiguration;
    
    int genStalled = 0; // number of generations without an increase in best fitness.
    int maxGenStalled = 5000;
    
// Store fitness
    double fitness[populationSize];
    
    // cout << "initialize " << endl;
// Initialize population
    vector < vector <int> > population[populationSize];
    // int numAdd = maxGenes - numInitialLoci;     // number of genes to add
    int maxGenes = numInitialLoci + numAdd;
    
//     cout << "Initial count = " << numInitialLoci << "; adding " << numAdd
//         << " sampled characters for total of " << maxGenes << " occupied cells of a possible "
//         << data.size() * data[0].size() << "." << endl;
    
    
// Function currently errors out when too many genes initially present.
// This should be changed to allow genes to be removed.
    cout << endl << "SETTING UP INITIAL POPULATION" << endl << endl;
    cout << "Population consists of " << populationSize << " individuals." << endl;
    if (numAdd < 0) { // subtract genes; not sure when this will be used...
        cout << numAdd << " populated cells will be removed." << endl;
        for (int i = 0; i < populationSize; i++) {
            population[i] = data;
            for (int j = numAdd; j < 0; j++) {
                mutate(population[i], legal, 0, 1);
            }
        }
    } else { // add genes
        cout << numAdd << " additional random taxon-character cells will be populated for each individual..." << endl;
        for (int i = 0; i < populationSize; i++) {
            population[i] = data;
            for (int j = 0; j < numAdd; j++) { // Randomly add numAdd genes in legal positions
                mutate(population[i], legal, 1, 0); // last 2 arguments are add, remove (boolean)
            }
        }
    }
    cout << "Done." << endl << endl;
    
    //cout << "Best solution from initial population = " << bestSolution << "." << endl;
    
// Calculate initial fitness. For each matrix, crunch through 1000 trees
    cout << "Calculating initial fitnesses for all individuals..." << endl;
    
    for (int i = 0; i < populationSize; i++) {
        fitness[i] = calcFit(referenceTaxonPresent, numTrees, population[i], numProcs);
    //    cout << "Fitness for individual " << i << " is: " << fitness[i] << "." << endl;
        if (fitness[i] > bestSolution) {
            bestSolution = fitness[i];
            bestConfiguration = population[i];
    //        cout << "Best solution from initial population = " << bestSolution << "." << endl;
        }
    }
    cout << "Done." << endl;
    
    cout << "Best solution from initial population = " << bestSolution << "." << endl;
    
    bool done = 0;        // flag to end evolution loop
    int gen = 0;         // counter for number of generations

// Check for out of bounds parameter values
    if (populationSize < 3) {
        cout << "Population size must be greater than 2" << endl;
        return(1);
    }
    if (tournamentSize < 2) {
        cout << "Tournament group size must be greater than 1" << endl;
        return(1);
    }

// Evolution loop
    cout << endl << "STARTING EVOLUTION" << endl << endl;
    cout << "Running for a total of " << maxgen << " generations." << endl << endl;
    // cout << "Progress:" << endl;
    while (!done) {
        gen = gen + 1;
        cout << "Generation " << gen;
// Choose tournament group
// initialize best and worst with random
        int random = (int)(rand() % (populationSize));
        int best = random;        // integer to hold position of best member of tournament group
        int worst = random;        // position of worst member of tournament group
        double bestFit = fitness[random];    // best fitness
        double worstFit = fitness[random];    // worst fitness
        
// selection step
        //for (int i = 1; i < tournamentSize; i++)
        for (int i = 0; i < tournamentSize; i++) {
            random = (int)(rand() % populationSize);
            if (fitness[random] > bestFit) {
                bestFit = fitness[random];
                best = random;
            }
            if (fitness[random] < worstFit) {
                worstFit = fitness[random];
                worst = random;
            }
        }
        
// make sure best and worst are not the same
// if they are the same, choose new best randomly from population
        while (best == worst) {
            best = (int)(rand() % (populationSize));
        }

// Copy best individual to worst individual
        population[worst] = population[best];
        
        
// enter stepwise function. too demanding.
//         if (gen % 1000 == 0) {
//             cout << ": Entering stepwise phase." << endl;
//             stepReplace (population[worst], legal, maxGenes, referenceTaxonPresent, numProcs);
//         }
//         else {
// Mutate with some probability, otherwise xover
            double randmu = (double) rand() / RAND_MAX;
            if (randmu < mutationProbability) {
// mutate; change one cell.
                mutate(population[worst], legal, 1, 1);
            } else {
// xover with random individual in population
                int randx = (int)(rand() % populationSize);
                while (randx == worst) {
                    randx = (int)(rand() % populationSize);
                }
                xover(population[worst], population[randx], legal, maxGenes);
            }
//        }
// Calculate fitness of new individual
        fitness[worst] = calcFit(referenceTaxonPresent, numTrees, population[worst], numProcs);
        
        if (fitness[worst] > bestSolution) {
            bestSolution = fitness[worst];
            bestConfiguration = population[worst];
            genStalled = 0;
        } else {
            genStalled += 1;
        }
        
// Check to see if max number of generations has been reached
        cout << ": best solution thus far = " << bestSolution << "." << endl;
        if (gen == maxgen || genStalled == maxGenStalled) {
            done = 1;
        }
    } 
// End of evolution loop

//     double bestFit = 0;
//     int bestInd = 0;
//     for (int i = 0; i < populationSize; i++)
//     {
//         if (fitness[i] > bestFit) {
//             bestFit = fitness[i];
//             bestInd = i;
//         }
//     }
    
//    cout << endl << "Best result found: " << bestFit << endl;
// copy best evolved matrix to original data matrix
//    data = population[bestInd];
    
    //printData(data);
    
    
    cout << endl << "Best result found: " << bestSolution << endl;
    data = bestConfiguration;
    
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
    for (int i = 0; i < numTaxa; i++) {
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

    int count0 = 0;        // number of legal absent positions
    int count1 = 0;        // number of legal present positions

// count number of legal positions
// probably a faster way to do this
    for (int i = 0; i < numTaxa; i++) {
        for (int j = 0; j < numLoci; j++) {
            if (legal[i][j]) {
                if (data[i][j]) {
                    count1++;
                } else {
                    count0++;
                }
            }
        }
    }    
    
// select random currently absent gene to add
    int radd = 0;
    if (add) {
        if (!count0) {
            cout << "All legal loci full, cannot add." << endl;
            return;
        } else {
            radd = rand() % count0;
            if (debugging) {cout << "Randomly adding taxon-gene to cell " << radd << "." << endl;}
        }
    }
// select random currently present gene to remove
    int rremove = 0;
    if (remove) {
        if (!count1) {
            cout << "All legal loci empty, cannot remove." << endl;
            return;
        } else {
            rremove = rand() % count1;
            if (debugging) {cout << "Randomly removing taxon-gene from cell " << rremove << "." << endl;}
        }
    }
        
// modify data, adding random gene, removing random gene
    count0 = -1;
    count1 = -1;

    for (int i = 0; i < numTaxa; i++) {
        for (int j = 0; j < numLoci; j++) {
            if (legal[i][j]) {
                if (data[i][j]) {
                    count1++;
                }
                if (!data[i][j]) {
                    count0++;
                }

// after finding correct position, change value
                if (remove) {
                    if (count1 == rremove) {
                        data[i][j] = 0;
                        rremove = -1;
                    }
                }
                if (add) {
                    if (count0 == radd) {
                        data[i][j] = 1;
                        radd = -1;
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
    int upper = lower + rand() % (numLoci-lower);
    for (int i = 0; i < numTaxa; i++) {
        for (int j = lower; j <= upper; j++) {
            recipient[i][j] = donor[i][j];
        }
    }
// randomly add or remove genes if too many or too few
    int numAdd = maxGenes - countLociPresent(recipient);    
    if (numAdd > 0) {
        for (int i = 0; i < numAdd; i++) {
            mutate(recipient, legal, 1, 0);
        }
    }
    if (numAdd < 0) {
        for (int i=numAdd; i<0; i++) {
            mutate(recipient, legal, 0, 1);
        }
    }
}

// Functions to add genes one by one to maximize decisivness. 
// 12/7/2011

// Stepwise addition function. Sequentially adds genes at best positions
/*
int stepAdd (vector < vector < int > > & data, vector < vector < int > > const& legal, int maxGenes)
{
    int numTaxa = data.size();
    int numLoci = data[0].size();
    int numInitialLoci = countLociPresent(data);
    
// use -1 to solve uninitializing problems
    int bestTax = -1;    // store best taxa
    int bestLoci = -1;    // store best loci
    double bestDecisive = calcFit(data); // store best decisiveness

// Add until reaching the maximum number of genes
    for (int l = 0; l < (maxGenes - numInitialLoci); l++)
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
        if (bestTax != -1 && bestLoci != -1)
        {
            data[bestTax][bestLoci] = 1;
        }
    }
    return(0);
}
*/


// int optimizeDecisivnessGA (vector < vector <int> > & data,
//     vector < vector <int> > const& legal, int const& numAdd, bool const& referenceTaxonPresent,
//     int const& numProcs)
// fitness[i] = calcFit(referenceTaxonPresent, numTrees, data, numProcs);

// 12/10/2011
/* Hill climbing approach similar to stepAdd. In this function genes are added randomly until
maxGenes reached. The genes are then sequentially moved from current location to best location
one by one. This process repeats until no moves increases fitness. Can be stochastic result,
depends on starting condition (initialized randomly). This could result in an infinite loop.
TODO: set maximum number of loops */

// This is pretty computationally demanding. Lowering number of trees evaluated.
int stepReplace (vector < vector < int > > & data, vector < vector < int > > const& legal,
    int const& maxGenes, bool const& referenceTaxonPresent, int const& numProcs)
{
    int maxLoops = 100;     // maximum number of iterations through while loop. This would be better as a user defined parameter.
    int numTaxa = data.size();
    int numLoci = data[0].size();
    int numInitialLoci = countLociPresent(data);
    int numTrees = 100;
    
    // fill matrix randomly 
    int numAdd = maxGenes - numInitialLoci; // number of genes to add

    if (numAdd < 0) { // subtract genes
        cout << "Genes will be removed" << endl;
        for (int j = numAdd; j < 0; j++) {
            mutate(data, legal, 0, 1);
        }

    } else { // add genes
        for (int j = 0; j < numAdd; j++) {
            mutate(data, legal, 1, 0);
        }
    }

    bool done = 0;     // flag for end of loop
    int loops = 0;  // counter for number of times through while loop
    
// probably better if taxa and loci selected randomly.
    
    while (!done) {
        loops++;
        int moves = 0;
        cout << "Step #" << loops << endl;
        vector <int> temp (numTaxa);
        for (int i = 0; i < numTaxa; i++) {
            temp[i] = i;
        }
        
        vector <int> taxaOrder (numTaxa);
        for (int i = 0; i < numTaxa; i++) {
            taxaOrder[i] = (int)(rand() % (int)temp.size());
        }
        
        for (int k = 0; k < numTaxa; k++) {
            int taxa = taxaOrder[k];
            for (int loci = 0; loci < numLoci; loci++) {
                if (data[taxa][loci] & legal[taxa][loci]) { // if loci present and movable
                    data[taxa][loci] = 0;    // remove datum

// use -1 to solve uninitializing problems
                    int bestTax = -1;    // store best taxa
                    int bestLoci = -1;    // store best loci
                    //int bestTax;    // store best taxa
                    //int bestLoci;    // store best loci
                    //double bestDecisive = calcFit(data); // store best decisiveness
                    double bestDecisive = calcFit(referenceTaxonPresent, numTrees, data, numProcs);
                    
                    // test all possible positions
                    for (int i = 0; i < numTaxa; i++) {
                        for (int j = 0; j < numLoci; j++) {
                            if (legal[i][j] & !data[i][j]) {
                                vector < vector < int > > temp = data;
                                temp[i][j] = 1;
                                //double tempDecisive = calcFit(temp);
                                double tempDecisive = calcFit(referenceTaxonPresent, numTrees, temp, numProcs);
                                if (tempDecisive >= bestDecisive) {
                                    bestDecisive = tempDecisive;
                                    bestTax = i;
                                    bestLoci = j;
                                }
                            }
                        }
                    }
                    if (bestTax != -1 && bestLoci != -1) {
                        data[bestTax][bestLoci] = 1;
                    }
                    data[bestTax][bestLoci] = 1; // add loci to best position
                    if ((bestTax != taxa) || (bestLoci != loci)) {
                        moves++;
                    }
                }
            }
        }
// done when no moves made
        if (moves < 1) done = 1; 
// check for maximum number of iterations
        if (loops >= maxLoops) {
            cout << "Max loops reached" << endl;
            done = 1;
        }
    }
    return(0);
}

vector < vector <int> > getLegalMatrix (vector < vector <int> > const& data)
{
    vector < vector <int> > legal;
    int numTaxa = data.size();
    int numLoci = data[0].size();
    
    legal.resize(numTaxa, vector <int> (numLoci, 1));
    
//    cout << "First:" << endl;
//    printData(legal);
    for (int i = 0; i < numTaxa; i++) {
        for (int j = 0; j < numLoci; j++) {
            legal[i][j] =! data[i][j]; // clever
        }
    }
//    cout << "New:" << endl;
//    printData(legal);
    return(legal);
}

// troubleshooting function
void printData (vector < vector <int> > const& data)
{
    for (int i = 0; i < (int)data.size(); i++) {
        for (int j = 0; j < (int)data[0].size(); j++) {
            cout << data[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void printGADataToFile (vector < vector <int> > const& data, vector <string> const& taxonNames,
    string const& nexusFileName, int & numAddGA)
{
    ofstream log;
    
    string outputFileName = "GA_matrix_" + nexusFileName + "_add-" + convertIntToString(numAddGA) + ".log";
    
    log.open(outputFileName.c_str());
    
    //log.open("GA_matrix.log", ios::app);
    
    int numTaxa = (int)data.size();
    int numLoci = (int)data[0].size();
    
    string maxString;
    int longestName = 0;
    
    for (int i = 0; i < numTaxa; i++) {
        string currentString = taxonNames[i];
        if (currentString.size() > maxString.size()) {
            maxString = currentString;
        }
    }
// Determine length of longest name
    for (string::const_iterator iterCharacters = maxString.begin(); iterCharacters < maxString.end(); ++iterCharacters) {
        longestName++;
    }
    
    for (int i = 0; i < numTaxa; i++) {
        log << " " ;
// Print out leading spaces
        string tempName = taxonNames[i];
        if (tempName.size() < maxString.size()) {
            string::size_type tempDiffSize;
            tempDiffSize = maxString.size() - tempName.size();
            for (string::size_type iterSpaces = 0; iterSpaces < tempDiffSize; iterSpaces++) {
                log << " ";
            }
        }
        log << tempName << " ";
        for (int j = 0; j < numLoci; j++) {
            log << data[i][j];
        }
    log << endl;
    }
}
