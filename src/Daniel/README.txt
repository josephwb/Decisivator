12/15/2011
The files optimize.cpp and optimize.h are the newest version of ga.cpp and ga.h. The name was changed to reflect the presence of multiple functions for the optimization of the character-by-taxon matrix.

Optimizing functions:
	The function optimizeDecisivnessGA() uses a genetic algorithm to optimally add positions in the matrix. 
	The function stepAdd() greedily adds positions to the matrix sequentially. The best position is added each time.
	The function stepReplace() first randomly adds positions to the matrix until the maximum number of positions is reached. It then sequentially moves each position to the optimal position in the matrix. It repeats this process until no moves increase the fitness of the matrix.


Acessory functions:
	The function calcFit() tests the fitness of each matrix. This should be a measure of decisiveness. This function accepts a vector of vectors and returns a number between 0 and 1 with 1 being the optimum.
	The function mutate() adds and/or subtracts a random position in the matrix.
	The function xover() copies a random part of one matrix into another. This is used by the genetic algorithm to produce new individuals.
	The function countLociPresent() counts the number of positions present in a matrix.

Note: the main function and the calcFit function are currentlu included only for testing purposes. The main function should be deleted and the calcFit function modified to return a measure of decisiveness.


