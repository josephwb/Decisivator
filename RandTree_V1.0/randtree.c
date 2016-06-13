/* CURRENT AS OF 24 Mar 2011 (14.37.49)
 *
 * randtree : generate a random binary tree with unbiased shape, keep track of labels for 
 * internal nodes (union of labels of children)
 *
 * Usage:
 * 	randtree <num>
 * 		<num> the max number of elements the stack.
 *
 * Modification notes
 * 18 Apr 2011 (16.21.02)
 *
 * Future features
 *
*/
#define DEBUG 0

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "stack.h"
#include "ptr_tree.h"

#define DEBUGSEED 257 

TYPE_PTRTREE makeRandomTree ( int numberOfLeaves ) ; 
TYPE_PTRTREE makeRandomTree ( int numberOfLeaves ) {
	TYPE_STACK randomForest ; 
	char ST_BUFFER [ MAX_LABEL ] ; 
	int sret = 0 ; 

	TYPE_PTRTREE leftTree = NULL ; 
	TYPE_PTRTREE rightTree = NULL ; 
	TYPE_PTRTREE newTree = NULL ; 	

	int i = 0 ;	
	int j = 0 ; 
	int top = 0 ; 
	int bottom = 0 ;

	int MaxErrCode = 0 ; 
	int retCode = 0 ; 
	int nextNodeID = 0 ; 
	unsigned LChildIndex , RChildIndex , newChildIndex = 0 ; 
	unsigned int iseed = (unsigned int)time(NULL);
					if ( DEBUG ) { iseed = DEBUGSEED ; }
	srand (iseed);

  	// set maximum array size, create the array
	MaxErrCode = STACK_INIT ( &randomForest , numberOfLeaves ) ;
	if ( MaxErrCode ) {
		printf ( "ERROR Forest could not be initialized (%d)\n", MaxErrCode ) ;
		exit ( MaxErrCode ) ; 
		} 
					if ( DEBUG ) {
						printf ( "..the stack randForrest before initialization\n" ) ; 
						STACK_DUMP ( &randomForest ) ; 
						}

	// add all leaves to the stack
	for ( i = 0 ; i < numberOfLeaves ; i++ ) {
		newTree = PTRTREE_NODE_CREATE ( ) ; 
		sret = sprintf ( ST_BUFFER, "\'%d\'" , nextNodeID ) ;
		if ( !sret ) { 
			printf ( 
				"ERROR (%d) creating label LABEL=%s\n" ,
				sret, ST_BUFFER
				) ;
			} 	
		MaxErrCode = PTRTREE_SET ( newTree , NULL , NULL, NULL, nextNodeID++ , ST_BUFFER ) ; 
		if ( MaxErrCode ) { 
			printf ( 
				"ERROR (%d) initializing new node [THIS=%p PARENT=%p CHILD=%p SIB=%p ID=%d LABEL=%s]\n" ,
				MaxErrCode, newTree , NULL , NULL, NULL, nextNodeID , ST_BUFFER
				) ;
			}
		MaxErrCode = STACK_PUSH ( &randomForest , newTree ) ; 
		}
					if ( DEBUG ) {
						printf ( "..the stack randForrest after initialization\n" ) ; 
						STACK_DUMP ( &randomForest ) ; 
						}

	// start merging random sub-trees
	for ( i = 0 ; i < numberOfLeaves  - 1; i++ ) {
		// choose two "live" trees at random, move left one to the basement, replace right one
		// invarients: trees in the basement are unreliable
		// 		trees in the house are all roots of trees, all are disjoint
		//  	there are numberOfLeaves  leaves, i internal nodes, and LEAVE+i total nodes
		if ( STACK_IS_EMPTY ( &randomForest ) ) { 
			printf ( "ERROR tried to choose a stack entry when there are none (i=%d top=%d bottom=%d)\n" 
				, i , top , bottom
				) ; 
			STACK_DUMP ( &randomForest ) ; 
			}
		else { // there are trees in this forest!

			int randInt = 0 ; 

			// choose random child in the house
			top = STACK_TOP_IS ( &randomForest ) + 1 ; // remember: _TOP_ is index of top element
			randInt = (int) rand () % ( top - STACK_BOTTOM_IS ( &randomForest ) ) ;
			LChildIndex = 
				STACK_BOTTOM_IS ( &randomForest ) + randInt ;
			leftTree = 
				STACK_GET ( &randomForest , LChildIndex ) ;

			// send child to the basement, unchanged
			MaxErrCode = 
				STACK_SWAP ( &randomForest , LChildIndex , STACK_BOTTOM_IS ( &randomForest ) ) ;
			LChildIndex = STACK_BOTTOM_IS ( &randomForest ) ;
			MaxErrCode = // send left child (unchanged) to the basement
				STACK_LIFT ( &randomForest ) ;

			// choose random child in the (now smaller) house
			randInt = (int) rand () % ( top - STACK_BOTTOM_IS ( &randomForest ) ) ;
			RChildIndex = 
				STACK_BOTTOM_IS ( &randomForest ) + randInt ;
			rightTree = 
				STACK_GET ( &randomForest , RChildIndex ) ;

			// merge the two children, replacing right one
			sret = sprintf ( 
			ST_BUFFER, 
			"\'%d-%d\'" , 
			 PTRTREE_NODE_ID ( leftTree ) ,
			 PTRTREE_NODE_ID ( rightTree )
			) ; // make the label
			if ( !sret ) { 
			printf ( 
				"ERROR (%d) creating label from leftID=%u and rightID=%u\n" ,
				sret, 
				PTRTREE_NODE_ID ( leftTree ) ,
				PTRTREE_NODE_ID ( rightTree ) 
				) ;
			MaxErrCode = sret ; 
			}
			newTree = PTRTREE_NODE_MERGE ( 
						leftTree , rightTree , nextNodeID++ , ST_BUFFER
						) ;
			MaxErrCode = STACK_PUT ( &randomForest , newTree , RChildIndex ) ; 

			}	
		if ( MaxErrCode ) {
			printf ( 
				"ERROR (%d) adding merged tree %p with ID %d and label %s\n" ,
				MaxErrCode, newTree , nextNodeID , ST_BUFFER
				) ;
			}
		} // else
		return ( STACK_GET ( &randomForest , STACK_TOP_IS ( &randomForest ) ) ) ; 
	}


/*
 * main() - 
 *
 * The entry point to the program.
 *
 */
int main(int argc, char *argv[]) {
	int MaxErrCode = 0 ; 
	unsigned LEAVES = (argc <= 1 ? 10 : atoi( argv[1] ) ) ;

	PTRTREE_DUMP_NEWICK ( makeRandomTree ( LEAVES ) ) ;

	return ( MaxErrCode ) ; 
	}  // main()