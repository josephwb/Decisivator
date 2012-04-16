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
 * jaf 07 Mar 2011 (11.33.42) created
 *
 * Future features
 *		add tests for error conditions
 *
*/
#define DEBUG 0

#include <stdio.h>
#include <stdlib.h>

#include "stack.h"
#include "ptr_tree.h"

int i;

TYPE_STACK randomForest ; 
char ST_BUFFER[255] ; 
int sret = 0 ;
unsigned LEAVES = 0 ;
unsigned MAX = 0 ;
int MaxErrCode = 0 ;

int i = 0 ;	
int j = 0 ; 
enum TreeType {LONGTREE, WIDETREE, FULLTREE} ; 
const char * TTLabel[] = { "LONGTREE" , "WIDETREE" , "FULLTREE" } ;
unsigned TTNUM = 3 ; 

TYPE_PTRTREE rootTree = NULL ; 
TYPE_PTRTREE thisTree = NULL ; 

void TreeHarnessTest ( enum TreeType treetype ) {
	int push_code = 0 ; 
	int nodeID = 0 ;

	for ( i = 0; i< LEAVES; i++ ) {
		thisTree = rootTree = PTRTREE_NODE_CREATE ( ) ;
		sret = sprintf ( ST_BUFFER, "\'%d (root)\'" , nodeID ) ;
		PTRTREE_SET ( rootTree, NULL , NULL, NULL, nodeID++, ST_BUFFER ) ; 
						if ( DEBUG ) {
					printf ( "\n**New tree for i=%d is: %p \n" , i , rootTree) ;
					PTRTREE_DUMP ( rootTree ) ; 
					} ;
		for ( j = 0; j <= i; j++ ) {	
			// build trees of depth j to put in ith slot
			sret = sprintf ( ST_BUFFER, "\'%d,%d\'" , 
				PTRTREE_NODE_ID ( rootTree ) , 
				nodeID
				) ;
			switch ( treetype ) {
				case WIDETREE : { 
						if ( DEBUG ) { printf ( "in WideTree\n" ) ; }
					PTRTREE_NODE_ADD_THISTREE ( rootTree , nodeID++ , ST_BUFFER ) ;
						if ( DEBUG ) { PTRTREE_DUMP ( rootTree ) ; } ; 
					break ; 
					}
				case LONGTREE : { 
						if ( DEBUG ) { printf ( "in LONGTREE\n" ) ; }
					PTRTREE_NODE_ADD_THISTREE ( thisTree , nodeID++ , ST_BUFFER ) ;
					PTRTREE_NODE_ADD_THISTREE ( thisTree , nodeID++ , ST_BUFFER ) ;
					thisTree = thisTree -> _PTRTREE_NODE_CHILD ; 
						if ( DEBUG ) { PTRTREE_DUMP ( rootTree ) ; } ; 
					break ;
					}
				case FULLTREE : { 
						if ( DEBUG ) { printf ( "in FULLTREE\n" ) ; }
					int k = 0;
					for ( k=0 ; k<=j ; k++ ) {	
						PTRTREE_NODE_ADD_THISTREE ( thisTree , nodeID++ , ST_BUFFER ) ;	
						}
					thisTree = thisTree -> _PTRTREE_NODE_CHILD ; 
					if ( DEBUG ) { PTRTREE_DUMP ( rootTree ) ; } ; 
					break ; 
					}
				}  // end of switch
			} ; // inner, j loop
					if ( DEBUG ) { STACK_DUMP ( &randomForest ) ; }
			push_code = STACK_PUSH ( &randomForest , rootTree ) ; 		
			if ( push_code ) {
				printf ( "\nERROR STACK_PUSH failed (%d)\n" , MaxErrCode ) ;
				exit ( MaxErrCode ) ;  
				} ; 
		} ; // outer, i loop
	} ; // end of tree harness

/*
 * main() - 
 *
 * The entry point to the program.
 *
 */
int main(int argc, char *argv[]) {
  	// set maximum array size
	LEAVES = (argc <= 1 ? 10 : atoi( argv[1] ) ) ;
  	MAX = (argc <= 1 ? 10 : 2 * LEAVES );

	printf ("\n**** beginning pointer tests ** \n" ) ; 
	enum TreeType treetype = 0 ;
	for ( treetype = 0 ; treetype < TTNUM ; treetype++ ) { 
		printf ( "** testing treetype %s\n" , TTLabel[treetype] ) ;
		MaxErrCode = STACK_INIT ( &randomForest , MAX ) ;
		if ( MaxErrCode ) {
		MaxErrCode = 99 ; 
		printf ( "ERROR Forest could not be initialized (%d)\n", MaxErrCode ) ;
		exit ( MaxErrCode ) ; 
		} ;
					if ( DEBUG ) { STACK_DUMP ( &randomForest ) ; }
		TreeHarnessTest ( treetype ) ; 
		STACK_DUMP ( &randomForest ) ; 
		}

	printf ("**** beginning STACK tests ** \n" ) ; 
	// *randomForrest is now a stack of full trees
	int retCode = 0 ; 
	TYPE_PTRTREE LeftTree, RightTree = NULL; 	

	printf ( "** beginning POP test ** \n" ) ; 
					if ( DEBUG ) { 
						STACK_DUMP ( &randomForest ) ;
						}
	LeftTree = STACK_POP ( &randomForest ) ; 
	STACK_DUMP ( &randomForest ) ; 
					if ( DEBUG ) {
						STACK_DUMP ( &randomForest ) ;
						printf ( "..after STACK_POP, STACK_DUMP returns %p\n" , LeftTree ) ; 
						}

	printf ("** beginning PUSH test ** \n" ) ; 
					if ( DEBUG ) { 
						STACK_DUMP ( &randomForest ) ;
						}
	retCode = STACK_PUSH ( &randomForest, LeftTree ) ; 
	STACK_DUMP ( &randomForest ) ; 
					if ( DEBUG ) {
						STACK_DUMP ( &randomForest ) ;
						printf ( "..after STACK_PUSH returns %d\n" , retCode ) ; 
						}

	printf ("** beginning LIFT test ** \n" ) ;
					if ( DEBUG ) { 
						STACK_DUMP ( &randomForest ) ;
						}
	retCode = STACK_LIFT ( &randomForest ) ; 
	STACK_DUMP ( &randomForest ) ; 
					if ( DEBUG ) {
						STACK_DUMP ( &randomForest ) ;
						printf ( "..after STACK_LIFT returns %d\n" , retCode ) ; 
						}

	printf ("** beginning SWAP test ** \n" ) ; 
					if ( DEBUG ) { 
						STACK_DUMP ( &randomForest ) ;
						}
	retCode = STACK_SWAP ( &randomForest , 
		STACK_TOP_IS ( &randomForest ) , 
		STACK_BOTTOM_IS ( &randomForest ) 
		) ;
					if ( DEBUG ) {
						STACK_DUMP ( &randomForest ) ;
						printf ( 
							"..after STACK_SWAP ( %u, %u ) returning %d\n" , 
							STACK_TOP_IS ( &randomForest ) , 
							STACK_BOTTOM_IS ( &randomForest ) ,
							retCode 
							) ; 
						}

	}  // main()
	