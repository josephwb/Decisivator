// from scratch 
// $Id: ptr_tree.c 139 2011-06-14 23:37:52Z jamesafoster $ rev $Rev: 139 $ on $Date: 2011-06-14 16:37:52 -0700 (Tue, 14 Jun 2011) $ by $Author: jamesafoster $
#include "ptr_tree.h"

/* CURRENT AS OF 26 Mar 2011 (15.05.38)
 *
 * pointer-based tree data structure 
 *
 * needed functionality:
 *	error checking in INIT() and KILL()
 *	error checking in DELETE() and ADD()
 *	add ability to specify number of nodes on init
 *	add taxaname lookup to NEWICK dump
 * expand from left and right child, to any number of children (using first, sib, last)
 *
*/

/*	TYPE_PTRTREE_NODE
		|---|---|---|---|---|
[i]	-->	^ * left child of node i
			^ * right child of node i
				^ * sibling of node i
					^ * parent of node i
						^ * label (any string) of node i
*/

/* ----------------------------------------------------------
 *	data structure and manipulation routines
 */

int PTRTREE_SET ( 
	TYPE_PTRTREE thisTree ,
	TYPE_PTRTREE parent ,
	TYPE_PTRTREE CHILD , 
	TYPE_PTRTREE sib , 
	int id ,
	char * label
	) {
					if ( DEBUG ) { 
						printf ( 
							"--in PTRTREE_NODE_SET (%p, %p, %p, %p, %d, %s)\n" ,
							thisTree, parent, CHILD, sib, id, label
							) ; 
						} 
	thisTree -> _PTRTREE_NODE_PARENT = parent ; 
	thisTree -> _PTRTREE_NODE_CHILD = CHILD ; 
	thisTree -> _PTRTREE_NODE_SIB = sib ; 
	thisTree -> _PTRTREE_NODE_ID = id ; 
	strcpy ( thisTree -> _PTRTREE_NODE_LABEL , label ) ; 
					if ( DEBUG ) { printf ( "--leaving PTRTREE_NODE_SET\n" ) ; }
	return 0 ;
	} 

// annihilate a tree
unsigned PTRTREE_KILL ( TYPE_PTRTREE thisTree ) { 
	free ( thisTree ) ;
	return ( 0 ) ; 
	}

unsigned PTRTREE_IS_EMPTY ( TYPE_PTRTREE thisTree ) {
if ( DEBUG ) { 
	printf ( "  in PTRTREE_IS_EMPTY ( %p)\treturning %u\n ", thisTree, thisTree == NULL ) ; } ;  
	return ( thisTree == NULL ) ;
	}

int PTRTREE_NODE_ID ( TYPE_PTRTREE thisTree ) {
if ( DEBUG ) { 
	printf ( "  in PTRTREE_NODE_ID ( %p)\treturning %d\n ", thisTree, thisTree -> _PTRTREE_NODE_ID ) ; } ;  
	return ( thisTree -> _PTRTREE_NODE_ID ) ;
	}

// create an empty node and return a pointer to it
TYPE_PTRTREE PTRTREE_NODE_CREATE ( void ) {
					if ( DEBUG ) { printf ( "--in PTRTREE_NODE_CREATE\n" ) ; }
	TYPE_PTRTREE newChild = malloc ( sizeof ( struct STRUCT_PTRTREE_NODE ) ) ;
	PTRTREE_SET ( newChild, NULL, NULL, NULL, -1, "" ) ; 
					if ( DEBUG ) { printf ( "--leaving PTRTREE_NODE_CREATE\n" ) ; } 
	return ( newChild ) ;
	}

/* add childTree as next child of parentTree, return pointer to new node. 
 * precondition: thisTree points to a node of the correct type (not NULL)
 *
 * WARNING: should return a non-zero code and a warning message if child already has parent
*/
unsigned PTRTREE_NODE_ADD ( TYPE_PTRTREE parentTree ,
						TYPE_PTRTREE childTree 
	) {
	parentTree -> _PTRTREE_NODE_CHILD = childTree ; 
	childTree -> _PTRTREE_NODE_PARENT = parentTree ;
	// WARNING: if childTree already has a parent, that node's children will still point to child
	return ( 0 ) ;
	}

// add a node to the tree pointed at, return pointer to new node. 
// precondition: thisTree points to a node of the correct type (not NULL)
unsigned PTRTREE_NODE_ADD_THISTREE ( 
		TYPE_PTRTREE thisTree ,
		const int id , 
		char * label 
		) {
					if ( DEBUG ) { 
						printf ( "-- in PTRTREE_NODE_ADD ( %p, %d, %s)\n" ,
							thisTree, id, label 
							) ; 
						} 
	unsigned errCode = 0 ; 
	TYPE_PTRTREE newChild = PTRTREE_NODE_CREATE ( ) ;
	errCode = PTRTREE_SET ( 
		newChild, thisTree , NULL , thisTree -> _PTRTREE_NODE_CHILD , id , label 
		) ; 
	thisTree -> _PTRTREE_NODE_CHILD = newChild ; 
					if ( DEBUG ) { 
						printf ( "-- leaving PTRTREE_NODE_ADD returning %d\n" , errCode ) ; 
						} 
	return ( errCode ) ; 
	}
	
/* merge two nodes in a new tree, return pointer to that tree
 * NOTE: parent of new tree is NULL and ID is the passed parameter
 * assumes right and left are roots of trees with disjoint sets of nodes
 * actually doesn't have to be disjoint!!
 * returns pointer to new tree with new root, children are right and left
 */
TYPE_PTRTREE PTRTREE_NODE_MERGE ( 
	TYPE_PTRTREE LeftTree ,
	TYPE_PTRTREE RightTree ,
	int newID , 
	char * newLabel
	) {
					if ( DEBUG ) { 
						printf ( "-- in PTRTREE_NODE_MERGE ( %p, %p, %d, %s)\n" ,
						LeftTree, RightTree, newID, newLabel
							) ; 
						} 
	unsigned errCode = 0; 
	TYPE_PTRTREE newRoot = PTRTREE_NODE_CREATE ( ) ; 
	errCode = PTRTREE_SET ( 
		newRoot , 
		NULL , 
		LeftTree , 
		NULL ,
		newID , 
		newLabel 
		) ; 
	LeftTree -> _PTRTREE_NODE_PARENT = newRoot ;
	RightTree -> _PTRTREE_NODE_PARENT = newRoot ;
	LeftTree -> _PTRTREE_NODE_SIB = RightTree ;
					if ( DEBUG ) { 
						printf ( "-- leaving PTRTREE_NODE_MERGE, returns %p\n" ,
							newRoot
							) ; 
						} 
	return ( newRoot ) ;
	} 

// delete a node from the tree, recursively deleting each child in a DFS
int PTRTREE_NODE_DELETE ( TYPE_PTRTREE thisTree ) {
					if ( DEBUG ) { 
						printf ( "-- in PTRTREE_NODE_DELETE ( %p)\n" , thisTree ) ; 
						}
	int errCode = 0;
	if ( PTRTREE_IS_EMPTY ( thisTree ) ) {
		printf ( "ERROR attempted to delete a non-existent tree\n" ) ; 
		errCode = 99 ;
		} 
	else {
		TYPE_PTRTREE nextChild = thisTree -> _PTRTREE_NODE_CHILD ; 
		while ( nextChild ) { 
		// is this a memory leak? Does the label get deallocated??
			PTRTREE_NODE_DELETE ( nextChild -> _PTRTREE_NODE_SIB ) ;
			free ( thisTree ) ;
			}
		}
					if ( DEBUG ) { 
						printf ( "-- Leaving PTRTREE_NODE_DELETE, returns %u\n" , errCode ) ; 
						} 
	return errCode ;
	}

// dump the tree in DFS order
unsigned PTRTREE_DUMP ( TYPE_PTRTREE thisTree ) {
				if ( DEBUG ) { printf ( "<<in PTRTREE_DUMP ( %p )\n" , thisTree ) ; } 
 	TYPE_PTRTREE nextChild = NULL ; 
	TYPE_PTRTREE nextSib = NULL ; 
	if ( thisTree ) { 
		if ( !PTRTREE_IS_EMPTY ( thisTree ) ) { 
			// ASSUMES that undefined values are printable and of the right type
			printf ( "\t&THIS=%p\t&PARENT=%p\tCHILD=%p\t&SIB=%p\tID=%d\tLabel=%s\n", 
				thisTree , 
				thisTree -> _PTRTREE_NODE_PARENT ,			
				thisTree -> _PTRTREE_NODE_CHILD , 
				thisTree -> _PTRTREE_NODE_SIB ,
				thisTree -> _PTRTREE_NODE_ID , 
				thisTree -> _PTRTREE_NODE_LABEL 
				) ;
			nextChild = thisTree -> _PTRTREE_NODE_CHILD ;
			nextSib = thisTree -> _PTRTREE_NODE_SIB ;
			if ( nextChild ) { PTRTREE_DUMP ( nextChild ) ; } ;
			if ( nextSib ) { PTRTREE_DUMP ( nextSib ) ; } ;
			} 
		else {
			printf ( "ERROR tried to dump an empty tree (1)\n" ) ; return ( 1 ) ; 
			}
		} 
	else { 	
		printf ( "ERROR tried to dump an empty tree (2)\n" ) ; return ( 2 ) ; 
		}
					if ( DEBUG ) { 
						printf ( ">>leaving PTRTREE_DUMP\n" ) ; } 
						return ( 0 ) ; 
						}

unsigned PTRTREE_DUMP_NEWICK ( TYPE_PTRTREE thisTree ) {
	unsigned retVal = _PTRTREE_DUMP_NEWICK_REC ( thisTree ) ;
	printf ( "\n" ) ; 
	return ( retVal ) ; 
	}

unsigned _PTRTREE_DUMP_NEWICK_REC ( TYPE_PTRTREE thisTree ) {
					if ( DEBUG ) { 
						printf ( "\n<<in _PTRTREE_DUMP_NEWICK_REC ( %p )\n" , thisTree ) ; 
						}
	int errCode = 0 ; 
	TYPE_PTRTREE nextChild = NULL ; 	
	if ( PTRTREE_IS_EMPTY ( thisTree ) ) { // nothing to see here. Move along.
		errCode = 1 ;
		printf ( "ERROR dumping empty tree %p (returning %d)\n" , thisTree , errCode ) ; 
		} 
	else { // there IS a tree here 	
		nextChild = thisTree -> _PTRTREE_NODE_CHILD ;
		if ( !PTRTREE_IS_EMPTY ( nextChild ) ) { // this is an internal node
			printf ( "(" ) ; 
			_PTRTREE_DUMP_NEWICK_REC ( nextChild ) ; // dump children, comma-separated
			nextChild = nextChild -> _PTRTREE_NODE_SIB ; // now dump the sibs, if any
			while ( !PTRTREE_IS_EMPTY ( nextChild ) ) { 
				printf ( "," ) ;
				_PTRTREE_DUMP_NEWICK_REC ( nextChild ) ;
				nextChild = nextChild -> _PTRTREE_NODE_SIB ;
				}  
			printf ( ")" ) ; 
			printf ( "ID%d" , thisTree -> _PTRTREE_NODE_ID ) ;
			}  // end internal node
		else { // this is a leaf
			printf ( "ID%d" , thisTree -> _PTRTREE_NODE_ID ) ;  
			}
					if ( DEBUG ) { 
						// ASSUMES that undefined values are printable and of the right type
						printf ( "[LABEL=%s Parent=%p Child=%p Sib=%p THIS=%p]\n", 
							thisTree -> _PTRTREE_NODE_LABEL ,
							( PTRTREE_IS_EMPTY ( thisTree ->  _PTRTREE_NODE_PARENT ) ?
								NULL :
								thisTree -> _PTRTREE_NODE_PARENT
								) ,
							( PTRTREE_IS_EMPTY ( thisTree -> _PTRTREE_NODE_CHILD ) ?
								NULL :
								thisTree -> _PTRTREE_NODE_CHILD
								) ,
							( PTRTREE_IS_EMPTY ( thisTree -> _PTRTREE_NODE_SIB ) ?
								NULL :
								thisTree -> _PTRTREE_NODE_SIB
								) ,
							thisTree
							) ;
						}
		}  // end of details when there is a tree
					if ( DEBUG ) { printf ( ">>leaving _PTRTREE_DUMP_NEWICK_REC\n" ) ; } 
	return ( errCode ) ; 
	} // end of dump

	
unsigned PTRTREE_DUMP_DOT ( TYPE_PTRTREE thisTree ) {
	printf ( "graph \"RandTree\" {\n" ) ; 
	int retVal = _PTRTREE_DUMP_DOT_REC ( thisTree ) ; 
	printf ( "}\n" ) ; 
	return ( retVal ) ; 
	}
	
unsigned _PTRTREE_DUMP_DOT_REC ( TYPE_PTRTREE thisTree ) {
					if ( DEBUG ) { 
						printf ( "\n<<in _PTRTREE_DUMP_DOT_REC ( %p )\n" , thisTree ) ; 
						}
	TYPE_PTRTREE nextTree = NULL ;
	TYPE_PTRTREE nextChild = NULL ;
	int errCode = 0 ; 	
	if ( PTRTREE_IS_EMPTY ( thisTree ) ) { // nothing to see here. Move along.
		errCode = 1 ;
		printf ( "ERROR dumping empty tree %p (returning %d)\n" , thisTree , errCode ) ; 
		} 
	else { // there is a tree here
		printf ( "\tID%d [label=\"LABEL=%s" , 
				thisTree -> _PTRTREE_NODE_ID ,
				thisTree -> _PTRTREE_NODE_LABEL 
				) ; 
					if ( DEBUG ) {
						printf ( " Parent=%p Child=%p Sib=%p THIS=%p" , 
							( PTRTREE_IS_EMPTY ( thisTree ->  _PTRTREE_NODE_PARENT ) ?
								NULL :
								thisTree -> _PTRTREE_NODE_PARENT
								) ,
							( PTRTREE_IS_EMPTY ( thisTree -> _PTRTREE_NODE_CHILD ) ?
								NULL :
								thisTree -> _PTRTREE_NODE_CHILD
								) ,
							( PTRTREE_IS_EMPTY ( thisTree -> _PTRTREE_NODE_SIB ) ?
								NULL :
								thisTree -> _PTRTREE_NODE_SIB
								) ,
							thisTree
							) ;
						}
		printf ( "\"]\n" ) ; 

		nextTree = nextChild = thisTree -> _PTRTREE_NODE_CHILD ;
		if ( !PTRTREE_IS_EMPTY ( nextTree ) ) { // not a leaf
			printf ( "\tID%d -- { " , thisTree -> _PTRTREE_NODE_ID ) ; // edges to each child
			while ( nextTree ) {
				printf ( "ID%d " , nextTree -> _PTRTREE_NODE_ID ) ;
				nextTree = nextTree -> _PTRTREE_NODE_SIB ; 
				}
			printf ( "}\n" ) ; 

			while ( nextChild ) { // BF recursive traversal
				_PTRTREE_DUMP_DOT_REC ( nextChild ) ;
				nextChild = nextChild -> _PTRTREE_NODE_SIB ; 
				}
			}
					if ( DEBUG ) { printf ( ">>leaving _PTRTREE_DUMP_DOT_REC\n" ) ; } 
		}
	return ( errCode ) ; 
	}
