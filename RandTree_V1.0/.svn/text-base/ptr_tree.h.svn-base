#ifndef _PTRTREE_H
#define _PTRTREE_H

#ifndef DEBUG
#define DEBUG 0
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define PT 255				// number of nodes per increment of size
#define MAX_LABEL 255		// longest label allowable

/* 
 * individual nodes
 * 08 Apr 2011 (14.07.10)
 *
*/

/*	TYPE_PTRTREE_NODE
		|---|---|---|-------|
[i]	-->	^ * first child of node i
			^ * sibling of node i
				^ * parent of node i
					^ * label (any string) of node i
*/

struct STRUCT_PTRTREE_NODE ; 

typedef struct STRUCT_PTRTREE_NODE {
	struct STRUCT_PTRTREE_NODE * _PTRTREE_NODE_PARENT ; 
	struct STRUCT_PTRTREE_NODE * _PTRTREE_NODE_CHILD ; 
	struct STRUCT_PTRTREE_NODE * _PTRTREE_NODE_SIB ; 
	int _PTRTREE_NODE_ID ;     
	char _PTRTREE_NODE_LABEL[MAX_LABEL] ;
	}  *TYPE_PTRTREE ; 

//typedef struct STRUCT_PTRTREE_NODE *TYPE_PTRTREE  ; 


/*************** manipulation routines functions ***************/
// create and initialize a new tree root
int PTRTREE_SET ( 	TYPE_PTRTREE thisTree ,
					TYPE_PTRTREE parent ,
					TYPE_PTRTREE firstChild , 					
					TYPE_PTRTREE sib , 
					int id ,
					char * label
					) ;

// annihilate a tree
unsigned PTRTREE_KILL ( TYPE_PTRTREE thisTree ) ;

unsigned PTRTREE_IS_EMPTY ( TYPE_PTRTREE thisTree );

int PTRTREE_NODE_ID ( TYPE_PTRTREE thisTree ) ;

TYPE_PTRTREE PTRTREE_NODE_CREATE ( void ) ;

// add a node as the designated next_node of this tree, or create if thistree==NULL 
// return pointer to the new node
unsigned PTRTREE_NODE_ADD_THISTREE ( TYPE_PTRTREE thisTree, 
					const int ID,
					char * label
					) ;

// add childTree to parentTree, or create if thistree==NULL 
// return pointer to the new node
unsigned PTRTREE_NODE_ADD ( TYPE_PTRTREE parentTree, 
					TYPE_PTRTREE childTree ) ;

// merge two nodes in the tree, creating and returning a new root 
TYPE_PTRTREE PTRTREE_NODE_MERGE ( 
						TYPE_PTRTREE LeftTree , 
						TYPE_PTRTREE RightTree ,
						int ID , 
						char * NewLabel
						) ;

// delete a node from the tree
int PTRTREE_NODE_DELETE ( TYPE_PTRTREE thisTree ) ;

/*************** random tree maker: the main reason for this code ***********/
TYPE_PTRTREE makeRandomTree ( int numberOfLeaves ) ; 

/*************** output functions ***************/
unsigned PTRTREE_DUMP ( TYPE_PTRTREE thisTree ) ;

unsigned PTRTREE_DUMP_NEWICK ( TYPE_PTRTREE thisTree ) ;
unsigned _PTRTREE_DUMP_NEWICK_REC ( TYPE_PTRTREE thisTree ) ;

unsigned PTRTREE_DUMP_DOT ( TYPE_PTRTREE thisTree ) ;
unsigned _PTRTREE_DUMP_DOT_REC ( TYPE_PTRTREE thisTree ) ;


#endif

