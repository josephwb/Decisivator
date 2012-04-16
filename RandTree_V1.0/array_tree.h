#ifndef _ARRAYTREE_H
#define _ARRAYTREE_H

/* CURRENT AS OF 24 Mar 2011 (14.40.53) */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define _ARRAYTREE_CHUNK_SIZE 100				// number of nodes per increment of size


#define TYPE_STACK_DATA TYPE_ARRAYTREE_NODE		// change this to change node type, here and stack.h

/* 
 * array-based tree data structure and methods
 *
 * needed functionality:
 * 	add sib pointers?
 *  add dynamic memory management
 *
 * 16 Mar 2011 (15.47.48) JAF
 *
*/

/*	TYPE_ARRAYTREE_NODE
		|---|---|---|---|
[i]		^ left child of node i
			^ right child of node i
				^ parent of node i
					^ label (any string) of node i
*/

/* ----------------------------------------------------------
 *	data structure and manipulation routines
 */

// typedef for the 
typedef struct {
	int _ARRAYTREE_NODE_LCHILD ; 
	int _ARRAYTREE_NODE_RCHILD ; 
	int _ARRAYTREE_NODE_PARENT ; 
	char * 	_ARRAYTREE_NODE_LABEL ;
	} TYPE_ARRAYTREE_NODE ; 

typedef enum { LC_INDEX, RC_INDEX, PARENT_INDEX } TYPE_CHILD_INDEX ;

typedef struct {
	int _ARRAYTREE_NEXTNODE ;
	int _ARRAYTREE_MAXNODE ;
	TYPE_ARRAYTREE_NODE * _ARRAYTREE_DATA ;
	} TYPE_ARRAYTREE ; 

int ARRAYTREE_DUMP ( TYPE_ARRAYTREE * thistree ) ;

int ARRAYTREE_NODE_DUMP ( TYPE_ARRAYTREE_NODE * thisnode ) ;

char * ARRAYTREE_DUMP_NEXUS ( TYPE_ARRAYTREE * thisTree, int startHere ) ;

char * ARRAYTREE_DUMP_NODE_NEXUS ( TYPE_ARRAYTREE_NODE * thisTreeNode ) ;

int ARRAYTREE_IS_FULL ( TYPE_ARRAYTREE * thistree ) ;

int ARRAYTREE_IS_EMPTY ( TYPE_ARRAYTREE * thistree );

// create and initialize a new tree
int ARRAYTREE_INIT ( TYPE_ARRAYTREE * thistree , unsigned sizeofarray ) ;

// annihilate a tree
int ARRAYTREE_KILL ( TYPE_ARRAYTREE * thistree );

// add a node to the tree
int ARRAYTREE_NODE_ADD ( TYPE_ARRAYTREE * thistree, 
						int parenindex, 
						TYPE_CHILD_INDEX add_here, 
						char * label ) ;

// merge two nodes in the tree, creating a new root
int ARRAYTREE_NODE_MERGE ( TYPE_ARRAYTREE_NODE * thistree, 
						int nodeToRemove, 
						int ResultNode 
						) ;

// delete a node from the tree
int ARRAYTREE_NODE_DELETE ( TYPE_ARRAYTREE * thistree, int node ) ;

/* ----------------------------------------------------------
 *	pointer into a tree, and pointer manipulation routines
 */

// Check validity of index into the tree
int ARRAYNODE_PTR_IS_VALID ( TYPE_ARRAYTREE * into_tree, int ptr) ;

#endif