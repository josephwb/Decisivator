#include "array_tree.h"

/* CURRENT AS OF 24 Mar 2011 (14.38.09)
 *
 * array-based tree data structure 
 *
 * needed functionality:
 *	error checking in INIT() and KILL()
 *	error checking in DELETE() and ADD()
 *	change ADD to not require left and right children 
 *	add ability to specify number of nodes on init
 *	add taxaname lookup to nexus dump
 *
 * 16 Mar 2011 (15.47.48) JAF
 *
*/

/*
		|---|---|---|---|
[i]		^ left child of node i
			^ right child of node i
				^ parent of node i
					^ label (any string) of node i
*/

/* ----------------------------------------------------------
 *	data structure and manipulation routines
 */

// create and initialize a new tree
int ARRAYTREE_INIT ( TYPE_ARRAYTREE * thistree , unsigned sizeOfArray ) {
	if ( sizeOfArray ) {
		_ARRAYTREE_MAXNODE = sizeOfArray ; 
		}
	else {
		_ARRAYTREE_MAXNODE = _ARRAYTREE_CHUNK_SIZE ; 
	]
	thistree -> _ARRAYTREE_NEXTNODE = 0 ;	

	// *** add ERROR CHECKING ****
	thistree -> _ARRAYTREE_DATA = ( TYPE_ARRAYTREE_NODE * ) malloc ( _ARRAYTREE_MAXNODE * sizeof (  TYPE_ARRAYTREE_NODE ) );
	} ;

// annihilate a tree
int ARRAYTREE_KILL ( TYPE_ARRAYTREE * thistree ) { 
	free ( thistree -> _ARRAYTREE_DATA ) ;
	} ;

int ARRAYTREE_IS_FULL ( TYPE_ARRAYTREE * thistree ) {
	return ( thistree -> _ARRAYTREE_NEXTNODE == thistree -> _ARRAYTREE_MAXNODE ? 0 : 1 ) ;
	} ;

int ARRAYNODE_PTR_IS_VALID ( TYPE_ARRAYTREE * into_tree, int ptr) {
	return	( ( ptr < into_tree -> _ARRAYTREE_MAXNODE ) ? 1 : 0 ) ; 
	} ;

int ARRAYTREE_IS_EMPTY ( TYPE_ARRAYTREE * thistree ) {
	return ( !thistree -> _ARRAYTREE_NEXTNODE ) ;
	} ;

// add a node to the tree at the leaf addnode
// thistree -> [ NEXT | MAX | [ node0 node1 ... parentNode ... ]
//                                              ^ parentIndex
int ARRAYTREE_NODE_ADD ( TYPE_ARRAYTREE * thistree , 
		int parentIndex , // ERROR CHECKING needed
		TYPE_CHILD_INDEX add_here , 
		char * label ) {
	
	if ( ARRAYTREE_IS_FULL ( thistree ) ) {
		return ( 1 ) ;
		} 
	else {
		TYPE_ARRAYTREE_NODE parentNode = thistree -> _ARRAYTREE_DATA [ parentIndex ] ;
		int newNodeIndex = thistree -> _ARRAYTREE_NEXTNODE++ ; // node with leaf to be changed ;	
		TYPE_ARRAYTREE_NODE thisnode = thistree -> _ARRAYTREE_DATA [ newNodeIndex ] ;
		thisnode._ARRAYTREE_NODE_LCHILD = 0 ;			
		thisnode._ARRAYTREE_NODE_RCHILD = 0 ;
		thisnode._ARRAYTREE_NODE_PARENT = parentIndex ;
		if ( add_here == LC_INDEX ) {
			parentNode._ARRAYTREE_NODE_LCHILD = parentIndex ;
			} 
		else {
			parentNode._ARRAYTREE_NODE_RCHILD = parentIndex ;
			} ; 
		if ( strcpy ( thisnode._ARRAYTREE_NODE_LABEL, label ) ) {
			printf ( "ERROR failed adding node lable %s\n" , label ) ;
			} ;
		} ;
		
	return 0;
	} ;
	
	// merge two nodes in the tree, changing the SECOND
int ARRAYTREE_NODE_MERGE ( TYPE_ARRAYTREE_NODE * thistree, 
						int nodeToRemove, 
						int nodeToPromote
						) {
	// put the nodes in place
	// 1. move nodeToRemove to left, under the floor
	// 2. make new parent for these two nodes
	// 3. push onto array
	_TYPE_STACK ArrayOfNodes = thistree ->  _ARRAYTREE_DATA ;
	int ceiling = ArrayOfNodes -> _ARRAYTREE_DATA._STACK_TOP ;
	int floor = ArrayOfNodes -> _ARRAYTREE_DATA._STACK_FLOOR ;
	TYPE_ARRAYTREE_NODE RemoveMe = ArrayOfNodes [ nodeToRemove ] ;
	TYPE_ARRAYTREE_NODE PromoteMe = ArrayOfNodes [ nodeToPromote ] ;
	TYPE_ARRAYTREE_NODE NewTreeRoot = ArrayOfNodes [ ceiling ] ;
	
	RemoveMe._ARRAYTREE_NODE_PARENT = ceiling ;
	PromoteMe._ARRAYTREE_NODE_PARENT = ceiling ;
	NewTreeRoot._ARRAYTREE_NODE_LCHILD = nodeToRemove ;
	NewTreeRoot._ARRAYTREE_NODE_RCHILD = nodeToPromote ;
	NewTreeRoot._ARRAYTREE_NODE_PARENT = ceiling ;
	NewTreeRoot._ARRAYTREE_NODE_LABLE = scanf ( "%s %s",
		RemoveMe._ARRAYTREE_NODE_LABLE ,
		PromoteMe._ARRAYTREE_NODE_LABLE
		) ;
	
	STACK_swap_elements ( ArrayOfNodes, int nodeToRemove, int floor );
	int STACK_lift ( ArrayOfNodes ) ;
	STACK_push ( ArrayOfNodes, NewTreeRoot );
	
	} ;

// delete a node from the tree
int ARRAYTREE_NODE_DELETE ( TYPE_ARRAYTREE * thistree, int thisnode ) {
	int next_node = 0;
	if ( ARRAYTREE_IS_EMPTY ( thistree ) ) {
		return ( 1 ); 
		} ;
	for ( next_node = thisnode; next_node < thistree -> _ARRAYTREE_MAXNODE; ) {
		thistree -> _ARRAYTREE_DATA [ next_node ] = thistree -> _ARRAYTREE_DATA [ next_node++ ] ;
		} ;
	return 0 ;
	} ;

int ARRAYTREE_DUMP ( TYPE_ARRAYTREE * thistree ) {
	printf ( "--Dumping tree\tNext %d\tMax %d\n", 
		thistree -> _ARRAYTREE_NEXTNODE,
		thistree -> _ARRAYTREE_MAXNODE 
		) ;
	printf ( "  Data\n" ) ; 
	int i = 0;
	TYPE_ARRAYTREE_NODE thisnode ; 
	for ( i = 0; i < thistree -> _ARRAYTREE_MAXNODE; i++ ) {
		thisnode = thistree -> _ARRAYTREE_DATA [ i ] ;
		ARRAYTREE_NODE_DUMP ( ( TYPE_ARRAYTREE_NODE *) &thistree ) ;
		} ;
	printf ( "\n" ) ; 
	} ;

char * ARRAYTREE_DUMP_NEXUS ( TYPE_ARRAYTREE * thisTree, int startHere ) {
	if ( !ARRAYTREE_IS_EMPTY ( thisTree ) {
		printf ( "%s\n", 
			ARRAYTREE_DUMP_NODE_NEXUS ( thisTree, 
			thisTree -> _ARRAYTREE_DATA [ starthere ] ) ;
		}
	else {
		printf ( "tree is empty/n" ) ;
		}
	} ;

char * ARRAYTREE_DUMP_NODE_NEXUS ( TYPE_ARRAYTREE_NODE * thisTreeNode ) {
	printf ( "(%d,%d:%s)",
		thisTreeNode -> _ARRAYTREE_NODE_LCHILD,
		thisTreeNode -> _ARRAYTREE_NODE_RCHILD,
		thisTreeNode -> _ARRAYTREE_NODE_LABLE,
		) ;
	if ARRAYNODE_PTR_IS_VALID ( thisTreeNode, thisTreeNode -> _ARRAYTREE_NODE_LCHILD ) {
		ARRAYTREE_DUMP_NODE_NEXUS ( thisTreeNode -> _ARRAYTREE_NODE_LCHILD ) ;
		} ;
	if ARRAYNODE_PTR_IS_VALID ( thisTreeNode, thisTreeNode -> _ARRAYTREE_NODE_LCHILD ) {
		ARRAYTREE_DUMP_NODE_NEXUS ( thisTreeNode -> _ARRAYTREE_NODE_RCHILD ) ;
		} ;
	printf ( "\n" ) ; 
	} ;

int ARRAYTREE_NODE_DUMP ( TYPE_ARRAYTREE_NODE * thisnode ) {
	printf (" (Left=%d Right=%d Parent=%d Label=%s)\n",
		thisnode -> _ARRAYTREE_NODE_LCHILD ,
		thisnode -> _ARRAYTREE_NODE_RCHILD ,
		thisnode -> _ARRAYTREE_NODE_PARENT ,
		thisnode -> _ARRAYTREE_NODE_LABEL ) ; 
	return 0;
	} ;