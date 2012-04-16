#ifndef _STACK_H
#define _STACK_H

#ifndef DEBUG
#define DEBUG 0
#endif

#include <stdio.h>
#include <stdlib.h>

/* copy and set these each compile
#include "--enterfilehere--"    	// h file that defines functions below
#define DDF --enter function here ; // extern DUMP_FUNCTION DDF ;
#define DDT --ebntertypehere-- ; 	// extern DUMP_FUNCTION DDF ; must be %p
*/

#include "ptr_tree.h"
#define DDT TYPE_PTRTREE 	// data type in array
#define DDF PTRTREE_DUMP_NEWICK // function to dump DDT type elements

typedef struct STRUCT_STACK { 
	unsigned _MAX_SIZE ;
	unsigned _STACK_TOP ;
	unsigned _STACK_FLOOR ;
	unsigned _STACK_FOUNDATION ;
	unsigned _STACK_BYTES ; 
	DDT * _STACK_DATA ;
	} TYPE_STACK ;

int STACK_INIT ( TYPE_STACK * stack, unsigned size ); // size in #/ptrs

int STACK_PUSH ( TYPE_STACK * stack, DDT element );

DDT STACK_POP ( TYPE_STACK * stack );

DDT STACK_GET ( TYPE_STACK * stack , unsigned here );

unsigned STACK_PUT ( TYPE_STACK * stack , DDT newElement , unsigned here );

int STACK_IS_FULL ( TYPE_STACK * stack );

int STACK_IS_EMPTY ( TYPE_STACK * stack );

int STACK_LIFT ( TYPE_STACK * stack );

int STACK_SWAP ( TYPE_STACK * stack, unsigned from, unsigned to );

int STACK_DUMP ( TYPE_STACK * stack ); 

unsigned STACK_TOP_IS ( TYPE_STACK * thisStack ) ;

unsigned STACK_BOTTOM_IS ( TYPE_STACK * thisStack ) ;

unsigned STACK_FOUNDATION_IS ( TYPE_STACK * thisStack ) ; 

#endif

