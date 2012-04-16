#include "stack.h"

int i;

/* 
 * stack-based data structure with a movable floor. NOT GENERAL. this is a STACK
 * of trees as implemented in ptr_tree
 *
 * needed functionality:
 * - change to stack of void pointers, so it can be used on anything
 *   use XOR in swap
 * - add error checking to STACK_INIT, STACK_PUSH, STACK_POP
 * - add dynamic resizing
 * - add dynamic allocation of data array on init()
 * - return size of stack, instead of void, on most functions
 * - add error checking to lift function
 * - add option to dump to declare format of the dump 
 *
 * JAF 07 Mar 2011 (12.07.30) : created
 *
*/

/*
	|-------|------------|-------|
	^ _STACK_FOUNDATION
	        ^ _STACK_FLOOR
	                     ^ _STACK_TOP
                                 ^ _MAX_SIZE 
    <------ _STACK_SIZE ---------> (in bytes)

	initialize stack pointer
*/

int STACK_INIT ( TYPE_STACK * thisStack, unsigned size ) { // size in #/ptrs
						if ( DEBUG ) { printf ( "in STACK_INIT ( %p, %u )\n" , thisStack, size ) ; }
 	thisStack -> _MAX_SIZE = size ; 
	thisStack -> _STACK_TOP = 0 ;
	thisStack -> _STACK_FLOOR = 0 ;
	thisStack -> _STACK_FOUNDATION = 0 ;
	thisStack -> _STACK_BYTES = size * sizeof( void * ) ;
	thisStack -> _STACK_DATA = malloc ( thisStack -> _STACK_BYTES ) ;
	return ( thisStack -> _STACK_DATA ) ? 0 : 666 ;
	}

/*  push an elPTR into stack
	precondition: the stack is not full
	add ERROR CHECKING
*/
int STACK_PUSH ( TYPE_STACK * thisStack, DDT elPTR ) {
					if ( DEBUG ) { 
						printf ( ">>in STACK_PUSH ( thisStack=%p, elPTR=%p )\n", thisStack, elPTR ) ; 
						} 
	int errCode = 0 ; 
	int nextElement = thisStack -> _STACK_TOP ;
	if ( !STACK_IS_FULL ( thisStack ) ) {
		thisStack -> _STACK_TOP = nextElement + 1 ; 
		thisStack -> _STACK_DATA [ nextElement ] = elPTR;
		} 
	else {
		printf ( "ERROR tried to push into a full stack at %p returning %d\n" , 
				thisStack , errCode 
				) ; 
		errCode = 1 ;
		} ;
					if ( DEBUG ) { 
									printf 
										( "<<leaving STACK_PUSH ( %p, %p )\tnextElement=%d TOP=%d\n", 
										thisStack, elPTR, nextElement, thisStack -> _STACK_TOP 
										) ;
						} 
	return ( errCode ) ;
	}

/*  pop an elPTR from stack
	precondition: stack is not empty
	add ERROR CHECKING
*/
DDT STACK_POP ( TYPE_STACK * thisStack ) {
					if ( DEBUG ) { printf ( ">>in STACK_POP ( thisStack=%p)\n", thisStack ) ; }
	DDT popElement = NULL ; 
	if ( !STACK_IS_EMPTY ( thisStack) ) {
		popElement = ( DDT ) thisStack -> _STACK_DATA [ --( thisStack -> _STACK_TOP ) ] ;
		} 
	else {
		printf ( "ERROR tried to pop and empty stack at %p\n" , thisStack ) ; 
		} ;
					if ( DEBUG ) { 
						printf ( "<<leaving STACK_POP, returning %p\n" , popElement ) ;
						} 	
	return ( popElement ) ; 
	}

/* get the element at index thisElement and return it (should be a pointer for randtree)
*/
DDT STACK_GET ( TYPE_STACK * stack , unsigned thisone ) {
							if ( DEBUG ) { 
						printf ( ">>entering STACK_GET (%p , %d)\n" , stack , thisone ) ;
						} 	
	DDT returnElement = NULL ;
	if ( thisone < stack -> _STACK_FLOOR ) { 
		printf ( "WARNING tried to get item %d from the basement (%d)\n" ,
			thisone , stack -> _STACK_FLOOR 
			) ;
		} 
	if ( (thisone < stack -> _STACK_TOP ) && ( thisone >= stack -> _STACK_FOUNDATION ) ) {
		returnElement = stack -> _STACK_DATA [ thisone ] ;
		}
	else {
		printf ( 
			"ERROR tried to get element %d from %p, when top is %d and foundation is %d\n" , 
			thisone , stack , stack -> _STACK_TOP , stack -> _STACK_FOUNDATION
			) ; 
			returnElement = NULL ;
		}
							if ( DEBUG ) { 
						printf ( "<<leaving STACK_GET, returning %p\n" , returnElement ) ;
						} 	
	return ( returnElement ) ; 
}

/* put the newElement at index "here" of stack replacing what's there
*/
unsigned STACK_PUT ( TYPE_STACK * stack , DDT newElement , unsigned here ) {
							if ( DEBUG ) { 
								printf ( 
									">>entering STACK_PUT(%p , %p , %d)\n" , 
									stack , newElement , here 
								) ;
								}
	int retCode = 0 ; 
	if ( here < stack -> _STACK_FLOOR ) { 
		printf ( "WARNING tried to get item %d from the basement (%d)\n" ,
			here , stack -> _STACK_FLOOR 
			) ;
		} 
	if ( 
				( here < ( stack -> _STACK_TOP ) )
				&&  ( here >= ( stack -> _STACK_FOUNDATION ) ) 
			) {
		stack -> _STACK_DATA [ here ] = newElement ;
		}
	else {
		printf ( 
			"ERROR tried to put element %p into location %d of %p, when top is %d and bottom is %d\n" , 
			newElement , here, stack , stack -> _STACK_TOP , stack -> _STACK_FOUNDATION
			) ; 
		retCode = 1 ; 
		}
							if ( DEBUG ) { 
								printf ( "<<leaving STACK_PUT, returning %d\n" , retCode ) ;
								}
	return ( retCode ) ; 
	}


/*  return 1 if stack is full, otherwise return 0
	add ERROR CHECKING
*/
int STACK_IS_FULL ( TYPE_STACK * thisStack )
{
	int returnValue = ( 
		( 
			(thisStack -> _STACK_TOP) == (thisStack -> _MAX_SIZE) 
		) ? 1 : 0 
	) ;
					if ( DEBUG ) { 
						printf (
							"..in STACK_IS_FULL ( thisStack=%p) TOP=%d FLOOR=%d, returning %d\n", 
							thisStack , thisStack -> _STACK_TOP, thisStack -> _STACK_FLOOR, returnValue
							) ; 
						}
	return ( returnValue );
}

/* return 1 if the stack is empty, otherwise return 0
	add ERROR CHECKING
*/
int STACK_IS_EMPTY ( TYPE_STACK * thisStack ) {
	int top =thisStack -> _STACK_TOP ;
	int bottom = thisStack -> _STACK_FLOOR ;
	int returnValue = ( top == bottom ) ;
					if ( DEBUG ) { 
						printf ( "..in STACK_IS_EMPTY ( thisStack=%p) TOP=%d FLOOR=%d, returning %d\n", 
						thisStack , top , bottom , returnValue 
						) ;
						}
	return ( returnValue ) ; 
	}

/* STACK_LIFT 
	lift the floor of the stack, if possible. else report an error
	add ERROR CHECKING
*/
int STACK_LIFT ( TYPE_STACK * thisStack ) {
					if ( DEBUG ) { 
						printf ( "..in STACK_LIFT ( thisStack=%p)\n", 
						thisStack , 
						thisStack -> _STACK_FLOOR + 1
						) ;
						}
	thisStack -> _STACK_FLOOR++;
	return 0;
	}

/* STACK_SWAP 
 	swap data elPTRs at offsets from floor of from and to
	NOTE: indexes ar absolute (from foundation, not floor)
	NOTE: change to XOR for speed
*/ 
int STACK_SWAP ( TYPE_STACK * thisStack, unsigned from, unsigned to ) {
					if ( DEBUG ) { 
						printf ( "..in STACK_SWAP ( thisStack=%p from=%u to=%u)\n", thisStack , from, to ) ; 
						}
	DDT temp_value = ( DDT ) thisStack -> _STACK_DATA[ to ];
	( thisStack -> _STACK_DATA[ to ] ) = thisStack -> _STACK_DATA[ from ];
	( thisStack -> _STACK_DATA[ from ] ) = temp_value;
	return 0;
	}

/* return index of the top element (the element itself, not the empty slot following it)
*/
unsigned STACK_TOP_IS ( TYPE_STACK * thisStack ) {
						if ( DEBUG ) { 
						printf ( "..in STACK_TOP_IS ( thisStack=%p) returning %u\n", 
							thisStack , thisStack -> _STACK_TOP
							) ;
							}
	return ( thisStack -> _STACK_TOP - 1) ; 
	}

/* return index of the bottom element, which is the first one in the stack relative to foundation	
*/
unsigned STACK_BOTTOM_IS ( TYPE_STACK * thisStack ) {
						if ( DEBUG ) { 
						printf ( "..in STACK_BOTTOM_IS ( thisStack=%p) returning %u\n", 
							thisStack , thisStack -> _STACK_FLOOR
							) ;
							}
	return ( thisStack -> _STACK_FLOOR ) ; 
	}

/* return index of the foudnation element, which is the first in the array, regardless of bottom	
*/
unsigned STACK_FOUNDATION_IS ( TYPE_STACK * thisStack ) {
						if ( DEBUG ) { 
						printf ( "..in STACK_FOUDATION_IS ( thisStack=%p) returning %u\n", 
							thisStack , thisStack -> _STACK_FOUNDATION
							) ;
							}
	return ( thisStack -> _STACK_FOUNDATION ) ; 
	}

/* display the contents of the STACK */
int STACK_DUMP ( TYPE_STACK * thisStack ) {
	int errCode = 0 ;
						if ( DEBUG ) { printf (">>in STACK_DUMP ( %p ) \n", thisStack ) ; }
	printf ( "_STACK_FOUNDATION %d\t_STACK_FLOOR %d\t_STACK_TOP %d\t_MAX_SIZE %d\tBytes %d\n",
		thisStack -> _STACK_FOUNDATION,
		thisStack -> _STACK_FLOOR,
		thisStack -> _STACK_TOP,
		thisStack -> _MAX_SIZE,
		thisStack -> _STACK_BYTES
		) ;
	for (i=thisStack -> _STACK_FOUNDATION; i < thisStack -> _STACK_TOP; i++) {
		printf ( "DATA[%d]\n" , i ) ;
		DDT nextElement = (DDT) ( thisStack -> _STACK_DATA[i] );
		if ( !DDF ( nextElement ) ) { errCode = 1; }
		printf ( ";\n" ) ; 
					if ( DEBUG ) { 
					printf ( " _STACK_DATA [ %u ] = %p -> %p \n" , i, thisStack , nextElement ) ; 
					} ;
		} ;
	printf("\n");
					if ( DEBUG ) { 
						printf ("<< leaving STACK_DUMP ( %p ) returning errCode %d\n", 
						thisStack , errCode
						 ) ; 
						}
	return errCode;
	}
	

