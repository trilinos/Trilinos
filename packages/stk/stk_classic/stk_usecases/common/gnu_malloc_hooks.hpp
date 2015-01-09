/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef _gnu_malloc_hooks_hpp_
#define _gnu_malloc_hooks_hpp_

//'gnu malloc hooks' are functions that give access to some statistics
//that are collected/logged during calls to malloc/realloc/free.
//These are only available when using g++/gcc.
//All prototypes and code for these is completely protected by the
//macro USE_GNU_MALLOC_HOOKS. If that macro is not defined, then the
//compiler doesn't see anything associated with gnu malloc hooks.

//if the macro USE_GNU_MALLOC_HOOKS is defined, then declare prototypes
//for the malloc-hook functions.
//Otherwise, #define them to be nothing so they vanish the same way
//'assert' vanishes when NDEBUG is defined.

#if defined(USE_GNU_MALLOC_HOOKS)

void disable_gnu_malloc_hooks();

void reset_malloc_stats();

/** Return how many MB of heap memory has been allocated. */
double alloc_MB();

/** Return how many separate blocks the heap memory is allocated in. */
unsigned alloc_blks();

/** Return how many MB of heap memory has been freed. */
double freed_MB();

/** Return how many blocks of heap memory has been freed. */
unsigned freed_blks();

#else
//USE_GNU_MALLOC_HOOKS is not defined, so #define
//the function prototypes to be nothing.

#define disable_gnu_malloc_hooks() ((void)0)
#define reset_malloc_stats() ((void)0)
#define alloc_MB() (0.0)
#define alloc_blks() (0)
#define freed_MB() (0.0)
#define freed_blks() (0)

#endif

#endif

