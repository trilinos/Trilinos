#ifndef __SUPERLU_UTIL /* allow multiple inclusions */
#define __SUPERLU_UTIL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <assert.h>

/*
 * Macros
 */
#ifndef USER_ABORT
#define USER_ABORT(msg) superlu_abort_and_exit_dist(msg)
#endif

#define ABORT(err_msg) \
 { char msg[256];\
   sprintf(msg,"%s at line %d in file %s\n",err_msg,__LINE__, __FILE__);\
   USER_ABORT(msg); }


#ifndef USER_MALLOC
#define USER_MALLOC(size) superlu_malloc_dist(size)
#endif

#define SUPERLU_MALLOC(size) USER_MALLOC(size)

#ifndef USER_FREE
#define USER_FREE(addr) superlu_free_dist(addr)
#endif

#define SUPERLU_FREE(addr) USER_FREE(addr)

#define CHECK_MALLOC(pnum, where) {                 \
    extern int_t superlu_malloc_total;        \
    printf("(%2d) %s: malloc_total %d Bytes\n",     \
	   pnum, where, superlu_malloc_total); \
}

#define MAX(x, y) 	( (x) > (y) ? (x) : (y) )
#define MIN(x, y) 	( (x) < (y) ? (x) : (y) )

/* 
 * Constants 
 */
#define EMPTY	(-1)
#ifndef FALSE
#define FALSE	(0)
#endif
#ifndef TRUE
#define TRUE	(1)
#endif

/*
 * Type definitions
 */
typedef float    flops_t;
typedef unsigned char Logical;
#ifdef _CRAY
#define int short
#endif

/* 
 * The following enumerate type is used by the statistics variable 
 * SuperLUStat, to keep track of flop count and time spent at various stages.
 *
 * Note that not all of the fields are disjoint.
 */
typedef enum {
    COLPERM, /* find a column ordering that minimizes fills */
    ROWPERM, /* find a row ordering maximizes diagonal. */
    RELAX,   /* find artificial supernodes */
    ETREE,   /* compute column etree */
    EQUIL,   /* equilibrate the original matrix */
    SYMBFAC, /* symbolic factorization. */
    DIST,    /* distribute matrix. */
    FACT,    /* perform LU factorization */
    COMM,    /* communication for factorization */
    SOL_COMM,/* communication for solve */
    RCOND,   /* estimate reciprocal condition number */
    SOLVE,   /* forward and back solves */
    REFINE,  /* perform iterative refinement */
    FLOAT,   /* time spent in floating-point operations */
    TRSV,    /* fraction of FACT spent in xTRSV */
    GEMV,    /* fraction of FACT spent in xGEMV */
    FERR,    /* estimate error bounds after iterative refinement */
    NPHASES  /* total number of phases */
} PhaseType;

typedef struct {
    int     *panel_histo; /* histogram of panel size distribution */
    double  *utime;       /* running time at various phases */
    flops_t *ops;         /* operation count at various phases */
    int     TinyPivots;   /* number of tiny pivots */
    int     RefineSteps;  /* number of iterative refinement steps */
} SuperLUStat_t;

/* Headers for 2 types of dynamatically managed memory */
typedef struct e_node {
    int size;      /* length of the memory that has been used */
    void *mem;     /* pointer to the new malloc'd store */
} ExpHeader;

typedef struct {
    int  size;
    int  used;
    int  top1;  /* grow upward, relative to &array[0] */
    int  top2;  /* grow downward */
    void *array;
} LU_stack_t;

/* Constants */
#define GluIntArray(n)   (5 * (n) + 5)
#define NO_MEMTYPE  4      /* 0: lusup;
			      1: ucol;
			      2: lsub;
			      3: usub */

#if 0
/* Macros to manipulate stack */
#define StackFull(x)         ( x + stack.used >= stack.size )
#define NotDoubleAlign(addr) ( (long)addr & 7 )
#define DoubleAlign(addr)    ( ((long)addr + 7) & ~7L )
#define TempSpace(n, w)      ( (2*w + 4 + NO_MARKER)*m*sizeof(int) + \
			      (w + 1)*n*sizeof(double) )
#define DSS_Reduce(alpha)        ((alpha + 1) / 2)  /* i.e. (alpha-1)/2 + 1 */

#define FIRSTCOL_OF_SNODE(i)	(xsup[i])
#endif


#if ( PROFlevel>=1 )
#define TIC(t)          t = SuperLU_timer_()
#define TOC(t2, t1)     t2 = SuperLU_timer_() - t1
#else
#define TIC(t)
#define TOC(t2, t1)
#endif

#endif /* __SUPERLU_UTIL */
