/*
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 */

/* 
 * This header file is to be included in source files z*.c
 */
#ifndef __SUPERLU_DCOMPLEX /* allow multiple inclusions */
#define __SUPERLU_DCOMPLEX

#include <mpi.h>

typedef struct { double r, i; } doublecomplex;

/*
 * These variables will be defined to be MPI datatypes for complex
 * and double complex. I'm too lazy to declare
 * these guys external in every file that needs them.
 */
extern MPI_Datatype SuperLU_MPI_DOUBLE_COMPLEX;


/* Macro definitions */

/* Complex Addition c = a + b */
#define z_add(c, a, b) { (c)->r = (a)->r + (b)->r; \
			 (c)->i = (a)->i + (b)->i; }

/* Complex Subtraction c = a - b */
#define z_sub(c, a, b) { (c)->r = (a)->r - (b)->r; \
			 (c)->i = (a)->i - (b)->i; }

/* Complex-Double Multiplication */
#define zd_mult(c, a, b) { (c)->r = (a)->r * (b); \
                           (c)->i = (a)->i * (b); }

/* Complex-Complex Multiplication */
#define zz_mult(c, a, b) { \
	double cr, ci; \
    	cr = (a)->r * (b)->r - (a)->i * (b)->i; \
    	ci = (a)->i * (b)->r + (a)->r * (b)->i; \
    	(c)->r = cr; \
    	(c)->i = ci; \
    }

/* Complex equality testing */
#define z_eq(a, b)  ( (a)->r == (b)->r && (a)->i == (b)->i )


#ifdef __cplusplus
extern "C" {
#endif

/* Prototypes for functions in dcomplex.c */
void   z_div(doublecomplex *, doublecomplex *, doublecomplex *);
double z_abs(doublecomplex *);     /* exact */
double z_abs1(doublecomplex *);    /* approximate */
void   z_exp(doublecomplex *, doublecomplex *);
void   d_cnjg(doublecomplex *r, doublecomplex *z);
double d_imag(doublecomplex *);


#ifdef __cplusplus
  }
#endif


#endif  /* __SUPERLU_DCOMPLEX */
