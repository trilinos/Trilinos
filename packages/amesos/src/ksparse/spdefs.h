/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

/*
 *  DATA STRUCTURE AND MACRO DEFINITIONS for Sparse.
 *
 *  Author:                     Advising professor:
 *      Kenneth S. Kundert          Alberto Sangiovanni-Vincentelli
 *      UC Berkeley
 *
 *  This file contains common type definitions and macros for the sparse
 *  matrix routines.  These definitions are of no interest to the user.
 */


/*
 *  Revision and copyright information.
 *
 *  Copyright (c) 1985,86,87,88,89,90
 *  by Kenneth S. Kundert and the University of California.
 *
 *  Permission to use, copy, modify, and distribute this software and
 *  its documentation for any purpose and without fee is hereby granted,
 *  provided that the copyright notices appear in all copies and
 *  supporting documentation and that the authors and the University of
 *  California are properly credited.  The authors and the University of
 *  California make no representations as to the suitability of this
 *  software for any purpose.  It is provided `as is', without express
 *  or implied warranty.
 *
 *  $Date$
 *  $Revision$
 */


/* Density of index arrays for fast access to arbitrary elements in matrix */
/* Tried 2 and 4, but the speed up was negligible, and memory use higher */
#define IND_DENSITY 1

static long padsize, padshift;

/*
 *  IMPORTS
 */

#include <stdio.h>
#include "misc.h"
#include "sppars.h"


/*
 *   MACRO DEFINITIONS
 *
 *   Macros are distinguished by using solely capital letters in their
 *   identifiers.  This contrasts with C defined identifiers which are strictly
 *   lower case, and program variable and procedure names which use both upper
 *   and lower case.
 */

/* Begin macros. */

/* Boolean data type */
#define  BOOLEAN        int
#define  NO             0
#define  YES            1
#define  NOT            !
#define  AND            &&
#define  OR             ||

/* NULL pointer */
#ifndef  NULL
#define  NULL           0
#endif

#define  SPARSE_ID      0x772773        /* Arbitrary (is Sparse on phone). */
#define  IS_SPARSE(matrix)      ((matrix) != NULL &&            \
                                 (matrix)->ID == SPARSE_ID)
#define  IS_VALID(matrix)       ((matrix) != NULL &&            \
                                 (matrix)->ID == SPARSE_ID &&   \
                                 (matrix)->Error >= spOKAY &&   \
                                 (matrix)->Error < spFATAL)
#define  IS_FACTORED(matrix)    ((matrix)->Factored && !(matrix)->NeedsOrdering)

/* Macro commands */
/* Macro functions that return the maximum or minimum independent of type. */
#define  MAX(a,b)           ((a) > (b) ? (a) : (b))
#define  MIN(a,b)           ((a) < (b) ? (a) : (b))

/* Macro function that returns the absolute value of a floating point number. */
#define  ABS(a)             ((a) < 0 ? -(a) : (a))

/* Macro function that returns the square of a number. */
#define  SQR(a)             ((a)*(a))

/* Macro procedure that swaps two entities. */
#define  SWAP(type, a, b)   {type swapx; swapx = a; a = b; b = swapx;}

/* Macro function that returns the approx absolute value of a complex number. */
#if spCOMPLEX
#define  ELEMENT_MAG(ptr)   (ABS((ptr)->Real) + ABS((ptr)->Imag))
#else
#define  ELEMENT_MAG(ptr)   ((ptr)->Real < 0.0 ? -(ptr)->Real : (ptr)->Real)
#endif

/* Complex assignment statements. */
#define  CMPLX_ASSIGN(to,from)  \
{   (to).Real = (from).Real;    \
    (to).Imag = (from).Imag;    \
}
#define  CMPLX_CONJ_ASSIGN(to,from)     \
{   (to).Real = (from).Real;            \
    (to).Imag = -(from).Imag;           \
}
#define  CMPLX_NEGATE_ASSIGN(to,from)   \
{   (to).Real = -(from).Real;           \
    (to).Imag = -(from).Imag;           \
}
#define  CMPLX_CONJ_NEGATE_ASSIGN(to,from)      \
{   (to).Real = -(from).Real;                   \
    (to).Imag = (from).Imag;                    \
}
#define  CMPLX_CONJ(a)  (a).Imag = -(a).Imag
#define  CMPLX_NEGATE(a)        \
{   (a).Real = -(a).Real;       \
    (a).Imag = -(a).Imag;       \
}

/* Macro that returns the approx magnitude (L-1 norm) of a complex number. */
#define  CMPLX_1_NORM(a)        (ABS((a).Real) + ABS((a).Imag))

/* Macro that returns the approx magnitude (L-infinity norm) of a complex. */
#define  CMPLX_INF_NORM(a)      (MAX (ABS((a).Real),ABS((a).Imag)))

/* Macro function that returns the magnitude (L-2 norm) of a complex number. */
#define  CMPLX_2_NORM(a)        (sqrt((a).Real*(a).Real + (a).Imag*(a).Imag))

/* Macro function that performs complex addition. */
#define  CMPLX_ADD(to,from_a,from_b)            \
{   (to).Real = (from_a).Real + (from_b).Real;  \
    (to).Imag = (from_a).Imag + (from_b).Imag;  \
}

/* Macro function that performs complex subtraction. */
#define  CMPLX_SUBT(to,from_a,from_b)           \
{   (to).Real = (from_a).Real - (from_b).Real;  \
    (to).Imag = (from_a).Imag - (from_b).Imag;  \
}

/* Macro function that is equivalent to += operator for complex numbers. */
#define  CMPLX_ADD_ASSIGN(to,from)      \
{   (to).Real += (from).Real;           \
    (to).Imag += (from).Imag;           \
}

/* Macro function that is equivalent to -= operator for complex numbers. */
#define  CMPLX_SUBT_ASSIGN(to,from)     \
{   (to).Real -= (from).Real;           \
    (to).Imag -= (from).Imag;           \
}

/* Macro function that multiplies a complex number by a scalar. */
#define  SCLR_MULT(to,sclr,cmplx)       \
{   (to).Real = (sclr) * (cmplx).Real;  \
    (to).Imag = (sclr) * (cmplx).Imag;  \
}

/* Macro function that multiply-assigns a complex number by a scalar. */
#define  SCLR_MULT_ASSIGN(to,sclr)      \
{   (to).Real *= (sclr);                \
    (to).Imag *= (sclr);                \
}

/* Macro function that multiplies two complex numbers. */
#define  CMPLX_MULT(to,from_a,from_b)           \
{   (to).Real = (from_a).Real * (from_b).Real - \
                (from_a).Imag * (from_b).Imag;  \
    (to).Imag = (from_a).Real * (from_b).Imag + \
                (from_a).Imag * (from_b).Real;  \
}

/* Macro function that implements to *= from for complex numbers. */
#define  CMPLX_MULT_ASSIGN(to,from)             \
{   RealNumber to_real_ = (to).Real;            \
    (to).Real = to_real_ * (from).Real -        \
                (to).Imag * (from).Imag;        \
    (to).Imag = to_real_ * (from).Imag +        \
                (to).Imag * (from).Real;        \
}

/* Macro function that multiplies two complex numbers, the first of which is
 * conjugated. */
#define  CMPLX_CONJ_MULT(to,from_a,from_b)      \
{   (to).Real = (from_a).Real * (from_b).Real + \
                (from_a).Imag * (from_b).Imag;  \
    (to).Imag = (from_a).Real * (from_b).Imag - \
                (from_a).Imag * (from_b).Real;  \
}

/* Macro function that multiplies two complex numbers and then adds them
 * to another. to = add + mult_a * mult_b */
#define  CMPLX_MULT_ADD(to,mult_a,mult_b,add)                   \
{   (to).Real = (mult_a).Real * (mult_b).Real -                 \
                (mult_a).Imag * (mult_b).Imag + (add).Real;     \
    (to).Imag = (mult_a).Real * (mult_b).Imag +                 \
                (mult_a).Imag * (mult_b).Real + (add).Imag;     \
}

/* Macro function that subtracts the product of two complex numbers from
 * another.  to = subt - mult_a * mult_b */
#define  CMPLX_MULT_SUBT(to,mult_a,mult_b,subt)                 \
{   (to).Real = (subt).Real - (mult_a).Real * (mult_b).Real +   \
                              (mult_a).Imag * (mult_b).Imag;    \
    (to).Imag = (subt).Imag - (mult_a).Real * (mult_b).Imag -   \
                              (mult_a).Imag * (mult_b).Real;    \
}

/* Macro function that multiplies two complex numbers and then adds them
 * to another. to = add + mult_a* * mult_b where mult_a* represents mult_a
 * conjugate. */
#define  CMPLX_CONJ_MULT_ADD(to,mult_a,mult_b,add)              \
{   (to).Real = (mult_a).Real * (mult_b).Real +                 \
                (mult_a).Imag * (mult_b).Imag + (add).Real;     \
    (to).Imag = (mult_a).Real * (mult_b).Imag -                 \
                (mult_a).Imag * (mult_b).Real + (add).Imag;     \
}

/* Macro function that multiplies two complex numbers and then adds them
 * to another. to += mult_a * mult_b */
#define  CMPLX_MULT_ADD_ASSIGN(to,from_a,from_b)        \
{   (to).Real += (from_a).Real * (from_b).Real -        \
                 (from_a).Imag * (from_b).Imag;         \
    (to).Imag += (from_a).Real * (from_b).Imag +        \
                 (from_a).Imag * (from_b).Real;         \
}

/* Macro function that multiplies two complex numbers and then subtracts them
 * from another. */
#define  CMPLX_MULT_SUBT_ASSIGN(to,from_a,from_b)       \
{   (to).Real -= (from_a).Real * (from_b).Real -        \
                 (from_a).Imag * (from_b).Imag;         \
    (to).Imag -= (from_a).Real * (from_b).Imag +        \
                 (from_a).Imag * (from_b).Real;         \
}

/* Macro function that multiplies two complex numbers and then adds them
 * to the destination. to += from_a* * from_b where from_a* represents from_a
 * conjugate. */
#define  CMPLX_CONJ_MULT_ADD_ASSIGN(to,from_a,from_b)   \
{   (to).Real += (from_a).Real * (from_b).Real +        \
                 (from_a).Imag * (from_b).Imag;         \
    (to).Imag += (from_a).Real * (from_b).Imag -        \
                 (from_a).Imag * (from_b).Real;         \
}

/* Macro function that multiplies two complex numbers and then subtracts them
 * from the destination. to -= from_a* * from_b where from_a* represents from_a
 * conjugate. */
#define  CMPLX_CONJ_MULT_SUBT_ASSIGN(to,from_a,from_b)  \
{   (to).Real -= (from_a).Real * (from_b).Real +        \
                 (from_a).Imag * (from_b).Imag;         \
    (to).Imag -= (from_a).Real * (from_b).Imag -        \
                 (from_a).Imag * (from_b).Real;         \
}

/*
 * Macro functions that provide complex division.
 */

/* Complex division:  to = num / den */
#define CMPLX_DIV(to,num,den)                                           \
{   RealNumber  r_, s_;                                                 \
    if (((den).Real >= (den).Imag AND (den).Real > -(den).Imag) OR      \
        ((den).Real < (den).Imag AND (den).Real <= -(den).Imag))        \
    {   r_ = (den).Imag / (den).Real;                                   \
        s_ = (den).Real + r_*(den).Imag;                                \
        (to).Real = ((num).Real + r_*(num).Imag)/s_;                    \
        (to).Imag = ((num).Imag - r_*(num).Real)/s_;                    \
    }                                                                   \
    else                                                                \
    {   r_ = (den).Real / (den).Imag;                                   \
        s_ = (den).Imag + r_*(den).Real;                                \
        (to).Real = (r_*(num).Real + (num).Imag)/s_;                    \
        (to).Imag = (r_*(num).Imag - (num).Real)/s_;                    \
    }                                                                   \
}

/* Complex division and assignment:  num /= den */
#define CMPLX_DIV_ASSIGN(num,den)                                       \
{   RealNumber  r_, s_, t_;                                             \
    if (((den).Real >= (den).Imag AND (den).Real > -(den).Imag) OR      \
        ((den).Real < (den).Imag AND (den).Real <= -(den).Imag))        \
    {   r_ = (den).Imag / (den).Real;                                   \
        s_ = (den).Real + r_*(den).Imag;                                \
        t_ = ((num).Real + r_*(num).Imag)/s_;                           \
        (num).Imag = ((num).Imag - r_*(num).Real)/s_;                   \
        (num).Real = t_;                                                \
    }                                                                   \
    else                                                                \
    {   r_ = (den).Real / (den).Imag;                                   \
        s_ = (den).Imag + r_*(den).Real;                                \
        t_ = (r_*(num).Real + (num).Imag)/s_;                           \
        (num).Imag = (r_*(num).Imag - (num).Real)/s_;                   \
        (num).Real = t_;                                                \
    }                                                                   \
}

/* Complex reciprocation:  to = 1.0 / den */
#define CMPLX_RECIPROCAL(to,den)                                        \
{   RealNumber  r_;                                                     \
    if (((den).Real >= (den).Imag AND (den).Real > -(den).Imag) OR      \
        ((den).Real < (den).Imag AND (den).Real <= -(den).Imag))        \
    {   r_ = (den).Imag / (den).Real;                                   \
        (to).Imag = -r_*((to).Real = 1.0/((den).Real + r_*(den).Imag)); \
    }                                                                   \
    else                                                                \
    {   r_ = (den).Real / (den).Imag;                                   \
        (to).Real = -r_*((to).Imag = -1.0/((den).Imag + r_*(den).Real));\
    }                                                                   \
}






/*
 *  ASSERT and ABORT
 *
 *  Macro used to assert that if the code is working correctly, then 
 *  a condition must be true.  If not, then execution is terminated
 *  and an error message is issued stating that there is an internal
 *  error and giving the file and line number.  These assertions are
 *  not evaluated unless the DEBUG flag is true.
 */

#if DEBUG
#define ASSERT(condition) if (NOT(condition)) ABORT()
#else
#define ASSERT(condition)
#endif

#ifndef ABORT
#if DEBUG
#define  ABORT()                                                        \
{   (void)fflush(stdout);                                               \
    (void)fprintf(stderr, "sparse: panic in file `%s' at line %d.\n",   \
            __FILE__, __LINE__);                                        \
    (void)fflush(stderr);                                               \
    abort();                                                            \
}
#else
#define  ABORT()
#endif
#endif





/*
 *  IMAGINARY VECTORS
 *
 *  The imaginary vectors iRHS and iSolution are only needed when the
 *  options spCOMPLEX and spSEPARATED_COMPLEX_VECTORS are set.  The following
 *  macro makes it easy to include or exclude these vectors as needed.
 */

#if spCOMPLEX AND spSEPARATED_COMPLEX_VECTORS
#define IMAG_VECTORS    , iRHS, iSolution
#define IMAG_RHS        , iRHS
#else
#define IMAG_VECTORS
#define IMAG_RHS
#endif

#ifdef SHARED_MEM
#include "shared_mem.h"
#define ALLOC(type,number)  ((type *)amalloc_SM((unsigned)(sizeof(type)*(number))))
#define PALLOC(type,number)  ((type *)amalloc_SM((unsigned)(padsize*(number+1))))
#ifndef REALLOC
#define REALLOC(ptr,type,number)  \
        ptr = (type *)arealloc_SM((void *)ptr,(unsigned)(sizeof(type)*(number)))
#endif
#ifndef FREE
#define FREE(ptr) { if ((ptr) != NULL) afree_SM((void *)(ptr)); (ptr) = NULL; }
#endif
#else
#define ALLOC(type,number)  ((type *)tmalloc((unsigned)(sizeof(type)*(number))))
#define PALLOC(type,number)  ((type *)tmalloc((unsigned)(padsize*(number+1))))
#ifndef REALLOC
#define REALLOC(ptr,type,number)  \
           ptr = (type *)trealloc((char *)ptr,(unsigned)(sizeof(type)*(number)))
#endif
#ifndef FREE
#define FREE(ptr) { if ((ptr) != NULL) txfree((char *)(ptr)); (ptr) = NULL; }
#endif
#endif /* SHARED_MEM */


/* Calloc that properly handles allocating a cleared vector. */
#define CALLOC(ptr,type,number)                         \
{   int i; ptr = ALLOC(type, number);                   \
    if (ptr != (type *)NULL)                            \
        for(i=(number)-1;i>=0; i--) ptr[i] = (type) 0;  \
}







/*
 *  REAL NUMBER
 */

/* Begin `RealNumber'. */

typedef  spREAL  RealNumber, *RealVector;








/*
 *  COMPLEX NUMBER DATA STRUCTURE
 *
 *  >>> Structure fields:
 *  Real  (RealNumber)
 *      The real portion of the number.  Real must be the first
 *      field in this structure.
 *  Imag  (RealNumber)
 *      The imaginary portion of the number. This field must follow
 *      immediately after Real.
 */

/* Begin `ComplexNumber'. */

typedef  struct
{   RealNumber  Real;
    RealNumber  Imag;
} ComplexNumber, *ComplexVector;


#include "spmat.h"







/*
 *  Function declarations
 */

#ifdef __STDC__
extern ElementPtr spcGetElement( MatrixPtr, int, int );
extern ElementPtr spcGetFillin( MatrixPtr, int, int );
extern ElementPtr spcFindElementInCol( MatrixPtr, ElementPtr*, int, int, int );
extern ElementPtr spcCreateElement( MatrixPtr, int, int, ElementPtr*, int );
extern void spcCreateInternalVectors( MatrixPtr );
extern void spcLinkRows( MatrixPtr );
extern void spcColExchange( MatrixPtr, int, int );
extern void spcRowExchange( MatrixPtr, int, int );
#else /* __STDC__ */
extern ElementPtr spcGetElement();
extern ElementPtr spcGetFillin();
extern ElementPtr spcFindElementInCol();
extern ElementPtr spcCreateElement();
extern void spcCreateInternalVectors();
extern void spcLinkRows();
extern void spcColExchange();
extern void spcRowExchange();
#endif /* __STDC__ */
