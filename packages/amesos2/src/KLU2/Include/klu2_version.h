// @HEADER
// *****************************************************************************
//                   KLU2: A Direct Linear Solver package
//
// Copyright 2011 NTESS and the KLU2 contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _TKLU_VERSION_H
#define _TKLU_VERSION_H

// TODO : Check SPLIT complex 

#ifdef DLONG
#define Int_id UF_long_id
#define Int_MAX UF_long_max
#else
#define Int_id "%d"
#define Int_MAX INT_MAX
#endif

#define NPRINT  

#define KLU2_BYTES(type,n) (sizeof (type) * (n))
#define KLU2_CEILING(b,u)  (((b)+(u)-1) / (u))
#define UNITS(type,n) (KLU2_CEILING (KLU2_BYTES (type,n), sizeof (Unit)))
#define DUNITS(type,n) (ceil (KLU2_BYTES (type, (double) n) / sizeof (Unit)))

#define GET_I_POINTER(LU, Xip, Xi, k) \
{ \
    Xi = (Int *) (LU + Xip [k]) ; \
}

#define GET_X_POINTER(LU, Xip, Xlen, Xx, k) \
{ \
    Xx = (Entry *) (LU + Xip [k] + UNITS (Int, Xlen [k])) ; \
}

#define GET_POINTER(LU, Xip, Xlen, Xi, Xx, k, xlen) \
{ \
    Unit *xp = LU + Xip [k] ; \
    xlen = Xlen [k] ; \
    Xi = (Int *) xp ; \
    Xx = (Entry *) (xp + UNITS (Int, xlen)) ; \
}

/* function names */
#define KLU_scale klu_scale
#define KLU_solve klu_solve
#define KLU_tsolve klu_tsolve
#define KLU_solve2 klu_solve2
#define KLU_tsolve2 klu_tsolve2
#define KLU_free_numeric klu_free_numeric
#define KLU_factor klu_factor
#define KLU_refactor klu_refactor
#define KLU_kernel_factor klu_kernel_factor 
#define KLU_lsolve klu_lsolve
#define KLU_ltsolve klu_ltsolve
#define KLU_usolve klu_usolve
#define KLU_utsolve klu_utsolve
#define KLU_kernel klu_kernel
#define KLU_valid klu_valid
#define KLU_valid_LU klu_valid_LU
#define KLU_sort klu_sort
#define KLU_rgrowth klu_rgrowth
#define KLU_rcond klu_rcond
#define KLU_extract klu_extract
#define KLU_condest klu_condest
#define KLU_flops klu_flops

#define KLU_analyze klu_analyze
#define KLU_analyze_given klu_analyze_given
#define KLU_alloc_symbolic klu_alloc_symbolic
#define KLU_free_symbolic klu_free_symbolic
#define KLU_defaults klu_defaults
#define KLU_free klu_free
#define KLU_malloc klu_malloc
#define KLU_realloc klu_realloc
#define KLU_add_size_t klu_add_size_t
#define KLU_mult_size_t klu_mult_size_t

#define KLU_symbolic klu_symbolic
#define KLU_numeric klu_numeric
#define KLU_common klu_common


/* -------------------------------------------------------------------------- */
/* Numerical relop macros for correctly handling the NaN case */
/* -------------------------------------------------------------------------- */

/*
SCALAR_IS_NAN(x):
    True if x is NaN.  False otherwise.  The commonly-existing isnan(x)
    function could be used, but it's not in Kernighan & Ritchie 2nd edition
    (ANSI C).  It may appear in <math.h>, but I'm not certain about
    portability.  The expression x != x is true if and only if x is NaN,
    according to the IEEE 754 floating-point standard.

SCALAR_IS_ZERO(x):
    True if x is zero.  False if x is nonzero, NaN, or +/- Inf.
    This is (x == 0) if the compiler is IEEE 754 compliant.

SCALAR_IS_NONZERO(x):
    True if x is nonzero, NaN, or +/- Inf.  False if x zero.
    This is (x != 0) if the compiler is IEEE 754 compliant.

SCALAR_IS_LTZERO(x):
    True if x is < zero or -Inf.  False if x is >= 0, NaN, or +Inf.
    This is (x < 0) if the compiler is IEEE 754 compliant.
*/

/* These all work properly, according to the IEEE 754 standard ... except on */
/* a PC with windows.  Works fine in Linux on the same PC... */
#define SCALAR_IS_NAN(x)        ((x) != (x))
#define SCALAR_IS_ZERO(x)       ((x) == 0.)
#define SCALAR_IS_NONZERO(x)    ((x) != 0.)
#define SCALAR_IS_LTZERO(x)     ((x) < 0.)


/* scalar absolute value macro. If x is NaN, the result is NaN: */
#define SCALAR_ABS(x) ((SCALAR_IS_LTZERO (x)) ? -(x) : (x))

/* print a scalar (avoid printing "-0" for negative zero).  */
#ifdef NPRINT
#define PRINT_SCALAR(a)
#else
#define PRINT_SCALAR(a) \
{ \
    if (SCALAR_IS_NONZERO (a)) \
    { \
        PRINTF ((" (%g)", (a))) ; \
    } \
    else \
    { \
        PRINTF ((" (0)")) ; \
    } \
}
#endif

/* -------------------------------------------------------------------------- */
/* Real floating-point arithmetic */
/* -------------------------------------------------------------------------- */

#ifndef COMPLEX

typedef double Unit ;
/*#define Entry double*/
/* TODO: Need to add namespace to these methods */

#define SPLIT(s)                    (1)
#define REAL(c)                     (Teuchos::ScalarTraits<Entry>::real(c))
#define IMAG(c)                     (Teuchos::ScalarTraits<Entry>::imag(c))
#define CLEAR(c)                    { (c) = 0. ; }
#define CLEAR_AND_INCREMENT(p)      { *p++ = 0. ; }
#define IS_NAN(a)                   SCALAR_IS_NAN (a) /* TODO : ???*/
#define IS_ZERO(a)                  SCALAR_IS_ZERO (a)
#define IS_NONZERO(a)               SCALAR_IS_NONZERO (a)
#define SCALE_DIV(c,s)              { (c) /= (s) ; }
#define SCALE_DIV_ASSIGN(a,c,s)     { a = c / s ; }
#define SCALE(c,s)                  { (c) *= (s) ; }
#define ASSEMBLE(c,a)               { (c) += (a) ; }
#define ASSEMBLE_AND_INCREMENT(c,p) { (c) += *p++ ; }
#define DECREMENT(c,a)              { (c) -= (a) ; }
#define MULT(c,a,b)                 { (c) = (a) * (b) ; }
#define MULT_CONJ(c,a,b)            { (c) = (a) * Teuchos::ScalarTraits<Entry>::conjugate(b) ; }
#define MULT_SUB(c,a,b)             { (c) -= (a) * (b) ; }
#define MULT_SUB_CONJ(c,a,b)        { (c) -= (a) * Teuchos::ScalarTraits<Entry>::conjugate(b) ; }
#define DIV(c,a,b)                  { (c) = KLU_ScalarTraits<Entry>::divide(a, b) ; }
#define RECIPROCAL(c)               { (c) = KLU_ScalarTraits<Entry>::reciprocal(c) ; }
#define DIV_CONJ(c,a,b)             { (c) = KLU_ScalarTraits<Entry>::divideConjugate(a, b) ; }
#define APPROX_ABS(s,a)             { (s) =  KLU_ScalarTraits<Entry>::approxABS(a) ; }
#define KLU2_ABS(s,a)                    { (s) =  KLU_ScalarTraits<Entry>::abs(a) ; }
#define PRINT_ENTRY(a)              PRINT_SCALAR (a)
#define CONJ(a,x)                   a = (Teuchos::ScalarTraits<Entry>::conjugate(x))

/* for flop counts */
#define MULTSUB_FLOPS   2.      /* c -= a*b */
#define DIV_FLOPS       1.      /* c = a/b */
#define ABS_FLOPS       0.      /* c = abs (a) */
#define ASSEMBLE_FLOPS  1.      /* c += a */
#define DECREMENT_FLOPS 1.      /* c -= a */
#define MULT_FLOPS      1.      /* c = a*b */
#define SCALE_FLOPS     1.      /* c = a/s */

#else

/* -------------------------------------------------------------------------- */
/* Complex floating-point arithmetic */
/* -------------------------------------------------------------------------- */

/*
    Note:  An alternative to this Double_Complex type would be to use a
    struct { double r ; double i ; }.  The problem with that method
    (used by the Sun Performance Library, for example) is that ANSI C provides
    no guarantee about the layout of a struct.  It is possible that the sizeof
    the struct above would be greater than 2 * sizeof (double).  This would
    mean that the complex BLAS could not be used.  The method used here avoids
    that possibility.  ANSI C *does* guarantee that an array of structs has
    the same size as n times the size of one struct.

    The ANSI C99 version of the C language includes a "double _Complex" type.
    It should be possible in that case to do the following:

    #define Entry double _Complex

    and remove the Double_Complex struct.  The macros, below, could then be
    replaced with instrinsic operators.  Note that the #define Real and
    #define Imag should also be removed (they only appear in this file).

    For the MULT, MULT_SUB, MULT_SUB_CONJ, and MULT_CONJ macros,
    the output argument c cannot be the same as any input argument.

*/

#if 0
typedef struct
{
    double component [2] ;      /* real and imaginary parts */

} Double_Complex ;

typedef Double_Complex Unit ;
/*#define Entry Double_Complex*/
#define Real component [0]
#define Imag component [1]

/* for flop counts */
#define MULTSUB_FLOPS   8.      /* c -= a*b */
#define DIV_FLOPS       9.      /* c = a/b */
#define ABS_FLOPS       6.      /* c = abs (a), count sqrt as one flop */
#define ASSEMBLE_FLOPS  2.      /* c += a */
#define DECREMENT_FLOPS 2.      /* c -= a */
#define MULT_FLOPS      6.      /* c = a*b */
#define SCALE_FLOPS     2.      /* c = a/s or c = a*s */

/* -------------------------------------------------------------------------- */

/* Return TRUE if a complex number is in split form, FALSE if in packed form */
#define SPLIT(sz) ((sz) != (double *) NULL)

/* c = (s1) + (s2)*i, if s2 is null, then X is in "packed" format (compatible
 * with Entry and ANSI C99 double _Complex type).  */
/*#define ASSIGN(c,s1,s2,p,split)       \
{ \
    if (split) \
    { \
        (c).Real = (s1)[p] ; \
        (c).Imag = (s2)[p] ; \
    }  \
    else \
    { \
        (c) = ((Entry *)(s1))[p] ; \
    }  \
}*/

/* -------------------------------------------------------------------------- */
#endif

/* -------------------------------------------------------------------------- */

/* print an entry (avoid printing "-0" for negative zero).  */
#ifdef NPRINT
#define PRINT_ENTRY(a)
#else
#define PRINT_ENTRY(a) \
{ \
    if (SCALAR_IS_NONZERO ((a).Real)) \
    { \
        PRINTF ((" (%g", (a).Real)) ; \
    } \
    else \
    { \
        PRINTF ((" (0")) ; \
    } \
    if (SCALAR_IS_LTZERO ((a).Imag)) \
    { \
        PRINTF ((" - %gi)", -(a).Imag)) ; \
    } \
    else if (SCALAR_IS_ZERO ((a).Imag)) \
    { \
        PRINTF ((" + 0i)")) ; \
    } \
    else \
    { \
        PRINTF ((" + %gi)", (a).Imag)) ; \
    } \
}
#endif

/* -------------------------------------------------------------------------- */

#endif  /* #ifndef COMPLEX */

#endif
