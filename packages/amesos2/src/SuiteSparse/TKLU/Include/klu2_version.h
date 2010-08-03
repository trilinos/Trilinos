#ifndef _TKLU_VERSION_H
#define _TKLU_VERSION_H

#ifdef DLONG
#define Int_id UF_long_id
#define Int_MAX UF_long_max
#else
#define Int_id "%d"
#define Int_MAX INT_MAX
#endif

#define NPRINT  

#define BYTES(type,n) (sizeof (type) * (n))
#define CEILING(b,u)  (((b)+(u)-1) / (u))
#define UNITS(type,n) (CEILING (BYTES (type,n), sizeof (Unit)))
#define DUNITS(type,n) (ceil (BYTES (type, (double) n) / sizeof (Unit)))

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

#ifdef DLONG

#define BTF_order amesos_btf_l_order
#define BTF_strongcomp amesos_btf_l_strongcomp

#define AMD_order amesos_amd_l_order
#define COLAMD amesos_colamd_l
#define COLAMD_recommended amesos_colamd_l_recommended

#else

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

#define BTF_order amesos_btf_order
#define BTF_strongcomp amesos_btf_strongcomp

#define AMD_order amesos_amd_order
#define COLAMD amesos_colamd
#define COLAMD_recommended amesos_colamd_recommended

#endif

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

template <typename T>
struct KLU_ScalarTraits
{
    typedef T magnitudeType ;
    static inline magnitudeType real(T c) {}
};

template <>
struct KLU_ScalarTraits<double>
{
    typedef double magnitudeType ;
    static inline double reciprocal (double c) { return 1.0/c ; }
    static inline double divide (double a, double b) { return a/b ; }
    static inline double divideConjugate (double a, double b) { return a/b ; }
    static inline magnitudeType approxABS (double a)
    {
        return (SCALAR_ABS (a));
    }
    static inline magnitudeType abs (double a)
    {
        return (SCALAR_ABS (a));
    }
};

template <typename T>
struct KLU_ScalarTraits<
std::complex<T>
>
{
    typedef std::complex<T> ComplexT ;
    typedef typename KLU_ScalarTraits<T>::magnitudeType magnitudeType ;

    static inline ComplexT reciprocal (ComplexT c) 
    {
        T r, den, cr, ci ;
        ComplexT ret ;
        cr = (Teuchos::ScalarTraits<ComplexT>::real(c)) ;
        ci = (Teuchos::ScalarTraits<ComplexT>::imag(c)) ;
        if (SCALAR_ABS (cr) >= SCALAR_ABS (ci))
        {
            r = ci / cr ;
            den = cr + r * ci ;
            ret = std::complex<T>(1.0 / den, -r / den) ;
        }
        else
        {
            r = cr / ci ;
            den = r * cr + ci ;
            ret = std::complex<T>(r / den, -1.0 / den) ;
        }
        return ret;
    }
    
    static inline ComplexT divide (ComplexT a, ComplexT b)
    {
        T r, den, ar, ai, br, bi ;
        ComplexT ret;

        br = (Teuchos::ScalarTraits<ComplexT>::real(b)) ;
        bi = (Teuchos::ScalarTraits<ComplexT>::imag(b)) ;
        ar = (Teuchos::ScalarTraits<ComplexT>::real(a)) ;
        ai = (Teuchos::ScalarTraits<ComplexT>::imag(a)) ;
        if (SCALAR_ABS (br) >= SCALAR_ABS (bi))
        {
            r = bi / br ;
            den = br + r * bi ;
            ret = std::complex<T>((ar + ai * r) / den, (ai - ar * r) / den) ;
        }
        else
        {
            r = br / bi ;
            den = r * br + bi ;
            ret = std::complex<T>((ar * r + ai) / den, (ai * r - ar) / den) ;
        }
        return ret;
    }
    
    static inline ComplexT divideConjugate (ComplexT a, ComplexT b)
    {
        T r, den, ar, ai, br, bi ;
        ComplexT ret;

        br = (Teuchos::ScalarTraits<ComplexT>::real(b)) ;
        bi = (Teuchos::ScalarTraits<ComplexT>::imag(b)) ;
        ar = (Teuchos::ScalarTraits<ComplexT>::real(a)) ;
        ai = (Teuchos::ScalarTraits<ComplexT>::imag(a)) ;
        if (SCALAR_ABS (br) >= SCALAR_ABS (bi))
        {
            r = (-bi) / br ;
            den = br - r * bi ;
            ret = std::complex<T>((ar + ai * r) / den, (ai - ar * r) / den) ;
        }
        else
        {
            r = br / (-bi) ;
            den =  r * br - bi;
            ret = std::complex<T>((ar * r + ai) / den, (ai * r - ar) / den) ;
        }
        return ret;
    }
    
    static inline magnitudeType approxABS (ComplexT a)
    {
        return ( SCALAR_ABS (Teuchos::ScalarTraits<ComplexT>::real(a)) + 
                    SCALAR_ABS (Teuchos::ScalarTraits<ComplexT>::imag(a)) ) ;
    }
    
    static inline magnitudeType abs (ComplexT a)
    {
        T r, ar, ai ;
        magnitudeType s;

        ar = SCALAR_ABS (Teuchos::ScalarTraits<ComplexT>::real(a)) ;
        ai = SCALAR_ABS (Teuchos::ScalarTraits<ComplexT>::imag(a)) ;
        if (ar >= ai)
        {
            if (ar + ai == ar)
            {
                (s) = ar ;
            }
            else
            {
                r = ai / ar ;
                (s) = ar * sqrt (1.0 + r*r) ;
            }
        }
        else
        {
            if (ai + ar == ai)
            {
                (s) = ai ;
            }
            else
            {
                r = ar / ai ;
                (s) = ai * sqrt (1.0 + r*r) ;
            }
        }
        return s;
    }
};

#ifndef COMPLEX

typedef double Unit ;
/*#define Entry double*/

#define SPLIT(s)                    (1)
/*#define REAL(c)                     (c)
#define IMAG(c)                     (0.)
#define ASSIGN(c,s1,s2,p,split)     { (c) = (s1)[p] ; }
#define MULT_CONJ(c,a,b)            { (c) = (a) * (b) ; }
#define MULT_SUB_CONJ(c,a,b)        { (c) -= (a) * (b) ; }*/
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
#define ABS(s,a)                    { (s) =  KLU_ScalarTraits<Entry>::abs(a) ; }
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

/* real part of c */
#define REAL(c) ((c).Real)

/* -------------------------------------------------------------------------- */

/* imag part of c */
#define IMAG(c) ((c).Imag)

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
#define CONJ(a, x) \
{ \
    a.Real = x.Real ; \
    a.Imag = -x.Imag ; \
}

/* c = 0 */
#define CLEAR(c) \
{ \
    (c).Real = 0. ; \
    (c).Imag = 0. ; \
}

/* -------------------------------------------------------------------------- */

/* *p++ = 0 */
#define CLEAR_AND_INCREMENT(p) \
{ \
    p->Real = 0. ; \
    p->Imag = 0. ; \
    p++ ; \
}

/* -------------------------------------------------------------------------- */

/* True if a == 0 */
#define IS_ZERO(a) \
    (SCALAR_IS_ZERO ((a).Real) && SCALAR_IS_ZERO ((a).Imag))

/* -------------------------------------------------------------------------- */

/* True if a is NaN */
#define IS_NAN(a) \
    (SCALAR_IS_NAN ((a).Real) || SCALAR_IS_NAN ((a).Imag))

/* -------------------------------------------------------------------------- */

/* True if a != 0 */
#define IS_NONZERO(a) \
    (SCALAR_IS_NONZERO ((a).Real) || SCALAR_IS_NONZERO ((a).Imag))

/* -------------------------------------------------------------------------- */

/* a = c/s */
#define SCALE_DIV_ASSIGN(a,c,s) \
{ \
    a.Real = c.Real / s ; \
    a.Imag = c.Imag / s ; \
}

/* c /= s */
#define SCALE_DIV(c,s) \
{ \
    (c).Real /= (s) ; \
    (c).Imag /= (s) ; \
}

/* -------------------------------------------------------------------------- */

/* c *= s */
#define SCALE(c,s) \
{ \
    (c).Real *= (s) ; \
    (c).Imag *= (s) ; \
}

/* -------------------------------------------------------------------------- */

/* c += a */
#define ASSEMBLE(c,a) \
{ \
    (c).Real += (a).Real ; \
    (c).Imag += (a).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c += *p++ */
#define ASSEMBLE_AND_INCREMENT(c,p) \
{ \
    (c).Real += p->Real ; \
    (c).Imag += p->Imag ; \
    p++ ; \
}

/* -------------------------------------------------------------------------- */

/* c -= a */
#define DECREMENT(c,a) \
{ \
    (c).Real -= (a).Real ; \
    (c).Imag -= (a).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c = a*b, assert because c cannot be the same as a or b */
#define MULT(c,a,b) \
{ \
    ASSERT (&(c) != &(a) && &(c) != &(b)) ; \
    (c).Real = (a).Real * (b).Real - (a).Imag * (b).Imag ; \
    (c).Imag = (a).Imag * (b).Real + (a).Real * (b).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c = a*conjugate(b), assert because c cannot be the same as a or b */
#define MULT_CONJ(c,a,b) \
{ \
    ASSERT (&(c) != &(a) && &(c) != &(b)) ; \
    (c).Real = (a).Real * (b).Real + (a).Imag * (b).Imag ; \
    (c).Imag = (a).Imag * (b).Real - (a).Real * (b).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c -= a*b, assert because c cannot be the same as a or b */
#define MULT_SUB(c,a,b) \
{ \
    ASSERT (&(c) != &(a) && &(c) != &(b)) ; \
    (c).Real -= (a).Real * (b).Real - (a).Imag * (b).Imag ; \
    (c).Imag -= (a).Imag * (b).Real + (a).Real * (b).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c -= a*conjugate(b), assert because c cannot be the same as a or b */
#define MULT_SUB_CONJ(c,a,b) \
{ \
    ASSERT (&(c) != &(a) && &(c) != &(b)) ; \
    (c).Real -= (a).Real * (b).Real + (a).Imag * (b).Imag ; \
    (c).Imag -= (a).Imag * (b).Real - (a).Real * (b).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c = a/b, be careful to avoid underflow and overflow */
#ifdef MATHWORKS
#define DIV(c,a,b) \
{ \
    (void) utDivideComplex ((a).Real, (a).Imag, (b).Real, (b).Imag, \
        &((c).Real), &((c).Imag)) ; \
}
#else
/* This uses ACM Algo 116, by R. L. Smith, 1962. */
/* c can be the same variable as a or b. */
/* Ignore NaN case for double relop br>=bi. */
#define DIV(c,a,b) \
{ \
    double r, den, ar, ai, br, bi ; \
    br = (b).Real ; \
    bi = (b).Imag ; \
    ar = (a).Real ; \
    ai = (a).Imag ; \
    if (SCALAR_ABS (br) >= SCALAR_ABS (bi)) \
    { \
        r = bi / br ; \
        den = br + r * bi ; \
        (c).Real = (ar + ai * r) / den ; \
        (c).Imag = (ai - ar * r) / den ; \
    } \
    else \
    { \
        r = br / bi ; \
        den = r * br + bi ; \
        (c).Real = (ar * r + ai) / den ; \
        (c).Imag = (ai * r - ar) / den ; \
    } \
}
#endif

/* -------------------------------------------------------------------------- */

/* c = 1/c, be careful to avoid underflow and overflow */
/* Not used if MATHWORKS is defined. */
/* This uses ACM Algo 116, by R. L. Smith, 1962. */
/* Ignore NaN case for double relop cr>=ci. */
#define RECIPROCAL(c) \
{ \
    double r, den, cr, ci ; \
    cr = (c).Real ; \
    ci = (c).Imag ; \
    if (SCALAR_ABS (cr) >= SCALAR_ABS (ci)) \
    { \
        r = ci / cr ; \
        den = cr + r * ci ; \
        (c).Real = 1.0 / den ; \
        (c).Imag = - r / den ; \
    } \
    else \
    { \
        r = cr / ci ; \
        den = r * cr + ci ; \
        (c).Real = r / den ; \
        (c).Imag = - 1.0 / den ; \
    } \
}


/* -------------------------------------------------------------------------- */

/* c = a/conjugate(b), be careful to avoid underflow and overflow */
#ifdef MATHWORKS
#define DIV_CONJ(c,a,b) \
{ \
    (void) utDivideComplex ((a).Real, (a).Imag, (b).Real, (-(b).Imag), \
        &((c).Real), &((c).Imag)) ; \
}
#else
/* This uses ACM Algo 116, by R. L. Smith, 1962. */
/* c can be the same variable as a or b. */
/* Ignore NaN case for double relop br>=bi. */
#define DIV_CONJ(c,a,b) \
{ \
    double r, den, ar, ai, br, bi ; \
    br = (b).Real ; \
    bi = (b).Imag ; \
    ar = (a).Real ; \
    ai = (a).Imag ; \
    if (SCALAR_ABS (br) >= SCALAR_ABS (bi)) \
    { \
        r = (-bi) / br ; \
        den = br - r * bi ; \
        (c).Real = (ar + ai * r) / den ; \
        (c).Imag = (ai - ar * r) / den ; \
    } \
    else \
    { \
        r = br / (-bi) ; \
        den =  r * br - bi; \
        (c).Real = (ar * r + ai) / den ; \
        (c).Imag = (ai * r - ar) / den ; \
    } \
}
#endif

/* -------------------------------------------------------------------------- */

/* approximate absolute value, s = |r|+|i| */
#define APPROX_ABS(s,a) \
{ \
    (s) = SCALAR_ABS ((a).Real) + SCALAR_ABS ((a).Imag) ; \
}

/* -------------------------------------------------------------------------- */

/* exact absolute value, s = sqrt (a.real^2 + amag^2) */
#ifdef MATHWORKS
#define ABS(s,a) \
{ \
    (s) = utFdlibm_hypot ((a).Real, (a).Imag) ; \
}
#else
/* Ignore NaN case for the double relops ar>=ai and ar+ai==ar. */
#define ABS(s,a) \
{ \
    double r, ar, ai ; \
    ar = SCALAR_ABS ((a).Real) ; \
    ai = SCALAR_ABS ((a).Imag) ; \
    if (ar >= ai) \
    { \
        if (ar + ai == ar) \
        { \
            (s) = ar ; \
        } \
        else \
        { \
            r = ai / ar ; \
            (s) = ar * sqrt (1.0 + r*r) ; \
        } \
    } \
    else \
    { \
        if (ai + ar == ai) \
        { \
            (s) = ai ; \
        } \
        else \
        { \
            r = ar / ai ; \
            (s) = ai * sqrt (1.0 + r*r) ; \
        } \
    } \
}
#endif
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
