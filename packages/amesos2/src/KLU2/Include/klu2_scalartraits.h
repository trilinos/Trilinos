// @HEADER
// *****************************************************************************
//                   KLU2: A Direct Linear Solver package
//
// Copyright 2011 NTESS and the KLU2 contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef KLU2_SCALARTRAITS_H
#define KLU2_SCALARTRAITS_H

template <typename T>
struct KLU_ScalarTraits
{
    typedef T magnitudeType ;
    static inline double reciprocal (double c) {return 0;}
    static inline double divide (double a, double b) {return 0.0;}
    static inline double divideConjugate (double a, double b) {return 0.0;}
    static inline magnitudeType approxABS (double a)
    {
    }
    static inline magnitudeType abs (double a)
    {
    }
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

template <>
struct KLU_ScalarTraits<float>
{
    typedef float magnitudeType ;
    static inline float reciprocal (float c) { return 1.0/c ; }
    static inline float divide (float a, float b) { return a/b ; }
    static inline float divideConjugate (float a, float b) { return a/b ; }
    static inline magnitudeType approxABS (float a)
    {
        return (SCALAR_ABS (a));
    }
    static inline magnitudeType abs (float a)
    {
        return (SCALAR_ABS (a));
    }
};

// mfh 13 Sep 2012: The Teuchos::ScalarTraits<std::complex<T> >
// specialization doesn't exist unless Teuchos was built with complex
// arithmetic support.  To enable complex arithmetic support in
// Teuchos, set the CMake Boolean option Teuchos_ENABLE_COMPLEX to ON
// at configure time.
#ifdef HAVE_TEUCHOS_COMPLEX

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

#endif // HAVE_TEUCHOS_COMPLEX

#endif
