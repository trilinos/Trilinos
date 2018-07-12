#ifndef SHYLUBASKER_SCALAR_TRAITS_HPP
#define SHYLUBASKER_SCALAR_TRAITS_HPP

#ifndef MY_SCALAR_ABS
#define MY_SCALAR_ABS(x) (((x) < 0) ? -(x) : (x))
#endif

namespace BaskerNS
{

  template <class T>
  struct Basker_ScalarTraits
  {
    typedef T magnitudeType;
    static inline double reciprocal(double c ){return 0;}
    static inline double divide(double a, double b){return 0;}
    static inline magnitudeType approxABS(double a){return 0;}
    static inline magnitudeType abs(double a){return 0;}
    static inline bool gt (double a, double b){return 0;}  
  };
  
  //double
  template <>
  struct Basker_ScalarTraits<double>
  {
    typedef double magnitudeType;
    static inline double reciprocal(double c){return 1.0/c;}
    static inline double divide(double a, double b){return a/b;}
    static inline magnitudeType approxABS(double a)
    { return (MY_SCALAR_ABS(a));}
    static inline magnitudeType abs(double a)
    { return (MY_SCALAR_ABS(a));}
    static inline bool gt (double a, double b){return (a>b);}
  };
  
  //float
  template <>
  struct Basker_ScalarTraits<float>
  {
    typedef float magnitudeType;
    static inline float reciprocal(float c){return 1.0/c;}
    static inline float divide(float a, float b){return a/b;}
    static inline magnitudeType approxABS(float a)
    { return (MY_SCALAR_ABS(a));}
    static inline magnitudeType abs(float a)
    { return (MY_SCALAR_ABS(a));}
    static inline bool gt(float a, float b){return (a>b);}
  };
  
  //complex
  #ifdef HAVE_TEUCHOS_COMPLEX
  template <class T>
  struct Basker_ScalarTraits< std::complex<T> >
  {
    typedef std::complex<T> ComplexT ;
    typedef typename Basker_ScalarTraits<T>::magnitudeType magnitudeType ;
    
    static inline ComplexT reciprocal (ComplexT c)
    {
      T r, den, cr, ci ;
      ComplexT ret ;
      cr = (Teuchos::ScalarTraits<ComplexT>::real(c)) ;
      ci = (Teuchos::ScalarTraits<ComplexT>::imag(c)) ;
      if (MY_SCALAR_ABS (cr) >= MY_SCALAR_ABS (ci))
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
      if (MY_SCALAR_ABS (br) >= MY_SCALAR_ABS (bi))
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
    
    static inline magnitudeType approxABS (ComplexT a)
    {
      return ( MY_SCALAR_ABS (Teuchos::ScalarTraits<ComplexT>::real(a)) +
               MY_SCALAR_ABS (Teuchos::ScalarTraits<ComplexT>::imag(a)) ) ;
    }
    
    static inline magnitudeType abs (ComplexT a)
    {
      T r, ar, ai ;
      magnitudeType s;
      
      ar = MY_SCALAR_ABS (Teuchos::ScalarTraits<ComplexT>::real(a)) ;
      ai = MY_SCALAR_ABS (Teuchos::ScalarTraits<ComplexT>::imag(a)) ;
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

    static inline bool gt(ComplexT a, ComplexT b)
    {
      return( (Teuchos::ScalarTraits<ComplexT>::real(a)+Teuchos::ScalarTraits<ComplexT>::imag(a)) > (Teuchos::ScalarTraits<ComplexT>::real(b) + Teuchos::ScalarTraits<ComplexT>::imag(b)));
    }
    
  };
  
  #else  //C++ complexx
#include <complex>

  template <class T>
  struct Basker_ScalarTraits< std::complex<T> >
  {
    typedef std::complex<T> ComplexT ;
    typedef typename Basker_ScalarTraits<T>::magnitudeType magnitudeType ;

    static inline ComplexT reciprocal (ComplexT c)
    {
      T r, den, cr, ci ;
      ComplexT ret ;
      cr = (std::real(c)) ;
      ci = (std::imag(c)) ;
      if (MY_SCALAR_ABS (cr) >= MY_SCALAR_ABS (ci))
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

      br = (std::real(b)) ;
      bi = (std::imag(b)) ;
      ar = (std::real(a)) ;
      ai = (std::imag(a)) ;
      if (MY_SCALAR_ABS (br) >= MY_SCALAR_ABS (bi))
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

    static inline magnitudeType approxABS (ComplexT a)
    {
      return ( MY_SCALAR_ABS (std::real(a)) +
          MY_SCALAR_ABS (std::imag(a)) ) ;
    }

    static inline magnitudeType abs (ComplexT a)
    {
      T r, ar, ai ;
      magnitudeType s;

      ar = MY_SCALAR_ABS (std::real(a)) ;
      ai = MY_SCALAR_ABS (std::imag(a)) ;
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

    static inline bool gt(ComplexT a, ComplexT b)
    {
      return ((std::real(a)+std::imag(a)) > (std::real(b) + std::imag(b)));
    }

  };
  #endif // HAVE _TEUCHOS_COMPLEX

}//end BaskerNS
#endif // BASKER_SCALAR_TRAIT_HPP
