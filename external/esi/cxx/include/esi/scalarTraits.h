#ifndef __ESI_scalarTraits_h
#define __ESI_scalarTraits_h

#include <math.h> //for the fabs function...

namespace esi {

/** The ESI scalarTraits file.

    For the most general type 'class T', we define aborting functions, 
    which should restrict implementations from using traits other than the
    specializations defined below.
*/

template<class T>
struct scalarTraits {

  typedef T magnitude_type;

  static inline magnitude_type magnitude(T a) { 
    cout << "esi::scalarTraits: unsupported scalar type." << endl; 
    abort(); 
    return(a);
  };

  static inline T random() {
    cout << "esi::scalarTraits: unsupported scalar type." << endl; 
    abort(); 
    return(a);
  };
  
  static inline const char * name() {
    cout << "esi::scalarTraits: unsupported scalar type." << endl; 
    abort(); 
    return(0); 
  };
};

template<>
struct scalarTraits<float> {

  typedef float magnitude_type;
  
  static inline magnitude_type magnitude(float a) { 
    return(fabs(a)); 
  };

  static inline float random() {
    float rnd = (float)rand()/RAND_MAX;
    return( (float)(-1.0 + 2.0*rnd) );
  };

  static inline const char * name() { 
    return("float"); 
  };
};


template<>
struct scalarTraits<double> {

  typedef double magnitude_type;

  static magnitude_type magnitude(double a) {
    return(fabs(a));
  };

  static inline double random() {
    double rnd = (double)rand()/RAND_MAX;
    return( (double)(-1.0 + 2.0*rnd) );
  };

  static inline const char * name() {
    return("double");
  };
};

#ifndef ESI_NO_COMPLEX
}; //close the esi namespace
#include <complex>
namespace esi {

template<>
struct scalarTraits< std::complex<float> > {

  typedef float magnitude_type;

  static magnitude_type magnitude(std::complex<float> a) {
    return(std::abs(a));
  };

  static inline std::complex<float> random() {
    float rnd1 = (float)rand()/RAND_MAX;
    float rnd2 = (float)rand()/RAND_MAX;
    return( std::complex<float>(-1.0+2.0*rnd1, -1.0+2.0*rnd2) );
  };

  static inline const char * name() {
    return("std::complex<float>");
  };
};

template<>
struct scalarTraits< std::complex<double> > {

  typedef double magnitude_type;

  static magnitude_type magnitude(std::complex<double> a) {
    return(std::abs(a));
  };
  
  static inline std::complex<double> random() {
    double rnd1 = (double)rand()/RAND_MAX;
    double rnd2 = (double)rand()/RAND_MAX;
    return( std::complex<double>(-1.0+2.0*rnd1, -1.0+2.0*rnd2) );
  };

  static inline const char * name() {
    return("std::complex<double>");
  };
};

#endif  /* ESI_NO_COMPLEX */


#ifdef ESI_NO_TYPENAME_KEYWORD
#define TYPENAME
#else
#define TYPENAME typename
#endif

};     // esi namespace
#endif

