#ifndef _TPETRA_SCALARTRAITS_H_
#define _TPETRA_SCALARTRAITS_H_

#include <math.h> //for the fabs function...

// Modified from the original ESI "ESI_scalarTraits.h" file. 
namespace TPetra {
  /** The TPetra ScalarTraits file.

      For the most general type 'class T', we define aborting functions, 
      which should restrict implementations from using traits other than the
      specializations defined below.
  */

  template<class T>
    struct ScalarTraits {

      typedef T magnitudeType;

      static inline magnitudeType magnitude(T a) { 
	cout << "TPetra::ScalarTraits: unsupported scalar type." << endl; 
	abort(); 
	return(a);
      };

      static inline T zero() {
	cout << "TPetra::ScalarTraits: unsupported scalar type." << endl; 
	abort();
	T a;
	return(a);
      };
	
      static inline T one() {
	cout << "TPetra::ScalarTraits: unsupported scalar type." << endl; 
	abort();
	T a;
	return(a);
      };
	
      static inline T random() {
	cout << "TPetra::ScalarTraits: unsupported scalar type." << endl; 
	abort();
	T a;
	return(a);
      };
      static inline const char * name() {
	cout << "TPetra::ScalarTraits: unsupported scalar type." << endl; 
	abort(); 
	return(0); 
      };
    };

  template<>
    struct ScalarTraits<float> {

      typedef float magnitudeType;
	
      static inline magnitudeType magnitude(float a) { 
	return(fabs(a)); 
      };

      static inline float zero()  {return(0.0);};
      static inline float one()   {return(1.0);};
/*
      static inline float eps()   {Petra_LAPACK lp; return(lp.SLAMCH('E'));};
      static inline float sfmin() {Petra_LAPACK lp; return(lp.SLAMCH('S'));};
      static inline float base()  {Petra_LAPACK lp; return(lp.SLAMCH('B'));};
      static inline float prec()  {Petra_LAPACK lp; return(lp.SLAMCH('P'));};
      static inline float t()     {Petra_LAPACK lp; return(lp.SLAMCH('N'));};
      static inline float rnd()   {Petra_LAPACK lp; return(lp.SLAMCH('R'));};
      static inline float emin()  {Petra_LAPACK lp; return(lp.SLAMCH('M'));};
      static inline float rmin()  {Petra_LAPACK lp; return(lp.SLAMCH('U'));};
      static inline float emax()  {Petra_LAPACK lp; return(lp.SLAMCH('L'));};
      static inline float rmax()  {Petra_LAPACK lp; return(lp.SLAMCH('O'));};
*/      

      static inline float random() {
	float rnd = (float)rand()/RAND_MAX;
	return( (float)(-1.0 + 2.0*rnd) );
      };

      static inline const char * name() { 
	return("float"); 
      };
    };


  template<>
    struct ScalarTraits<double> {

      typedef double magnitudeType;
	
      static inline magnitudeType magnitude(double a) {return(fabs(a));};

      static inline double zero()  {return(0.0);};
      static inline double one()   {return(1.0);};
/*
      static inline double eps()   {Petra_LAPACK lp; return(lp.DLAMCH('E'));};
      static inline double sfmin() {Petra_LAPACK lp; return(lp.DLAMCH('S'));};
      static inline double base()  {Petra_LAPACK lp; return(lp.DLAMCH('B'));};
      static inline double prec()  {Petra_LAPACK lp; return(lp.DLAMCH('P'));};
      static inline double t()     {Petra_LAPACK lp; return(lp.DLAMCH('N'));};
      static inline double rnd()   {Petra_LAPACK lp; return(lp.DLAMCH('R'));};
      static inline double emin()  {Petra_LAPACK lp; return(lp.DLAMCH('M'));};
      static inline double rmin()  {Petra_LAPACK lp; return(lp.DLAMCH('U'));};
      static inline double emax()  {Petra_LAPACK lp; return(lp.DLAMCH('L'));};
      static inline double rmax()  {Petra_LAPACK lp; return(lp.DLAMCH('O'));};
*/     

      static inline double random() {
	double rnd = (double)rand()/RAND_MAX;
	return( (double)(-1.0 + 2.0*rnd) );
      };

      static inline const char * name() { 
	return("double"); 
      };
    };

  /*
    If we're using a Sun compiler, version earlier than 5.0,
    then complex isn't available.
  */
#if defined(__SUNPRO_CC) && __SUNPRO_CC < 0x500
#define NO_COMPLEX
#endif

  /*
    If we're using the tflops Portland Group compiler, then complex isn't
    available. (As of July 21, 2000. abw)
  */
#if defined(__PGI) && defined(__i386)
#define NO_COMPLEX
#endif

//#define NO_COMPLEX
#ifndef NO_COMPLEX

#include <complex>

  template<>
    struct ScalarTraits< std::complex<float> > {

      typedef float magnitudeType;
      typedef std::complex<float> scalarType;

      static magnitudeType magnitude(std::complex<float> a) {
	return(std::abs(a));
      };

      static inline std::complex<float> zero()  {return(std::complex<float>(0.0,0.0));};
      static inline std::complex<float> one()   {return(std::complex<float>(1.0,0.0));};

     static inline std::complex<float> random() {
	//float rnd1 = (float)rand()/RAND_MAX;
	//float rnd2 = (float)rand()/RAND_MAX;
	float rnd1 = ScalarTraits<magnitudeType>::random();
	float rnd2 = ScalarTraits<magnitudeType>::random();
	return( std::complex<float>(rnd1, rnd2) );
      };

      static inline const char * name() {
	return("std::complex<float>");
      };
    };

  template<>
    struct ScalarTraits< std::complex<double> > {

      typedef double magnitudeType;

      static magnitudeType magnitude(std::complex<double> a) {
	return(std::abs(a));
      };

      static inline std::complex<double> zero()  {return(std::complex<double>(0.0,0.0));};
      static inline std::complex<double> one()   {return(std::complex<double>(1.0,0.0));};

      static inline std::complex<double> random() {
	//double rnd1 = (double)rand()/RAND_MAX;
	//double rnd2 = (double)rand()/RAND_MAX;
	double rnd1 = ScalarTraits<magnitudeType>::random();
	double rnd2 = ScalarTraits<magnitudeType>::random();
	return( std::complex<double>(rnd1, rnd2) );
      };

      static inline const char * name() {
	return("std::complex<double>");
      };
    };

#endif  /* NO_COMPLEX */


  /*
    If we're using a Sun compiler, version earlier than 5.0,
    then 'typename' isn't available.
  */
#if defined(__SUNPRO_CC) && __SUNPRO_CC < 0x500
#define TYPENAME
#else
#define TYPENAME typename
#endif

}     // TPetra namespace
#endif // _TPETRA_SCALARTRAITS_H_

