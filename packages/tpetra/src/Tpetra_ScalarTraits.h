// tpetra/src/Tpetra_ScalarTraits.h
// 16-May-2002 - Changed to use Epetra_LAPACK instead of Petra_LAPACK
// 16-May-2002 - Switched names from TPetra to Tpetra
// 28-May-2002 - Heroux fixes things.

#ifndef _TPETRA_SCALARTRAITS_H_
#define _TPETRA_SCALARTRAITS_H_

#include <cmath> //for the fabs function...
#include <complex>
#include "Tpetra_LAPACK.h"

// Modified from the original ESI "ESI_scalarTraits.h" file. 
namespace Tpetra {
  /** The Tpetra ScalarTraits file.

      For the most general type 'class T', we define aborting functions, 
      which should restrict implementations from using traits other than the
      specializations defined below.
  */

  template<class T>
    struct ScalarTraits {

      typedef T magnitudeType;

      static inline bool haveMachineParameters() {return(false);}; // Allows testing to see if scalar traits machine parameters defined 
      static inline T basicErrorMessage() {
	cout << "Tpetra::ScalarTraits: Machine parameters are undefined for this scalar type." << endl; 
	std::abort();
	T a;
	return(a);
      };
      static inline magnitudeType eps()   {return(basicErrorMessage());};
      static inline magnitudeType sfmin() {return(basicErrorMessage());};
      static inline magnitudeType base()  {return(basicErrorMessage());};
      static inline magnitudeType prec()  {return(basicErrorMessage());};
      static inline magnitudeType t()     {return(basicErrorMessage());};
      static inline magnitudeType rnd()   {return(basicErrorMessage());};
      static inline magnitudeType emin()  {return(basicErrorMessage());};
      static inline magnitudeType rmin()  {return(basicErrorMessage());};
      static inline magnitudeType emax()  {return(basicErrorMessage());};
      static inline magnitudeType rmax()  {return(basicErrorMessage());};
      static inline magnitudeType magnitude(T a) { 
	cout << "Tpetra::ScalarTraits: unsupported scalar type." << endl; 
	std::abort(); 
	return(a);
      };

      static inline T zero() {
	cout << "Tpetra::ScalarTraits: unsupported scalar type." << endl; 
	std::abort();
	T a;
	return(a);
      };
	
      static inline T one() {
	cout << "Tpetra::ScalarTraits: unsupported scalar type." << endl; 
	std::abort();
	T a;
	return(a);
      };
	
      static inline T random() {
	cout << "Tpetra::ScalarTraits: unsupported scalar type." << endl; 
	std::abort();
	T a;
	return(a);
      };
      static inline const char * name() {
	cout << "Tpetra::ScalarTraits: unsupported scalar type." << endl; 
	std::abort(); 
	return(0); 
      };
    };

  template<>
    struct ScalarTraits<int> {

      typedef long long magnitudeType;
	
      static inline magnitudeType magnitude(int a) { 
	return(abs(a)); 
      };
      static inline bool haveMachineParameters() {return(false);}; // Allows testing to see if scalar traits machine parameters defined 
      static inline magnitudeType basicErrorMessage() {
	cout << "Tpetra::ScalarTraits: Machine parameters are undefined for this scalar type." << endl; 
	std::abort();
	magnitudeType a;
	return(a);
      };
      static inline magnitudeType eps()   {return(basicErrorMessage());};
      static inline magnitudeType sfmin() {return(basicErrorMessage());};
      static inline magnitudeType base()  {return(basicErrorMessage());};
      static inline magnitudeType prec()  {return(basicErrorMessage());};
      static inline magnitudeType t()     {return(basicErrorMessage());};
      static inline magnitudeType rnd()   {return(basicErrorMessage());};
      static inline magnitudeType emin()  {return(basicErrorMessage());};
      static inline magnitudeType rmin()  {return(basicErrorMessage());};
      static inline magnitudeType emax()  {return(basicErrorMessage());};
      static inline magnitudeType rmax()  {return(basicErrorMessage());};

      static inline int zero()  {return(0);};
      static inline int one()   {return(1);};
      

      static inline int random() {return(rand()/RAND_MAX);};

      static inline const char * name() { 
	return("int"); 
      };
    };


  template<>
    struct ScalarTraits<float> {

      static inline bool haveMachineParameters() {return(true);}; // Allows testing to see if scalar traits machine parameters defined 
      typedef float magnitudeType;
	
      static inline magnitudeType magnitude(float a) { 
	return(fabs(a)); 
      };

      static inline float zero()  {return(0.0);};
      static inline float one()   {return(1.0);};
      static inline float eps()   {LAPACK lp; return(lp.SLAMCH('E'));};
      static inline float sfmin() {LAPACK lp; return(lp.SLAMCH('S'));};
      static inline float base()  {LAPACK lp; return(lp.SLAMCH('B'));};
      static inline float prec()  {LAPACK lp; return(lp.SLAMCH('P'));};
      static inline float t()     {LAPACK lp; return(lp.SLAMCH('N'));};
      static inline float rnd()   {LAPACK lp; return(lp.SLAMCH('R'));};
      static inline float emin()  {LAPACK lp; return(lp.SLAMCH('M'));};
      static inline float rmin()  {LAPACK lp; return(lp.SLAMCH('U'));};
      static inline float emax()  {LAPACK lp; return(lp.SLAMCH('L'));};
      static inline float rmax()  {LAPACK lp; return(lp.SLAMCH('O'));};
      

      static inline float random() {
	float rnd = (float) std::rand()/RAND_MAX;
	return( (float)(-1.0 + 2.0*rnd) );
      };

      static inline const char * name() { 
	return("float"); 
      };
    };


  template<>
    struct ScalarTraits<double> {

      static inline bool haveMachineParameters() {return(true);}; // Allows testing to see if scalar traits machine parameters defined 
      typedef double magnitudeType;
	
      static inline magnitudeType magnitude(double a) {return(fabs(a));};

      static inline double zero()  {return(0.0);};
      static inline double one()   {return(1.0);};
      static inline double eps()   {LAPACK lp; return(lp.DLAMCH('E'));};
      static inline double sfmin() {LAPACK lp; return(lp.DLAMCH('S'));};
      static inline double base()  {LAPACK lp; return(lp.DLAMCH('B'));};
      static inline double prec()  {LAPACK lp; return(lp.DLAMCH('P'));};
      static inline double t()     {LAPACK lp; return(lp.DLAMCH('N'));};
      static inline double rnd()   {LAPACK lp; return(lp.DLAMCH('R'));};
      static inline double emin()  {LAPACK lp; return(lp.DLAMCH('M'));};
      static inline double rmin()  {LAPACK lp; return(lp.DLAMCH('U'));};
      static inline double emax()  {LAPACK lp; return(lp.DLAMCH('L'));};
      static inline double rmax()  {LAPACK lp; return(lp.DLAMCH('O'));};
      

      static inline double random() {
	double rnd = (double) std::rand()/RAND_MAX;
	return( (double)(-1.0 + 2.0*rnd) );
      };

      static inline const char * name() { 
	return("double");
      };
    };


//  If we're using a Sun compiler, version earlier than 5.0,
//  then complex isn't available.

#if defined(__SUNPRO_CC) && __SUNPRO_CC < 0x500
#define NO_COMPLEX
#endif


//  If we're using the tflops Portland Group compiler, then complex isn't
//  available. (As of July 21, 2000. abw)

#if defined(__PGI) && defined(__i386)
#define NO_COMPLEX
#endif

#ifndef NO_COMPLEX

#include <complex>

  template<> struct ScalarTraits< std::complex<float> > {

      static inline bool haveMachineParameters() {return(true);}; // Allows testing to see if scalar traits machine parameters defined 
      typedef float magnitudeType;

      static magnitudeType magnitude(std::complex<float> a) {
	return(std::abs(a));
     };

      static inline std::complex<float> zero()  {return(std::complex<float>(0.0,0.0));};
      static inline std::complex<float> one()   {return(std::complex<float>(1.0,0.0));};
      static inline float eps()   {LAPACK lp; return(lp.SLAMCH('E'));};
      static inline float sfmin() {LAPACK lp; return(lp.SLAMCH('S'));};
      static inline float base()  {LAPACK lp; return(lp.SLAMCH('B'));};
      static inline float prec()  {LAPACK lp; return(lp.SLAMCH('P'));};
      static inline float t()     {LAPACK lp; return(lp.SLAMCH('N'));};
      static inline float rnd()   {LAPACK lp; return(lp.SLAMCH('R'));};
      static inline float emin()  {LAPACK lp; return(lp.SLAMCH('M'));};
      static inline float rmin()  {LAPACK lp; return(lp.SLAMCH('U'));};
      static inline float emax()  {LAPACK lp; return(lp.SLAMCH('L'));};
      static inline float rmax()  {LAPACK lp; return(lp.SLAMCH('O'));};

     static inline std::complex<float> random() {
	//float rnd1 = (float)std::rand()/RAND_MAX;
	//float rnd2 = (float)std::rand()/RAND_MAX;
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

      static inline bool haveMachineParameters() {return(true);}; // Allows testing to see if scalar traits machine parameters defined 
      typedef double magnitudeType;

      static magnitudeType magnitude(std::complex<double> a) {
	return(std::abs(a));
      };

      static inline std::complex<double> zero()  {return(std::complex<double>(0.0,0.0));};
      static inline std::complex<double> one()   {return(std::complex<double>(1.0,0.0));};
      static inline double eps()   {LAPACK lp; return(lp.DLAMCH('E'));};
      static inline double sfmin() {LAPACK lp; return(lp.DLAMCH('S'));};
      static inline double base()  {LAPACK lp; return(lp.DLAMCH('B'));};
      static inline double prec()  {LAPACK lp; return(lp.DLAMCH('P'));};
      static inline double t()     {LAPACK lp; return(lp.DLAMCH('N'));};
      static inline double rnd()   {LAPACK lp; return(lp.DLAMCH('R'));};
      static inline double emin()  {LAPACK lp; return(lp.DLAMCH('M'));};
      static inline double rmin()  {LAPACK lp; return(lp.DLAMCH('U'));};
      static inline double emax()  {LAPACK lp; return(lp.DLAMCH('L'));};
      static inline double rmax()  {LAPACK lp; return(lp.DLAMCH('O'));};

      static inline std::complex<double> random() {
	//double rnd1 = (double)std::rand()/RAND_MAX;
	//double rnd2 = (double)std::rand()/RAND_MAX;
	double rnd1 = ScalarTraits<magnitudeType>::random();
	double rnd2 = ScalarTraits<magnitudeType>::random();
	return( std::complex<double>(rnd1, rnd2) );
      };

      static inline const char * name() {
	return("std::complex<double>");
      };
    };

#endif  // NO_COMPLEX


//  If we're using a Sun compiler, version earlier than 5.0,
//  then 'typename' isn't available.

#if defined(__SUNPRO_CC) && __SUNPRO_CC < 0x500
#define TYPENAME
#else
#define TYPENAME typename
#endif

}     // Tpetra namespace
#endif // _TPETRA_SCALARTRAITS_H_

