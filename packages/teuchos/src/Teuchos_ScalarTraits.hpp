/*Paul
16-May-2002 - Changed to use Epetra_LAPACK instead of Petra_LAPACK
16-May-2002 - Switched names from TPetra to Tpetra
28-May-2002 - Heroux fixes things.
06-August-2002 Switched to images (nothing changed). Cleaned up a bit.
19-Nov-2002 ScalarTraits throws exceptions now instead of calling std::abort()
04-Dec-2002 Moved configs to Tpetra_ConfigDefs.h.
*/

// Kris
// 06.18.03 -- Minor formatting changes
//          -- Changed calls to LAPACK objects to use new <OType, SType> templates
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_SCALARTRAITS_HPP_
#define _TEUCHOS_SCALARTRAITS_HPP_

#include <cmath>
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_ConfigDefs.hpp"
#ifndef NO_COMPLEX
#include <complex>
#define ComplexFloat std::complex<float>
#define ComplexDouble std::complex<double>
#endif

namespace Teuchos {
  /** The Teuchos ScalarTraits file.
      
      For the most general type 'class T', we define aborting functions, 
      which should restrict implementations from using traits other than the
      specializations defined below.
  */
  
  template<class T>
  struct ScalarTraits 
  {
    typedef T magnitudeType;
    
    static inline int undefinedParameters()
    {
#ifndef TEUCHOS_NO_ERROR_REPORTS
      cerr << endl << "Teuchos::ScalarTraits: Machine parameters are undefined for this scalar type." << endl;
#endif
      return(-1);
    }

    static inline int unsupportedType()
    {
#ifndef TEUCHOS_NO_ERROR_REPORTS
      cerr << endl << "Teuchos::ScalarTraits: unsupported scalar type." << endl;
#endif
      return(-2);
    }
    
    static inline bool haveMachineParameters() { return false; };
    static inline magnitudeType eps()   { throw(undefinedParameters()); };
    static inline magnitudeType sfmin() { throw(undefinedParameters()); };
    static inline magnitudeType base()  { throw(undefinedParameters()); };
    static inline magnitudeType prec()  { throw(undefinedParameters()); };
    static inline magnitudeType t()     { throw(undefinedParameters()); };
    static inline magnitudeType rnd()   { throw(undefinedParameters()); };
    static inline magnitudeType emin()  { throw(undefinedParameters()); };
    static inline magnitudeType rmin()  { throw(undefinedParameters()); };
    static inline magnitudeType emax()  { throw(undefinedParameters()); };
    static inline magnitudeType rmax()  { throw(undefinedParameters()); };
    static inline magnitudeType magnitude(T a) { throw(unsupportedType()); };
    static inline T zero()                     { throw(unsupportedType()); };
    static inline T one()                      { throw(unsupportedType()); };
    static inline T random()                   { throw(unsupportedType()); };
    static inline const char* name()           { throw(unsupportedType()); };
  };
  
  template<>
  struct ScalarTraits<int>
  {
    typedef long long magnitudeType;
    
    static inline int undefinedParameters()
    {
#ifndef TEUCHOS_NO_ERROR_REPORTS
      cerr << endl << "Teuchos::ScalarTraits: Machine parameters are undefined for this scalar type." << endl;
#endif
      return(-1);
    }

    static inline bool haveMachineParameters() { return false; };
    static inline magnitudeType eps()   { throw(undefinedParameters()); };
    static inline magnitudeType sfmin() { throw(undefinedParameters()); };
    static inline magnitudeType base()  { throw(undefinedParameters()); };
    static inline magnitudeType prec()  { throw(undefinedParameters()); };
    static inline magnitudeType t()     { throw(undefinedParameters()); };
    static inline magnitudeType rnd()   { throw(undefinedParameters()); };
    static inline magnitudeType emin()  { throw(undefinedParameters()); };
    static inline magnitudeType rmin()  { throw(undefinedParameters()); };
    static inline magnitudeType emax()  { throw(undefinedParameters()); };
    static inline magnitudeType rmax()  { throw(undefinedParameters()); };
    static inline magnitudeType magnitude(int a) { return abs(a); };
    static inline int zero()  { return 0; };
    static inline int one()   { return 1; };
    static inline int random() { return(rand() / RAND_MAX); };
    static inline const char * name() { return "int"; };
  };
  
  
  template<>
  struct ScalarTraits<float>
  {
    typedef float magnitudeType;
    static inline bool haveMachineParameters() { return true; };
    static inline float eps()   { LAPACK<int, float> lp; return lp.LAMCH('E'); };
    static inline float sfmin() { LAPACK<int, float> lp; return lp.LAMCH('S'); };
    static inline float base()  { LAPACK<int, float> lp; return lp.LAMCH('B'); };
    static inline float prec()  { LAPACK<int, float> lp; return lp.LAMCH('P'); };
    static inline float t()     { LAPACK<int, float> lp; return lp.LAMCH('N'); };
    static inline float rnd()   { LAPACK<int, float> lp; return lp.LAMCH('R'); };
    static inline float emin()  { LAPACK<int, float> lp; return lp.LAMCH('M'); };
    static inline float rmin()  { LAPACK<int, float> lp; return lp.LAMCH('U'); };
    static inline float emax()  { LAPACK<int, float> lp; return lp.LAMCH('L'); };
    static inline float rmax()  { LAPACK<int, float> lp; return lp.LAMCH('O'); };
    static inline magnitudeType magnitude(float a) { return fabs(a); };    
    static inline float zero()  { return(0.0); };
    static inline float one()   { return(1.0); };    
    static inline float random() { float rnd = (float)std::rand() / RAND_MAX; return (float)(-1.0 + 2.0 * rnd); };
    static inline const char* name() { return "float"; };
  };
  
  
  template<>
  struct ScalarTraits<double>
  {
    typedef double magnitudeType;
    static inline bool haveMachineParameters() { return true; };
    static inline magnitudeType magnitude(double a) { return fabs(a); };
    static inline double zero()  { return 0.0; };
    static inline double one()   { return 1.0; };
    static inline double eps()   { LAPACK<int, double> lp; return lp.LAMCH('E'); };
    static inline double sfmin() { LAPACK<int, double> lp; return lp.LAMCH('S'); };
    static inline double base()  { LAPACK<int, double> lp; return lp.LAMCH('B'); };
    static inline double prec()  { LAPACK<int, double> lp; return lp.LAMCH('P'); };
    static inline double t()     { LAPACK<int, double> lp; return lp.LAMCH('N'); };
    static inline double rnd()   { LAPACK<int, double> lp; return lp.LAMCH('R'); };
    static inline double emin()  { LAPACK<int, double> lp; return lp.LAMCH('M'); };
    static inline double rmin()  { LAPACK<int, double> lp; return lp.LAMCH('U'); };
    static inline double emax()  { LAPACK<int, double> lp; return lp.LAMCH('L'); };
    static inline double rmax()  { LAPACK<int, double> lp; return lp.LAMCH('O'); };
    static inline double random() { double rnd = (double)std::rand() / RAND_MAX; return (double)(-1.0 + 2.0 * rnd); };
    static inline const char* name() { return "double"; };
  };
  
#ifndef NO_COMPLEX
  
  template<> 
  struct ScalarTraits<ComplexFloat>
  {
    typedef float magnitudeType;
    static inline bool haveMachineParameters() { return true; };
    static inline float eps()   { LAPACK<int, float> lp; return lp.LAMCH('E'); };
    static inline float sfmin() { LAPACK<int, float> lp; return lp.LAMCH('S'); };
    static inline float base()  { LAPACK<int, float> lp; return lp.LAMCH('B'); };
    static inline float prec()  { LAPACK<int, float> lp; return lp.LAMCH('P'); };
    static inline float t()     { LAPACK<int, float> lp; return lp.LAMCH('N'); };
    static inline float rnd()   { LAPACK<int, float> lp; return lp.LAMCH('R'); };
    static inline float emin()  { LAPACK<int, float> lp; return lp.LAMCH('M'); };
    static inline float rmin()  { LAPACK<int, float> lp; return lp.LAMCH('U'); };
    static inline float emax()  { LAPACK<int, float> lp; return lp.LAMCH('L'); };
    static inline float rmax()  { LAPACK<int, float> lp; return lp.LAMCH('O'); };
    static magnitudeType magnitude(ComplexFloat a) { return std::abs(a); };
    static inline ComplexFloat zero()  { return ComplexFloat(0.0, 0.0); };
    static inline ComplexFloat one()   { return ComplexFloat(1.0, 0.0); };
    static inline ComplexFloat random()
    {
      float rnd1 = ScalarTraits<magnitudeType>::random();
      float rnd2 = ScalarTraits<magnitudeType>::random();
      return ComplexFloat(rnd1, rnd2);
    };
    static inline const char* name() { return "std::complex<float>"; };
  };
  
  template<>
  struct ScalarTraits<ComplexDouble>
  {
    typedef double magnitudeType;
    static inline bool haveMachineParameters() { return true; };
    static inline double eps()   { LAPACK<int, double> lp; return lp.LAMCH('E'); };
    static inline double sfmin() { LAPACK<int, double> lp; return lp.LAMCH('S'); };
    static inline double base()  { LAPACK<int, double> lp; return lp.LAMCH('B'); };
    static inline double prec()  { LAPACK<int, double> lp; return lp.LAMCH('P'); };
    static inline double t()     { LAPACK<int, double> lp; return lp.LAMCH('N'); };
    static inline double rnd()   { LAPACK<int, double> lp; return lp.LAMCH('R'); };
    static inline double emin()  { LAPACK<int, double> lp; return lp.LAMCH('M'); };
    static inline double rmin()  { LAPACK<int, double> lp; return lp.LAMCH('U'); };
    static inline double emax()  { LAPACK<int, double> lp; return lp.LAMCH('L'); };
    static inline double rmax()  { LAPACK<int, double> lp; return lp.LAMCH('O'); };
    static magnitudeType magnitude(ComplexDouble a) { return std::abs(a); };
    static inline ComplexDouble zero()  {return ComplexDouble(0.0,0.0); };
    static inline ComplexDouble one()   {return ComplexDouble(1.0,0.0); };    
    static inline ComplexDouble random()
    {
      double rnd1 = ScalarTraits<magnitudeType>::random();
      double rnd2 = ScalarTraits<magnitudeType>::random();
      return ComplexDouble(rnd1, rnd2);
    };
    static inline const char* name() { return "std::complex<double>"; };
  };
  
#endif  // NO_COMPLEX
  
} // Teuchos namespace

#endif // _TEUCHOS_SCALARTRAITS_HPP_
