// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// Kris
// 06.18.03 -- Minor formatting changes
//          -- Changed calls to LAPACK objects to use new <OType, SType> templates
// 07.08.03 -- Move into Teuchos package/namespace
// 07.11.03 -- Added ScalarTraits for ARPREC::mp_real
// 07.14.03 -- Fixed int rand() function (was set up to return a floating-point style random number)
// 07.17.03 -- Added squareroot() function

#ifndef _TEUCHOS_SCALARTRAITS_HPP_
#define _TEUCHOS_SCALARTRAITS_HPP_

/*! \file Teuchos_ScalarTraits.hpp
    \brief Defines basic traits for the scalar field type
*/

#include "Teuchos_LAPACK.hpp"

#ifdef HAVE_TEUCHOS_ARPREC
#include "mp/mpreal.h"
#endif

/*! \struct Teuchos::ScalarTraits
    \brief This structure defines some basic traits for a scalar field type.

    Scalar traits are an essential part of templated codes.  This structure offers
    the basic traits of the templated scalar type, like defining zero and one,
    and basic functions on the templated scalar type, like performing a square root.

    The functions in the templated base unspecialized struct are designed not to
    compile (giving a nice compile-time error message) and therefore specializations
    must be written for Scalar types actually used.

    \note 
     <ol>
       <li> The default defined specializations are provided for \c int, \c float, and \c double.
     	 <li> ScalarTraits can be used with the Arbitrary Precision Library ( \c http://crd.lbl.gov/~dhbailey/mpdist/ )
            by configuring Teuchos with \c --enable-teuchos-arprec and giving the appropriate paths to ARPREC.
            Then ScalarTraits has the specialization: \c mp_real.
     	 <li> If Teuchos is configured with \c --enable-teuchos-complex then ScalarTraits also has
            a parital specialization for all complex numbers of the form <tt>complex<T></tt>.
     </ol>
*/

namespace Teuchos {

  template<class T>
  struct UndefinedScalarTraits
  {
    // This function should not compile if there is an attempt to instantiate!
    static inline T notDefined() { return T::this_type_is_missing_a_specialization(); };
  };

  template<class T>
  struct ScalarTraits
  {
    //! Madatory typedef for result of magnitude
    typedef T magnitudeType;
    //! Does this scalar type have machine-specific parameters.
    static inline bool haveMachineParameters() { return false; };
    //! Returns relative machine precision.
    static inline magnitudeType eps()   { return UndefinedScalarTraits<T>::notDefined(); };
    //! Returns safe minimum (sfmin), such that 1/sfmin does not overflow.
    static inline magnitudeType sfmin() { return UndefinedScalarTraits<T>::notDefined(); };
    //! Returns the base of the machine.
    static inline magnitudeType base()  { return UndefinedScalarTraits<T>::notDefined(); };
    //! Returns \c eps*base.
    static inline magnitudeType prec()  { return UndefinedScalarTraits<T>::notDefined(); };
    //! Returns the number of (base) digits in the mantissa.
    static inline magnitudeType t()     { return UndefinedScalarTraits<T>::notDefined(); };
    //! Returns 1.0 when rounding occurs in addition, 0.0 otherwise
    static inline magnitudeType rnd()   { return UndefinedScalarTraits<T>::notDefined(); };
    //! Returns the minimum exponent before (gradual) underflow.
    static inline magnitudeType emin()  { return UndefinedScalarTraits<T>::notDefined(); };
    //! Returns the underflow threshold - \c base^(emin-1)
    static inline magnitudeType rmin()  { return UndefinedScalarTraits<T>::notDefined(); };
    //! Returns the largest exponent before overflow.
    static inline magnitudeType emax()  { return UndefinedScalarTraits<T>::notDefined(); };
    //! Overflow theshold - \c (base^emax)*(1-eps)
    static inline magnitudeType rmax()  { return UndefinedScalarTraits<T>::notDefined(); };
    //! Returns the magnitudeType of the scalar type \c a.
    static inline magnitudeType magnitude(T a) { return UndefinedScalarTraits<T>::notDefined(); };
    //! Returns representation of zero for this scalar type.
    static inline T zero()                     { return UndefinedScalarTraits<T>::notDefined(); };
    //! Returns representation of one for this scalar type.
    static inline T one()                      { return UndefinedScalarTraits<T>::notDefined(); };
    //! Seed the random number generator returned by <tt>random()</tt>.
    static inline void seedrandom(unsigned int s) { int i; T t = &i; };
    //! Returns a random number (between -one() and +one()) of this scalar type.
    static inline T random()                   { return UndefinedScalarTraits<T>::notDefined(); };
    //! Returns the name of this scalar type.
    static inline std::string name()           { (void)UndefinedScalarTraits<T>::notDefined(); return 0; };
    //! Returns a number of magnitudeType that is the square root of this scalar type \c x. 
    static inline magnitudeType squareroot(T x) { return UndefinedScalarTraits<T>::notDefined(); };
  };
  
#ifndef DOXYGEN_SHOULD_SKIP_THIS

  template<>
  struct ScalarTraits<int>
  {
    typedef int magnitudeType;
    static inline bool haveMachineParameters() { return false; };
    // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
    static inline magnitudeType magnitude(int a) { return abs(a); };
    static inline int zero()  { return 0; };
    static inline int one()   { return 1; };
    static inline void seedrandom(unsigned int s) { srand(s); };
    //static inline int random() { return (-1 + 2*rand()); };  // RAB: This version should be used to be consistent with others
    static inline int random() { return rand(); };             // RAB: This version should be used for an unsigned int, not int
    static inline std::string name() { return "int"; };
    static inline int squareroot(int x) { return (int) sqrt((double) x); };
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
    static inline void seedrandom(unsigned int s) { srand(s); };
    static inline float random() { float rnd = (float) rand() / RAND_MAX; return (float)(-1.0 + 2.0 * rnd); };
    static inline std::string name() { return "float"; };
    static inline float squareroot(float x) { return sqrt(x); };
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
    static inline void seedrandom(unsigned int s) { srand(s); };
    static inline double random() { double rnd = (double) rand() / RAND_MAX; return (double)(-1.0 + 2.0 * rnd); };
    static inline std::string name() { return "double"; };
    static inline double squareroot(double x) { return sqrt(x); };
  };

#ifdef HAVE_TEUCHOS_ARPREC

  template<>
  struct ScalarTraits<mp_real>
  {
    typedef mp_real magnitudeType;
    static inline bool haveMachineParameters() { return false; };
    // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
    static magnitudeType magnitude(mp_real a) { return abs(a); };
    static inline mp_real zero() { mp_real zero = 0.0; return zero; };
    static inline mp_real one() { mp_real one = 1.0; return one; };    
    static inline void seedrandom(unsigned int s) { 
      long int seedVal = static_cast<long int>(s);
      srand48(seedVal);
    };
    static inline mp_real random() { return mp_rand(); };
    static inline std::string name() { return "mp_real"; };
    static inline mp_real squareroot(mp_real x) { return sqrt(x); };
  };
  
#endif // HAVE_TEUCHOS_ARPREC
 
#if ( defined(HAVE_COMPLEX) || defined(HAVE_COMPLEX_H) ) && defined(HAVE_TEUCHOS_COMPLEX)

  // Partial specialization for complex numbers templated on real type T
  template<class T> 
  struct ScalarTraits<
#if defined(HAVE_COMPLEX)
    std::complex<T>
#elif  defined(HAVE_COMPLEX_H)
    ::complex<T>
#endif
    >
  {
#if defined(HAVE_COMPLEX)
    typedef std::complex<T>  ComplexT;
#elif  defined(HAVE_COMPLEX_H)
    typedef ::complex<T>     ComplexT;
#endif
    typedef typename ScalarTraits<T>::magnitudeType magnitudeType;
    static inline bool haveMachineParameters() { return true; };
    static inline magnitudeType eps()          { return ScalarTraits<magnitudeType>::eps(); };
    static inline magnitudeType sfmin()        { return ScalarTraits<magnitudeType>::sfmin(); };
    static inline magnitudeType base()         { return ScalarTraits<magnitudeType>::base(); };
    static inline magnitudeType prec()         { return ScalarTraits<magnitudeType>::prec(); };
    static inline magnitudeType t()            { return ScalarTraits<magnitudeType>::t(); };
    static inline magnitudeType rnd()          { return ScalarTraits<magnitudeType>::rnd(); };
    static inline magnitudeType emin()         { return ScalarTraits<magnitudeType>::emin(); };
    static inline magnitudeType rmin()         { return ScalarTraits<magnitudeType>::rmin(); };
    static inline magnitudeType emax()         { return ScalarTraits<magnitudeType>::emax(); };
    static inline magnitudeType rmax()         { return ScalarTraits<magnitudeType>::rmax(); };
    static magnitudeType magnitude(ComplexT a) { return std::fabs(a); };
    static inline ComplexT zero()              { return ComplexT(ScalarTraits<magnitudeType>::zero(),ScalarTraits<magnitudeType>::zero()); };
    static inline ComplexT one()               { return ComplexT(ScalarTraits<magnitudeType>::one(),ScalarTraits<magnitudeType>::zero());; };
    static inline void seedrandom(unsigned int s) { ScalarTraits<magnitudeType>::seedrandom(s); };
    static inline ComplexT random()
    {
      const T rnd1 = ScalarTraits<magnitudeType>::random();
      const T rnd2 = ScalarTraits<magnitudeType>::random();
      return ComplexT(rnd1,rnd2);
    };
    static inline std::string name() { return std::string("std::complex<")+std::string(ScalarTraits<magnitudeType>::name())+std::string(">"); };
    // This will only return one of the square roots of x, the other can be obtained by taking its conjugate
    static inline ComplexT squareroot(ComplexT x)
    {
      typedef ScalarTraits<magnitudeType>  STMT;
      const T r  = x.real(), i = x.imag();
      const T a  = STMT::squareroot((r*r)+(i*i));
      const T nr = STMT::squareroot((a+r)/2);
      const T ni = STMT::squareroot((a-r)/2);
      return ComplexT(nr,ni);
    };
  };

#endif //  HAVE_COMPLEX || HAVE_COMPLEX_H

#endif // DOXYGEN_SHOULD_SKIP_THIS

} // Teuchos namespace

#endif // _TEUCHOS_SCALARTRAITS_HPP_
