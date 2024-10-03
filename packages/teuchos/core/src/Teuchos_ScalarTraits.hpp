// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// NOTE: Before adding specializations of ScalarTraits, make sure that they do
// not duplicate specializations already present in PyTrilinos (see
// packages/PyTrilinos/src/Teuchos_Traits.i)

// NOTE: halfPrecision and doublePrecision are not currently implemented for
// ARPREC, GMP or the ordinal types (e.g., int, char)

#ifndef _TEUCHOS_SCALARTRAITS_HPP_
#define _TEUCHOS_SCALARTRAITS_HPP_

/*! \file Teuchos_ScalarTraits.hpp
    \brief Defines basic traits for the scalar field type.
*/

#include "Teuchos_ConfigDefs.hpp"

#ifdef HAVE_TEUCHOSCORE_KOKKOS
#include "Kokkos_Complex.hpp"
#endif // HAVE_TEUCHOSCORE_KOKKOS

#ifdef HAVE_TEUCHOS_ARPREC
#include <arprec/mp_real.h>
#endif

#ifdef HAVE_TEUCHOSCORE_QUADMATH
#include <quadmath.h>

// Teuchos_ConfigDefs.hpp includes <iostream>, which includes
// <ostream>.  If this ever changes, include <ostream> here.

namespace std {

/// \brief Overload of operator<<(std::ostream&, const T&) for T =
///   __float128.
///
/// GCC / libquadmath doesn't implement an std::ostream operator<< for
/// __float128, so we have to write our own.  At least libquadmath
/// provides a printing function specifically for __float128.
///
/// FIXME (mfh 19 Mar 2015) This will break if users have already
/// defined their own operator<< in this namespace.
ostream&
operator<< (std::ostream& out, const __float128& x);

/// \brief Overload of operator>>(std::istream&, const T&) for T =
///   __float128.
///
/// GCC / libquadmath doesn't implement an std::istream operator>> for
/// __float128, so we have to write our own.  At least libquadmath
/// provides an input function specifically for __float128.
///
/// FIXME (mfh 10 Sep 2015) This will break if users have already
/// defined their own operator>> in this namespace.
istream&
operator>> (std::istream& in, __float128& x);

} // namespace std

#endif // HAVE_TEUCHOSCORE_QUADMATH

#ifdef HAVE_TEUCHOS_QD
#include <qd/qd_real.h>
#include <qd/dd_real.h>
#endif

#ifdef HAVE_TEUCHOS_GNU_MP
#include <gmp.h>
#include <gmpxx.h>
#endif


#include "Teuchos_ScalarTraitsDecl.hpp"


namespace Teuchos {


#ifndef DOXYGEN_SHOULD_SKIP_THIS


TEUCHOSCORE_LIB_DLL_EXPORT
void throwScalarTraitsNanInfError( const std::string &errMsg );


template<class Scalar>
bool generic_real_isnaninf(const Scalar &x)
{
#ifdef HAVE_TEUCHOSCORE_CXX11
  if (std::isnan(x)) return true;
  if (std::isinf(x)) return true;
  return false;
#else
  typedef std::numeric_limits<Scalar> STD_NL;
  // IEEE says this should fail for NaN (not all compilers do not obey IEEE)
  const Scalar tol = 1.0; // Any (bounded) number should do!
  if (!(x <= tol) && !(x > tol)) return true;
  // Use fact that Inf*0 = NaN (again, all compilers do not obey IEEE)
  Scalar z = static_cast<Scalar>(0.0) * x;
  if (!(z <= tol) && !(z > tol)) return true;
  // As a last result use comparisons
  if (x == STD_NL::infinity() || x == -STD_NL::infinity()) return true;
  // We give up and assume the number is finite
  return false;
#endif
}


#define TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR( VALUE, MSG ) \
  if (isnaninf(VALUE)) { \
    std::ostringstream omsg; \
    omsg << MSG; \
    Teuchos::throwScalarTraitsNanInfError(omsg.str());  \
  }


template<>
struct ScalarTraits<char>
{
  typedef char magnitudeType;
  typedef char halfPrecision;
  typedef char doublePrecision;
  typedef char coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static inline magnitudeType magnitude(char a) { return static_cast<char>(std::fabs(static_cast<double>(a))); }
  static inline char zero()  { return 0; }
  static inline char one()   { return 1; }
  static inline char conjugate(char x) { return x; }
  static inline char real(char x) { return x; }
  static inline char imag(char) { return 0; }
  static inline bool isnaninf(char ) { return false; }
  static inline void seedrandom(unsigned int s) {
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  //static inline char random() { return (-1 + 2*rand()); } // RAB: This version should be used to be consistent with others
  static inline char random() { return std::rand(); } // RAB: This version should be used for an unsigned char, not char
  static inline std::string name() { return "char"; }
  static inline char squareroot(char x) { return (char) std::sqrt((double) x); }
  static inline char pow(char x, char y) { return (char) std::pow((double)x,(double)y); }
  static inline char log(char x) { return static_cast<char> (std::log (static_cast<double> (x))); }
  static inline char log10(char x) { return static_cast<char> (std::log10 (static_cast<double> (x))); }
};


template<>
struct ScalarTraits<short int>
{
  typedef short int magnitudeType;
  typedef short int halfPrecision;
  typedef short int doublePrecision;
  typedef short int coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static inline magnitudeType magnitude(short int a) { return static_cast<short int>(std::fabs(static_cast<double>(a))); }
  static inline short int zero()  { return 0; }
  static inline short int one()   { return 1; }
  static inline short int conjugate(short int x) { return x; }
  static inline short int real(short int x) { return x; }
  static inline short int imag(short int) { return 0; }
  static inline bool isnaninf(short int) { return false; }
  static inline void seedrandom(unsigned int s) {
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  //static inline int random() { return (-1 + 2*rand()); }  // RAB: This version should be used to be consistent with others
  static inline short int random() { return std::rand(); }             // RAB: This version should be used for an unsigned int, not int
  static inline std::string name() { return "short int"; }
  static inline short int squareroot(short int x) { return (short int) std::sqrt((double) x); }
  static inline short int pow(short int x, short int y) { return (short int) std::pow((double)x,(double)y); }
  static inline short int log(short int x) { return static_cast<short int> (std::log (static_cast<double> (x))); }
  static inline short int log10(short int x) { return static_cast<short int> (std::log10 (static_cast<double> (x))); }
};

template<>
struct ScalarTraits<unsigned short int>
{
  typedef unsigned short int magnitudeType;
  typedef unsigned short int halfPrecision;
  typedef unsigned short int doublePrecision;
  typedef unsigned short int coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static inline magnitudeType magnitude(unsigned short int a) { return static_cast<unsigned short int>(std::fabs(static_cast<double>(a))); }
  static inline unsigned short int zero()  { return 0; }
  static inline unsigned short int one()   { return 1; }
  static inline unsigned short int conjugate(unsigned short int x) { return x; }
  static inline unsigned short int real(unsigned short int x) { return x; }
  static inline unsigned short int imag(unsigned short int) { return 0; }
  static inline bool isnaninf(unsigned short int) { return false; }
  static inline void seedrandom(unsigned int s) {
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  //static inline int random() { return (-1 + 2*rand()); }  // RAB: This version should be used to be consistent with others
  static inline unsigned short int random() { return std::rand(); }             // RAB: This version should be used for an unsigned int, not int
  static inline std::string name() { return "unsigned short int"; }
  static inline unsigned short int squareroot(unsigned short int x) { return (unsigned short int) std::sqrt((double) x); }
  static inline unsigned short int pow(unsigned short int x, unsigned short int y) { return (unsigned short int) std::pow((double)x,(double)y); }
  static inline unsigned short int log(unsigned short int x) { return static_cast<unsigned short int> (std::log (static_cast<double> (x))); }
  static inline unsigned short int log10(unsigned short int x) { return static_cast<unsigned short int> (std::log10 (static_cast<double> (x))); }
};


template<>
struct ScalarTraits<int>
{
  typedef int magnitudeType;
  typedef int halfPrecision;
  typedef int doublePrecision;
  typedef int coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static inline magnitudeType magnitude(int a) { return static_cast<int>(std::fabs(static_cast<double>(a))); }
  static inline int zero()  { return 0; }
  static inline int one()   { return 1; }
  static inline int conjugate(int x) { return x; }
  static inline int real(int x) { return x; }
  static inline int imag(int) { return 0; }
  static inline bool isnaninf(int) { return false; }
  static inline void seedrandom(unsigned int s) {
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  //static inline int random() { return (-1 + 2*rand()); }  // RAB: This version should be used to be consistent with others
  static inline int random() { return std::rand(); }             // RAB: This version should be used for an unsigned int, not int
  static inline std::string name() { return "int"; }
  static inline int squareroot(int x) { return (int) std::sqrt((double) x); }
  static inline int pow(int x, int y) { return (int) std::pow((double)x,(double)y); }
  static inline int log(int x) { return static_cast<int> (std::log (static_cast<double> (x))); }
  static inline int log10(int x) { return static_cast<int> (std::log10 (static_cast<double> (x))); }
};


template<>
struct ScalarTraits<unsigned int>
{
  typedef unsigned int magnitudeType;
  typedef unsigned int halfPrecision;
  typedef unsigned int doublePrecision;
  typedef unsigned int coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static inline magnitudeType magnitude(unsigned int a) { return static_cast<unsigned int>(std::fabs(static_cast<double>(a))); }
  static inline unsigned int zero()  { return 0; }
  static inline unsigned int one()   { return 1; }
  static inline unsigned int conjugate(unsigned int x) { return x; }
  static inline unsigned int real(unsigned int x) { return x; }
  static inline unsigned int imag(unsigned int) { return 0; }
  static inline bool isnaninf(unsigned int) { return false; }
  static inline void seedrandom(unsigned int s) {
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  //static inline int random() { return (-1 + 2*rand()); }  // RAB: This version should be used to be consistent with others
  static inline unsigned int random() { return std::rand(); }             // RAB: This version should be used for an unsigned int, not int
  static inline std::string name() { return "unsigned int"; }
  static inline unsigned int squareroot(unsigned int x) { return (unsigned int) std::sqrt((double) x); }
  static inline unsigned int pow(unsigned int x, unsigned int y) { return (unsigned int) std::pow((double)x,(double)y); }
  static inline unsigned int log(unsigned int x) { return static_cast<unsigned int> (std::log (static_cast<double> (x))); }
  static inline unsigned int log10(unsigned int x) { return static_cast<unsigned int> (std::log10 (static_cast<double> (x))); }
};


template<>
struct ScalarTraits<long int>
{
  typedef long int magnitudeType;
  typedef long int halfPrecision;
  typedef long int doublePrecision;
  typedef long int coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static inline magnitudeType magnitude(long int a) { return static_cast<long int>(std::fabs(static_cast<double>(a))); }
  static inline long int zero()  { return 0; }
  static inline long int one()   { return 1; }
  static inline long int conjugate(long int x) { return x; }
  static inline long int real(long int x) { return x; }
  static inline long int imag(long int) { return 0; }
  static inline bool isnaninf(long int) { return false; }
  static inline void seedrandom(unsigned int s) {
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  //static inline int random() { return (-1 + 2*rand()); }  // RAB: This version should be used to be consistent with others
  static inline long int random() { return std::rand(); }             // RAB: This version should be used for an unsigned int, not int
  static inline std::string name() { return "long int"; }
  static inline long int squareroot(long int x) { return (long int) std::sqrt((double) x); }
  static inline long int pow(long int x, long int y) { return (long int) std::pow((double)x,(double)y); }
  // Note: Depending on the number of bits in long int, the cast from
  // long int to double may not be exact.
  static inline long int log(long int x) { return static_cast<long int> (std::log (static_cast<double> (x))); }
  static inline long int log10(long int x) { return static_cast<long int> (std::log10 (static_cast<double> (x))); }
};


template<>
struct ScalarTraits<long unsigned int>
{
  typedef long unsigned int magnitudeType;
  typedef long unsigned int halfPrecision;
  typedef long unsigned int doublePrecision;
  typedef long unsigned int coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static inline magnitudeType magnitude(long unsigned int a) { return static_cast<long unsigned int>(std::fabs(static_cast<double>(a))); }
  static inline long unsigned int zero()  { return 0; }
  static inline long unsigned int one()   { return 1; }
  static inline long unsigned int conjugate(long unsigned int x) { return x; }
  static inline long unsigned int real(long unsigned int x) { return x; }
  static inline long unsigned int imag(long unsigned int) { return 0; }
  static inline bool isnaninf(long unsigned int) { return false; }
  static inline void seedrandom(unsigned int s) {
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  //static inline int random() { return (-1 + 2*rand()); }  // RAB: This version should be used to be consistent with others
  static inline long unsigned int random() { return std::rand(); }             // RAB: This version should be used for an unsigned int, not int
  static inline std::string name() { return "long unsigned int"; }
  static inline long unsigned int squareroot(long unsigned int x) { return (long unsigned int) std::sqrt((double) x); }
  static inline long unsigned int pow(long unsigned int x, long unsigned int y) { return (long unsigned int) std::pow((double)x,(double)y); }
  // Note: Depending on the number of bits in long unsigned int, the
  // cast from long unsigned int to double may not be exact.
  static inline long unsigned int log(long unsigned int x) { return static_cast<long unsigned int> (std::log (static_cast<double> (x))); }
  static inline long unsigned int log10(long unsigned int x) { return static_cast<long unsigned int> (std::log10 (static_cast<double> (x))); }
};


template<>
struct ScalarTraits<long long int>
{
  typedef long long int magnitudeType;
  typedef long long int halfPrecision;
  typedef long long int doublePrecision;
  typedef long long int coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static inline magnitudeType magnitude(long long int a) { return static_cast<long long int>(std::fabs(static_cast<double>(a))); }
  static inline long long int zero()  { return 0; }
  static inline long long int one()   { return 1; }
  static inline long long int conjugate(long long int x) { return x; }
  static inline long long int real(long long int x) { return x; }
  static inline long long int imag(long long int) { return 0; }
  static inline bool isnaninf(long long int) { return false; }
  static inline void seedrandom(unsigned int s) {
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  //static inline int random() { return (-1 + 2*rand()); }  // RAB: This version should be used to be consistent with others
  static inline long long int random() { return std::rand(); }             // RAB: This version should be used for an unsigned int, not int
  static inline std::string name() { return "long long int"; }
  static inline long long int squareroot(long long int x) { return (long long int) std::sqrt((double) x); }
  static inline long long int pow(long long int x, long long int y) { return (long long int) std::pow((double)x,(double)y); }
  // Note: Depending on the number of bits in long long int, the cast
  // from long long int to double may not be exact.
  static inline long long int log(long long int x) { return static_cast<long long int> (std::log (static_cast<double> (x))); }
  static inline long long int log10(long long int x) { return static_cast<long long int> (std::log10 (static_cast<double> (x))); }
};

template<>
struct ScalarTraits<unsigned long long int>
{
  typedef unsigned long long int magnitudeType;
  typedef unsigned long long int halfPrecision;
  typedef unsigned long long int doublePrecision;
  typedef unsigned long long int coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static inline magnitudeType magnitude(unsigned long long int a) { return static_cast<unsigned long long int>(std::fabs(static_cast<double>(a))); }
  static inline unsigned long long int zero()  { return 0; }
  static inline unsigned long long int one()   { return 1; }
  static inline unsigned long long int conjugate(unsigned long long int x) { return x; }
  static inline unsigned long long int real(unsigned long long int x) { return x; }
  static inline unsigned long long int imag(unsigned long long int) { return 0; }
  static inline bool isnaninf(unsigned long long int) { return false; }
  static inline void seedrandom(unsigned int s) {
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  //static inline int random() { return (-1 + 2*rand()); }  // RAB: This version should be used to be consistent with others
  static inline unsigned long long int random() { return std::rand(); }             // RAB: This version should be used for an unsigned int, not int
  static inline std::string name() { return "unsigned long long int"; }
  static inline unsigned long long int squareroot(unsigned long long int x) { return (unsigned long long int) std::sqrt((double) x); }
  static inline unsigned long long int pow(unsigned long long int x, unsigned long long int y) { return (unsigned long long int) std::pow((double)x,(double)y); }
  // Note: Depending on the number of bits in unsigned long long int,
  // the cast from unsigned long long int to double may not be exact.
  static inline unsigned long long int log(unsigned long long int x) { return static_cast<unsigned long long int> (std::log (static_cast<double> (x))); }
  static inline unsigned long long int log10(unsigned long long int x) { return static_cast<unsigned long long int> (std::log10 (static_cast<double> (x))); }
};


#ifdef HAVE_TEUCHOS___INT64

template<>
struct ScalarTraits<__int64>
{
  typedef __int64 magnitudeType;
  typedef __int64 halfPrecision;
  typedef __int64 doublePrecision;
  typedef __int64 coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static inline magnitudeType magnitude(__int64 a) { return static_cast<__int64>(std::fabs(static_cast<double>(a))); }
  static inline __int64 zero()  { return 0; }
  static inline __int64 one()   { return 1; }
  static inline __int64 conjugate(__int64 x) { return x; }
  static inline __int64 real(__int64 x) { return x; }
  static inline __int64 imag(__int64) { return 0; }
  static inline void seedrandom(unsigned int s) {
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  //static inline int random() { return (-1 + 2*rand()); }  // RAB: This version should be used to be consistent with others
  static inline __int64 random() { return std::rand(); }             // RAB: This version should be used for an unsigned int, not int
  static inline std::string name() { return "__int64"; }
  static inline __int64 squareroot(__int64 x) { return (__int64) std::sqrt((double) x); }
  static inline __int64 pow(__int64 x, __int64 y) { return (__int64) std::pow((double)x,(double)y); }
  // Note: Depending on the number of bits in __int64, the cast
  // from __int64 to double may not be exact.
  static inline __int64 log(__int64 x) { return static_cast<__int64> (std::log (static_cast<double> (x))); }
  static inline __int64 log10(__int64 x) { return static_cast<__int64> (std::log10 (static_cast<double> (x))); }
};

template<>
struct ScalarTraits<unsigned __int64>
{
  typedef unsigned __int64 magnitudeType;
  typedef unsigned __int64 halfPrecision;
  typedef unsigned __int64 doublePrecision;
  typedef unsigned __int64 coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static inline magnitudeType magnitude(unsigned __int64 a) { return static_cast<unsigned __int64>(std::fabs(static_cast<double>(a))); }
  static inline unsigned __int64 zero()  { return 0; }
  static inline unsigned __int64 one()   { return 1; }
  static inline unsigned __int64 conjugate(unsigned __int64 x) { return x; }
  static inline unsigned __int64 real(unsigned __int64 x) { return x; }
  static inline unsigned __int64 imag(unsigned __int64) { return 0; }
  static inline void seedrandom(unsigned int s) {
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  //static inline int random() { return (-1 + 2*rand()); }  // RAB: This version should be used to be consistent with others
  static inline unsigned __int64 random() { return std::rand(); }             // RAB: This version should be used for an unsigned int, not int
  static inline std::string name() { return "unsigned __int64"; }
  static inline unsigned __int64 squareroot(unsigned __int64 x) { return (unsigned __int64) std::sqrt((double) x); }
  static inline unsigned __int64 pow(unsigned __int64 x, unsigned __int64 y) { return (unsigned __int64) std::pow((double)x,(double)y); }
  // Note: Depending on the number of bits in unsigned __int64,
  // the cast from unsigned __int64 to double may not be exact.
  static inline unsigned __int64 log(unsigned __int64 x) { return static_cast<unsigned __int64> (std::log (static_cast<double> (x))); }
  static inline unsigned __int64 log10(unsigned __int64 x) { return static_cast<unsigned __int64> (std::log10 (static_cast<double> (x))); }
};

#endif // HAVE_TEUCHOS___INT64


#ifndef __sun
extern TEUCHOSCORE_LIB_DLL_EXPORT const float flt_nan;
#endif


template<>
struct ScalarTraits<float>
{
  typedef float magnitudeType;
  typedef float halfPrecision; // should become IEEE754-2008 binary16 or fp16 later, perhaps specified at configure according to architectural support
  typedef double doublePrecision;
  typedef float coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;
  static inline float eps()   {
    return std::numeric_limits<float>::epsilon();
  }
  static inline float sfmin() {
    return std::numeric_limits<float>::min();
  }
  static inline float base()  {
    return static_cast<float>(std::numeric_limits<float>::radix);
  }
  static inline float prec()  {
    return eps()*base();
  }
  static inline float t()     {
    return static_cast<float>(std::numeric_limits<float>::digits);
  }
  static inline float rnd()   {
    return ( std::numeric_limits<float>::round_style == std::round_to_nearest ? one() : zero() );
  }
  static inline float emin()  {
    return static_cast<float>(std::numeric_limits<float>::min_exponent);
  }
  static inline float rmin()  {
    return std::numeric_limits<float>::min();
  }
  static inline float emax()  {
    return static_cast<float>(std::numeric_limits<float>::max_exponent);
  }
  static inline float rmax()  {
    return std::numeric_limits<float>::max();
  }
  static inline magnitudeType magnitude(float a)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        a, "Error, the input value to magnitude(...) a = " << a << " can not be NaN!" );
#endif
      return std::fabs(a);
    }
  static inline float zero()  { return(0.0f); }
  static inline float one()   { return(1.0f); }
  static inline float conjugate(float x)   { return(x); }
  static inline float real(float x) { return x; }
  static inline float imag(float) { return zero(); }
  static inline float nan() {
#ifdef __sun
    return 0.0f/std::sin(0.0f);
#else
    return std::numeric_limits<float>::quiet_NaN();
#endif
  }
  static inline bool isnaninf(float x) {
    return generic_real_isnaninf<float>(x);
  }
  static inline void seedrandom(unsigned int s) {
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  static inline float random() { float rnd = (float) std::rand() / static_cast<float>(RAND_MAX); return (-1.0f + 2.0f * rnd); }
  static inline std::string name() { return "float"; }
  static inline float squareroot(float x)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        x, "Error, the input value to squareroot(...) x = " << x << " can not be NaN!" );
#endif
      errno = 0;
      const float rtn = std::sqrt(x);
      if (errno)
        return nan();
      return rtn;
    }
  static inline float pow(float x, float y) { return std::pow(x,y); }
  static inline float pi() { return 3.14159265358979323846f; }
  static inline float log(float x) { return std::log(x); }
  static inline float log10(float x) { return std::log10(x); }
};


#ifndef __sun
extern TEUCHOSCORE_LIB_DLL_EXPORT const double dbl_nan;
#endif


template<>
struct ScalarTraits<double>
{
  typedef double magnitudeType;
  typedef float halfPrecision;
  /* there are different options as to how to double "double"
     - QD's DD (if available)
     - ARPREC
     - GNU MP
     - a true hardware quad

     in the shortterm, this should be specified at configure time. I have inserted a configure-time option (--enable-teuchos-double-to-dd)
     which uses QD's DD when available. This must be used alongside --enable-teuchos-qd.
   */
#if defined(HAVE_TEUCHOS_DOUBLE_TO_QD)
  typedef dd_real doublePrecision;
#elif defined(HAVE_TEUCHOS_DOUBLE_TO_ARPREC)
  typedef mp_real doublePrecision;
#else
  typedef double doublePrecision;     // don't double "double" in this case
#endif
  typedef double coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;
  static inline double eps()   {
    return std::numeric_limits<double>::epsilon();
  }
  static inline double sfmin() {
    return std::numeric_limits<double>::min();
  }
  static inline double base()  {
    return std::numeric_limits<double>::radix;
  }
  static inline double prec()  {
    return eps()*base();
  }
  static inline double t()     {
    return std::numeric_limits<double>::digits;
  }
  static inline double rnd()   {
    return ( std::numeric_limits<double>::round_style == std::round_to_nearest ? double(1.0) : double(0.0) );
  }
  static inline double emin()  {
    return std::numeric_limits<double>::min_exponent;
  }
  static inline double rmin()  {
    return std::numeric_limits<double>::min();
  }
  static inline double emax()  {
    return std::numeric_limits<double>::max_exponent;
  }
  static inline double rmax()  {
    return std::numeric_limits<double>::max();
  }
  static inline magnitudeType magnitude(double a)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        a, "Error, the input value to magnitude(...) a = " << a << " can not be NaN!" );
#endif
      return std::fabs(a);
    }
  static inline double zero()  { return 0.0; }
  static inline double one()   { return 1.0; }
  static inline double conjugate(double x)   { return(x); }
  static inline double real(double x) { return(x); }
  static inline double imag(double) { return(0); }
  static inline double nan() {
#ifdef __sun
    return 0.0/std::sin(0.0);
#else
    return std::numeric_limits<double>::quiet_NaN();
#endif
  }
  static inline bool isnaninf(double x) {
    return generic_real_isnaninf<double>(x);
  }
  static inline void seedrandom(unsigned int s) {
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  static inline double random() { double rnd = (double) std::rand() / RAND_MAX; return (double)(-1.0 + 2.0 * rnd); }
  static inline std::string name() { return "double"; }
  static inline double squareroot(double x)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        x, "Error, the input value to squareroot(...) x = " << x << " can not be NaN!" );
#endif
      errno = 0;
      const double rtn = std::sqrt(x);
      if (errno)
        return nan();
      return rtn;
    }
  static inline double pow(double x, double y) { return std::pow(x,y); }
  static inline double pi() { return 3.14159265358979323846; }
  static inline double log(double x) { return std::log(x); }
  static inline double log10(double x) { return std::log10(x); }
};


#ifdef HAVE_TEUCHOS_LONG_DOUBLE
template<>
struct ScalarTraits<long double>
{
  typedef long double magnitudeType;
  typedef double halfPrecision;
  typedef double doublePrecision;    
  typedef long double coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;
  static inline long double eps()   {
    return std::numeric_limits<long double>::epsilon();
  }
  static inline long double sfmin() {
    return std::numeric_limits<long double>::min();
  }
  static inline long double base()  {
    return std::numeric_limits<long double>::radix;
  }
  static inline long double prec()  {
    return eps()*base();
  }
  static inline long double t()     {
    return std::numeric_limits<long double>::digits;
  }
  static inline long double rnd()   {
    return ( std::numeric_limits<long double>::round_style == std::round_to_nearest ? (long double)(1.0) : (long double)(0.0) );
  }
  static inline long double emin()  {
    return std::numeric_limits<long double>::min_exponent;
  }
  static inline long double rmin()  {
    return std::numeric_limits<long double>::min();
  }
  static inline long double emax()  {
    return std::numeric_limits<long double>::max_exponent;
  }
  static inline long double rmax()  {
    return std::numeric_limits<long double>::max();
  }
  static inline magnitudeType magnitude(long double a)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        a, "Error, the input value to magnitude(...) a = " << a << " can not be NaN!" );
#endif
      return std::fabs(a);
    }
  static inline long double zero()  { return 0.0; }
  static inline long double one()   { return 1.0; }
  static inline long double conjugate(long double x)   { return(x); }
  static inline long double real(long double x) { return(x); }
  static inline long double imag(long double) { return(0); }
  static inline long double nan() {
#ifdef __sun
    return 0.0/std::sin(0.0);
#else
    return std::numeric_limits<long double>::quiet_NaN();
#endif
  }
  static inline bool isnaninf(long double x) {
    return generic_real_isnaninf<long double>(x);
  }
  static inline void seedrandom(unsigned int s) {
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  static inline long double random() { long double rnd = (long double) std::rand() / RAND_MAX; return (long double)(-1.0 + 2.0 * rnd); }
  static inline std::string name() { return "long double"; }
  static inline long double squareroot(long double x)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        x, "Error, the input value to squareroot(...) x = " << x << " can not be NaN!" );
#endif
      errno = 0;
      const long double rtn = std::sqrt(x);
      if (errno)
        return nan();
      return rtn;
    }
  static inline long double pow(long double x, long double y) { return std::pow(x,y); }
  static inline long double pi() { return 3.14159265358979323846264338327950288419716939937510; }
  static inline long double log(double x) { return std::log(x); }
  static inline long double log10(double x) { return std::log10(x); }
};
#endif

#ifdef HAVE_TEUCHOSCORE_QUADMATH

template<>
struct ScalarTraits<__float128> {
  typedef __float128 magnitudeType;
  // Unfortunately, we can't rely on a standard __float256 type.  It
  // might make sense to use qd_real, but mixing quadmath and QD might
  // cause unforeseen issues.
  typedef __float128 doublePrecision;
  typedef double halfPrecision;
  typedef __float128 coordinateType;

  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;

  static __float128 eps () {
    return FLT128_EPSILON;
  }
  static __float128 sfmin () {
    return FLT128_MIN; // ???
  }
  static __float128 base () {
    return 2;
  }
  static __float128 prec () {
    return eps () * base ();
  }
  static __float128 t () {
    return FLT128_MANT_DIG;
  }
  static __float128 rnd () {
    return 1.0;
  }
  static __float128 emin () {
    return FLT128_MIN_EXP;
  }
  static __float128 rmin () {
    return FLT128_MIN; // ??? // should be base^(emin-1)
  }
  static __float128 emax () {
    return FLT128_MAX_EXP;
  }
  static __float128 rmax () {
    return FLT128_MAX; // ??? // should be (base^emax)*(1-eps)
  }
  static magnitudeType magnitude (const __float128& x) {
    return fabsq (x);
  }
  static __float128 zero () {
    return 0.0;
  }
  static __float128 one () {
    return 1.0;
  }
  static __float128 conjugate (const __float128& x) {
    return x;
  }
  static __float128 real (const __float128& x) {
    return x;
  }
  static __float128 imag (const __float128& /* x */) {
    return 0.0;
  }
  static __float128 nan () {
    return strtoflt128 ("NAN()", NULL); // ???
  }
  static bool isnaninf (const __float128& x) {
    return isinfq (x) || isnanq (x);
  }
  static inline void seedrandom (unsigned int s) {
    std::srand (s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random ();
#endif
  }
  static __float128 random () {
    // Half the smallest normalized double, is the scaling factor of
    // the lower-order term in the double-double representation.
    const __float128 scalingFactor =
      static_cast<__float128> (std::numeric_limits<double>::min ()) /
      static_cast<__float128> (2.0);
    const __float128 higherOrderTerm =
      static_cast<__float128> (ScalarTraits<double>::random ());
    const __float128 lowerOrderTerm =
      static_cast<__float128> (ScalarTraits<double>::random ()) *
      scalingFactor;
    return higherOrderTerm + lowerOrderTerm;
  }
  static std::string name () {
    return "__float128";
  }
  static __float128 squareroot (const __float128& x) {
    return sqrtq (x);
  }
  static __float128 pow (const __float128& x, const __float128& y) {
    return powq (x, y);
  }
  static __float128 pi() { return 3.14159265358979323846; }
  static __float128 log (const __float128& x) {
    return logq (x);
  }
  static __float128 log10 (const __float128& x) {
    return log10q (x);
  }
};
#endif // HAVE_TEUCHOSCORE_QUADMATH



#ifdef HAVE_TEUCHOS_QD

bool operator&&(const dd_real &a, const dd_real &b);
bool operator&&(const qd_real &a, const qd_real &b);

template<>
struct ScalarTraits<dd_real>
{
  typedef dd_real magnitudeType;
  typedef double halfPrecision;
  typedef qd_real doublePrecision;
  typedef dd_real coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;
  static inline dd_real eps()   { return std::numeric_limits<dd_real>::epsilon(); }
  static inline dd_real sfmin() { return std::numeric_limits<dd_real>::min(); }
  static inline dd_real base()  { return std::numeric_limits<dd_real>::radix; }
  static inline dd_real prec()  { return eps()*base(); }
  static inline dd_real t()     { return std::numeric_limits<dd_real>::digits; }
  static inline dd_real rnd()   { return ( std::numeric_limits<dd_real>::round_style == std::round_to_nearest ? dd_real(1.0) : dd_real(0.0) ); }
  static inline dd_real emin()  { return std::numeric_limits<dd_real>::min_exponent; }
  static inline dd_real rmin()  { return std::numeric_limits<dd_real>::min(); }
  static inline dd_real emax()  { return std::numeric_limits<dd_real>::max_exponent; }
  static inline dd_real rmax()  { return std::numeric_limits<dd_real>::max(); }
  static inline magnitudeType magnitude(dd_real a)
  {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        a, "Error, the input value to magnitude(...) a = " << a << " can not be NaN!" );
#endif
      return ::abs(a);
  }
  static inline dd_real zero()  { return dd_real(0.0); }
  static inline dd_real one()   { return dd_real(1.0); }
  static inline dd_real conjugate(dd_real x)   { return(x); }
  static inline dd_real real(dd_real x) { return x ; }
  static inline dd_real imag(dd_real) { return zero(); }
  static inline dd_real nan() { return dd_real::_nan; }
  static inline bool isnaninf(dd_real x) { return isnan(x) || isinf(x); }
  static inline void seedrandom(unsigned int s) {
    // ddrand() uses std::rand(), so the std::srand() is our seed
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  static inline dd_real random() { return ddrand(); }
  static inline std::string name() { return "dd_real"; }
  static inline dd_real squareroot(dd_real x)
  {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        x, "Error, the input value to squareroot(...) x = " << x << " can not be NaN!" );
#endif
      return ::sqrt(x);
  }
  static inline dd_real pow(dd_real x, dd_real y) { return ::pow(x,y); }
  static inline dd_real pi() { return 3.14159265358979323846; }
  // dd_real puts its transcendental functions in the global namespace.
  static inline dd_real log(dd_real x) { return ::log(x); }
  static inline dd_real log10(dd_real x) { return ::log10(x); }
};


template<>
struct ScalarTraits<qd_real>
{
  typedef qd_real magnitudeType;
  typedef dd_real halfPrecision;
  typedef qd_real doublePrecision;
  typedef qd_real coordinateType;
  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;
  static inline qd_real eps()   { return std::numeric_limits<qd_real>::epsilon(); }
  static inline qd_real sfmin() { return std::numeric_limits<qd_real>::min(); }
  static inline qd_real base()  { return std::numeric_limits<qd_real>::radix; }
  static inline qd_real prec()  { return eps()*base(); }
  static inline qd_real t()     { return std::numeric_limits<qd_real>::digits; }
  static inline qd_real rnd()   { return ( std::numeric_limits<qd_real>::round_style == std::round_to_nearest ? qd_real(1.0) : qd_real(0.0) ); }
  static inline qd_real emin()  { return std::numeric_limits<qd_real>::min_exponent; }
  static inline qd_real rmin()  { return std::numeric_limits<qd_real>::min(); }
  static inline qd_real emax()  { return std::numeric_limits<qd_real>::max_exponent; }
  static inline qd_real rmax()  { return std::numeric_limits<qd_real>::max(); }
  static inline magnitudeType magnitude(qd_real a)
  {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        a, "Error, the input value to magnitude(...) a = " << a << " can not be NaN!" );
#endif
      return ::abs(a);
  }
  static inline qd_real zero()  { return qd_real(0.0); }
  static inline qd_real one()   { return qd_real(1.0); }
  static inline qd_real conjugate(qd_real x)   { return(x); }
  static inline qd_real real(qd_real x) { return x ; }
  static inline qd_real imag(qd_real) { return zero(); }
  static inline qd_real nan() { return qd_real::_nan; }
  static inline bool isnaninf(qd_real x) { return isnan(x) || isinf(x); }
  static inline void seedrandom(unsigned int s) {
    // qdrand() uses std::rand(), so the std::srand() is our seed
    std::srand(s);
#ifdef __APPLE__
    // throw away first random number to address bug 3655
    // http://software.sandia.gov/bugzilla/show_bug.cgi?id=3655
    random();
#endif
  }
  static inline qd_real random() { return qdrand(); }
  static inline std::string name() { return "qd_real"; }
  static inline qd_real squareroot(qd_real x)
  {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        x, "Error, the input value to squareroot(...) x = " << x << " can not be NaN!" );
#endif
      return ::sqrt(x);
  }
  static inline qd_real pow(qd_real x, qd_real y) { return ::pow(x,y); }
  static inline qd_real pi() { return 3.14159265358979323846; }
  // qd_real puts its transcendental functions in the global namespace.
  static inline qd_real log(qd_real x) { return ::log(x); }
  static inline qd_real log10(qd_real x) { return ::log10(x); }
};


#endif  // HAVE_TEUCHOS_QD


#ifdef HAVE_TEUCHOS_GNU_MP


extern gmp_randclass gmp_rng;


/* Regarding halfPrecision, doublePrecision and mpf_class:
   Because the precision of an mpf_class float is not determined by the data type,
   there is no way to fill the typedefs for this object.

   Instead, we could create new data classes (e.g., Teuchos::MPF128, Teuchos::MPF256) for
   commonly used levels of precision, and fill out ScalarTraits for these. This would allow us
   to typedef the promotions and demotions in the appropriate way. These classes would serve to
   wrap an mpf_class object, calling the constructor for the appropriate precision, exposing the
   arithmetic routines but hiding the precision-altering routines.

   Alternatively (perhaps, preferably), would could create a single class templated on the precision (e.g., Teuchos::MPF<N>).
   Then we have a single (partially-specialized) implementation of ScalarTraits. This class, as above, must expose all of the
   operations expected of a scalar type; however, most of these can be trivially stolen from the gmpcxx.h header file

   CGB/RAB, 01/05/2009
*/
template<>
struct ScalarTraits<mpf_class>
{
  typedef mpf_class magnitudeType;
  typedef mpf_class halfPrecision;
  typedef mpf_class doublePrecision;
  typedef mpf_class coordinateType;
  static const bool isComplex = false;
  static const bool hasMachineParameters = false;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static magnitudeType magnitude(mpf_class a) { return std::abs(a); }
  static inline mpf_class zero() { mpf_class zero = 0.0; return zero; }
  static inline mpf_class one() { mpf_class one = 1.0; return one; }
  static inline mpf_class conjugate(mpf_class x) { return x; }
  static inline mpf_class real(mpf_class x) { return(x); }
  static inline mpf_class imag(mpf_class x) { return(0); }
  static inline bool isnaninf(mpf_class x) { return false; } // mpf_class currently can't handle nan or inf!
  static inline void seedrandom(unsigned int s) {
    unsigned long int seedVal = static_cast<unsigned long int>(s);
    gmp_rng.seed( seedVal );
  }
  static inline mpf_class random() {
    return gmp_rng.get_f();
  }
  static inline std::string name() { return "mpf_class"; }
  static inline mpf_class squareroot(mpf_class x) { return std::sqrt(x); }
  static inline mpf_class pow(mpf_class x, mpf_class y) { return pow(x,y); }
  // Todo: RAB: 2004/05/28: Add nan() and isnaninf() functions when needed!
};

#endif  // HAVE_TEUCHOS_GNU_MP

#ifdef HAVE_TEUCHOS_ARPREC

/* See discussion above for mpf_class, regarding halfPrecision and doublePrecision. Something similar will need to be done
   for ARPREC. */
template<>
struct ScalarTraits<mp_real>
{
  typedef mp_real magnitudeType;
  typedef double halfPrecision;
  typedef mp_real doublePrecision;
  typedef mp_real coordinateType;
  static const bool isComplex = false;
  static const bool isComparable = true;
  static const bool isOrdinal = false;
  static const bool hasMachineParameters = false;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static magnitudeType magnitude(mp_real a) { return abs(a); }
  static inline mp_real zero() { mp_real zero = 0.0; return zero; }
  static inline mp_real one() { mp_real one = 1.0; return one; }
  static inline mp_real conjugate(mp_real x) { return x; }
  static inline mp_real real(mp_real x) { return(x); }
  static inline mp_real imag(mp_real x) { return zero(); }
  static inline bool isnaninf(mp_real x) { return false; } // ToDo: Change this?
  static inline void seedrandom(unsigned int s) {
    long int seedVal = static_cast<long int>(s);
    srand48(seedVal);
  }
  static inline mp_real random() { return mp_rand(); }
  static inline std::string name() { return "mp_real"; }
  static inline mp_real squareroot(mp_real x) { return sqrt(x); }
  static inline mp_real pow(mp_real x, mp_real y) { return pow(x,y); }
  static inline mp_real pi() { return 3.14159265358979323846; }
  // Todo: RAB: 2004/05/28: Add nan() and isnaninf() functions when needed!
};


#endif // HAVE_TEUCHOS_ARPREC


#ifdef HAVE_TEUCHOS_COMPLEX


// Partial specialization for std::complex numbers templated on real type T
template<class T>
struct ScalarTraits<
  std::complex<T>
>
{
  typedef std::complex<T>  ComplexT;
  typedef std::complex<typename ScalarTraits<T>::halfPrecision> halfPrecision;
  typedef std::complex<typename ScalarTraits<T>::doublePrecision> doublePrecision;
  typedef typename ScalarTraits<T>::magnitudeType magnitudeType;
  typedef typename ScalarTraits<T>::coordinateType coordinateType;
  static const bool isComplex = true;
  static const bool isOrdinal = ScalarTraits<T>::isOrdinal;
  static const bool isComparable = false;
  static const bool hasMachineParameters = true;
  static inline magnitudeType eps()          { return ScalarTraits<magnitudeType>::eps(); }
  static inline magnitudeType sfmin()        { return ScalarTraits<magnitudeType>::sfmin(); }
  static inline magnitudeType base()         { return ScalarTraits<magnitudeType>::base(); }
  static inline magnitudeType prec()         { return ScalarTraits<magnitudeType>::prec(); }
  static inline magnitudeType t()            { return ScalarTraits<magnitudeType>::t(); }
  static inline magnitudeType rnd()          { return ScalarTraits<magnitudeType>::rnd(); }
  static inline magnitudeType emin()         { return ScalarTraits<magnitudeType>::emin(); }
  static inline magnitudeType rmin()         { return ScalarTraits<magnitudeType>::rmin(); }
  static inline magnitudeType emax()         { return ScalarTraits<magnitudeType>::emax(); }
  static inline magnitudeType rmax()         { return ScalarTraits<magnitudeType>::rmax(); }
  static magnitudeType magnitude(ComplexT a)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        a, "Error, the input value to magnitude(...) a = " << a << " can not be NaN!" );
#endif
      return std::abs(a);
    }
  static inline ComplexT zero()              { return ComplexT(ScalarTraits<magnitudeType>::zero(),ScalarTraits<magnitudeType>::zero()); }
  static inline ComplexT one()               { return ComplexT(ScalarTraits<magnitudeType>::one(),ScalarTraits<magnitudeType>::zero()); }
  static inline ComplexT conjugate(ComplexT a){ return ComplexT(a.real(),-a.imag()); }
  static inline magnitudeType real(ComplexT a) { return a.real(); }
  static inline magnitudeType imag(ComplexT a) { return a.imag(); }
  static inline ComplexT nan()               { return ComplexT(ScalarTraits<magnitudeType>::nan(),ScalarTraits<magnitudeType>::nan()); }
  static inline bool isnaninf(ComplexT x)    { return ScalarTraits<magnitudeType>::isnaninf(x.real()) || ScalarTraits<magnitudeType>::isnaninf(x.imag()); }
  static inline void seedrandom(unsigned int s) { ScalarTraits<magnitudeType>::seedrandom(s); }
  static inline ComplexT random()
    {
      const T rnd1 = ScalarTraits<magnitudeType>::random();
      const T rnd2 = ScalarTraits<magnitudeType>::random();
      return ComplexT(rnd1,rnd2);
    }
  static inline std::string name() { return std::string("std::complex<")+std::string(ScalarTraits<magnitudeType>::name())+std::string(">"); }
  // This will only return one of the square roots of x, the other can be obtained by taking its conjugate
  static inline ComplexT squareroot(ComplexT x)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        x, "Error, the input value to squareroot(...) x = " << x << " can not be NaN!" );
#endif
      typedef ScalarTraits<magnitudeType>  STMT;
      const T r  = x.real(), i = x.imag(), zero = STMT::zero(), two = 2.0;
      const T a  = STMT::squareroot((r*r)+(i*i));
      const T nr = STMT::squareroot((a+r)/two);
      const T ni = ( i == zero ? zero : STMT::squareroot((a-r)/two) );
      return ComplexT(nr,ni);
    }
    // 2010/03/19: rabartl: Above, I had to add the check for i == zero
    // to avoid a returned NaN on the Intel 10.1 compiler.  For some
    // reason, having these two squareroot calls in a row produce a NaN
    // in an optimized build with this compiler.  Amazingly, when I put
    // in print statements (i.e. std::cout << ...) the NaN went away and
    // the second squareroot((a-r)/two) returned zero correctly.  I
    // guess that calling the output routine flushed the registers or
    // something and restarted the squareroot rountine for this compiler
    // or something similar.  Actually, due to roundoff, it is possible that a-r
    // might be slightly less than zero (i.e. -1e-16) and that would cause
    // a possbile NaN return.  THe above if test is the right thing to do
    // I think and is very cheap.
  static inline ComplexT pow(ComplexT x, ComplexT y) { return pow(x,y); }
  static inline ComplexT pi() { return ScalarTraits<T>::pi(); }
};
#endif //  HAVE_TEUCHOS_COMPLEX

#ifdef HAVE_TEUCHOSCORE_KOKKOS
// Partial specialization for Kokkos::complex<T>
template<class T>
struct ScalarTraits<
  Kokkos::complex<T>
>
{
  typedef Kokkos::complex<T>  ComplexT;
  typedef Kokkos::complex<typename ScalarTraits<T>::halfPrecision> halfPrecision;
  typedef Kokkos::complex<typename ScalarTraits<T>::doublePrecision> doublePrecision;
  typedef typename ScalarTraits<T>::magnitudeType magnitudeType;
  typedef typename ScalarTraits<T>::coordinateType coordinateType;
  static const bool isComplex = true;
  static const bool isOrdinal = ScalarTraits<T>::isOrdinal;
  static const bool isComparable = false;
  static const bool hasMachineParameters = true;
  static inline magnitudeType eps()          { return ScalarTraits<magnitudeType>::eps(); }
  static inline magnitudeType sfmin()        { return ScalarTraits<magnitudeType>::sfmin(); }
  static inline magnitudeType base()         { return ScalarTraits<magnitudeType>::base(); }
  static inline magnitudeType prec()         { return ScalarTraits<magnitudeType>::prec(); }
  static inline magnitudeType t()            { return ScalarTraits<magnitudeType>::t(); }
  static inline magnitudeType rnd()          { return ScalarTraits<magnitudeType>::rnd(); }
  static inline magnitudeType emin()         { return ScalarTraits<magnitudeType>::emin(); }
  static inline magnitudeType rmin()         { return ScalarTraits<magnitudeType>::rmin(); }
  static inline magnitudeType emax()         { return ScalarTraits<magnitudeType>::emax(); }
  static inline magnitudeType rmax()         { return ScalarTraits<magnitudeType>::rmax(); }
  static magnitudeType magnitude(ComplexT a)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        a, "Error, the input value to magnitude(...) a = " << a << " can not be NaN!" );
#endif
      return std::abs(std::complex<T>(a));
    }
  static inline ComplexT zero()              { return ComplexT(ScalarTraits<magnitudeType>::zero(),ScalarTraits<magnitudeType>::zero()); }
  static inline ComplexT one()               { return ComplexT(ScalarTraits<magnitudeType>::one(),ScalarTraits<magnitudeType>::zero()); }
  static inline ComplexT conjugate(ComplexT a){ return ComplexT(a.real(),-a.imag()); }
  static inline magnitudeType real(ComplexT a) { return a.real(); }
  static inline magnitudeType imag(ComplexT a) { return a.imag(); }
  static inline ComplexT nan()               { return ComplexT(ScalarTraits<magnitudeType>::nan(),ScalarTraits<magnitudeType>::nan()); }
  static inline bool isnaninf(ComplexT x)    { return ScalarTraits<magnitudeType>::isnaninf(x.real()) || ScalarTraits<magnitudeType>::isnaninf(x.imag()); }
  static inline void seedrandom(unsigned int s) { ScalarTraits<magnitudeType>::seedrandom(s); }
  static inline ComplexT random()
    {
      const T rnd1 = ScalarTraits<magnitudeType>::random();
      const T rnd2 = ScalarTraits<magnitudeType>::random();
      return ComplexT(rnd1,rnd2);
    }
  static inline std::string name() { return std::string("Kokkos::complex<")+std::string(ScalarTraits<magnitudeType>::name())+std::string(">"); }
  // This will only return one of the square roots of x, the other can be obtained by taking its conjugate
  static inline ComplexT squareroot(ComplexT x)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        x, "Error, the input value to squareroot(...) x = " << x << " can not be NaN!" );
#endif
      typedef ScalarTraits<magnitudeType>  STMT;
      const T r  = x.real(), i = x.imag(), zero = STMT::zero(), two = 2.0;
      const T a  = STMT::squareroot((r*r)+(i*i));
      const T nr = STMT::squareroot((a+r)/two);
      const T ni = ( i == zero ? zero : STMT::squareroot((a-r)/two) );
      return ComplexT(nr,ni);
    }
  static inline ComplexT pow(ComplexT x, ComplexT y) { return pow(std::complex<T>(x), std::complex<T>(y)); }
  static inline ComplexT pi() { return ScalarTraits<T>::pi(); }
};
#endif // HAVE_TEUCHOSCORE_KOKKOS

#endif // DOXYGEN_SHOULD_SKIP_THIS

} // Teuchos namespace

#endif // _TEUCHOS_SCALARTRAITS_HPP_
