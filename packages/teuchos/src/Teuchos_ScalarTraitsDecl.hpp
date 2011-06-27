// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef _TEUCHOS_SCALARTRAITS_DECL_HPP_
#define _TEUCHOS_SCALARTRAITS_DECL_HPP_

/*! \file Teuchos_ScalarTraitsDecl.hpp
    \brief Declaration and default implementation for basic traits for the scalar field type.
*/
 
#include "Teuchos_ConfigDefs.hpp"


namespace Teuchos {


template <typename T>
struct UndefinedScalarTraits
{
  //! This function should not compile if there is an attempt to instantiate!
  static inline T notDefined() { return T::this_type_is_missing_a_specialization(); }
};


/* This is the default structure used by ScalarTraits<T> to produce a compile time
	error when the specialization does not exist for type <tt>T</tt>.
*/


/*! \brief This structure defines some basic traits for a scalar field type.
 *
 * Scalar traits are an essential part of templated codes.  This structure offers
 * the basic traits of the templated scalar type, like defining zero and one,
 * and basic functions on the templated scalar type, like performing a square root.
 * 
 * The functions in the templated base unspecialized struct are designed not to
 * compile (giving a nice compile-time error message) and therefore specializations
 * must be written for Scalar types actually used.
 * 
 * \note <ol>
 * 
 * <li> The default defined specializations are provided for \c int, \c float, and \c double.
 * 
 * <li> If Teuchos is configured with </tt>Teuchos_ENABLE_COMPLEX=ON</tt> then
 * ScalarTraits also has a partial specialization for all
 * <tt>std::complex</tt> numbers of the form <tt>std::complex<T></tt>.
 * 
 * </ol>
*/
template <typename T>
struct ScalarTraits
{
  //! Mandatory typedef for result of magnitude
  typedef T magnitudeType;
  //! Typedef for half precision
  typedef T halfPrecision;
  //! Typedef for double precision
  typedef T doublePrecision;
  //! Determines if scalar type is std::complex
  static const bool isComplex = false;
  //! Determines if scalar type is an ordinal type
  static const bool isOrdinal = false;
  //! Determines if scalar type supports relational operators such as <, >, <=, >=.
  static const bool isComparable = false;
  /** \brief Determines if scalar type have machine-specific parameters
   * (i.e. eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(),
   * rmax() are supported).
   */
  static const bool hasMachineParameters = false;
  //! Returns relative machine precision.
  static inline magnitudeType eps()   { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns safe minimum (sfmin), such that 1/sfmin does not overflow.
  static inline magnitudeType sfmin() { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns the base of the machine.
  static inline magnitudeType base()  { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns \c eps*base.
  static inline magnitudeType prec()  { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns the number of (base) digits in the mantissa.
  static inline magnitudeType t()     { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns 1.0 when rounding occurs in addition, 0.0 otherwise
  static inline magnitudeType rnd()   { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns the minimum exponent before (gradual) underflow.
  static inline magnitudeType emin()  { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns the underflow threshold - \c base^(emin-1)
  static inline magnitudeType rmin()  { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns the largest exponent before overflow.
  static inline magnitudeType emax()  { return UndefinedScalarTraits<T>::notDefined(); }
  //! Overflow theshold - \c (base^emax)*(1-eps)
  static inline magnitudeType rmax()  { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns the magnitudeType of the scalar type \c a.
  static inline magnitudeType magnitude(T a) { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns representation of zero for this scalar type.
  static inline T zero()                     { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns representation of one for this scalar type.
  static inline T one()                      { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns the real part of the scalar type \c a.
  static inline magnitudeType real(T a) { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns the imaginary part of the scalar type \c a.
  static inline magnitudeType imag(T a) { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns the conjugate of the scalar type \c a.
  static inline T conjugate(T a) { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns a number that represents NaN.
  static inline T nan()                      { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns <tt>true</tt> if <tt>x</tt> is NaN or Inf.
  static inline bool isnaninf(const T& x)     { return UndefinedScalarTraits<T>::notDefined(); }
  //! Seed the random number generator returned by <tt>random()</tt>.
  static inline void seedrandom(unsigned int s) { int i; T t = &i; }
  //! Returns a random number (between -one() and +one()) of this scalar type.
  static inline T random()                   { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns the name of this scalar type.
  static inline std::string name()           { (void)UndefinedScalarTraits<T>::notDefined(); return 0; }
  //! Returns a number of magnitudeType that is the square root of this scalar type \c x. 
  static inline T squareroot(T x) { return UndefinedScalarTraits<T>::notDefined(); }
  //! Returns the result of raising one scalar \c x to the power \c y.
  static inline T pow(T x, T y) { return UndefinedScalarTraits<T>::notDefined(); }
};

  
} // Teuchos namespace


#endif // _TEUCHOS_SCALARTRAITS_DECL_HPP_
