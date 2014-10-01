// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
//
// The forward-mode AD classes in Sacado are a derivative work of the
// expression template classes in the Fad package by Nicolas Di Cesare.
// The following banner is included in the original Fad source code:
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses,
//         templates : new C++ techniques
//            for scientific computing
//
//********************************************************
//
//  NumericalTraits class to illustrate TRAITS
//
//********************************************************
// @HEADER

#ifndef SACADO_TRAITS_HPP
#define SACADO_TRAITS_HPP

#include "Sacado_ConfigDefs.h"
#include <string>

#ifdef HAVE_SACADO_COMPLEX
#include <complex>
#endif

namespace Sacado {

  //! Base template specification for %Promote
  /*!
   * The %Promote classes provide a mechanism for computing the
   * promoted type of a binary operation.
   */
  template <typename A, typename B> struct Promote {};

  //! Specialization of %Promote for a single type
  template <typename A> struct Promote<A,A> {
    typedef A type;
  };

  //! Specialization of %Promote to builtin types
#define SACADO_PROMOTE_SPECIALIZATION(type1,type2,type3) \
  template <> struct Promote< type1, type2 > {             \
    typedef type3 type;                                    \
  };                                                       \
  template <> struct Promote< type2, type1 > {             \
    typedef type3 type;                                    \
  };

  SACADO_PROMOTE_SPECIALIZATION(double,float,double)
  SACADO_PROMOTE_SPECIALIZATION(double,long,double)
  SACADO_PROMOTE_SPECIALIZATION(double,int,double)
  SACADO_PROMOTE_SPECIALIZATION(float,long,float)
  SACADO_PROMOTE_SPECIALIZATION(float,int,float)
#ifdef HAVE_SACADO_COMPLEX
  SACADO_PROMOTE_SPECIALIZATION(std::complex<double>,std::complex<float>,std::complex<double>)
  SACADO_PROMOTE_SPECIALIZATION(std::complex<double>,double,std::complex<double>)
  SACADO_PROMOTE_SPECIALIZATION(std::complex<double>,float,std::complex<double>)
  SACADO_PROMOTE_SPECIALIZATION(std::complex<double>,long,std::complex<double>)
  SACADO_PROMOTE_SPECIALIZATION(std::complex<double>,int,std::complex<double>)
  SACADO_PROMOTE_SPECIALIZATION(std::complex<float>,float,std::complex<float>)
  SACADO_PROMOTE_SPECIALIZATION(std::complex<float>,long,std::complex<float>)
  SACADO_PROMOTE_SPECIALIZATION(std::complex<float>,int,std::complex<float>)
#endif // HAVE_SACADO_COMPLEX

#undef SACADO_PROMOTE_SPECIALIZATION

  //
  // We define defaults for all of the traits to make Sacado easier to use.
  // The default choices are based on what appears to be the "safest" choice
  // for any scalar type.  They may not work in all cases, in which case a
  // specialization should be provided.
  //

  //! Base template specification for %ScalarType
  /*!
   * The %ScalarType classes provide a mechanism for computing the
   * base underlying type of nested AD classes
   */
  template <typename T> struct ScalarType {
    typedef T type;
  };

  //! Specialization of %ScalarType for const types
  /*!
   * This should work for most types
   */
  template <typename T> struct ScalarType<const T> {
    typedef const typename ScalarType<T>::type type;
  };

  //! Base template specification for %ValueType
  /*!
   * The %ValueType classes provide a mechanism for computing the
   * the type stored in AD classes
   */
  template <typename T> struct ValueType {
    typedef T type;
  };

  //! Specialization of %ValueType for const types
  /*!
   * This should work for most types
   */
  template <typename T> struct ValueType<const T> {
    typedef const typename ValueType<T>::type type;
  };

  //! Base template specification for %IsADType
  /*!
   * The %IsADType classes provide a mechanism for computing the
   * determining whether a type is an AD type
   */
  template <typename T> struct IsADType {
    static const bool value = false;
  };

  //! Base template specification for %IsScalarType
  /*!
   * The %IsScalarType classes provide a mechanism for computing the
   * determining whether a type is a scalar type (float, double, etc...)
   */
  template <typename T> struct IsScalarType {
    static const bool value = false;
  };

  //! Base template specification for %Value
  /*!
   * The %Value functor returns the value of an AD type.
   */
  template <typename T> struct Value {
    KOKKOS_INLINE_FUNCTION
    static const T& eval(const T& x) { return x; }
  };

  //! Base template specification for %ScalarValue
  /*!
   * The %ScalarValue functor returns the base scalar value of an AD type,
   * i.e., something that isn't an AD type.
   */
  template <typename T> struct ScalarValue {
    KOKKOS_INLINE_FUNCTION
    static const T& eval(const T& x) { return x; }
  };

  //! Base template specification for marking constants
  template <typename T> struct MarkConstant {
    KOKKOS_INLINE_FUNCTION
    static void eval(T& x) {}
  };

  //! Base template specification for string names of types
  template <typename T> struct StringName {
    KOKKOS_INLINE_FUNCTION
    static std::string eval() { return ""; }
  };

  //! Base template specification for testing equivalence
  template <typename T> struct IsEqual {
    KOKKOS_INLINE_FUNCTION
    static bool eval(const T& x, const T& y) { return x == y; }
  };

  //! Base template specification for testing whether type is statically sized
  template <typename T> struct IsStaticallySized {
    static const bool value = false;
  };

  //! Base template specification for static size
  template <typename T> struct StaticSize {
    static const unsigned value = 0;
  };

  //! Specialization of above classes to builtin types
#define SACADO_BUILTIN_SPECIALIZATION(t,NAME)             \
  template <> struct ScalarType< t > {                    \
    typedef t type;                                       \
  };                                                      \
  template <> struct ValueType< t > {                     \
    typedef t type;                                       \
  };                                                      \
  template <> struct IsADType< t > {                      \
    static const bool value = false;                      \
  };                                                      \
  template <> struct IsScalarType< t > {                  \
    static const bool value = true;                       \
  };                                                      \
  template <> struct Value< t > {                         \
    KOKKOS_INLINE_FUNCTION                                \
    static const t& eval(const t& x) { return x; }        \
  };                                                      \
  template <> struct ScalarValue< t > {                   \
    KOKKOS_INLINE_FUNCTION                                \
    static const t& eval(const t& x) { return x; }        \
  };                                                      \
  template <> struct StringName< t > {                    \
    KOKKOS_INLINE_FUNCTION                                \
    static std::string eval() { return NAME; }            \
  };                                                      \
  template <> struct IsEqual< t > {                       \
    KOKKOS_INLINE_FUNCTION                                \
    static bool eval(const t& x, const t& y) {            \
      return x == y; }                                    \
  };                                                      \
  template <> struct IsStaticallySized< t > {             \
    static const bool value = true;                       \
  };

  SACADO_BUILTIN_SPECIALIZATION(char,"char")
  SACADO_BUILTIN_SPECIALIZATION(float,"float")
  SACADO_BUILTIN_SPECIALIZATION(double,"double")
  SACADO_BUILTIN_SPECIALIZATION(int,"int")
  SACADO_BUILTIN_SPECIALIZATION(unsigned int,"unsigned int")
  SACADO_BUILTIN_SPECIALIZATION(long,"long")
  SACADO_BUILTIN_SPECIALIZATION(unsigned long,"unsigned long")
  SACADO_BUILTIN_SPECIALIZATION(bool,"bool")
#ifdef HAVE_SACADO_COMPLEX
  SACADO_BUILTIN_SPECIALIZATION(std::complex<double>,"std::complex<double>")
  SACADO_BUILTIN_SPECIALIZATION(std::complex<float>,"std::complex<float>")
#endif

#undef SACADO_BUILTIN_SPECIALIZATION

} // namespace Sacado

#endif // SACADO_TRAITS_HPP
