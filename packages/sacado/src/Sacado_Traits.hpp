// $Id$ 
// $Source$ 
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
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
  template <> struct Promote< type1, type2 > {		   \
    typedef type3 type;			                   \
  };							   \
  template <> struct Promote< type2, type1 > {		   \
    typedef type3 type;				           \
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

  //! Base template specification for %ScalarType
  /*!
   * The %ScalarType classes provide a mechanism for computing the 
   * base underlying type of nested AD classes
   */
  template <typename T> struct ScalarType {};

  //! Base template specification for %ValueType
  /*!
   * The %ValueType classes provide a mechanism for computing the 
   * the type stored in AD classes
   */
  template <typename T> struct ValueType {};

  //! Base template specification for %ScalarValueType
  /*!
   * The %ScalarValueType classes provide a mechanism for computing the 
   * base underlying type of the value type of nested AD classes
   */
  template <typename T> struct ScalarValueType {};

  //! Base template specification for %IsADType
  /*!
   * The %IsADType classes provide a mechanism for computing the 
   * determining whether a type is an AD type
   */
  template <typename T> struct IsADType {};

  //! Base template specification for %IsScalarType
  /*!
   * The %IsScalarType classes provide a mechanism for computing the 
   * determining whether a type is a scalar type (float, double, etc...)
   */
  template <typename T> struct IsScalarType {};

  //! Base template specification for %Value
  /*!
   * The %Value functor returns the value of an AD type
   */
  template <typename T> struct Value {};

  //! Base template specification for marking constants
  template <typename T> struct MarkConstant {
    static void eval(T& x) {}
  };

  //! Specialization of above classes to builtin types
#define SACADO_BUILTIN_SPECIALIZATION(t)                  \
  template <> struct ScalarType< t > {		          \
    typedef t type;				          \
  };                                                      \
  template <> struct ValueType< t > {		          \
    typedef t type;				          \
  };                                                      \
  template <> struct ScalarValueType< t > {		  \
    typedef t type;				          \
  };                                                      \
  template <> struct IsADType< t > {		          \
    static const bool value = false;	       		  \
  };                                                      \
  template <> struct IsScalarType< t > {	          \
    static const bool value = true;	       		  \
  };                                                      \
  template <> struct Value< t > {		          \
    static const t& eval(const t& x) { return x; }        \
  };

  SACADO_BUILTIN_SPECIALIZATION(char)
  SACADO_BUILTIN_SPECIALIZATION(float)
  SACADO_BUILTIN_SPECIALIZATION(double)
  SACADO_BUILTIN_SPECIALIZATION(int)
  SACADO_BUILTIN_SPECIALIZATION(unsigned int)
  SACADO_BUILTIN_SPECIALIZATION(long)
  SACADO_BUILTIN_SPECIALIZATION(unsigned long)
  SACADO_BUILTIN_SPECIALIZATION(bool)
#ifdef HAVE_SACADO_COMPLEX
  SACADO_BUILTIN_SPECIALIZATION(std::complex<double>)
  SACADO_BUILTIN_SPECIALIZATION(std::complex<float>)
#endif

#undef SACADO_BUILTIN_SPECIALIZATION

} // namespace Sacado

// Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#endif

#endif // SACADO_TRAITS_HPP
