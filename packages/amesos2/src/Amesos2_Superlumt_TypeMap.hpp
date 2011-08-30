// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//
// @HEADER

/**
   \file   Amesos2_Superlumt_TypeMap.hpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Mon May 31 23:12:32 2010

   \brief Provides definition of SuperLU_MT types as well as
          conversions and type traits.

*/

#ifndef AMESOS2_SUPERLUMT_TYPEMAP_HPP
#define AMESOS2_SUPERLUMT_TYPEMAP_HPP

#include <functional>

#include <Teuchos_as.hpp>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <Teuchos_SerializationTraits.hpp>
#endif

#include "Amesos2_TypeMap.hpp"

namespace SLUMT {

typedef int int_t;

extern "C" {

#undef __SUPERLU_SUPERMATRIX
#include "supermatrix.h"	// for Dtype_t declaration

#ifdef HAVE_TEUCHOS_COMPLEX
  
namespace C {
#undef __SUPERLU_SCOMPLEX
#undef SCOMPLEX_INCLUDE
#include "slu_scomplex.h"     // single-precision complex data type definitions
}

namespace Z {
#undef __SUPERLU_DCOMPLEX
#undef DCOMPLEX_INCLUDE
#include "slu_dcomplex.h"     // double-precision complex data type definitions
}

#endif	// HAVE_TEUCHOS_COMPLEX

} // end extern "C"

  // Declare and specialize a std::binary_funtion class for
  // multiplication of SLUMT types
  template <typename slu_scalar_t, typename slu_mag_t>
  struct slu_mt_mult {};

  // This specialization handles the generic case were the scalar and
  // magnitude types are double or float.
  template <typename T>
  struct slu_mt_mult<T,T> : std::multiplies<T> {};

#ifdef HAVE_TEUCHOS_COMPLEX
  
  // For namespace/macro reasons, we prefix our variables with amesos_*
  template <>
  struct slu_mt_mult<C::complex,float>
    : std::binary_function<C::complex,float,C::complex> {
    C::complex operator()(C::complex amesos_c, float amesos_f) {
      C::complex amesos_cr;
      cs_mult(&amesos_cr, &amesos_c, amesos_f);	// cs_mult is a macro, so no namespacing
      return( amesos_cr );
    }
  };

  template <>
  struct slu_mt_mult<C::complex,C::complex>
    : std::binary_function<C::complex,C::complex,C::complex> {
    C::complex operator()(C::complex amesos_c1, C::complex amesos_c2) {
      C::complex amesos_cr;
      cc_mult(&amesos_cr, &amesos_c1, &amesos_c2); // cc_mult is a macro, so no namespacing
      return( amesos_cr );
    }
  };
    
  template <>
  struct slu_mt_mult<Z::doublecomplex,double>
    : std::binary_function<Z::doublecomplex,double,Z::doublecomplex> {
    Z::doublecomplex operator()(Z::doublecomplex amesos_z, double amesos_d) {
      Z::doublecomplex amesos_zr;
      zd_mult(&amesos_zr, &amesos_z, amesos_d);	// zd_mult is a macro, so no namespacing
      return( amesos_zr );
    }
  };

  template <>
  struct slu_mt_mult<Z::doublecomplex,Z::doublecomplex>
    : std::binary_function<Z::doublecomplex,Z::doublecomplex,Z::doublecomplex> {
    Z::doublecomplex operator()(Z::doublecomplex amesos_z1, Z::doublecomplex amesos_z2) {
      Z::doublecomplex amesos_zr;
      zz_mult(&amesos_zr, &amesos_z1, &amesos_z2);    // zz_mult is a macro, so no namespacing
      return( amesos_zr );
    }
  };

#endif	// HAVE_TEUCHOS_COMPLEX
} // end namespace SLUMt
#ifdef HAVE_TEUCHOS_COMPLEX

/* ==================== Conversion ==================== */
namespace Teuchos {

/**
 * \defgroup slu_conversion Conversion definitions for SLU types.
 *
 * Define specializations of Teuchos::as<> for the SLU types.
 *
 * These specializations are meant to work with any complex data type that
 * implements the same interface as the STL complex type.
 *
 * @{
 */
template <typename TypeFrom>
class ValueTypeConversionTraits<SLUMT::C::complex, TypeFrom>
{
public:
  static SLUMT::C::complex convert( const TypeFrom t )
    {
      SLUMT::C::complex ret;
      ret.r = Teuchos::as<float>(t.real());
      ret.i = Teuchos::as<float>(t.imag());
      return( ret );
    }

  static SLUMT::C::complex safeConvert( const TypeFrom t )
    {
      SLUMT::C::complex ret;
      ret.r = Teuchos::as<float>(t.real());
      ret.i = Teuchos::as<float>(t.imag());
      return( ret );
    }
};


template <typename TypeFrom>
class ValueTypeConversionTraits<SLUMT::Z::doublecomplex, TypeFrom>
{
public:
  static SLUMT::Z::doublecomplex convert( const TypeFrom t )
    {
      SLUMT::Z::doublecomplex ret;
      ret.r = Teuchos::as<double>(t.real());
      ret.i = Teuchos::as<double>(t.imag());
      return( ret );
    }

  static SLUMT::Z::doublecomplex safeConvert( const TypeFrom t )
    {
      SLUMT::Z::doublecomplex ret;
      ret.r = Teuchos::as<double>(t.real());
      ret.i = Teuchos::as<double>(t.imag());
      return( ret );
    }
};


// Also convert from SLU types
template <typename TypeTo>
class ValueTypeConversionTraits<TypeTo, SLUMT::C::complex>
{
public:
  static TypeTo convert( const SLUMT::C::complex t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }

  // No special checks for safe Convert
  static TypeTo safeConvert( const SLUMT::C::complex t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }
};


template <typename TypeTo>
class ValueTypeConversionTraits<TypeTo, SLUMT::Z::doublecomplex>
{
public:
  static TypeTo convert( const SLUMT::Z::doublecomplex t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }

  // No special checks for safe Convert
  static TypeTo safeConvert( const SLUMT::Z::doublecomplex t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }
};

template <typename Ordinal>
class SerializationTraits<Ordinal,SLUMT::C::complex>
  : public DirectSerializationTraits<Ordinal,SLUMT::C::complex>
{};

template <typename Ordinal>
class SerializationTraits<Ordinal,SLUMT::Z::doublecomplex>
  : public DirectSerializationTraits<Ordinal,SLUMT::Z::doublecomplex>
{};

//@}  End Conversion group

} // end namespace Teuchos



/**
 * \defgroup slu_std_operators
 *
 * @{
 */
namespace std {
  // C++-style output functions for Superlumt complex types
  ostream& operator<<(ostream& out, const SLUMT::C::complex c);

  ostream& operator<<(ostream& out, const SLUMT::Z::doublecomplex z);
  
  //@} End std operators group
}

#endif	// HAVE_TEUCHOS_COMPLEX


namespace Amesos2 {

template <class, class> class Superlumt;

/* Specialize the Amesos2::TypeMap struct for Superlumt types
 *
 * \cond Superlumt_type_specializations 
 */
template <>
struct TypeMap<Superlumt,float>
{
  static SLUMT::Dtype_t dtype;
  typedef float type;
  typedef float magnitude_type;
};


template <>
struct TypeMap<Superlumt,double>
{
  static SLUMT::Dtype_t dtype;
  typedef double type;
  typedef double magnitude_type;
};


#ifdef HAVE_TEUCHOS_COMPLEX
template <>
struct TypeMap<Superlumt,std::complex<float> >
{
  static SLUMT::Dtype_t dtype;
  typedef SLUMT::C::complex type;
  typedef float magnitude_type;
};


template <>
struct TypeMap<Superlumt,std::complex<double> >
{
  static SLUMT::Dtype_t dtype;
  typedef SLUMT::Z::doublecomplex type;
  typedef double magnitude_type;
};


template <>
struct TypeMap<Superlumt,SLUMT::C::complex>
{
  static SLUMT::Dtype_t dtype;
  typedef SLUMT::C::complex type;
  typedef float magnitude_type;
};


template <>
struct TypeMap<Superlumt,SLUMT::Z::doublecomplex>
{
  static SLUMT::Dtype_t dtype;
  typedef SLUMT::Z::doublecomplex type;
  typedef double magnitude_type;
};

#endif	// HAVE_TEUCHOS_COMPLEX

/* \endcond Superlumt_type_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_SUPERLUMT_TYPEMAP_HPP
