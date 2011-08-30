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
   \file   Amesos2_Superludist_TypeMap.hpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Tue Jun 21 13:37:41 MDT 2011

   \brief Provides definition of SuperLU_DIST types as well as
          conversions and type traits.

          SuperLU_DIST provides definitions for real and complex
          double data-types.

*/

#ifndef AMESOS2_SUPERLUDIST_TYPEMAP_HPP
#define AMESOS2_SUPERLUDIST_TYPEMAP_HPP

#include <functional>

#include <Teuchos_as.hpp>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <Teuchos_SerializationTraits.hpp>
#endif

#include "Amesos2_TypeMap.hpp"

namespace SLUD {

extern "C" {

  // undefine compiler guard in case we also have the sequential
  // SuperLU enabled
#undef __SUPERLU_SUPERMATRIX
#include "superlu_defs.h"

  namespace D {
#include "superlu_ddefs.h"	// double-precision real definitions
  }

#ifdef HAVE_TEUCHOS_COMPLEX
  namespace Z {
#include "superlu_zdefs.h"     // double-precision complex definitions
  }
#endif  // HAVE_TEUCHOS_COMPLEX

} // end extern "C"
#ifdef HAVE_TEUCHOS_COMPLEX

  // Declare and specialize a std::binary_funtion class for
  // multiplication of SLUD types
  template <typename slu_scalar_t, typename slu_mag_t>
  struct slu_mt_mult {};

  // This specialization handles the generic case were the scalar and
  // magnitude types are double or float.
  template <typename T>
  struct slu_mt_mult<T,T> : std::multiplies<T> {};

  // For namespace/macro reasons, we prefix our variables with amesos_*
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
} // end namespace SLUD
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
class ValueTypeConversionTraits<SLUD::Z::doublecomplex, TypeFrom>
{
public:
  static SLUD::Z::doublecomplex convert( const TypeFrom t )
    {
      SLUD::Z::doublecomplex ret;
      ret.r = Teuchos::as<double>(t.real());
      ret.i = Teuchos::as<double>(t.imag());
      return( ret );
    }

  static SLUD::Z::doublecomplex safeConvert( const TypeFrom t )
    {
      SLUD::Z::doublecomplex ret;
      ret.r = Teuchos::as<double>(t.real());
      ret.i = Teuchos::as<double>(t.imag());
      return( ret );
    }
};


// Also convert from SLU types
template <typename TypeTo>
class ValueTypeConversionTraits<TypeTo, SLUD::Z::doublecomplex>
{
public:
  static TypeTo convert( const SLUD::Z::doublecomplex t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }

  // No special checks for safe Convert
  static TypeTo safeConvert( const SLUD::Z::doublecomplex t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }
};

template <typename Ordinal>
class SerializationTraits<Ordinal,SLUD::Z::doublecomplex>
  : public DirectSerializationTraits<Ordinal,SLUD::Z::doublecomplex>
{};

//@}  End Conversion group

} // end namespace Teuchos



/**
 * \defgroup slu_std_operators
 *
 * @{
 */
namespace std {
  // C++-style output functions for Superludist complex types
  ostream& operator<<(ostream& out, const SLUD::Z::doublecomplex z);
  
  //@} End std operators group
}
#endif	// HAVE_TEUCHOS_COMPLEX



namespace Amesos2 {

template <class, class> class Superludist;

/* Specialize the Amesos2::TypeMap struct for SuperLU_DIST types
 *
 * \cond Superludist_type_specializations 
 */
template <>
struct TypeMap<Superludist,double>
{
  static const SLUD::Dtype_t dtype = SLUD::SLU_D;
  typedef double type;
  typedef double magnitude_type;
  typedef SLUD::D::LUstruct_t LUstruct_t;
  typedef SLUD::D::SOLVEstruct_t SOLVEstruct_t;
};

#ifdef HAVE_TEUCHOS_COMPLEX
template <>
struct TypeMap<Superludist,std::complex<double> >
{
  static const SLUD::Dtype_t dtype = SLUD::SLU_Z;
  typedef SLUD::Z::doublecomplex type;
  typedef double magnitude_type;
  typedef SLUD::Z::LUstruct_t LUstruct_t;
  typedef SLUD::Z::SOLVEstruct_t SOLVEstruct_t;
};

  // It probably won't happen, but what if someone does create a
  // matrix or multivector with the SuperLU_DIST doublecomplex type
  // directly?
template <>
struct TypeMap<Superludist,SLUD::Z::doublecomplex>
{
  static const SLUD::Dtype_t dtype = SLUD::SLU_Z;
  typedef SLUD::Z::doublecomplex type;
  typedef double magnitude_type;
  typedef SLUD::Z::LUstruct_t LUstruct_t;
  typedef SLUD::Z::SOLVEstruct_t SOLVEstruct_t;
};

#endif	// HAVE_TEUCHOS_COMPLEX

/* \endcond Superludist_type_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_SUPERLUDIST_TYPEMAP_HPP
