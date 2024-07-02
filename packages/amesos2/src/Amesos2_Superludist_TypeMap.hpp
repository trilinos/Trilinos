// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

//#if SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
//#endif

#include <functional>

#include <Teuchos_as.hpp>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <Teuchos_SerializationTraits.hpp>
#endif

#include "Amesos2_TypeMap.hpp"

#ifdef KOKKOS_ENABLE_CUDA
  #include <cublas_v2.h>
  #include <cuda_runtime_api.h>
#endif


namespace SLUD {

#if SUPERLU_DIST_MAJOR_VERSION > 4
// SuperLU_Dist before major version 5 does not contain the config file
#include "superlu_dist_config.h" // provides define for size 32 or 64 int_t
#endif

  /// use the same function with name space in the macro
#define USER_FREE(addr) SLUD::superlu_free_dist(addr)

  // undefine compiler guard in case we also have the sequential
  // SuperLU enabled
#undef __SUPERLU_SUPERMATRIX
#include "superlu_defs.h"
//

#if SUPERLU_DIST_MAJOR_VERSION > 4
  typedef superlu_dist_options_t   amesos2_superlu_dist_options_t;
  typedef superlu_dist_mem_usage_t amesos2_superlu_dist_mem_usage_t;
#define AMESOS2_ENABLES_SUPERLUDIST_VERSION5_AND_HIGHER 1
#else
  typedef superlu_options_t        amesos2_superlu_dist_options_t;
  typedef mem_usage_t              amesos2_superlu_dist_mem_usage_t;
#endif


  namespace D {
#include "superlu_ddefs.h"	// double-precision real definitions
  }

#if defined(HAVE_TEUCHOS_COMPLEX)  && !defined(__clang__)
  namespace Z {
#include "superlu_zdefs.h"     // double-precision complex definitions
  }
#endif  // HAVE_TEUCHOS_COMPLEX

// multiplication of SLUD types
template <typename slu_scalar_t, typename slu_mag_t>
struct slu_dist_mult {};

// This specialization handles the generic case were the scalar and
// magnitude types are double or float.
template <typename T>
struct slu_dist_mult<T,T> : std::multiplies<T> {};

// For namespace/macro reasons, we prefix our variables with amesos_*
template <>
struct slu_dist_mult<double,double>
{
  double operator()(double a, double b) {
    return( a*b );
  }
};

#if defined(HAVE_TEUCHOS_COMPLEX)  && !defined(__clang__)

  template <>
  struct slu_dist_mult<Z::doublecomplex,double>
  {
    Z::doublecomplex operator()(Z::doublecomplex amesos_z, double amesos_d) {
      Z::doublecomplex amesos_zr;
      zd_mult(&amesos_zr, &amesos_z, amesos_d);	// zd_mult is a macro, so no namespacing
      return( amesos_zr );
    }
  };

  template <>
  struct slu_dist_mult<Z::doublecomplex,Z::doublecomplex>
  {
    Z::doublecomplex operator()(Z::doublecomplex amesos_z1, Z::doublecomplex amesos_z2) {
      Z::doublecomplex amesos_zr;
      zz_mult(&amesos_zr, &amesos_z1, &amesos_z2);    // zz_mult is a macro, so no namespacing
      return( amesos_zr );
    }
  };
#endif	// HAVE_TEUCHOS_COMPLEX
} // end namespace SLUD
#if defined(HAVE_TEUCHOS_COMPLEX)  && !defined(__clang__)


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
#if SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
  typedef SLUD::D::dLUstruct_t LUstruct_t;
  typedef SLUD::D::dSOLVEstruct_t SOLVEstruct_t;
  typedef SLUD::D::dScalePermstruct_t ScalePermstruct_t;
#else
  typedef SLUD::D::LUstruct_t LUstruct_t;
  typedef SLUD::D::SOLVEstruct_t SOLVEstruct_t;
  typedef SLUD::ScalePermstruct_t ScalePermstruct_t;
#endif
};

#if defined(HAVE_TEUCHOS_COMPLEX) && !defined(__clang__)
template <>
struct TypeMap<Superludist,std::complex<double> >
{
  static const SLUD::Dtype_t dtype = SLUD::SLU_Z;
  typedef SLUD::Z::doublecomplex type;
  typedef double magnitude_type;
#if SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
  typedef SLUD::Z::zLUstruct_t LUstruct_t;
  typedef SLUD::Z::zSOLVEstruct_t SOLVEstruct_t;
  typedef SLUD::Z::zScalePermstruct_t ScalePermstruct_t;
#else
  typedef SLUD::Z::LUstruct_t LUstruct_t;
  typedef SLUD::Z::SOLVEstruct_t SOLVEstruct_t;
  typedef SLUD::ScalePermstruct_t ScalePermstruct_t;
#endif
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
#if SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
  typedef SLUD::Z::zLUstruct_t LUstruct_t;
  typedef SLUD::Z::zSOLVEstruct_t SOLVEstruct_t;
  typedef SLUD::Z::zScalePermstruct_t ScalePermstruct_t;
#else
  typedef SLUD::Z::LUstruct_t LUstruct_t;
  typedef SLUD::Z::SOLVEstruct_t SOLVEstruct_t;
  typedef SLUD::ScalePermstruct_t ScalePermstruct_t;
#endif
};

#endif	// HAVE_TEUCHOS_COMPLEX

/* \endcond Superludist_type_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_SUPERLUDIST_TYPEMAP_HPP
