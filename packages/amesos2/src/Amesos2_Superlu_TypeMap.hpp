/**
   \file   Amesos2_Superlu_TypeMap.hpp
   \author Eric Bavier <etbavier@etbavier@sandia.gov>
   \date   Mon May 31 23:12:32 2010

   \brief Provides definition of SuperLU types as well as conversions and type
          traits.

*/

#ifndef AMESOS2_SUPERLU_TYPEMAP_HPP
#define AMESOS2_SUPERLU_TYPEMAP_HPP

#include <complex>

#include <Teuchos_as.hpp>

#include "Amesos2_TypeMap.hpp"

namespace SLU {

typedef int int_t;

extern "C" {

#include "slu_Cnames.h"
#include "supermatrix.h"
#include "slu_util.h"

namespace C {
#include "slu_scomplex.h"     // single-precision complex data type definitions
}

namespace Z {
#include "slu_dcomplex.h"     // double-precision complex data type definitions
}

} // end extern "C"
} // end namespace SLU


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
class ValueTypeConversionTraits<SLU::C::complex, TypeFrom>
{
public:
  static SLU::C::complex convert( const TypeFrom t )
    {
      SLU::C::complex ret;
      ret.r = Teuchos::as<float>(t.real());
      ret.i = Teuchos::as<float>(t.imag());
      return( ret );
    }

  static SLU::C::complex safeConvert( const TypeFrom t )
    {
      SLU::C::complex ret;
      ret.r = Teuchos::as<float>(t.real());
      ret.i = Teuchos::as<float>(t.imag());
      return( ret );
    }
};


template <typename TypeFrom>
class ValueTypeConversionTraits<SLU::Z::doublecomplex, TypeFrom>
{
public:
  static SLU::Z::doublecomplex convert( const TypeFrom t )
    {
      SLU::Z::doublecomplex ret;
      ret.r = Teuchos::as<double>(t.real());
      ret.i = Teuchos::as<double>(t.imag());
      return( ret );
    }

  static SLU::Z::doublecomplex safeConvert( const TypeFrom t )
    {
      SLU::Z::doublecomplex ret;
      ret.r = Teuchos::as<double>(t.real());
      ret.i = Teuchos::as<double>(t.imag());
      return( ret );
    }
};


// Also convert from SLU types
template <typename TypeTo>
class ValueTypeConversionTraits<TypeTo, SLU::C::complex>
{
public:
  static TypeTo convert( const SLU::C::complex t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }

  // No special checks for safe Convert
  static TypeTo safeConvert( const SLU::C::complex t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }
};


template <typename TypeTo>
class ValueTypeConversionTraits<TypeTo, SLU::Z::doublecomplex>
{
public:
  static TypeTo convert( const SLU::Z::doublecomplex t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }

  // No special checks for safe Convert
  static TypeTo safeConvert( const SLU::Z::doublecomplex t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }
};

template <typename Ordinal>
class SerializationTraits<Ordinal,SLU::C::complex>
  : public DirectSerializationTraits<Ordinal,SLU::C::complex>
{};

template <typename Ordinal>
class SerializationTraits<Ordinal,SLU::Z::doublecomplex>
  : public DirectSerializationTraits<Ordinal,SLU::Z::doublecomplex>
{};

//@}  End Conversion group

} // end namespace Teuchos


namespace Amesos {

template <class, class> class Superlu;

/* Specialize the Amesos::TypeMap struct for Superlu types
 *
 * \cond Superlu_type_specializations 
 */
template <>
struct TypeMap<Superlu,float>
{
  static SLU::Dtype_t dtype;
  typedef float type;
  typedef float magnitude_type;
};

SLU::Dtype_t TypeMap<Superlu,float>::dtype = SLU::SLU_S;


template <>
struct TypeMap<Superlu,double>
{
  static SLU::Dtype_t dtype;
  typedef double type;
  typedef double magnitude_type;
};

SLU::Dtype_t TypeMap<Superlu,double>::dtype = SLU::SLU_D;


template <>
struct TypeMap<Superlu,std::complex<float> >
{
  static SLU::Dtype_t dtype;
  typedef SLU::C::complex type;
  typedef float magnitude_type;
};

SLU::Dtype_t TypeMap<Superlu,std::complex<float> >::dtype = SLU::SLU_C;


template <>
struct TypeMap<Superlu,std::complex<double> >
{
  static SLU::Dtype_t dtype;
  typedef SLU::Z::doublecomplex type;
  typedef double magnitude_type;
};

SLU::Dtype_t TypeMap<Superlu,std::complex<double> >::dtype = SLU::SLU_Z;


template <>
struct TypeMap<Superlu,SLU::C::complex>
{
  static SLU::Dtype_t dtype;
  typedef SLU::C::complex type;
  typedef float magnitude_type;
};

SLU::Dtype_t TypeMap<Superlu,SLU::C::complex>::dtype = SLU::SLU_C;


template <>
struct TypeMap<Superlu,SLU::Z::doublecomplex>
{
  static SLU::Dtype_t dtype;
  typedef SLU::Z::doublecomplex type;
  typedef double magnitude_type;
};

SLU::Dtype_t TypeMap<Superlu,SLU::Z::doublecomplex>::dtype = SLU::SLU_Z;

/* \endcond Superlu_type_specializations */


} // end namespace Amesos

#endif  // AMESOS2_SUPERLU_TYPEMAP_HPP
