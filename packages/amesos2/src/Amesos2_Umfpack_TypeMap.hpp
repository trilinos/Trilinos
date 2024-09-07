// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_UMFPACK_TYPEMAP_HPP
#define AMESOS2_UMFPACK_TYPEMAP_HPP

#include <functional>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include <Teuchos_as.hpp>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <Teuchos_SerializationTraits.hpp>
#endif

#include "Amesos2_TypeMap.hpp"

/* ==================== Conversion ==================== */
namespace Teuchos {

/**
 * \defgroup umfpack_conversion Conversion definitions for UMFPACK types.
 *
 * Define specializations of Teuchos::as<> for the UMFPACK types.
 *
 * These specializations are meant to work with any complex data type that
 * implements the same interface as the STL complex type.
 *
 * @{
 */

#ifdef HAVE_TEUCHOS_COMPLEX

// Provide conversion from std::complex<float> to std::complex<double>
template <typename TypeFrom>
class ValueTypeConversionTraits<std::complex<double>, TypeFrom>
{
public:
  static std::complex<double> convert( const TypeFrom t )
    {
      return std::complex<double>(
        Teuchos::as<double>(t.real()),
        Teuchos::as<double>(t.imag()));
    }

  static std::complex<double> safeConvert( const TypeFrom t )
    {
      return std::complex<double>(
        Teuchos::as<double>(t.real()),
        Teuchos::as<double>(t.imag()));
    }
};


// Also convert from UMFPACK types - convert back to std::complex<float> to std::complex<double>
template <typename TypeTo>
class ValueTypeConversionTraits<TypeTo, std::complex<double>>
{
public:
  static TypeTo convert( const std::complex<double> t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.real() );
      value_type ret_i = Teuchos::as<value_type>( t.imag() );
      return ( TypeTo( ret_r, ret_i ) );
    }

  // No special checks for safe Convert
  static TypeTo safeConvert( const std::complex<double> t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.real() );
      value_type ret_i = Teuchos::as<value_type>( t.imag() );
      return ( TypeTo( ret_r, ret_i ) );
    }
};

#endif  // HAVE_TEUCHOS_COMPLEX

//@}  End Conversion group

} // end namespace Teuchos

namespace Amesos2 {

template <class, class> class Umfpack;

/* Specialize the Amesos2::TypeMap struct for Umfpack types
 *
 * \cond Umfpack_type_specializations
 */

template <>
struct TypeMap<Umfpack,float> // provide conversion from float to double
{
  typedef double type;
  typedef double magnitude_type;
};

template <>
struct TypeMap<Umfpack,double>
{
  typedef double type;
  typedef double magnitude_type;
};


#ifdef HAVE_TEUCHOS_COMPLEX

template <>
struct TypeMap<Umfpack,std::complex<float> > // provide conversion from std::complex<float> to std::complex<double>
{
  typedef Kokkos::complex<double> type;
  typedef double magnitude_type;
};

template <>
struct TypeMap<Umfpack,std::complex<double> >
{
  typedef Kokkos::complex<double> type;
  typedef double magnitude_type;
};


#endif  // HAVE_TEUCHOS_COMPLEX

/* \endcond Umfpack_type_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_UMFPACK_TYPEMAP_HPP
