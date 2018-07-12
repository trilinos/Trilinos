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
   \file   Amesos2_ShyLUBasker_TypeMap.hpp
   \author Joshua Dennis Booth <jdbooth@sandia.gov>
           Siva Rajamanickam <srajama@sandia.gov>
           Nathan Ellingwood <ndellin@sandia.gov>

   \brief Provides definition of ShyLUBasker types as well as conversions and type
          traits.

*/

#ifndef AMESOS2_SHYLUBASKER_TYPEMAP_HPP
#define AMESOS2_SHYLUBASKER_TYPEMAP_HPP

#include <functional>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include <Teuchos_as.hpp>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <Teuchos_SerializationTraits.hpp>
#endif

#include "Amesos2_TypeMap.hpp"


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

#ifndef HAVE_AMESOS2_KLU2

template <>
class ValueTypeConversionTraits<std::complex<double>, std::complex<float> >
{
public:
  static std::complex<double> convert( const std::complex<float>  t )
    {
      std::complex<double> ret(Teuchos::as<double>(t.real()),
                                 Teuchos::as<double>(t.imag()));
      return( ret );
    }

  static std::complex<double> safeConvert( const std::complex<float>  t )
    {
      std::complex<double> ret(Teuchos::as<double>(t.real()),
                                 Teuchos::as<double>(t.imag()));
      return( ret );
    }
};


template <>
class ValueTypeConversionTraits<std::complex<float> , std::complex<double> >
{
public:
  static std::complex<float>  convert( const std::complex<double> t )
    {
      float ret_r = Teuchos::as<float>( t.real() );
      float ret_i = Teuchos::as<float>( t.imag() );
      std::complex<float> ret (ret_r,  ret_i);
      return (ret);
    }

  // No special checks for safe Convert
  static std::complex<float>  safeConvert( const std::complex<double> t )
    {
      float ret_r = Teuchos::as<float>( t.real() );
      float ret_i = Teuchos::as<float>( t.imag() );
      std::complex<float> ret (ret_r,  ret_i);
      return (ret);
    }
};


#endif
//@}  End Conversion group


} // end namespace Teuchos

#endif	// HAVE_TEUCHOS_COMPLEX


namespace Amesos2 {

template <class, class> class ShyLUBasker;

/* Specialize the Amesos2::TypeMap struct for ShyLUBasker types
 * TODO: Mostly dummy assignments as ShyLUBasker is templated. Remove if possible.
 *
 * \cond ShyLUBasker_type_specializations
 */

template <>
struct TypeMap<ShyLUBasker,float>
{
  static float dtype;
  typedef float type;
  typedef float magnitude_type;
};


template <>
struct TypeMap<ShyLUBasker,double>
{
  static double dtype;
  typedef double type;
  typedef double magnitude_type;
};


#ifdef HAVE_TEUCHOS_COMPLEX

template <>
struct TypeMap<ShyLUBasker,std::complex<float> >
{
  static std::complex<double> dtype;
  typedef std::complex<double> type;
  typedef double magnitude_type;
};


template <>
struct TypeMap<ShyLUBasker,std::complex<double> >
{
  static std::complex<double> dtype;
  typedef std::complex<double> type;
  typedef double magnitude_type;
};


#endif  // HAVE_TEUCHOS_COMPLEX

/* \endcond ShyLUBasker_type_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_SHYLUBASKER_TYPEMAP_HPP
