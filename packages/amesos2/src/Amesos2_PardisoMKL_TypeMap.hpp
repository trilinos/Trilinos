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
   \file   Amesos2_PardisoMKL_TypeMap.hpp
   \author John Doe <jd@sandia.gov>
   \date

   \brief Provides definition of PardisoMKL types as well as
          conversions and type traits.  For the purpose of
          demonstration, we assume that PardisoMKL has defined its own
          complex data-types called `complex' and `doublecomplex'.
*/

#ifndef AMESOS2_PARDISOMKL_TYPEMAP_HPP
#define AMESOS2_PARDISOMKL_TYPEMAP_HPP

#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include <Teuchos_as.hpp>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <Teuchos_SerializationTraits.hpp>
#endif

#include "Amesos2_TypeMap.hpp"

namespace Amesos2{
  namespace PMKL {
#   include <mkl_dss.h>
#   include <mkl_types.h>
  } // end namespace PMKL
} // end namespace Amesos2


/* ==================== Conversion ====================
 *
 * Define here, in the Teuchos namespace, any conversions between
 * commonly used date types and the solver-specific data types.  Use
 * template specializations of the Teuchos::ValueTypeConversionTraits
 * class.
 */
namespace Teuchos {

  template <typename TypeFrom>
  class ValueTypeConversionTraits<Amesos2::PMKL::_MKL_Complex8, TypeFrom>
  {
  public:
    static Amesos2::PMKL::_MKL_Complex8 convert( const TypeFrom t )
    {                           // adapt conversion as necessary
      Amesos2::PMKL::_MKL_Complex8 ret;
      ret.real = Teuchos::as<float>(t.real());
      ret.imag = Teuchos::as<float>(t.imag());
      return( ret );
    }

    static Amesos2::PMKL::_MKL_Complex8 safeConvert( const TypeFrom t )
    {                           // adapt conversion as necessary
      Amesos2::PMKL::_MKL_Complex8 ret;
      ret.real = Teuchos::as<float>(t.real());
      ret.imag = Teuchos::as<float>(t.imag());
      return( ret );
    }
  };


  template <typename TypeFrom>
  class ValueTypeConversionTraits<Amesos2::PMKL::_DOUBLE_COMPLEX_t, TypeFrom>
  {
  public:
    static Amesos2::PMKL::_DOUBLE_COMPLEX_t convert( const TypeFrom t )
    {                           // adapt conversion as necessary
      Amesos2::PMKL::_DOUBLE_COMPLEX_t ret;
      ret.r = Teuchos::as<double>(t.real());
      ret.i = Teuchos::as<double>(t.imag());
      return( ret );
    }

    static Amesos2::PMKL::_DOUBLE_COMPLEX_t safeConvert( const TypeFrom t )
    {                           // adapt conversion as necessary
      Amesos2::PMKL::_DOUBLE_COMPLEX_t ret;
      ret.r = Teuchos::as<double>(t.real());
      ret.i = Teuchos::as<double>(t.imag());
      return( ret );
    }
  };


  // Also convert *from* New_Solver types
  template <typename TypeTo>
  class ValueTypeConversionTraits<TypeTo, Amesos2::PMKL::_MKL_Complex8>
  {
  public:
    static TypeTo convert( const Amesos2::PMKL::_MKL_Complex8 t )
    {                           // adapt conversion as necessary
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.real );
      value_type ret_i = Teuchos::as<value_type>( t.imag );
      return ( TypeTo( ret_r, ret_i ) );
    }

    static TypeTo safeConvert( const Amesos2::PMKL::_MKL_Complex8 t )
    {                           // adapt conversion as necessary
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.real );
      value_type ret_i = Teuchos::as<value_type>( t.imag );
      return ( TypeTo( ret_r, ret_i ) );
    }
  };


  template <typename TypeTo>
  class ValueTypeConversionTraits<TypeTo, Amesos2::PMKL::_DOUBLE_COMPLEX_t>
  {
  public:
    static TypeTo convert( const Amesos2::PMKL::_DOUBLE_COMPLEX_t t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }

    // No special checks for safe Convert
    static TypeTo safeConvert( const Amesos2::PMKL::_DOUBLE_COMPLEX_t t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }
  };

  //@}  End Conversion group

} // end namespace Teuchos


namespace Amesos2 {

  // forward declaration due to circular reference
  template <class, class> class PardisoMKL;

  /* Specialize the Amesos::TypeMap struct for PardisoMKL types.
   *
   * Additional nested types may be added without harm.  For an example, look at
   * Amesos2_Superlu_TypeMap.hpp
   */
  template <>
  struct TypeMap<PardisoMKL,float>
  {
    typedef PMKL::_REAL_t type;
    typedef PMKL::_REAL_t magnitude_type;
  };


  template <>
  struct TypeMap<PardisoMKL,double>
  {
    typedef PMKL::_DOUBLE_PRECISION_t type;
    typedef PMKL::_DOUBLE_PRECISION_t magnitude_type;
  };

#ifdef HAVE_TEUCHOS_COMPLEX

  /*
   * We map the std complex types to the appropriate PardisoMKL complex
   * types.
   */

  template <>
  struct TypeMap<PardisoMKL,std::complex<float> >
  {
    typedef PMKL::_MKL_Complex8 type;
    typedef PMKL::_REAL_t magnitude_type;
  };


  template <>
  struct TypeMap<PardisoMKL,std::complex<double> >
  {
    typedef PMKL::_DOUBLE_COMPLEX_t type;
    typedef PMKL::_DOUBLE_PRECISION_t magnitude_type;
  };


  template <>
  struct TypeMap<PardisoMKL,PMKL::_MKL_Complex8>
  {
    typedef PMKL::_MKL_Complex8 type;
    typedef PMKL::_REAL_t magnitude_type;
  };


  template <>
  struct TypeMap<PardisoMKL,PMKL::_DOUBLE_COMPLEX_t>
  {
    typedef PMKL::_DOUBLE_COMPLEX_t type;
    typedef PMKL::_DOUBLE_PRECISION_t magnitude_type;
  };
#endif  // HAVE_TEUCHOS_COMPLEX

  template <>
  struct TypeMap<PardisoMKL,int>
  {
    typedef PMKL::_INTEGER_t type;
  };

  template <>
  struct TypeMap<PardisoMKL,long long int>
  {
    typedef long long int type;
  };

  /*
   * We check whether the size of long int is bigger than an int.  If
   * it is, then long int should be the same size as a long long int,
   * so we can safely promote.  Otherwise, long int will probably be
   * the same size as int, and we can safely treat it as such.
   */
  template <>
  struct TypeMap<PardisoMKL,long int>
  {
    typedef Meta::if_then_else<
      sizeof(int) < sizeof(long int),
      TypeMap<PardisoMKL,long long int>::type,
      TypeMap<PardisoMKL,int>::type >::type type;
  };

} // end namespace Amesos

#endif  // AMESOS2_PARDISOMKL_TYPEMAP_HPP
