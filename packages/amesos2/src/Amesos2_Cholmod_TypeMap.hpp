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

#ifndef AMESOS2_CHOLMOD_TYPEMAP_HPP
#define AMESOS2_CHOLMOD_TYPEMAP_HPP

#include <functional>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include <Teuchos_as.hpp>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <Teuchos_SerializationTraits.hpp>
#endif

#include "Amesos2_TypeMap.hpp"

namespace Amesos2{
  namespace CHOL {
    
    struct complex
    {
      double complexpair[2];
    };

    #   include "cholmod.h"
    //#   include <mkl_types.h>
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
  class ValueTypeConversionTraits<Amesos2::CHOL::complex, TypeFrom>
  {
  public:
    static Amesos2::CHOL::complex convert( const TypeFrom t )
    {                           // adapt conversion as necessary
      Amesos2::CHOL::complex ret;
      //ret.r = Teuchos::as<float>(t.real());
      //ret.i = Teuchos::as<float>(t.imag());
      ret.complexpair[0] = Teuchos::as<float>(t.real());
      ret.complexpair[1] = Teuchos::as<float>(t.imag());
      return( ret );
    }

    static Amesos2::CHOL::complex safeConvert( const TypeFrom t )
    {                           // adapt conversion as necessary
      Amesos2::CHOL::complex ret;
      //ret.r = Teuchos::as<float>(t.real());
      //ret.i = Teuchos::as<float>(t.imag());
      ret.complexpair[0] = Teuchos::as<float>(t.real());
      ret.complexpair[1] = Teuchos::as<float>(t.imag());
      return( ret );
    }
  };

  // Also convert *from* New_Solver types
  template <typename TypeTo>
  class ValueTypeConversionTraits<TypeTo, Amesos2::CHOL::complex>
  {
  public:
    static TypeTo convert( const Amesos2::CHOL::complex t )
    {                           // adapt conversion as necessary
      typedef typename TypeTo::value_type value_type;
      //value_type ret_r = Teuchos::as<value_type>( t.real );
      //value_type ret_i = Teuchos::as<value_type>( t.imag );
      value_type ret_r = Teuchos::as<value_type>(t.complexpair[0]);
      value_type ret_i = Teuchos::as<value_type>(t.complexpair[1]);
      return ( TypeTo( ret_r, ret_i ) );
    }

    static TypeTo safeConvert( const Amesos2::CHOL::complex t )
    {                           // adapt conversion as necessary
      typedef typename TypeTo::value_type value_type;
      //value_type ret_r = Teuchos::as<value_type>( t.real );
      //value_type ret_i = Teuchos::as<value_type>( t.imag );
      value_type ret_r = Teuchos::as<value_type>(t.complexpair[0]);
      value_type ret_i = Teuchos::as<value_type>(t.complexpair[1]);
      return ( TypeTo( ret_r, ret_i ) );
    }
  };

  //@}  End Conversion group

} // end namespace Teuchos


namespace Amesos2 {

  // forward declaration due to circular reference
  template <class, class> class Cholmod;

  template <>
  struct TypeMap<Cholmod,float>
  {
    typedef float type;
    typedef float magnitude_type;
  };

  template <>
  struct TypeMap<Cholmod,double>
  {
    typedef double type;
    typedef double magnitude_type;
  };

#ifdef HAVE_TEUCHOS_COMPLEX

  template <>
  struct TypeMap<Cholmod,std::complex<double> >
  {
    typedef CHOL::complex type;
    typedef double magnitude_type;
  };

#endif  // HAVE_TEUCHOS_COMPLEX

  /* \endcond Choldmod_type_specializations */

} // end namespace Amesos

#endif  // AMESOS2_CHOLMOD_TYPEMAP_HPP
