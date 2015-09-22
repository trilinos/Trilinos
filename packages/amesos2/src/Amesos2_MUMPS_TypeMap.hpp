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
   \file   Amesos2_MUMPS_TypeMap.hpp
   \author Joshua Dennis Booth <jdbooth@sandia.gov>

   \brief Provides definition of MUMPS types as well as conversions and type
          traits.

*/

#ifndef AMESOS2_MUMPS_TYPEMAP_HPP
#define AMESOS2_MUMPS_TYPEMAP_HPP

//#include <functional>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include <Teuchos_as.hpp>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <Teuchos_SerializationTraits.hpp>
#endif

#include "Amesos2_TypeMap.hpp"

namespace Amesos2
{
  namespace MUMPST
  {
    extern "C"
    {
    #include "smumps_c.h"
    #include "dmumps_c.h"
    #ifdef HAVE_TEUCHOS_COMPLEX
    #include "cmumps_c.h"
    #include "zmumps_c.h"
    #endif
    }
  }
}


namespace Amesos2 {

template <class, class> class MUMPS;

/* Specialize the Amesos2::TypeMap struct for Mumps types
 * TODO: Mostly dummy assignments as MUMPS is templated. Remove if possible.
 *
 * \cond Mumps_type_specializations
 */


  template <>
  struct TypeMap<MUMPS,float>
  {
    typedef float type;
    typedef float magnitude_type;
    typedef MUMPST::SMUMPS_STRUC_C MUMPS_STRUC_C; 
  };
 

  template <>
  struct TypeMap<MUMPS,double>
  {
    typedef double type;
    typedef double magnitude_type;
    typedef MUMPST::DMUMPS_STRUC_C MUMPS_STRUC_C; 
  };
  
#ifdef HAVE_TEUCHOS_COMPLEX
  
  template <>
  struct TypeMap<MUMPS,std::complex<float> >
  {
    typedef MUMPST::CMUMPS_COMPLEX type;
    typedef MUMPST::CMUMPS_COMPLEX  magnitude_type;
    typedef MUMPST::CMUMPS_STRUC_C MUMPS_STRUC_C;
  };
    
  template <>
  struct TypeMap<MUMPS,std::complex<double> >
  {
    typedef MUMPST::ZMUMPS_COMPLEX type;
    typedef MUMPST::ZMUMPS_COMPLEX   magnitude_type;
    typedef MUMPST::ZMUMPS_STRUC_C MUMPS_STRUC_C; 
  };
  
  template <>
  struct TypeMap<MUMPS, MUMPST::CMUMPS_COMPLEX>
  {
    typedef MUMPST::CMUMPS_COMPLEX type;
    typedef MUMPST::CMUMPS_COMPLEX magnitude_type;
    typedef MUMPST::CMUMPS_STRUC_C MUMPS_STRUC_C;
  };
  
  template <>
  struct TypeMap<MUMPS, MUMPST::ZMUMPS_COMPLEX>
  {
    typedef MUMPST::ZMUMPS_COMPLEX type;
    typedef MUMPST::ZMUMPS_COMPLEX magnitude_type;
    typedef MUMPST::ZMUMPS_STRUC_C MUMPS_STRUC_C;
  };
  
#endif  // HAVE_TEUCHOS_COMPLEX

/* \endcond MUMPS_type_specializations */

} // end namespace Amesos2


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

  
  template <typename TypeFrom> //float based complex
  class ValueTypeConversionTraits<Amesos2::MUMPST::CMUMPS_COMPLEX, TypeFrom >
  {
  public:
    static Amesos2::MUMPST::CMUMPS_COMPLEX convert( const TypeFrom t )
    {
      Amesos2::MUMPST::CMUMPS_COMPLEX ret;
      ret.r = Teuchos::as<float>(t.real());
      ret.i = Teuchos::as<float>(t.imag());
      return( ret );
    }

    static Amesos2::MUMPST::CMUMPS_COMPLEX safeConvert( const TypeFrom t )
    {
      Amesos2::MUMPST::CMUMPS_COMPLEX ret;
      ret.r = Teuchos::as<float>(t.real());
      ret.i = Teuchos::as<float>(t.imag());
      return( ret );
    }
  };
  

  
  template <typename TypeFrom> //double based
  class ValueTypeConversionTraits<Amesos2::MUMPST::ZMUMPS_COMPLEX , TypeFrom >
  {
  public:
    static Amesos2::MUMPST::ZMUMPS_COMPLEX  convert( const TypeFrom t )
    {
      Amesos2::MUMPST::ZMUMPS_COMPLEX ret;
      ret.r = Teuchos::as<double>( t.real() );
      ret.i = Teuchos::as<double>( t.imag() );
      return (ret);
    }

    // No special checks for safe Convert
    static Amesos2::MUMPST::ZMUMPS_COMPLEX  safeConvert( const TypeFrom t )
    {
      Amesos2::MUMPST::ZMUMPS_COMPLEX ret;
      ret.r = Teuchos::as<double>( t.real() );
      ret.i = Teuchos::as<double>( t.imag() );
      return (ret);
    }
  };
  
  template <typename TypeTo>
  class ValueTypeConversionTraits<TypeTo, Amesos2::MUMPST::CMUMPS_COMPLEX>
  {
  public:
    static TypeTo convert(const Amesos2::MUMPST::CMUMPS_COMPLEX t)
    {
      typedef typename TypeTo::value_type   value_type;
      value_type ret_r = Teuchos::as<value_type>(t.r);
      value_type ret_i = Teuchos::as<value_type>(t.i);
      return (TypeTo(ret_r, ret_i));
    }
    //No special checks for safe convert
    static TypeTo safeConvert(const Amesos2::MUMPST::CMUMPS_COMPLEX t)
    {
      typedef typename TypeTo::value_type   value_type;
      value_type ret_r = Teuchos::as<value_type>(t.r);
      value_type ret_i = Teuchos::as<value_type>(t.i);
      return (TypeTo(ret_r, ret_i));
    }
   
  };
  
  template <typename TypeTo>
  class ValueTypeConversionTraits<TypeTo, Amesos2::MUMPST::ZMUMPS_COMPLEX>
  {
  public:
    static TypeTo convert(const Amesos2::MUMPST::ZMUMPS_COMPLEX t)
    {
      typedef typename TypeTo::value_type   value_type;
      value_type ret_r = Teuchos::as<value_type>(t.r);
      value_type ret_i = Teuchos::as<value_type>(t.i);
      return (TypeTo(ret_r, ret_i));
    }
    //No special checks for safe convert
    static TypeTo safeConvert(const Amesos2::MUMPST::ZMUMPS_COMPLEX t)
    {
      typedef typename TypeTo::value_type   value_type;
      value_type ret_r = Teuchos::as<value_type>(t.r);
      value_type ret_i = Teuchos::as<value_type>(t.i);
      return (TypeTo(ret_r, ret_i));
    }
   
  };


//#endif
//@}  End Conversion group


} // end namespace Teuchos

#endif	// HAVE_TEUCHOS_COMPLEX



#endif  // AMESOS2_MUMPS_TYPEMAP_HPP
