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
   \file   Amesos2_Superlu_TypeMap.hpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Mon May 31 23:12:32 2010

   \brief Provides definition of SuperLU types as well as conversions and type
	  traits.

*/

#ifndef AMESOS2_SUPERLU_TYPEMAP_HPP
#define AMESOS2_SUPERLU_TYPEMAP_HPP

#include <functional>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include <Teuchos_as.hpp>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <Teuchos_SerializationTraits.hpp>
#endif

#include "Amesos2_TypeMap.hpp"

/* The SuperLU comples headers file only need to be included if
   complex has been enabled in Teuchos.  In addition we only need to
   define the conversion and printing functions if complex has been
   enabled. */
namespace SLU {

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
#endif // HAVE_TEUCHOS_COMPLEX

} // end extern "C"

} // end namespace SLU

#ifdef HAVE_TEUCHOS_COMPLEX

/* ==================== Conversion ==================== */
namespace Teuchos {

/**
 * \defgroup slu_conversion Conversion definitions for SLU types.
 *
 * Define specializations of Teuchos::as<> for the SLU types.
 * @{
 */

template <>
class ValueTypeConversionTraits<SLU::C::complex, Kokkos::complex<float>>
{
public:
  static SLU::C::complex convert( const Kokkos::complex<float> t ) {
    SLU::C::complex ret;
    ret.r = t.real();
    ret.i = t.imag();
    return( ret );
  }

  static SLU::C::complex safeConvert( const Kokkos::complex<float> t ) {
    SLU::C::complex ret;
    ret.r = t.real();
    ret.i = t.imag();
    return( ret );
  }
};

template <>
class ValueTypeConversionTraits<SLU::Z::doublecomplex, Kokkos::complex<double>>
{
public:
  static SLU::Z::doublecomplex convert( const Kokkos::complex<double> t ) {
    SLU::Z::doublecomplex ret;
    ret.r = t.real();
    ret.i = t.imag();
    return( ret );
  }

  static SLU::Z::doublecomplex safeConvert( const Kokkos::complex<double> t ) {
    SLU::Z::doublecomplex ret;
    ret.r = t.real();
    ret.i = t.imag();
    return( ret );
  }
};

// Also convert from SLU types

template <>
class ValueTypeConversionTraits<Kokkos::complex<float>, SLU::C::complex>
{
public:
  static Kokkos::complex<float> convert( const SLU::C::complex t ) {
    return ( Kokkos::complex<float>( t.r, t.i ) );
  }

  static Kokkos::complex<float> safeConvert( const SLU::C::complex t ) {
    return ( Kokkos::complex<float>( t.r, t.i ) );
  }
};

template <>
class ValueTypeConversionTraits<Kokkos::complex<double>, SLU::Z::doublecomplex>
{
public:
  static Kokkos::complex<double> convert( const SLU::Z::doublecomplex t ) {
    return ( Kokkos::complex<double>( t.r, t.i ) );
  }

  static Kokkos::complex<double> safeConvert( const SLU::Z::doublecomplex t ) {
    return ( Kokkos::complex<double>( t.r, t.i ) );
  }
};

//@}  End Conversion group

} // end namespace Teuchos

#endif  // HAVE_TEUCHOS_COMPLEX

namespace Amesos2 {

template <class, class> class Superlu;

/* Specialize the Amesos2::TypeMap struct for Superlu types
 *
 * \cond Superlu_type_specializations
 */
template <>
struct TypeMap<Superlu,float>
{
  static SLU::Dtype_t dtype;
  typedef float convert_type;
  typedef float type;
  typedef float magnitude_type;
};


template <>
struct TypeMap<Superlu,double>
{
  static SLU::Dtype_t dtype;
  typedef double convert_type;
  typedef double type;
  typedef double magnitude_type;
};


#ifdef HAVE_TEUCHOS_COMPLEX

template <>
struct TypeMap<Superlu,std::complex<float> >
{
  static SLU::Dtype_t dtype;
  typedef SLU::C::complex convert_type; // to create array before calling superlu
  typedef Kokkos::complex<float> type;
  typedef float magnitude_type;
};


template <>
struct TypeMap<Superlu,std::complex<double> >
{
  static SLU::Dtype_t dtype;
  typedef SLU::Z::doublecomplex convert_type; // to create array before calling superlu
  typedef Kokkos::complex<double> type;
  typedef double magnitude_type;
};


template <>
struct TypeMap<Superlu,Kokkos::complex<float> >
{
  static SLU::Dtype_t dtype;
  typedef SLU::C::complex convert_type; // to create array before calling superlu
  typedef Kokkos::complex<float> type;
  typedef float magnitude_type;
};


template <>
struct TypeMap<Superlu,Kokkos::complex<double> >
{
  static SLU::Dtype_t dtype;
  typedef SLU::Z::doublecomplex convert_type; // to create array before calling superlu
  typedef Kokkos::complex<double> type;
  typedef double magnitude_type;
};


#endif  // HAVE_TEUCHOS_COMPLEX

/* \endcond Superlu_type_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_SUPERLU_TYPEMAP_HPP
