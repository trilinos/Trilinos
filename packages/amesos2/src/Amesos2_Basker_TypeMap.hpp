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
   \file   Amesos2_Basker_TypeMap.hpp
   \author Joshua Dennis Booth <jdbooth@sandia.gov>
           Siva Rajamanickam <srajama@sandia.gov>

   \brief Provides definition of Basker types

*/

#ifndef AMESOS2_BASKER_TYPEMAP_HPP
#define AMESOS2_BASKER_TYPEMAP_HPP

#include <functional>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include <Teuchos_as.hpp>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <Teuchos_SerializationTraits.hpp>
#endif

#include "Amesos2_TypeMap.hpp"

namespace Amesos2 {

template <class, class> class Basker;

/* Specialize the Amesos2::TypeMap struct for Basker types
 * \cond Basker_type_specializations
 */

template <>
struct TypeMap<Basker,float>
{
  static double dtype;
  typedef double type;
};

template <>
struct TypeMap<Basker,double>
{
  static double dtype;
  typedef double type;
};


#ifdef HAVE_TEUCHOS_COMPLEX

template <>
struct TypeMap<Basker,std::complex<float> >
{
  static std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};

template <>
struct TypeMap<Basker,std::complex<double> >
{
  static std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};

template <>
struct TypeMap<Basker,Kokkos::complex<float> >
{
  static std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};

template <>
struct TypeMap<Basker,Kokkos::complex<double> >
{
  static std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};


#endif  // HAVE_TEUCHOS_COMPLEX

/* \endcond Basker_type_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_BASKER_TYPEMAP_HPP
