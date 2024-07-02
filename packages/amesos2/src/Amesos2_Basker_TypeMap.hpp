// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  typedef double dtype;
  typedef double type;
};

template <>
struct TypeMap<Basker,double>
{
  typedef double dtype;
  typedef double type;
};


#ifdef HAVE_TEUCHOS_COMPLEX

template <>
struct TypeMap<Basker,std::complex<float> >
{
  typedef std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};

template <>
struct TypeMap<Basker,std::complex<double> >
{
  typedef std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};

template <>
struct TypeMap<Basker,Kokkos::complex<float> >
{
  typedef std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};

template <>
struct TypeMap<Basker,Kokkos::complex<double> >
{
  typedef std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};


#endif  // HAVE_TEUCHOS_COMPLEX

/* \endcond Basker_type_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_BASKER_TYPEMAP_HPP
