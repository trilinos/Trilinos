// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  typedef float dtype;
  typedef float type;
};

template <>
struct TypeMap<ShyLUBasker,double>
{
  typedef double dtype;
  typedef double type;
};


#ifdef HAVE_TEUCHOS_COMPLEX

template <>
struct TypeMap<ShyLUBasker,std::complex<float> >
{
  typedef std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};

template <>
struct TypeMap<ShyLUBasker,std::complex<double> >
{
  typedef std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};

template <>
struct TypeMap<ShyLUBasker,Kokkos::complex<float> >
{
  typedef std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};

template <>
struct TypeMap<ShyLUBasker,Kokkos::complex<double> >
{
  typedef std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};

#endif  // HAVE_TEUCHOS_COMPLEX

/* \endcond ShyLUBasker_type_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_SHYLUBASKER_TYPEMAP_HPP
