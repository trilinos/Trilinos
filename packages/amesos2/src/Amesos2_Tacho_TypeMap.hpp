// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_TACHO_TYPEMAP_HPP
#define AMESOS2_TACHO_TYPEMAP_HPP

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

template <class, class> class TachoSolver;

/* Specialize the Amesos2::TypeMap struct for Tacho types
 *
 * \cond Tacho_type_specializations
 */

template <>
struct TypeMap<TachoSolver,float>
{
  typedef float type;
  typedef float magnitude_type;
};

template <>
struct TypeMap<TachoSolver,double>
{
  typedef double type;
  typedef double magnitude_type;
};


#ifdef HAVE_TEUCHOS_COMPLEX


template <>
struct TypeMap<TachoSolver,std::complex<float> >
{
  typedef Kokkos::complex<float> type;
  typedef float magnitude_type;
};

template <>
struct TypeMap<TachoSolver,std::complex<double> >
{
  typedef Kokkos::complex<double> type;
  typedef double magnitude_type;
};

template <>
struct TypeMap<TachoSolver,Kokkos::complex<float> >
{
  typedef Kokkos::complex<float> type;
  typedef float magnitude_type;
};

template <>
struct TypeMap<TachoSolver,Kokkos::complex<double> >
{
  typedef Kokkos::complex<double> type;
  typedef double magnitude_type;
};

#endif  // HAVE_TEUCHOS_COMPLEX

/* \endcond Tacho_type_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_TACHO_TYPEMAP_HPP
