// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \file   Amesos2_KLU2_TypeMap.hpp
   \author Siva Rajamanickam <srajama@sandia.gov>

   \brief Provides definition of KLU2 types as well as conversions and type
	  traits.

*/

#ifndef AMESOS2_KLU2_TYPEMAP_HPP
#define AMESOS2_KLU2_TYPEMAP_HPP

#include <functional>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include <Teuchos_as.hpp>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <Teuchos_SerializationTraits.hpp>
#endif

#include "Amesos2_TypeMap.hpp"

namespace KLU2 {

#include "klu2_ext.hpp"	// for Dtype_t declaration

} // end namespace KLU

namespace Amesos2 {

template <class, class> class KLU2;

/* Specialize the Amesos2::TypeMap struct for KLU2 types
 * TODO: Mostly dummy assignments as KLU2 is templated. Remove if possible.
 *
 * \cond KLU2_type_specializations
 */
template <>
struct TypeMap<KLU2,float>
{
  typedef float dtype;
  typedef float type;
};

template <>
struct TypeMap<KLU2,double>
{
  typedef double dtype;
  typedef double type;
};

#ifdef HAVE_TEUCHOS_COMPLEX

template <>
struct TypeMap<KLU2,std::complex<float> >
{
  typedef std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};

template <>
struct TypeMap<KLU2,std::complex<double> >
{
  typedef std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};

template <>
struct TypeMap<KLU2,Kokkos::complex<float> >
{
  typedef std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};

template <>
struct TypeMap<KLU2,Kokkos::complex<double> >
{
  typedef std::complex<double> dtype;
  typedef Kokkos::complex<double> type;
};

#endif  // HAVE_TEUCHOS_COMPLEX

/* \endcond KLU2_type_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_SUPERLU_TYPEMAP_HPP
