// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \file   Amesos2_Cholmod_TypeMap.hpp
   \author John Doe <jd@sandia.gov>
   \date

   \brief Provides definition of Cholmod types as well as
          conversions and type traits.
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

namespace Amesos2 {

  // forward declaration due to circular reference
  template <class, class> class Cholmod;

  template <>
  struct TypeMap<Cholmod,float> // Cholmod does not support float yet
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
   typedef Kokkos::complex<double> type;
   typedef double magnitude_type;
 };

 template <>
 struct TypeMap<Cholmod,Kokkos::complex<double> >
 {
   typedef Kokkos::complex<double> type;
   typedef double magnitude_type;
 };

 template <>
 struct TypeMap<Cholmod,std::complex<float> > // Cholmod does not support float yet
 {
   typedef Kokkos::complex<float> type;
   typedef float magnitude_type;
 };

 template <>
 struct TypeMap<Cholmod,Kokkos::complex<float> > // Cholmod does not support float yet
 {
   typedef Kokkos::complex<float> type;
   typedef float magnitude_type;
 };

#endif  // HAVE_TEUCHOS_COMPLEX

  /* \endcond Choldmod_type_specializations */

} // end namespace Amesos

#endif  // AMESOS2_CHOLMOD_TYPEMAP_HPP
