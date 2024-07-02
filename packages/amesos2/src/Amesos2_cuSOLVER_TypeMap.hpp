// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_CUSOLVER_TYPEMAP_HPP
#define AMESOS2_CUSOLVER_TYPEMAP_HPP

#include <functional>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include "Amesos2_TypeMap.hpp"

namespace Amesos2 {

  // forward declaration due to circular reference
  template <class, class> class cuSOLVER;

  template <>
  struct TypeMap<cuSOLVER,float>
  {
    typedef float type;
    typedef float magnitude_type;
  };

  template <>
  struct TypeMap<cuSOLVER,double>
  {
    typedef double type;
    typedef double magnitude_type;
  };

#ifdef HAVE_TEUCHOS_COMPLEX

  template <>
  struct TypeMap<cuSOLVER,std::complex<float> >
  {
    typedef Kokkos::complex<float> type;
    typedef float magnitude_type;
  };

  template <>
  struct TypeMap<cuSOLVER,Kokkos::complex<float> >
  {
    typedef Kokkos::complex<float> type;
    typedef float magnitude_type;
  };

  template <>
  struct TypeMap<cuSOLVER,std::complex<double> >
  {
    typedef Kokkos::complex<double> type;
    typedef double magnitude_type;
  };

  template <>
  struct TypeMap<cuSOLVER,Kokkos::complex<double> >
  {
    typedef Kokkos::complex<double> type;
    typedef double magnitude_type;
  };

#endif  // HAVE_TEUCHOS_COMPLEX

} // end namespace Amesos

#endif  // AMESOS2_CUSOLVER_TYPEMAP_HPP
