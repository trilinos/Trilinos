// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_IFPACK2_MP_VECTOR_HPP
#define STOKHOS_IFPACK2_MP_VECTOR_HPP

// This header file should be included whenever compiling any Ifpack2
// code with Stokhos scalar types

// MP includes and specializations
#include "Stokhos_Tpetra_MP_Vector.hpp"

// Specialization of LocalReciprocalThreshold functor in
// Ifpack2_Details_Chebyshev_def.hpp
namespace Ifpack2 {
namespace Details {

template<class XV, class SizeType>
struct LocalReciprocalThreshold;

template<class S, class ... P, class SizeType>
struct LocalReciprocalThreshold<
  Kokkos::View< Sacado::MP::Vector<S>*,P...>, SizeType > {
  typedef Kokkos::View<Sacado::MP::Vector<S>*,P...> XV;

  static void
  compute (const XV& X,
           const typename XV::non_const_value_type& minVal)
  {
    if (!Sacado::is_constant(minVal)) {
      Kokkos::Impl::raise_error(
        "LocalReciprocalThreshold not implemented for non-constant minVal");
    }

    typedef typename Kokkos::FlatArrayType<XV>::type Flat_XV;
    Flat_XV flat_X = X;
    LocalReciprocalThreshold< Flat_XV, SizeType >::compute( flat_X,
                                                            minVal.coeff(0) );
  }
};

}
}

#endif // STOKHOS_IFPACK2_MP_VECTOR_HPP
