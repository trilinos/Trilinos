// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_ZERO_SENSITIVITIES_HPP
#define PANZER_ZERO_SENSITIVITIES_HPP

#include <iostream>
#include "Panzer_Traits.hpp" // for scalar types
#include "Kokkos_Core.hpp"

namespace panzer {

  //! Zero out AD types so that we can drop sensitivity contributions
  // template<typename ScalarT> 
  // void zeroSensitivities(ScalarT& s);

  //! Specialization for Residual
  KOKKOS_INLINE_FUNCTION
  void zeroSensitivities(panzer::Traits::RealType& /* s */) {}

  //! Specialization for Fad type Jacobian
  KOKKOS_INLINE_FUNCTION
  void zeroSensitivities(panzer::Traits::FadType& s) 
  {
    s.zero();
  }

  //! Specialization for Fad type Hessian
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  KOKKOS_INLINE_FUNCTION
  void zeroSensitivities(panzer::Traits::HessianType& s) 
  {
    s.val().zero(); // is this right...? What does zero sensitivities mean here?
    s.zero();
  }
#endif

  //! Specialization for Kokkos View reference type Jacobian
  /*
  template<typename S, unsigned l, unsigned s, typename b>
  inline 
  void 
  zeroSensitivities(Sacado::Fad::ViewFad<S,l,s,b> v) 
  {
    v.zero();
  }
  */
  template <typename T>
  KOKKOS_INLINE_FUNCTION
  typename Sacado::mpl::enable_if< Sacado::IsView<T> >::type
  zeroSensitivities(T x) {
    //std::cout << "zeroSensitivities - ViewFad" << std::endl;
      x.zero();
  }

}

#endif
