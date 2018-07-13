// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
