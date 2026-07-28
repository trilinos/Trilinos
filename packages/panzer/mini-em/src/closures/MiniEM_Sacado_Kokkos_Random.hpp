// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_SACADO_KOKKOS_RANDOM_HPP
#define MINIEM_SACADO_KOKKOS_RANDOM_HPP


#include "Kokkos_Random.hpp"


namespace Kokkos {

  /** \brief Specializes Kokkos::rand for panzer::Traits::FadType,
    * drawing random values from the underlying base scalar's generator
    * and wrapping them as (non-differentiated) Fad values. Needed so
    * Kokkos random number generation can be used directly in AD
    * (Jacobian/Tangent) evaluation types.
    */
  template<class Generator>
  struct rand<Generator,panzer::Traits::FadType> {
    typedef panzer::Traits::FadType Scalar;
    typedef typename Scalar::value_type BaseScalar;
    typedef rand<Generator,BaseScalar> BaseRand;

    KOKKOS_INLINE_FUNCTION
    static Scalar max() { return BaseRand::max(); }

    KOKKOS_INLINE_FUNCTION
    static Scalar draw(Generator& gen) {
      return BaseRand::draw(gen);
    }

    KOKKOS_INLINE_FUNCTION
    static Scalar draw(Generator& gen, const Scalar& range) {
      return BaseRand::draw(gen, range.val());
    }

    KOKKOS_INLINE_FUNCTION
    static Scalar draw(Generator& gen, const Scalar& start, const Scalar& end) {
      return BaseRand::draw(gen, start.val(), end.val());
    }
  };

}



#endif
