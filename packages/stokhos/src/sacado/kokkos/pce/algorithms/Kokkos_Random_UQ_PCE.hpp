// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_RANDOM_UQ_PCE_HPP
#define KOKKOS_RANDOM_UQ_PCE_HPP

#include "Stokhos_ConfigDefs.h"
#if defined(HAVE_STOKHOS_KOKKOS)

#include "Sacado_UQ_PCE.hpp"
#include "Kokkos_View_UQ_PCE.hpp"
#include "Kokkos_Random.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos random functions for Sacado::UQ::PCE scalar type
//----------------------------------------------------------------------------

namespace Kokkos {

  template<class Generator, class Storage>
  struct rand<Generator,Sacado::UQ::PCE<Storage> > {
    typedef Sacado::UQ::PCE<Storage> Scalar;
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
      return BaseRand::draw(gen, range.coeff(0));
    }

    KOKKOS_INLINE_FUNCTION
    static Scalar draw(Generator& gen, const Scalar& start, const Scalar& end) {
      return BaseRand::draw(gen, start.coeff(0), end.coeff(0));
    }
  };

  template<class S, class ... P, class RandomPool>
  void fill_random( const View<Sacado::UQ::PCE<S>**,P...>& a,
                    RandomPool g,
                    const Sacado::UQ::PCE<S>& begin,
                    const Sacado::UQ::PCE<S>& end ) {
    typedef View<Sacado::UQ::PCE<S>**,P...> Vector;
    typename Kokkos::FlatArrayType<Vector>::type a_flat = a;
    fill_random( a_flat, g, begin.fastAccessCoeff(0), end.fastAccessCoeff(0) );
  }

}

#endif

#endif /* #ifndef KOKKOS_RANDOM_UQ_PCE_HPP */
