// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_RANDOM_MP_VECTOR_HPP
#define KOKKOS_RANDOM_MP_VECTOR_HPP

#include "Stokhos_ConfigDefs.h"
#if defined(HAVE_STOKHOS_KOKKOS)

#include "Sacado_MP_Vector.hpp"
#include "Kokkos_View_MP_Vector.hpp"
#include "Kokkos_Random.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos random functions for Sacado::MP::Vector scalar type
//----------------------------------------------------------------------------

namespace Kokkos {

  template<class Generator, class Storage>
  struct rand<Generator,Sacado::MP::Vector<Storage> > {
    typedef Sacado::MP::Vector<Storage> Scalar;
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

}

#endif

#endif /* #ifndef KOKKOS_RANDOM_MP_VECTOR_HPP */
