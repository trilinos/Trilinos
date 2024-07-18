// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
//
// The forward-mode AD classes in Sacado are a derivative work of the
// expression template classes in the Fad package by Nicolas Di Cesare.
// The following banner is included in the original Fad source code:
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses,
//         templates : new C++ techniques
//            for scientific computing
//
//********************************************************
//
//  A short implementation ( not all operators and
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
// @HEADER

#ifndef SACADO_FAD_SLFAD_MP_VECTOR_HPP
#define SACADO_FAD_SLFAD_MP_VECTOR_HPP

#include "Sacado_Fad_SLFad.hpp"
#include "Sacado_Fad_ExprSpec_MP_Vector.hpp"
#include "Sacado_Fad_GeneralFad_MP_Vector.hpp"
#include "Sacado_Fad_Ops_MP_Vector.hpp"

namespace Stokhos {
  template <typename Ord, typename Val, int Num, typename Dev>
  class StaticFixedStorage;
}

namespace Sacado {

  namespace MP {
    template <typename S> class Vector;
  }

  namespace Fad {

    template <typename Ord, typename Val, int VecNum, typename Dev, int Num>
    struct ExprSpec< SLFad< Sacado::MP::Vector< Stokhos::StaticFixedStorage<Ord,Val,VecNum,Dev> >, Num > > {
      typedef ExprSpecMPVector type;
    };

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_SLFAD_MP_VECTOR_HPP
