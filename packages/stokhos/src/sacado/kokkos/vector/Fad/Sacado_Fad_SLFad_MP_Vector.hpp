// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
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
