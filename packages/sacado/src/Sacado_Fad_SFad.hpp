// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
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

#ifndef SACADO_FAD_SFAD_HPP
#define SACADO_FAD_SFAD_HPP

#include "Sacado_ConfigDefs.h"

#include "Sacado_Fad_Exp_SFad.hpp"

namespace Sacado {
  namespace Fad {
    template <typename T, int Num>
    using SFad = Exp::GeneralFad< Exp::StaticFixedStorage<T,Num> >;
  }
}

#include "Sacado_Fad_ViewFad.hpp"

#endif // SACADO_FAD_SFAD_HPP
