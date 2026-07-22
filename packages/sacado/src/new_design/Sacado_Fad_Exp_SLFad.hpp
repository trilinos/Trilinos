// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_SLFAD_HPP
#define SACADO_FAD_EXP_SLFAD_HPP

#include "Sacado_Fad_Exp_GeneralFad.hpp"
#include "Sacado_Fad_Exp_StaticStorage.hpp"

namespace Sacado {

  namespace Fad {
  namespace Exp {

    template <typename T, int Num>
    using SLFad = GeneralFad< StaticStorage<T,Num> >;

  } // namespace Exp
  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_EXP_SLFAD_HPP
