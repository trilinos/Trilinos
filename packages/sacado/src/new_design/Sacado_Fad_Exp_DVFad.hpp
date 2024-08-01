// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_DVFAD_HPP
#define SACADO_FAD_EXP_DVFAD_HPP

#include "Sacado_Fad_Exp_GeneralFad.hpp"
#include "Sacado_Fad_Exp_VectorDynamicStorage.hpp"

namespace Sacado {

  namespace Fad {
  namespace Exp {

    template <typename T>
    using DVFad = GeneralFad< VectorDynamicStorage<T> >;

  } // namespace Exp
  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_EXP_DVFAD_HPP
