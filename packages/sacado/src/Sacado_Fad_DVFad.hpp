// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_DVFAD_HPP
#define SACADO_FAD_DVFAD_HPP

#include "Sacado_ConfigDefs.h"

#include "Sacado_Fad_Exp_DVFad.hpp"

namespace Sacado {
  namespace Fad {
    template <typename T>
    using DVFad = Exp::GeneralFad< Exp::VectorDynamicStorage<T> >;
  }
}

#endif // SACADO_FAD_DFAD_HPP
