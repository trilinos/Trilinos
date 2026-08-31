// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_DFAD_HPP
#define SACADO_FAD_DFAD_HPP

#include "Sacado_Fad_GeneralFad.hpp"
#include "Sacado_Fad_DynamicStorage.hpp"

namespace Sacado {
  namespace Fad {

    template <typename T>
    using DFad = GeneralFad< DynamicStorage<T> >;

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_DFAD_HPP
