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

#include "Sacado_ConfigDefs.h"

#ifdef SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#include "Sacado_Fad_Exp_DFad.hpp"

namespace Sacado {
  namespace Fad {
    template <typename T>
    using DFad = Exp::GeneralFad< Exp::DynamicStorage<T> >;
  }
}

#else

#include "Sacado_Fad_GeneralFadExpr.hpp"
#include "Sacado_Fad_DFadTraits.hpp"
#include "Sacado_Fad_DynamicStorage.hpp"

#define FAD_NS Fad
#include "Sacado_Fad_DFad_tmpl.hpp"
#undef FAD_NS

#endif // SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#include "Sacado_Fad_ViewFad.hpp"

#endif // SACADO_FAD_DFAD_HPP
