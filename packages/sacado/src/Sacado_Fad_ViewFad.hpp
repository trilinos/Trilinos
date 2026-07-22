// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_VIEWFAD_HPP
#define SACADO_FAD_VIEWFAD_HPP

#include "Sacado_ConfigDefs.h"

#ifdef SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#include "Sacado_Fad_Exp_ViewFad.hpp"

namespace Sacado {
  namespace Fad {
    template <typename T, unsigned static_length, unsigned static_stride,
              typename U>
    using ViewFad =
      Exp::GeneralFad< Exp::ViewStorage<T,static_length,static_stride,U> >;
  }
}

#else

#include "Sacado_Fad_GeneralFadExpr.hpp"
#include "Sacado_Fad_ViewFadTraits.hpp"
#include "Sacado_Fad_ViewStorage.hpp"

#define FAD_NS Fad
#include "Sacado_Fad_ViewFad_tmpl.hpp"
#undef FAD_NS

#endif // SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#endif // SACADO_FAD_VIEWFAD_HPP
