// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_VIEWSTORAGE_HPP
#define SACADO_FAD_VIEWSTORAGE_HPP

#include "Sacado_ConfigDefs.h"

#include "Sacado_Fad_Exp_ViewStorage.hpp"

namespace Sacado {
  namespace Fad {

    template <typename T, unsigned static_length, unsigned static_stride,
              typename U>
    using ViewStorage = Exp::ViewStorage<T,static_length,static_stride,U>;

  }
}


#endif // SACADO_FAD_VIEWSTORAGE_HPP
