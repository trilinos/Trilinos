// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_DYNAMICSTORAGE_HPP
#define SACADO_FAD_DYNAMICSTORAGE_HPP

#include "Sacado_ConfigDefs.h"

#include "Sacado_Fad_Exp_DynamicStorage.hpp"

namespace Sacado {
  namespace Fad {

    template <typename T, typename U = T>
    using DynamicStorage = Exp::DynamicStorage<T,U>;

  }
}


#endif // SACADO_FAD_DYNAMICSTORAGE_HPP
