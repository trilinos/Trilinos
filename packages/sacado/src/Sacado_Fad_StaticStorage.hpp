// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_STATICSTORAGE_HPP
#define SACADO_FAD_STATICSTORAGE_HPP

#include "Sacado_ConfigDefs.h"

#include "Sacado_Fad_Exp_StaticStorage.hpp"

namespace Sacado {
  namespace Fad {

    template <typename T, int N>
    using StaticStorage = Exp::StaticStorage<T,N>;

  }
}


#endif // SACADO_FAD_STATICSTORAGE_HPP
