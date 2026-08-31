// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_VECTORDYNAMICSTORAGE_HPP
#define SACADO_FAD_VECTORDYNAMICSTORAGE_HPP

#include "Sacado_ConfigDefs.h"

#include "Sacado_Fad_Exp_VectorDynamicStorage.hpp"

namespace Sacado {
  namespace Fad {

    template <typename T, typename U = T>
    using VectorDynamicStorage = Exp::VectorDynamicStorage<T,U>;

  }
}


#endif // SACADO_FAD_VECTORDYNAMICSTORAGE_HPP
