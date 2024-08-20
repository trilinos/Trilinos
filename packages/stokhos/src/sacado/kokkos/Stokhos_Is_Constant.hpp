// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_IS_CONSTANT_HPP
#define STOKHOS_IS_CONSTANT_HPP

#include "Kokkos_Macros.hpp"

namespace Sacado {

  // Simple function to determine whether a UQ type (Ensemble, PCE) is
  // constant.  Defaults to true for all types.
  template <typename T>
  KOKKOS_INLINE_FUNCTION
  bool is_constant(const T& x) {
    return true;
  }

}

#endif // STOKHOS_IS_CONSTANT_HPP
