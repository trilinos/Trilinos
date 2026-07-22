// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_VIEW_UQ_PCE_HPP
#define KOKKOS_VIEW_UQ_PCE_HPP

namespace Kokkos {

  // Get global sparse 3 tensor
  template <typename cijk_type>
  cijk_type& getGlobalCijkTensor() {
    static cijk_type cijk;
    return cijk;
  }

  // Set global sparse 3 tensor
  template <typename cijk_type>
  void setGlobalCijkTensor(const cijk_type& cijk) {
    cijk_type& global_cijk = getGlobalCijkTensor<cijk_type>();
    global_cijk = cijk;
  }

}

#include "KokkosExp_View_UQ_PCE_Contiguous.hpp"

#endif /* #ifndef KOKKOS_VIEW_UQ_PCE_HPP */
