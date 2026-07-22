// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_VIEW_MP_VECTOR_HPP
#define KOKKOS_VIEW_MP_VECTOR_HPP

namespace Kokkos {

  // Global vector size to be used in code that hasn't been refactored
  // into passing the vector size (e.g., Tpetra)
  extern unsigned global_sacado_mp_vector_size;

}

#include "KokkosExp_View_MP_Vector_Contiguous.hpp"

#endif /* #ifndef KOKKOS_VIEW_MP_VECTOR_HPP */
