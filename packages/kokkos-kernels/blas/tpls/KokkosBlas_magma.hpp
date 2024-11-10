//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSBLAS_MAGMA_HPP_
#define KOKKOSBLAS_MAGMA_HPP_

// If LAPACK TPL is enabled, it is preferred over magma's LAPACK
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
#include "magma_v2.h"

namespace KokkosBlas {
namespace Impl {

struct MagmaSingleton {
  MagmaSingleton();

  static MagmaSingleton& singleton();
};

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_MAGMA

#endif  // KOKKOSBLAS_MAGMA_HPP_
