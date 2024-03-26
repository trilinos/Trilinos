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

#ifndef KOKKOSLAPACK_MAGMA_HPP_
#define KOKKOSLAPACK_MAGMA_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
#include "magma_v2.h"

namespace KokkosLapack {
namespace Impl {

// Declaration of the singleton for cusolver
// this is the only header that needs to be
// included when using cusolverDn.
struct MagmaSingleton {
  MagmaSingleton();

  static MagmaSingleton& singleton();
};

}  // namespace Impl
}  // namespace KokkosLapack
#endif

#endif  // KOKKOSLAPACK_MAGMA_HPP_
