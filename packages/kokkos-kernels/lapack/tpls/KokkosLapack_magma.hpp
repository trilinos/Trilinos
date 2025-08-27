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

#include <KokkosKernels_config.h>
#include "KokkosKernels_Singleton.hpp"

#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
#include "magma_v2.h"

#ifdef KOKKOSKERNELS_ENABLE_COMPONENT_BLAS

// If both BLAS and LAPACK components are enabled,
// don't use two separate singletons for MAGMA
// (both would try to call the same magma_init/magma_finalize).
//
// Instead use the KokkosBlas::Impl::MagmaSingleton for both.
#include "KokkosBlas_magma.hpp"

namespace KokkosLapack {
namespace Impl {
using KokkosBlas::Impl::MagmaSingleton;
}
}  // namespace KokkosLapack
#else

namespace KokkosLapack {
namespace Impl {

// Declaration of the singleton for cusolver
// this is the only header that needs to be
// included when using cusolverDn.
struct MagmaSingleton {
  MagmaSingleton();
  ~MagmaSingleton();

  static MagmaSingleton& singleton();

  static bool is_initialized();

 private:
  static KokkosKernels::Impl::Singleton<MagmaSingleton>& get_instance();
};

}  // namespace Impl
}  // namespace KokkosLapack
#endif
#endif

#endif  // KOKKOSLAPACK_MAGMA_HPP_
