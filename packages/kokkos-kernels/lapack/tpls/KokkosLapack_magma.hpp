// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
