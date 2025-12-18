// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS_MAGMA_HPP_
#define KOKKOSBLAS_MAGMA_HPP_

#include "KokkosKernels_Singleton.hpp"

// If LAPACK TPL is enabled, it is preferred over magma's LAPACK
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
#include "magma_v2.h"

namespace KokkosBlas {
namespace Impl {

struct MagmaSingleton {
  MagmaSingleton();
  ~MagmaSingleton();

  static bool is_initialized();
  static MagmaSingleton& singleton();

 private:
  static KokkosKernels::Impl::Singleton<MagmaSingleton>& get_instance();
};

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_MAGMA

#endif  // KOKKOSBLAS_MAGMA_HPP_
