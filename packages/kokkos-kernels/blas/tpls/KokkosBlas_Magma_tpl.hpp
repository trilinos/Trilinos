// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBLAS_MAGMA_TPL_HPP_
#define KOKKOSBLAS_MAGMA_TPL_HPP_

#if defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)
#include <KokkosBlas_magma.hpp>

namespace KokkosBlas {
namespace Impl {

MagmaSingleton::MagmaSingleton() {
  magma_int_t stat = magma_init();
  if (stat != MAGMA_SUCCESS) Kokkos::abort("MAGMA initialization failed\n");
}

MagmaSingleton::~MagmaSingleton() { magma_finalize(); }

MagmaSingleton& MagmaSingleton::singleton() { return get_instance().get(); }

bool MagmaSingleton::is_initialized() { return get_instance().is_initialized(); }

KokkosKernels::Impl::Singleton<MagmaSingleton>& MagmaSingleton::get_instance() {
  static KokkosKernels::Impl::Singleton<MagmaSingleton> s;
  return s;
}

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)

#endif  // KOKKOSBLAS_MAGMA_TPL_HPP_
