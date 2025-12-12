// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSLAPACK_MAGMA_TPL_HPP_
#define KOKKOSLAPACK_MAGMA_TPL_HPP_

#if defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)

#include <KokkosLapack_magma.hpp>

// If BLAS is enabled, then we will only use the
// KokkosBlas::Impl::MagmaSingleton. Don't redefine the KokkosLapack::
// version here.
#ifndef KOKKOSKERNELS_ENABLE_COMPONENT_BLAS
namespace KokkosLapack {
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
}  // namespace KokkosLapack
#endif

#endif  // defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)

#endif  // KOKKOSLAPACK_MAGMA_TPL_HPP_
