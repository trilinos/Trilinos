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
