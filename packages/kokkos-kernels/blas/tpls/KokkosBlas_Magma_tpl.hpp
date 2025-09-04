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
