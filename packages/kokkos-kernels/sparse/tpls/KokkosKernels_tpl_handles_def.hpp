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

#ifndef KOKKOSKERNELS_TPL_HANDLES_DEF_HPP_
#define KOKKOSKERNELS_TPL_HANDLES_DEF_HPP_

#include "KokkosKernels_tpl_handles_decl.hpp"

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"

namespace KokkosKernels {
namespace Impl {

CusparseSingleton::CusparseSingleton() {
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreate(&cusparseHandle));

  Kokkos::push_finalize_hook([&]() { cusparseDestroy(cusparseHandle); });
}

CusparseSingleton& CusparseSingleton::singleton() {
  static CusparseSingleton s;
  return s;
}

}  // namespace Impl
}  // namespace KokkosKernels
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
#include "KokkosSparse_Utils_rocsparse.hpp"

namespace KokkosKernels {
namespace Impl {

RocsparseSingleton::RocsparseSingleton() {
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_handle(&rocsparseHandle));

  Kokkos::push_finalize_hook([&]() { KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_destroy_handle(rocsparseHandle)); });
}

RocsparseSingleton& RocsparseSingleton::singleton() {
  static RocsparseSingleton s;
  return s;
}

}  // namespace Impl
}  // namespace KokkosKernels
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE

#endif  // KOKKOSKERNELS_TPL_HANDLES_DEF_HPP_
