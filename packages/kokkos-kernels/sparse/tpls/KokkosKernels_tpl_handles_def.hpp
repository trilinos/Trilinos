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

CusparseSingleton::CusparseSingleton() { KOKKOSSPARSE_IMPL_CUSPARSE_SAFE_CALL(cusparseCreate(&cusparseHandle)); }

CusparseSingleton& CusparseSingleton::singleton() {
  std::unique_ptr<CusparseSingleton>& instance = get_instance();
  if (!instance) {
    instance = std::make_unique<CusparseSingleton>();
    Kokkos::push_finalize_hook([&]() {
      KOKKOSSPARSE_IMPL_CUSPARSE_SAFE_CALL(cusparseDestroy(instance->cusparseHandle));
      instance.reset();
    });
  }
  return *instance;
}

bool CusparseSingleton::is_initialized() { return get_instance() != nullptr; }

std::unique_ptr<CusparseSingleton>& CusparseSingleton::get_instance() {
  static std::unique_ptr<CusparseSingleton> s;
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
  KOKKOSSPARSE_IMPL_ROCSPARSE_SAFE_CALL(rocsparse_create_handle(&rocsparseHandle));
}

RocsparseSingleton& RocsparseSingleton::singleton() {
  std::unique_ptr<RocsparseSingleton>& instance = get_instance();
  if (!instance) {
    instance = std::make_unique<RocsparseSingleton>();
    Kokkos::push_finalize_hook([&]() {
      KOKKOSSPARSE_IMPL_ROCSPARSE_SAFE_CALL(rocsparse_destroy_handle(instance->rocsparseHandle));
      instance.reset();
    });
  }
  return *instance;
}

bool RocsparseSingleton::is_initialized() { return get_instance() != nullptr; }

std::unique_ptr<RocsparseSingleton>& RocsparseSingleton::get_instance() {
  static std::unique_ptr<RocsparseSingleton> s;
  return s;
}

}  // namespace Impl
}  // namespace KokkosKernels
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE

#endif  // KOKKOSKERNELS_TPL_HANDLES_DEF_HPP_
