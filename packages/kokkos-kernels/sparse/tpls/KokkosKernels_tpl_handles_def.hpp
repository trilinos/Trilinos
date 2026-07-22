// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSKERNELS_TPL_HANDLES_DEF_HPP_
#define KOKKOSKERNELS_TPL_HANDLES_DEF_HPP_

#include "KokkosKernels_tpl_handles_decl.hpp"

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"

namespace KokkosKernels {
namespace Impl {

CusparseSingleton::CusparseSingleton() { KOKKOSSPARSE_IMPL_CUSPARSE_SAFE_CALL(cusparseCreate(&cusparseHandle)); }
CusparseSingleton::~CusparseSingleton() { KOKKOSSPARSE_IMPL_CUSPARSE_SAFE_CALL(cusparseDestroy(cusparseHandle)); }

CusparseSingleton& CusparseSingleton::singleton() { return get_instance().get(); }

bool CusparseSingleton::is_initialized() { return get_instance().is_initialized(); }

KokkosKernels::Impl::Singleton<CusparseSingleton>& CusparseSingleton::get_instance() {
  static KokkosKernels::Impl::Singleton<CusparseSingleton> s;
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

RocsparseSingleton::~RocsparseSingleton() {
  KOKKOSSPARSE_IMPL_ROCSPARSE_SAFE_CALL(rocsparse_destroy_handle(rocsparseHandle));
}

RocsparseSingleton& RocsparseSingleton::singleton() { return get_instance().get(); }

bool RocsparseSingleton::is_initialized() { return get_instance().is_initialized(); }

KokkosKernels::Impl::Singleton<RocsparseSingleton>& RocsparseSingleton::get_instance() {
  static KokkosKernels::Impl::Singleton<RocsparseSingleton> s;
  return s;
}

}  // namespace Impl
}  // namespace KokkosKernels
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE

#endif  // KOKKOSKERNELS_TPL_HANDLES_DEF_HPP_
