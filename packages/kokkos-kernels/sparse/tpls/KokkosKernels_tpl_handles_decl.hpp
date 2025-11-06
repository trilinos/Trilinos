// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSKERNELS_TPL_HANDLES_DECL_HPP_
#define KOKKOSKERNELS_TPL_HANDLES_DECL_HPP_

#include "KokkosBlas_tpl_spec.hpp"
#include "KokkosKernels_Singleton.hpp"

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "KokkosSparse_Utils_cusparse.hpp"

namespace KokkosKernels {
namespace Impl {

struct CusparseSingleton {
  cusparseHandle_t cusparseHandle;

  CusparseSingleton();
  ~CusparseSingleton();

  static bool is_initialized();
  static CusparseSingleton& singleton();

 private:
  static KokkosKernels::Impl::Singleton<CusparseSingleton>& get_instance();
};

}  // namespace Impl
}  // namespace KokkosKernels
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
#include <rocsparse/rocsparse.h>

namespace KokkosKernels {
namespace Impl {

struct RocsparseSingleton {
  rocsparse_handle rocsparseHandle;

  RocsparseSingleton();
  ~RocsparseSingleton();

  static bool is_initialized();
  static RocsparseSingleton& singleton();

 private:
  static KokkosKernels::Impl::Singleton<RocsparseSingleton>& get_instance();
};

}  // namespace Impl
}  // namespace KokkosKernels
#endif

#endif  // KOKKOSKERNELS_TPL_HANDLES_DECL_HPP_
