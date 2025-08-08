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
