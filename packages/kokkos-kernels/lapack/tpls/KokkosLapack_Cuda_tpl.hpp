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
#ifndef KOKKOSLAPACK_CUDA_TPL_HPP_
#define KOKKOSLAPACK_CUDA_TPL_HPP_

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)
#include "KokkosLapack_cusolver.hpp"

namespace KokkosLapack {
namespace Impl {

CudaLapackSingleton::CudaLapackSingleton() {
  cusolverStatus_t stat = cusolverDnCreate(&handle);
  if (stat != CUSOLVER_STATUS_SUCCESS) Kokkos::abort("CUSOLVER initialization failed\n");

  Kokkos::push_finalize_hook([&]() { cusolverDnDestroy(handle); });
}

CudaLapackSingleton& CudaLapackSingleton::singleton() {
  static CudaLapackSingleton s;
  return s;
}

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // defined (KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)

#endif  // KOKKOSLAPACK_CUDA_TPL_HPP_
