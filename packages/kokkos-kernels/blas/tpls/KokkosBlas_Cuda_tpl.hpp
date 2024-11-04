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
#ifndef KOKKOSBLAS_CUDA_TPL_HPP_
#define KOKKOSBLAS_CUDA_TPL_HPP_

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

CudaBlasSingleton::CudaBlasSingleton() {
  cublasStatus_t stat = cublasCreate(&handle);
  if (stat != CUBLAS_STATUS_SUCCESS) Kokkos::abort("CUBLAS initialization failed\n");
}

CudaBlasSingleton& CudaBlasSingleton::singleton() {
  std::unique_ptr<CudaBlasSingleton>& instance = get_instance();
  if (!instance) {
    instance = std::make_unique<CudaBlasSingleton>();
    Kokkos::push_finalize_hook([&]() {
      cublasDestroy(instance->handle);
      instance.reset();
    });
  }
  return *instance;
}

bool CudaBlasSingleton::is_initialized() { return get_instance() != nullptr; }

std::unique_ptr<CudaBlasSingleton>& CudaBlasSingleton::get_instance() {
  static std::unique_ptr<CudaBlasSingleton> s;
  return s;
}

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // defined (KOKKOSKERNELS_ENABLE_TPL_CUBLAS)

#endif  // KOKKOSBLAS_CUDA_TPL_HPP_
