// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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

CudaBlasSingleton::~CudaBlasSingleton() { cublasDestroy(handle); }

CudaBlasSingleton& CudaBlasSingleton::singleton() { return get_instance().get(); }

bool CudaBlasSingleton::is_initialized() { return get_instance().is_initialized(); }

KokkosKernels::Impl::Singleton<CudaBlasSingleton>& CudaBlasSingleton::get_instance() {
  static KokkosKernels::Impl::Singleton<CudaBlasSingleton> s;
  return s;
}

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // defined (KOKKOSKERNELS_ENABLE_TPL_CUBLAS)

#endif  // KOKKOSBLAS_CUDA_TPL_HPP_
