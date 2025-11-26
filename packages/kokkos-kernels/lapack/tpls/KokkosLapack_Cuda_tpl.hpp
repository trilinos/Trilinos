// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSLAPACK_CUDA_TPL_HPP_
#define KOKKOSLAPACK_CUDA_TPL_HPP_

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)
#include "KokkosLapack_cusolver.hpp"

namespace KokkosLapack {
namespace Impl {

CudaLapackSingleton::CudaLapackSingleton() {
  cusolverStatus_t stat = cusolverDnCreate(&handle);
  if (stat != CUSOLVER_STATUS_SUCCESS) Kokkos::abort("CUSOLVER initialization failed\n");
}

CudaLapackSingleton::~CudaLapackSingleton() { cusolverDnDestroy(handle); }

CudaLapackSingleton& CudaLapackSingleton::singleton() { return get_instance().get(); }

bool CudaLapackSingleton::is_initialized() { return get_instance().is_initialized(); }

KokkosKernels::Impl::Singleton<CudaLapackSingleton>& CudaLapackSingleton::get_instance() {
  static KokkosKernels::Impl::Singleton<CudaLapackSingleton> s;
  return s;
}

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // defined (KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)

#endif  // KOKKOSLAPACK_CUDA_TPL_HPP_
