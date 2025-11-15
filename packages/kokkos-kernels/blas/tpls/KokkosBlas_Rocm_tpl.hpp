// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBLAS_ROCM_TPL_HPP_
#define KOKKOSBLAS_ROCM_TPL_HPP_

#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

RocBlasSingleton::RocBlasSingleton() { KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_create_handle(&handle)); }
RocBlasSingleton::~RocBlasSingleton() { KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_destroy_handle(handle)); }

RocBlasSingleton& RocBlasSingleton::singleton() { return get_instance().get(); }

bool RocBlasSingleton::is_initialized() { return get_instance().is_initialized(); }

KokkosKernels::Impl::Singleton<RocBlasSingleton>& RocBlasSingleton::get_instance() {
  static KokkosKernels::Impl::Singleton<RocBlasSingleton> s;
  return s;
}

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // defined (KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)

#endif  // KOKKOSBLAS_ROCM_TPL_HPP_
