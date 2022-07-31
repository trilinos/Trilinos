#ifndef KOKKOSBLAS_ROCM_TPL_HPP_
#define KOKKOSBLAS_ROCM_TPL_HPP_

#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

RocBlasSingleton::RocBlasSingleton() {
  KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_create_handle(&handle));

  Kokkos::push_finalize_hook(
      [&]() { KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_destroy_handle(handle)); });
}

RocBlasSingleton& RocBlasSingleton::singleton() {
  static RocBlasSingleton s;
  return s;
}

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // defined (KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)

#endif  // KOKKOSBLAS_ROCM_TPL_HPP_
