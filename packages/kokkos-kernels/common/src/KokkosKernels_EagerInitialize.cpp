// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "KokkosKernels_EagerInitialize.hpp"
#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"

// Include the minimal set of headers that declare all TPL singletons
#ifdef KOKKOSKERNELS_ENABLE_COMPONENT_BLAS
#include "KokkosBlas_tpl_spec.hpp"  //cuBLAS, rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
#include "KokkosBlas_magma.hpp"
#endif
#endif

#ifdef KOKKOSKERNELS_ENABLE_COMPONENT_SPARSE
// note: this file declares both cuSPARSE and rocSPARSE singletons
#include "KokkosKernels_tpl_handles_decl.hpp"
#endif

#ifdef KOKKOSKERNELS_ENABLE_COMPONENT_LAPACK
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSOLVER
#include "KokkosLapack_cusolver.hpp"
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
#include "KokkosLapack_magma.hpp"
#endif
#endif

namespace KokkosKernels {
void eager_initialize() {
  if (!Kokkos::is_initialized()) {
    throw std::runtime_error("Kokkos::intialize must be called before KokkosKernels::eager_initialize");
  }
#ifdef KOKKOSKERNELS_ENABLE_COMPONENT_BLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
  (void)KokkosBlas::Impl::CudaBlasSingleton::singleton();
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
  (void)KokkosBlas::Impl::RocBlasSingleton::singleton();
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  (void)KokkosBlas::Impl::MagmaSingleton::singleton();
#endif
#endif

#ifdef KOKKOSKERNELS_ENABLE_COMPONENT_SPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  (void)KokkosKernels::Impl::CusparseSingleton::singleton();
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
  (void)KokkosKernels::Impl::RocsparseSingleton::singleton();
#endif
#endif

#ifdef KOKKOSKERNELS_ENABLE_COMPONENT_LAPACK
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSOLVER
  (void)KokkosLapack::Impl::CudaLapackSingleton::singleton();
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  (void)KokkosLapack::Impl::MagmaSingleton::singleton();
#endif
#endif
}
}  // namespace KokkosKernels
