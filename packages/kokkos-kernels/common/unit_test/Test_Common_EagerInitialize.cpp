// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include "Kokkos_Core.hpp"
#include "KokkosKernels_config.h"
#include "KokkosKernels_EagerInitialize.hpp"

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

// Count the number of singletons which are currently initialized,
// and the numInitialized number of singleton classes that are currently enabled
// (based on which TPLs and components were enabled at configure-time)
void countSingletons(int& numInitialized, int& numEnabled) {
  numInitialized = 0;
  numEnabled     = 0;
#ifdef KOKKOSKERNELS_ENABLE_COMPONENT_BLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
  numEnabled++;
  if (KokkosBlas::Impl::CudaBlasSingleton::is_initialized()) numInitialized++;
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
  numEnabled++;
  if (KokkosBlas::Impl::RocBlasSingleton::is_initialized()) numInitialized++;
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  numEnabled++;
  if (KokkosBlas::Impl::MagmaSingleton::is_initialized()) numInitialized++;
#endif
#endif

#ifdef KOKKOSKERNELS_ENABLE_COMPONENT_SPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  numEnabled++;
  if (KokkosKernels::Impl::CusparseSingleton::is_initialized()) numInitialized++;
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
  numEnabled++;
  if (KokkosKernels::Impl::RocsparseSingleton::is_initialized()) numInitialized++;
#endif
#endif

#ifdef KOKKOSKERNELS_ENABLE_COMPONENT_LAPACK
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSOLVER
  numEnabled++;
  if (KokkosLapack::Impl::CudaLapackSingleton::is_initialized()) numInitialized++;
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  numEnabled++;
  if (KokkosLapack::Impl::MagmaSingleton::is_initialized()) numInitialized++;
#endif
#endif
}

int main() {
  int numInitialized, numEnabled;
  Kokkos::initialize();
  {
    // Check that no singletons are already initialized.
    countSingletons(numInitialized, numEnabled);
    if (numInitialized != 0)
      throw std::runtime_error("At least one singleton was initialized before it should have been");
    KokkosKernels::eager_initialize();
    // Check that all singletons are now initialized.
    countSingletons(numInitialized, numEnabled);
    std::cout << "Kokkos::eager_initialize() set up " << numInitialized << " of " << numEnabled << " TPL singletons.\n";
    if (numInitialized != numEnabled)
      throw std::runtime_error("At least one singleton was not initialized by eager_initialize()");
  }
  Kokkos::finalize();
  // Finally, make sure that all singletons were finalized during Kokkos::finalize().
  countSingletons(numInitialized, numEnabled);
  if (numInitialized != 0)
    throw std::runtime_error("At least one singleton was not correctly finalized by Kokkos::finalize()");
  return 0;
}
