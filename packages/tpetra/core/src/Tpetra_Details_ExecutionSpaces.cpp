// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_ExecutionSpaces.hpp"

#include <sstream>
#include <vector>

namespace Tpetra {
namespace Details {
namespace Spaces {

#if defined(KOKKOS_ENABLE_CUDA)
// cuda has default stream priority 0
CudaInfo::CudaInfo() : initialized_(false), mediumPrio_(0) {}
#endif

void lazy_init() {
#if defined(KOKKOS_ENABLE_CUDA)
  if (!cudaInfo.initialized_) {
    cudaInfo.initialized_ = true;
    TPETRA_DETAILS_SPACES_CUDA_RUNTIME(cudaEventCreateWithFlags(
        &cudaInfo.execSpaceWaitEvent_, cudaEventDisableTiming));
    TPETRA_DETAILS_SPACES_CUDA_RUNTIME(cudaDeviceGetStreamPriorityRange(
        &cudaInfo.lowPrio_, &cudaInfo.highPrio_));

    // We expect
    //   medium priority should be 0
    //   lower numbers to be higher priorities
    //   low is at least as good as medium
    //   medium is at least as good as high
    if (!(cudaInfo.lowPrio_ >= cudaInfo.mediumPrio_ &&
          cudaInfo.mediumPrio_ >= cudaInfo.highPrio_)) {
      std::stringstream ss;
      ss << "CUDA stream priority does not follow assumptions."
         << " low=" << cudaInfo.lowPrio_ << " medium=" << cudaInfo.mediumPrio_
         << " high=" << cudaInfo.highPrio_
         << " Please report this to the Tpetra developers.";
      throw std::runtime_error(ss.str());
    }
  }
#endif
}

/* -----------------------------
    Space Management Singletons
   -----------------------------*/

#if defined(KOKKOS_ENABLE_CUDA)
/*extern*/ InstanceLifetimeManager<Kokkos::Cuda> cudaSpaces;
/*extern*/ CudaInfo cudaInfo;
#endif
#ifdef KOKKOS_ENABLE_SERIAL
/*extern*/ InstanceLifetimeManager<Kokkos::Serial> serialSpaces;
#endif
#ifdef KOKKOS_ENABLE_OPENMP
/*extern*/ InstanceLifetimeManager<Kokkos::OpenMP> openMPSpaces;
#endif
#ifdef KOKKOS_ENABLE_HIP
/*extern*/ InstanceLifetimeManager<Kokkos::HIP> HIPSpaces;
#endif
#ifdef KOKKOS_ENABLE_SYCL
/*extern*/ InstanceLifetimeManager<Kokkos::Experimental::SYCL> SYCLSpaces;
#endif

} // namespace Spaces
} // namespace Details
} // namespace Tpetra
