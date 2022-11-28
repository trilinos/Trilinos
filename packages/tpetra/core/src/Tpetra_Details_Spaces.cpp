#include "Tpetra_Details_Spaces.hpp"

#include <vector>
#include <sstream>

namespace Tpetra {
namespace Details {
namespace Spaces {

#ifdef KOKKOS_ENABLE_CUDA

// cuda has default stream priority 0
CudaInfo::CudaInfo() : initialized_(false), mediumPrio_(0) {}

/*extern*/ CudaInfo cudaInfo;
#endif

void lazy_init() {
#ifdef KOKKOS_ENABLE_CUDA
    if (!cudaInfo.initialized_) {
        cudaInfo.initialized_ = true;
        CUDA_RUNTIME(cudaEventCreateWithFlags(&cudaInfo.execSpaceWaitEvent_, cudaEventDisableTiming));
        CUDA_RUNTIME(cudaDeviceGetStreamPriorityRange(&cudaInfo.lowPrio_, &cudaInfo.highPrio_));

        // We expect 
        //   medium priority should be 0
        //   lower numbers to be higher priorities
        //   low is at least as good as medium
        //   medium is at least as good as high
        if (!(cudaInfo.lowPrio_ >= cudaInfo.mediumPrio_ && cudaInfo.mediumPrio_ >= cudaInfo.highPrio_)) {
            std::stringstream ss;
            ss << "CUDA stream priority does not follow assumptions."
               << " low=" << cudaInfo.lowPrio_
               << " medium=" << cudaInfo.mediumPrio_
               << " high=" << cudaInfo.highPrio_
               << " Please report this to the Tptera developers.";
            throw std::runtime_error(ss.str());
        }
    }
#endif
}

#ifdef KOKKOS_ENABLE_CUDA
/*extern*/ InstanceLifetimeManager<Kokkos::Cuda> cudaSpaces;
#endif
#ifdef KOKKOS_ENABLE_SERIAL
/*extern*/ InstanceLifetimeManager<Kokkos::Serial> serialSpaces;
#endif
#ifdef KOKKOS_ENABLE_OPENMP
/*extern*/ InstanceLifetimeManager<Kokkos::OpenMP> openMPSpaces;
#endif
#ifdef KOKKOS_ENABLE_HIP
/*extern*/ InstanceLifetimeManager<Kokkos::OpenMP> HIPSpaces;
#endif
#ifdef KOKKOS_ENABLE_SYCL
/*extern*/ InstanceLifetimeManager<Kokkos::OpenMP> SYCLSpaces;
#endif

} // namespace Spaces
} // namespace Details
} // namespace Tpetra
