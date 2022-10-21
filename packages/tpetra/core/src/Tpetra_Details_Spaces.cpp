#include "Tpetra_Details_Spaces.hpp"

#include <vector>
#include <iostream>

namespace Tpetra {
namespace Details {
namespace Spaces {

#ifdef KOKKOS_ENABLE_CUDA
/*extern*/ CudaPriorityRange cudaPriorityRange;
/*extern*/ cudaEvent_t execSpaceWaitEvent; // see exec_space_wait
#endif

/*extern*/ bool initialized = false;

void lazy_init() {
    if (initialized) { return; }
#ifdef KOKKOS_ENABLE_CUDA
    CUDA_RUNTIME(cudaEventCreateWithFlags(&execSpaceWaitEvent, cudaEventDisableTiming));
    CUDA_RUNTIME(cudaDeviceGetStreamPriorityRange(&cudaPriorityRange.low, &cudaPriorityRange.high));
#endif
    initialized = true;
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

} // namespace Spaces
} // namespace Details
} // namespace Tpetra
