#include "Tpetra_Spaces.hpp"

#include <vector>
#include <iostream>

#define LOG_WARN(x) std::cerr << __FILE__ << ":" << __LINE__ << ": [WARN] " << x << std::endl;
#define LOG_INFO(x) std::cerr << __FILE__ << ":" << __LINE__ << ": [INFO] " << x << std::endl;

namespace Tpetra {

namespace Spaces {

namespace detail {

// Tpetra's managed spaces
#ifdef KOKKOS_ENABLE_CUDA
/*extern*/ std::vector<Kokkos::Cuda> cudaSpaces[static_cast<int>(Priority::NUM_LEVELS)];
/*extern*/ std::unordered_map<Kokkos::Cuda, cudaStream_t> cudaStreams; // track for optimized inter-space sync
#endif
#ifdef KOKKOS_ENABLE_SERIAL
/*extern*/ std::vector<Kokkos::Serial> serialSpaces[static_cast<int>(Priority::NUM_LEVELS)];
#endif
#ifdef KOKKOS_ENABLE_OPENMP
/*extern*/ std::vector<Kokkos::OpenMP> openMPSpaces[static_cast<int>(Priority::NUM_LEVELS)];
#endif

#ifdef KOKKOS_ENABLE_CUDA
/*extern*/ CudaPriorityRange cudaPriorityRange;
/*extern*/ cudaEvent_t execSpaceWaitEvent; // see exec_space_wait
#endif

void initialize() {
#ifdef KOKKOS_ENABLE_CUDA
    CUDA_RUNTIME(cudaEventCreateWithFlags(&execSpaceWaitEvent, cudaEventDisableTiming));
    CUDA_RUNTIME(cudaDeviceGetStreamPriorityRange(&cudaPriorityRange.low, &cudaPriorityRange.high));
    LOG_INFO("CUDA Priority: low:    " << cudaPriorityRange.low);
    LOG_INFO("CUDA Priority: medium: " << cudaPriorityRange.medium);
    LOG_INFO("CUDA Priority: high:   " << cudaPriorityRange.high);
#endif
}

void finalize() {
#ifdef KOKKOS_ENABLE_CUDA
    CUDA_RUNTIME(cudaEventDestroy(detail::execSpaceWaitEvent));
    for (int i = 0; i < static_cast<int>(Priority::NUM_LEVELS); ++i) {
        detail::cudaSpaces[i].clear();
    }
    detail::cudaStreams.clear();
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    for (int i = 0; i < static_cast<int>(Priority::NUM_LEVELS); ++i) {
        detail::openMpSpaces[i].clear();
    }
#endif
#ifdef KOKKOS_ENABLE_SERIAL
    for (int i = 0; i < static_cast<int>(Priority::NUM_LEVELS); ++i) {
        detail::serialSpaces[i].clear();
    }
#endif
}

} // namespace detail
} // namespace Spaces
} // namespace Tpetra

#undef LOG_WARN
#undef LOG_INFO