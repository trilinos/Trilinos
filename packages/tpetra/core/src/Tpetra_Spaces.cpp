#include "Tpetra_Spaces.hpp"

#include <vector>
#include <iostream>

#define LOG_WARN(x) std::cerr << __FILE__ << ":" << __LINE__ << ": [WARN] " << x << std::endl;
#define LOG_INFO(x) std::cerr << __FILE__ << ":" << __LINE__ << ": [INFO] " << x << std::endl;

namespace Tpetra {
namespace Spaces {
namespace detail {

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
    // LOG_INFO("CUDA Priority: low:    " << cudaPriorityRange.low);
    // LOG_INFO("CUDA Priority: medium: " << cudaPriorityRange.medium);
    // LOG_INFO("CUDA Priority: high:   " << cudaPriorityRange.high);
#endif
    // LOG_INFO("Tpetra::Spaces::detail::lazy_init() done");
    initialized = true;
}

} // namespace detail
} // namespace Spaces
} // namespace Tpetra

#undef LOG_WARN
#undef LOG_INFO