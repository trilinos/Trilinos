#include "Tpetra_Space.hpp"

#include <vector>
#include <iostream>

#define LOG_WARN(x) std::cerr << __FILE__ << ":" << __LINE__ << ": [WARN] " << x << std::endl;
#define LOG_INFO(x) std::cerr << __FILE__ << ":" << __LINE__ << ": [INFO] " << x << std::endl;

namespace Tpetra {



namespace detail {

/*extern*/ CudaPriorityRange cudaPriorityRange;

void lazy_init() {
    #ifdef KOKKOS_ENABLE_CUDA
    if (!cudaPriorityRange.isSet)  {
        cudaDeviceGetStreamPriorityRange(&cudaPriorityRange.low, &cudaPriorityRange.high);
        LOG_INFO("LOPRIO: " << cudaPriorityRange.low);
        LOG_INFO("HIPRIO: " << cudaPriorityRange.high);
        cudaPriorityRange.isSet = true;
    }
    #endif
}

} // namespace detail




void drop_exec_spaces() {
    LOG_INFO("drop_exec_spaces");
#ifdef KOKKOS_ENABLE_CUDA
    for (int i = 0; i < static_cast<int>(Priority::NUM_LEVELS); ++i) {
        detail::cudaSpaces[i].clear();
    }
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

}

#undef LOG_WARN
#undef LOG_INFO