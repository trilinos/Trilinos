#include "Tpetra_SpaceManager.hpp"

namespace Tpetra {

SpaceManager::~SpaceManager() {
#ifdef KOKKOS_ENABLE_CUDA
    for (int i = 0; i < static_cast<int>(Spaces::Priority::NUM_LEVELS); ++i) {
        cudaSpaces[i].clear();
    }
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    for (int i = 0; i < static_cast<int>(Spaces::Priority::NUM_LEVELS); ++i) {
        openMPSpaces[i].clear();
    }
#endif
#ifdef KOKKOS_ENABLE_SERIAL
    for (int i = 0; i < static_cast<int>(Spaces::Priority::NUM_LEVELS); ++i) {
        serialSpaces[i].clear();
    }
#endif
}


} // namespace Tpetra