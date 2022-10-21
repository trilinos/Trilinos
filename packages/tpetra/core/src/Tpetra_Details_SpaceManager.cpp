#include "Tpetra_Details_SpaceManager.hpp"

namespace Tpetra {
namespace Details {



SpaceManager2::SpaceManager2() {
    // priority only implemented for CUDA, other instances are all equivalent
#ifdef KOKKOS_ENABLE_CUDA
    for (int i = 0; i < static_cast<int>(Spaces::Priority::NUM_LEVELS); ++i) {
        cudaSpace[i] = Spaces::make_instance<Kokkos::Cuda>(static_cast<Spaces::Priority>(i));
    }
#endif
}

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

} // namespace Details
} // namespace Tpetra