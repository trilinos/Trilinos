#include "Tpetra_Space.hpp"

#include <vector>

namespace Tpetra {


/// Tpetra's internal spaces
std::vector<Kokkos::DefaultExecutionSpace> spaces;

Kokkos::DefaultExecutionSpace &get_exec_space(int i) {
    if (i < 0) {
        throw std::runtime_error("requested exec space < 0 from get_exec_space");
    }

    if (spaces.size() <= size_t(i)) {
#ifdef KOKKOS_ENABLE_CUDA
        cudaStream_t stream;
        int loPrio, hiPrio;
        cudaDeviceGetStreamPriorityRange(&loPrio, &hiPrio);
        std::cerr << "LOPRIO: " << loPrio << "\n";
        std::cerr << "HIPRIO: " << hiPrio << "\n";
        cudaStreamCreateWithPriority(&stream, cudaStreamNonBlocking, loPrio);
        spaces.push_back(Kokkos::Cuda(stream, true /*Kokkos will manage this stream*/));
#else
        spaces.push_back(Kokkos::DefaultExecSpace());
#endif
    }
    return spaces[i];
}

void drop_exec_spaces() {
    spaces.clear();
}

}