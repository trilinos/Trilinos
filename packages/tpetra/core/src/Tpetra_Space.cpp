#include "Tpetra_Space.hpp"

#include <vector>

namespace Tpetra {


/// Tpetra's internal spaces
std::vector<Kokkos::Cuda> cudaSpaces;
std::vector<Kokkos::Serial> serialSpaces;
std::vector<Kokkos::OpenMP> openMPSpaces;

template<>
Kokkos::Serial &get_space(int i) {
    if (i < 0) {
        throw std::runtime_error("requested exec space < 0 from get_space");
    }

    while (serialSpaces.size() <= size_t(i)) {
        serialSpaces.push_back(Kokkos::Serial());
    }
    return serialSpaces[i];
}

template<>
Kokkos::Cuda &get_space(int i) {
    if (i < 0) {
        throw std::runtime_error("requested exec space < 0 from get_space");
    }

    while (cudaSpaces.size() <= size_t(i)) {
        cudaStream_t stream;
        int loPrio, hiPrio;
        cudaDeviceGetStreamPriorityRange(&loPrio, &hiPrio);
        std::cerr << "LOPRIO: " << loPrio << "\n";
        std::cerr << "HIPRIO: " << hiPrio << "\n";
        cudaStreamCreateWithPriority(&stream, cudaStreamNonBlocking, loPrio);
        cudaSpaces.push_back(Kokkos::Cuda(stream, true /*Kokkos will manage this stream*/));
        // cudaSpaces.push_back(Kokkos::DefaultExecSpace());
    }
    return cudaSpaces[i];
}

template<>
Kokkos::OpenMP &get_space(int i) {
    if (i < 0) {
        throw std::runtime_error("requested exec space < 0 from get_space");
    }

    while (openMPSpaces.size() <= size_t(i)) {
        openMPSpaces.push_back(Kokkos::OpenMP());
    }
    return openMPSpaces[i];
}


void drop_exec_spaces() {
    cudaSpaces.clear();
}

}