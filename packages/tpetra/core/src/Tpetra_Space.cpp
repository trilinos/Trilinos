#include "Tpetra_Space.hpp"

#include <vector>
#include <iostream>

#define LOG_WARN(x) std::cerr << __FILE__ << ":" << __LINE__ << ": [WARN] " << x << std::endl;
#define LOG_INFO(x) std::cerr << __FILE__ << ":" << __LINE__ << ": [INFO] " << x << std::endl;

namespace Tpetra {

constexpr int NUM_SPACES_WARN_LIMIT = 128;

/// Tpetra's internal spaces
#ifdef KOKKOS_ENABLE_CUDA
std::vector<Kokkos::Cuda> cudaSpaces;
#endif
#ifdef KOKKOS_ENABLE_SERIAL
std::vector<Kokkos::Serial> serialSpaces;
#endif
#ifdef KOKKOS_ENABLE_OPENMP
std::vector<Kokkos::OpenMP> openMPSpaces;
#endif

#ifdef KOKKOS_ENABLE_SERIAL
template<>
Kokkos::Serial &get_space(int i) {
    if (i < 0) {
        throw std::runtime_error("requested exec space < 0 from get_space");
    }
    if (i >= NUM_SPACES_WARN_LIMIT) {
        std::cerr << "WARNING: requested space " << i << std::endl;
    }

    while (serialSpaces.size() <= size_t(i)) {
        serialSpaces.push_back(Kokkos::Serial());
    }
    return serialSpaces[i];
}
#endif

#ifdef KOKKOS_ENABLE_CUDA
template<>
Kokkos::Cuda &get_space(int i) {
    if (i < 0) {
        throw std::runtime_error("requested exec space < 0 from get_space");
    }
    if (i >= NUM_SPACES_WARN_LIMIT) {
        std::cerr << "WARNING: requested space " << i << std::endl;
    }

    while (cudaSpaces.size() <= size_t(i)) {
        cudaStream_t stream;
        int loPrio, hiPrio;
        cudaDeviceGetStreamPriorityRange(&loPrio, &hiPrio);
        std::cerr << "LOPRIO: " << loPrio << "\n";
        std::cerr << "HIPRIO: " << hiPrio << "\n";
        cudaStreamCreateWithPriority(&stream, cudaStreamNonBlocking, loPrio);
        cudaSpaces.push_back(Kokkos::Cuda(stream, true /*Kokkos will manage this stream*/));
    }
    return cudaSpaces[i];
}
#endif

#ifdef KOKKOS_ENABLE_OPENMP
template<>
Kokkos::OpenMP &get_space(int i) {
    if (i < 0) {
        throw std::runtime_error("requested exec space < 0 from get_space");
    }
    if (i >= NUM_SPACES_WARN_LIMIT) {
        std::cerr << "WARNING: requested space " << i << std::endl;
    }

    while (openMPSpaces.size() <= size_t(i)) {
        openMPSpaces.push_back(Kokkos::OpenMP());
    }
    return openMPSpaces[i];
}
#endif


void drop_exec_spaces() {
    LOG_INFO("drop_exec_spaces");
#ifdef KOKKOS_ENABLE_CUDA
    cudaSpaces.clear();
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    openMPSpaces.clear();
#endif
#ifdef KOKKOS_ENABLE_SERIAL
    serialSpaces.clear();
#endif
}

}

#undef LOG_WARN
#undef LOG_INFO