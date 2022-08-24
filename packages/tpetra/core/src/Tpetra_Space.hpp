#ifndef TPETRA_SPACE_HPP
#define TPETRA_SPACE_HPP

#include <vector>
#include <iostream>

#include <Kokkos_Core.hpp>

/*! \file

Interface for Tpetra's managed Kokkos execution spaces

A space is addressed by a <priority, index> tuple.
Even if a space does not support priorities (e.g. Kokkos::Serial), 
spaces <Priority::high, 0> and <Priority::low, 0> are different.
*/

namespace Tpetra {

    enum class Priority {
        low = 0,
        medium = 1,
        high = 2,
        NUM_LEVELS = 3 // not to be used as a priority
    };

namespace detail {

    // print a warning to stderr if more than this many spaces are used
    constexpr int NUM_SPACES_WARN_LIMIT = 128;

    // query the runtime to map Tpetra::Priorities to the implementation priority
    void lazy_init();

#ifdef KOKKOS_ENABLE_CUDA
struct CudaPriorityRange {
    bool isSet = false;
    int low;
    int medium = 0; // cudaDeviceGetStreamPriorityRange has 0 for the default priority
    int high;
};

extern CudaPriorityRange cudaPriorityRange;
#endif // KOKKOS_ENABLE_CUDA

// Tpetra's managed spaces
#ifdef KOKKOS_ENABLE_CUDA
std::vector<Kokkos::Cuda> cudaSpaces[static_cast<int>(Priority::NUM_LEVELS)];
#endif
#ifdef KOKKOS_ENABLE_SERIAL
std::vector<Kokkos::Serial> serialSpaces[static_cast<int>(Priority::NUM_LEVELS)];
#endif
#ifdef KOKKOS_ENABLE_OPENMP
std::vector<Kokkos::OpenMP> openMPSpaces[static_cast<int>(Priority::NUM_LEVELS)];
#endif


// Tpetra's managed spaces
#ifdef KOKKOS_ENABLE_CUDA
template <typename Space> 
using IsCuda = std::enable_if_t<std::is_same<Space, Kokkos::Cuda>::value, bool>;
#endif
#ifdef KOKKOS_ENABLE_SERIAL
template <typename Space> 
using IsSerial = std::enable_if_t<std::is_same<Space, Kokkos::Serial>::value, bool>;
#endif
#ifdef KOKKOS_ENABLE_OPENMP
template <typename Space> 
using IsOpenMP = std::enable_if_t<std::is_same<Space, Kokkos::OpenMP>::value, bool>;
#endif

#ifdef KOKKOS_ENABLE_CUDA
    template <typename Space, Priority priority, 
    IsCuda<Space> = true >
    std::vector<Space> &spaces() {
        return cudaSpaces[static_cast<int>(priority)];
    }
#endif
#ifdef KOKKOS_ENABLE_SERIAL
    template <typename Space, Priority priority, 
    IsSerial<Space> = true >
    std::vector<Space> &spaces() {
        return serialSpaces[static_cast<int>(priority)];
    }
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    template <typename Space, Priority priority, 
    IsOpenMP<Space> = true >
    std::vector<Space> &spaces() {
        return openMPSpaces[static_cast<int>(priority)];
    }
#endif


} // namespace detail

/* Get execution space i for a given priority

    In some algorithms, the desired priority of independent operations
    may be statically known

    Catch-all when we don't implement priority for spaces (non-CUDA)
*/
template <typename Space, Priority priority = Priority::medium>
Space &get_space(int i) {
    if (i < 0) {
        throw std::runtime_error("requested exec space < 0 from get_space");
    }
    if (i >= detail::NUM_SPACES_WARN_LIMIT) {
        std::cerr << "WARNING: requested space " << i << std::endl;
    }

    while (detail::spaces<Space, priority>().size() <= size_t(i)) {
        detail::spaces<Space, priority>().push_back(Space());
    }
    return detail::spaces<Space, priority>()[i];
}

/* Implement priority for CUDA spaces
*/
#ifdef KOKKOS_ENABLE_CUDA
template <typename Space, Priority priority, 
detail::IsCuda<Space> = true >
Kokkos::Cuda &get_space(int i) {

    detail::lazy_init();

    if (i < 0) {
        throw std::runtime_error("requested exec space < 0 from get_space");
    }
    if (i >= detail::NUM_SPACES_WARN_LIMIT) {
        std::cerr << "WARNING: requested space " << i << std::endl;
    }

    while (detail::spaces<Space, priority>().size() <= size_t(i)) {
        cudaStream_t stream;
        int prio;
        switch (priority) {
            case Priority::high: prio = detail::cudaPriorityRange.high; break;
            case Priority::medium: prio = detail::cudaPriorityRange.medium; break;
            case Priority::low: prio = detail::cudaPriorityRange.low; break;
            default: throw std::runtime_error("unexpected Tpetra Space priority");
        }
        cudaStreamCreateWithPriority(&stream, cudaStreamNonBlocking, prio);
        detail::spaces<Space, priority>().push_back(
            Kokkos::Cuda(stream, true /*Kokkos will manage this stream*/)
        );
    }
    return detail::spaces<Space, priority>()[i];

}
#endif

/* get Space i with priority prio
*/
template <typename Space>
Space &get_space(int i, const Priority &prio) {
    switch(prio) {
        case Priority::high: return get_space<Space, Priority::high>(i);
        case Priority::medium: return get_space<Space, Priority::medium>(i);
        case Priority::low: return get_space<Space, Priority::low>(i);
        default: throw std::runtime_error("unexpected Tpetra Space priority");
    }
}

void drop_exec_spaces();

}

#endif // TPETRA_SPACE_HPP
