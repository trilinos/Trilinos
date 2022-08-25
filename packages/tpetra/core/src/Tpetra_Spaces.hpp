#ifndef TPETRA_SPACE_HPP
#define TPETRA_SPACE_HPP

#include <vector>
#include <iostream>
#include <unordered_map>
#include <functional>

#include <Kokkos_Core.hpp>

/*! \file

Interface for Tpetra's managed Kokkos execution spaces

A space is addressed by a <priority, index> tuple.
Even if a space does not support priorities (e.g. Kokkos::Serial), 
spaces <Priority::high, 0> and <Priority::low, 0> are different.

*/

/* define a hash of Kokkos::Cuda so it can be used in an std::unordered map
   Assuming the instance ID uniquely distinguishes execution spaces
*/
template<>
struct std::hash<Kokkos::Cuda>
{
    std::size_t operator()(Kokkos::Cuda const& s) const noexcept
    {
        // Invoke expression for non-static member function needs an eplicit instance to operate on
        // This std::result_of (invoke_result in C++20) just is the type that Kokkos uses as an instance ID in case they change it
        return std::hash<std::result_of<decltype(&Kokkos::Cuda::impl_instance_id)(Kokkos::Cuda)>::type>{}(s.impl_instance_id());
    }
};

/* define equality for Kokkos::Cuda so it can be used in an std::unordered map
   Assuming the instance ID uniquely distinguishes execution spaces
*/
namespace Kokkos {
    bool operator==(const Kokkos::Cuda &lhs, const Kokkos::Cuda &rhs) {
        return lhs.impl_instance_id() == rhs.impl_instance_id();
    }
}


namespace Tpetra {

    namespace Spaces {

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
    void finalize();

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
std::unordered_map<Kokkos::Cuda, cudaStream_t> cudaStreams; // track for optimized inter-space sync
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
template <typename Space> 
using NotCuda = std::enable_if_t<!std::is_same<Space, Kokkos::Cuda>::value, bool>;
template <typename S1, typename S2>
using BothCuda = std::enable_if_t<
    std::is_same<S1, Kokkos::Cuda>::value
    && std::is_same<S2, Kokkos::Cuda>::value
, bool>;
template <typename S1, typename S2>
using NotBothCuda = std::enable_if_t<
    !std::is_same<S1, Kokkos::Cuda>::value
    || !std::is_same<S2, Kokkos::Cuda>::value
, bool>;
#endif // KOKKOS_ENABLE_CUDA
#ifdef KOKKOS_ENABLE_SERIAL
template <typename Space> 
using IsSerial = std::enable_if_t<std::is_same<Space, Kokkos::Serial>::value, bool>;
#endif // KOKKOS_ENABLE_SERIAL
#ifdef KOKKOS_ENABLE_OPENMP
template <typename Space> 
using IsOpenMP = std::enable_if_t<std::is_same<Space, Kokkos::OpenMP>::value, bool>;
#endif // KOKKOS_ENABLE_OPENMP

#ifdef KOKKOS_ENABLE_CUDA
    template <typename Space, Priority priority, 
    IsCuda<Space> = true >
    std::vector<Space> &spaces() {
        return cudaSpaces[static_cast<int>(priority)];
    }
#endif // KOKKOS_ENABLE_CUDA
#ifdef KOKKOS_ENABLE_SERIAL
    template <typename Space, Priority priority, 
    IsSerial<Space> = true >
    std::vector<Space> &spaces() {
        return serialSpaces[static_cast<int>(priority)];
    }
#endif // KOKKOS_ENABLE_SERIAL
#ifdef KOKKOS_ENABLE_OPENMP
    template <typename Space, Priority priority, 
    IsOpenMP<Space> = true >
    std::vector<Space> &spaces() {
        return openMPSpaces[static_cast<int>(priority)];
    }
#endif // KOKKOS_ENABLE_OPENMP


} // namespace detail

/* Get execution space i for a given priority

    In some algorithms, the desired priority of independent operations
    may be statically known

    Catch-all when we don't implement priority for spaces (non-CUDA)
*/
template <typename Space, Priority priority = Priority::medium
#ifdef KOKKOS_ENABLE_CUDA
, detail::NotCuda<Space> = true
#endif
>
Space &get(int i) {
    if (i < 0) {
        throw std::runtime_error("requested exec space < 0 from Spaces::get");
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
template <typename Space, Priority priority = Priority::medium, 
detail::IsCuda<Space> = true >
Kokkos::Cuda &get(int i) {

    detail::lazy_init();

    if (i < 0) {
        throw std::runtime_error("requested exec space < 0 from Spaces::get");
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
        Kokkos::Cuda space (stream, true /*Kokkos will manage this stream*/);
        detail::spaces<Space, priority>().push_back(space);
        detail::cudaStreams[space] = stream;
    }
    return detail::spaces<Space, priority>()[i];

}
#endif

/* get Space i with priority prio
*/
template <typename Space>
Space &get(int i, const Priority &prio) {
    switch(prio) {
        case Priority::high: return get<Space, Priority::high>(i);
        case Priority::medium: return get<Space, Priority::medium>(i);
        case Priority::low: return get<Space, Priority::low>(i);
        default: throw std::runtime_error("unexpected Tpetra Space priority");
    }
}


/* cause future work submitted to waiter to wait for the current work in waitee to finish
  may return immediately (e.g., before waitee's work is done)
*/
template <typename S1, typename S2
#ifdef KOKKOS_ENABLE_CUDA
, detail::NotBothCuda<S1, S2> = true
#endif
>
void exec_space_wait(S1 &waitee, S2 &waiter) {
    (void) waiter;
    waitee.fence();
}

#ifdef KOKKOS_ENABLE_CUDA
template <typename S1, typename S2,
detail::BothCuda<S1, S2> = true>
void exec_space_wait(S1 &waitee, S2 &waiter) {

    // TODO:move to init
    cudaEvent_t e1;
    cudaEventCreateWithFlags(&e1, cudaEventDisableTiming);

    cudaEventRecord(e1, detail::cudaStreams[waitee]);
    cudaStreamWaitEvent(detail::cudaStreams[waiter], e1, 0 /*flags*/);

    // TODO: move to cleanup
    // can't do this right here
    cudaEventDestroy(e1);
}
#endif

} // namespace Spaces

} // namespace Tpetra

#endif // TPETRA_SPACE_HPP
