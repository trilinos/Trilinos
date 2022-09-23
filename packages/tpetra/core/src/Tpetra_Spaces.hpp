#ifndef TPETRA_SPACE_HPP
#define TPETRA_SPACE_HPP

#include <vector>
#include <iostream>
#include <sstream>

#include <Kokkos_Core.hpp>
#include "Tpetra_Details_Behavior.hpp"

/*! \file

Interface for Tpetra's managed Kokkos execution spaces

A space is addressed by a <priority, index> tuple.
Even if a space does not support priorities (e.g. Kokkos::Serial), 
spaces <Priority::high, 0> and <Priority::low, 0> are different.

*/


namespace Tpetra {
namespace Spaces {

enum class Priority {
    low = 0,
    medium = 1,
    high = 2,
    NUM_LEVELS = 3 // not to be used as a priority
};

namespace detail {

/* This tracks whether spaces have been initialized.
   Not all unit-tests call Tpetra::initialize, so we
   have to do our own lazy-init on each call
*/
extern bool initialized;

#ifdef KOKKOS_ENABLE_SERIAL
inline void print_space(const Kokkos::Serial &space) {
    std::cerr << "[serial]";
}
#endif
#ifdef KOKKOS_ENABLE_CUDA
inline void print_space(const Kokkos::Cuda &space) {
    std::cerr << uintptr_t(space.cuda_stream());
}
#endif
#ifdef KOKKOS_ENABLE_OPENMP
inline void print_space(const Kokkos::OpenMP &space) {
    std::cerr << "[OpenMP]";
}
#endif

#ifdef KOKKOS_ENABLE_CUDA
inline void success_or_throw(cudaError_t err, const char *file, const int line) {
    if (err != cudaSuccess) {
        std::stringstream ss;
        ss << file << ":" << line << ": ";
        ss << cudaGetErrorString(err);
        throw std::runtime_error(ss.str());
    }
}
#define CUDA_RUNTIME(x) Tpetra::Spaces::detail::success_or_throw((x), __FILE__, __LINE__)
#endif // KOKKOS_ENABLE_CUDA

/* query the runtime to map Tpetra::Priorities to the implementation priority
   set up `cudaEvent_t`s 
*/
void lazy_init();
void initialize();
void finalize();


#ifdef KOKKOS_ENABLE_CUDA
struct CudaPriorityRange {
    bool isSet = false;
    int low;
    int medium = 0; // cudaDeviceGetStreamPriorityRange has 0 for the default priority
    int high;
};
extern CudaPriorityRange cudaPriorityRange;
extern cudaEvent_t execSpaceWaitEvent; // see exec_space_wait
#endif // KOKKOS_ENABLE_CUDA

// Tpetra's managed spaces
#ifdef KOKKOS_ENABLE_CUDA
extern std::vector<Kokkos::Cuda> cudaSpaces[static_cast<int>(Priority::NUM_LEVELS)];
#endif
#ifdef KOKKOS_ENABLE_SERIAL
extern std::vector<Kokkos::Serial> serialSpaces[static_cast<int>(Priority::NUM_LEVELS)];
#endif
#ifdef KOKKOS_ENABLE_OPENMP
extern std::vector<Kokkos::OpenMP> openMPSpaces[static_cast<int>(Priority::NUM_LEVELS)];
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
        lazy_init();
        return cudaSpaces[static_cast<int>(priority)];
    }
#endif // KOKKOS_ENABLE_CUDA
#ifdef KOKKOS_ENABLE_SERIAL
    template <typename Space, Priority priority, 
    IsSerial<Space> = true >
    std::vector<Space> &spaces() {
        lazy_init();
        return serialSpaces[static_cast<int>(priority)];
    }
#endif // KOKKOS_ENABLE_SERIAL
#ifdef KOKKOS_ENABLE_OPENMP
    template <typename Space, Priority priority, 
    IsOpenMP<Space> = true >
    std::vector<Space> &spaces() {
        lazy_init();
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
Space &get(int i = 0) {
    detail::lazy_init();
    if (i < 0) {
        throw std::runtime_error("requested exec space < 0 from Spaces::get");
    }
    if (i >= Tpetra::Details::Behavior::spacesIdWarnLimit()) {
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
Kokkos::Cuda &get(int i = 0) {
    detail::lazy_init();

    if (i < 0) {
        throw std::runtime_error("requested exec space < 0 from Spaces::get");
    }
    if (i >= Tpetra::Details::Behavior::spacesIdWarnLimit()) {
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
        CUDA_RUNTIME(cudaStreamCreateWithPriority(&stream, cudaStreamNonBlocking, prio));
        std::cerr << __FILE__ << ":" << __LINE__ << ": stream " << uintptr_t(stream) << " with prio " << prio << "\n";
        Kokkos::Cuda space (stream, true /*Kokkos will manage this stream*/);
        detail::spaces<Space, priority>().push_back(space);
    }
    return detail::spaces<Space, priority>()[i];

}
#endif

/* get Space i with priority prio
*/
template <typename Space>
Space &get(int i, const Priority &prio) {
    detail::lazy_init();
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
void exec_space_wait(const char *msg, const S1 &waitee, const S2 &waiter) {
    detail::lazy_init();
    (void) waiter;
    waitee.fence(msg);
}

#ifdef KOKKOS_ENABLE_CUDA
template <typename S1, typename S2,
detail::BothCuda<S1, S2> = true>
void exec_space_wait(const char */*msg*/, const S1 &waitee, const S2 &waiter) {
    detail::lazy_init();
    /* cudaStreamWaitEvent is not affected by later calls to cudaEventRecord, even if it overwrites
       the state of a shared event
       this means we only need one event even if many exec_space_waits are in flight at the same time
    */

    // TODO: add profiling hooks
    CUDA_RUNTIME(cudaEventRecord(detail::execSpaceWaitEvent, waitee.cuda_stream()));
    CUDA_RUNTIME(cudaStreamWaitEvent(waiter.cuda_stream(), detail::execSpaceWaitEvent, 0 /*flags*/));
}
#endif

/* cause future work submitted to waiter to wait for the current work in waitee to finish
  may return immediately (e.g., before waitee's work is done)
*/
template <typename S1, typename S2>
void exec_space_wait(const S1 &waitee, const S2 &waiter) {
    detail::lazy_init();
    exec_space_wait("anonymous", waitee, waiter);
}

template <typename ExecutionSpace>
constexpr KOKKOS_INLINE_FUNCTION bool is_gpu_exec_space() {
  return false;
}

#ifdef KOKKOS_ENABLE_CUDA
template <>
constexpr KOKKOS_INLINE_FUNCTION bool is_gpu_exec_space<Kokkos::Cuda>() {
  return true;
}
#endif

#ifdef KOKKOS_ENABLE_HIP
template <>
constexpr KOKKOS_INLINE_FUNCTION bool
is_gpu_exec_space<Kokkos::Experimental::HIP>() {
  return true;
}
#endif

#ifdef KOKKOS_ENABLE_SYCL
template <>
constexpr KOKKOS_INLINE_FUNCTION bool
is_gpu_exec_space<Kokkos::Experimental::SYCL>() {
  return true;
}
#endif


} // namespace Spaces

} // namespace Tpetra

#endif // TPETRA_SPACE_HPP
