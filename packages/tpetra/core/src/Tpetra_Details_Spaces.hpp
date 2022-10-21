#ifndef TPETRA_DETAILS_SPACEs_HPP
#define TPETRA_DETAILS_SPACEs_HPP

#include <vector>
#include <iostream>
#include <sstream>

#include <Kokkos_Core.hpp>
#include "Tpetra_Details_Behavior.hpp"
#include "Teuchos_RCP.hpp"

/*! \file

Interface for Tpetra's managed Kokkos execution spaces

A space is addressed by a <priority, index> tuple.
Even if a space does not support priorities (e.g. Kokkos::Serial), 
spaces <Priority::high, 0> and <Priority::low, 0> are different.

*/

#define THROW_RUNTIME(x) { \
    std::stringstream ss; \
    ss << __FILE__ << ":" << __LINE__ << ": " << x; \
    throw std::runtime_error(ss.str()); \
}

namespace Tpetra {
namespace Details {
namespace Spaces {

enum class Priority {
    low = 0,
    medium = 1,
    high = 2,
    NUM_LEVELS = 3 // not to be used as a priority
};

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
#define CUDA_RUNTIME(x) Tpetra::Details::Spaces::success_or_throw((x), __FILE__, __LINE__)
#endif // KOKKOS_ENABLE_CUDA

/*! \brief Automaticallyed called by functions in the Tpetra::Spaces namespace

    * Prepares resources for Kokkos::CUDA exec space instance sync
    * Tpetra::Priority to CUDA stream priorities
*/
void lazy_init();


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

/*! \brief Construct a Kokkos execution space instance with the requested priority
*/
template <typename ExecSpace, Priority priority = Priority::medium
#ifdef KOKKOS_ENABLE_CUDA
, NotCuda<ExecSpace> = true
#endif // KOKKOS_ENABLE_CUDA
>
ExecSpace make_instance() {
    return ExecSpace();
}

#ifdef KOKKOS_ENABLE_CUDA
template <typename ExecSpace, Priority priority = Priority::medium, 
IsCuda<ExecSpace> = true >
Kokkos::Cuda make_instance() {
    lazy_init(); // CUDA priorities
    cudaStream_t stream;
    int prio;
    switch (priority) {
        case Priority::high: prio = cudaPriorityRange.high; break;
        case Priority::medium: prio = cudaPriorityRange.medium; break;
        case Priority::low: prio = cudaPriorityRange.low; break;
        default: throw std::runtime_error("unexpected static Tpetra Space priority");
    }
    CUDA_RUNTIME(cudaStreamCreateWithPriority(&stream, cudaStreamNonBlocking, prio));
    // std::cerr << __FILE__ << ":" << __LINE__ << ": stream " << uintptr_t(stream) << " with prio " << prio << "\n";
    return Kokkos::Cuda(stream, true /*Kokkos will manage this stream*/);
}
#endif // KOKKOS_ENABLE_CUDA

/*! \brief Construct a Kokkos execution space instance with the requested priority
*/
template <typename ExecSpace>
ExecSpace make_instance(const Priority &prio) {
    switch (prio) {
        case Priority::high: return make_instance<ExecSpace, Priority::high>();
        case Priority::medium: return make_instance<ExecSpace, Priority::medium>();
        case Priority::low: return make_instance<ExecSpace, Priority::low>();
        default: throw std::runtime_error("unexpected dynamic Tpetra Space priority");
    }
}

/*! \brief Provides reusable Kokkos execution space instances

    Holds a weak RCP to exec space instances, but provides strong RCPs to users.
    When all strong RCP holders go away, the instance will also go away, restricting
    the lifetime of the underlying instance to Tpetra objects.
    This prevents the instance outliving Kokkos;
*/
template <typename ExecSpace>
class InstanceLifetimeManager {
public:
    using execution_space = ExecSpace;
    using rcp_type = Teuchos::RCP<const execution_space>;

    /*! \brief Retrieve a strong `Teuchos::RCP<const ExecSpace>` to instance `i`

        \tparam priority the Spaces::Details::Priority of the instance
    */
    template <
      Priority priority = Priority::medium
    >
    rcp_type space_instance(int i = 0) const {
        constexpr int p = static_cast<int>(priority);
        static_assert(p < sizeof(instances), "Spaces::Priority enum error");

        if (i < 0) {
            THROW_RUNTIME("requested instance id " << i << " (< 0)");
        }
        if (i > Tpetra::Details::Behavior::spacesIdWarnLimit()) {
            THROW_RUNTIME(
                "requested instance id " << i << " (> " 
                << Tpetra::Details::Behavior::spacesIdWarnLimit()
                << ") set by TPETRA_SPACES_ID_WARN_LIMT");
        }
        

        // make sure we can st

        if (i <= instances[p].size()) {
            instances[p].resize(i+1);
        }

        if (!instances[p][i]) {
            rcp_type r = Teuchos::RCP<const execution_space>(
                new ExecSpace(make_instance<ExecSpace, priority>())
            );
            instances[p][i] = r.create_weak();
            return r;
        }

        return instances[p][i].create_strong();
    }
private:
    // one vector of instances for each priority level
    std::vector<rcp_type> instances[static_cast<int>(Spaces::Priority::NUM_LEVELS)];
};

#ifdef KOKKOS_ENABLE_CUDA
extern InstanceLifetimeManager<Kokkos::Cuda> cudaSpaces;
#endif
#ifdef KOKKOS_ENABLE_SERIAL
extern InstanceLifetimeManager<Kokkos::Serial> serialSpaces;
#endif
#ifdef KOKKOS_ENABLE_OPENMP
extern InstanceLifetimeManager<Kokkos::OpenMP> openMPSpaces;
#endif

template <typename ExecSpace, Priority priority = Priority::medium>
Teuchos::RCP<const ExecSpace> space_instance(int i = 0) {

#ifdef KOKKOS_ENABLE_CUDA
    if (std::is_same<ExecSpace, Kokkos::Cuda>::value) {
        return cudaSpaces.space_instance(i);
    }
#endif

#ifdef KOKKOS_ENABLE_SERIAL
    if (std::is_same<ExecSpace, Kokkos::Serial>::value) {
        return serialSpaces.space_instance(i);
    }
#endif

#ifdef KOKKOS_ENABLE_OPENMP
    if (std::is_same<ExecSpace, Kokkos::OpenMP>::value) {
        return openMPSpaces.space_instance(i);
    }
#endif

    throw std::runtime_error("not implemented for space");

}


/* cause future work submitted to waiter to wait for the current work in waitee to finish
  may return immediately (e.g., before waitee's work is done)
*/
template <typename S1, typename S2
#ifdef KOKKOS_ENABLE_CUDA
, NotBothCuda<S1, S2> = true
#endif
>
void exec_space_wait(const char *msg, const S1 &waitee, const S2 &waiter) {
    lazy_init();
    (void) waiter;
    waitee.fence(msg);
}

#ifdef KOKKOS_ENABLE_CUDA
template <typename S1, typename S2,
BothCuda<S1, S2> = true>
void exec_space_wait(const char */*msg*/, const S1 &waitee, const S2 &waiter) {
    lazy_init();
    /* cudaStreamWaitEvent is not affected by later calls to cudaEventRecord, even if it overwrites
       the state of a shared event
       this means we only need one event even if many exec_space_waits are in flight at the same time
    */

    // TODO: add profiling hooks
    CUDA_RUNTIME(cudaEventRecord(execSpaceWaitEvent, waitee.cuda_stream()));
    CUDA_RUNTIME(cudaStreamWaitEvent(waiter.cuda_stream(), execSpaceWaitEvent, 0 /*flags*/));
}
#endif

/* cause future work submitted to waiter to wait for the current work in waitee to finish
  may return immediately (e.g., before waitee's work is done)
*/
template <typename S1, typename S2>
void exec_space_wait(const S1 &waitee, const S2 &waiter) {
    lazy_init();
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
} // namespace Details
} // namespace Tpetra

#undef THROW_RUNTIME

#endif // TPETRA_DETAILS_SPACEs_HPP

