#ifndef TPETRA_DETAILS_SPACEs_HPP
#define TPETRA_DETAILS_SPACEs_HPP

#include <vector>
#include <iostream>
#include <sstream>

#include <Kokkos_Core.hpp>
#include "Tpetra_Details_Behavior.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_Details_nvtx.hpp"

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

/*! \brief Priority interface for Tpetra's managed execution spaces

    Priority is best-effort. low <= medium <= high, but it may be the case that priority levels are equivalent.
*/
enum class Priority {
    low = 0,
    medium = 1,
    high = 2,
    NUM_LEVELS = 3 // not to be used as a priority
};

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

/*! \brief Should be called by all functions in the Tpetra::Details::Spaces namespace

    * Prepares resources for Kokkos::CUDA exec space instance sync
    * Maps Tpetra::Priority to CUDA stream priorities
*/
void lazy_init();


#ifdef KOKKOS_ENABLE_CUDA
struct CudaInfo {
    bool initialized_;
    int lowPrio_;
    int mediumPrio_; // same as CUDA default
    int highPrio_;
    cudaEvent_t execSpaceWaitEvent_;// see exec_space_wait

    CudaInfo();
    ~CudaInfo() = default; // execSpaceWaitEvent_ cleaned up by CUDA deinit
    CudaInfo(const CudaInfo &other) = delete;
    CudaInfo(CudaInfo &&other) = delete;
};
extern CudaInfo cudaInfo;
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
///\brief IsSerial<Space> is a type if Space is Kokkos::Serial
template <typename Space> 
using IsSerial = std::enable_if_t<std::is_same<Space, Kokkos::Serial>::value, bool>;
#endif // KOKKOS_ENABLE_SERIAL

#ifdef KOKKOS_ENABLE_OPENMP
///\brief IsOpenMP<Space> is a type if Space is Kokkos::OpenMP
template <typename Space> 
using IsOpenMP = std::enable_if_t<std::is_same<Space, Kokkos::OpenMP>::value, bool>;
#endif // KOKKOS_ENABLE_OPENMP

#ifdef KOKKOS_ENABLE_HIP
///\brief IsHIP<Space> is a type if Space is Kokkos::HIP
template <typename Space> 
using IsHIP = std::enable_if_t<std::is_same<Space, Kokkos::HIP>::value, bool>;
#endif // KOKKOS_ENABLE_HIP

#ifdef KOKKOS_ENABLE_SYCL
///\brief IsSYCL<Space> is a type if Space is Kokkos::SYCL
template <typename Space> 
using IsSYCL = std::enable_if_t<std::is_same<Space, Kokkos::SYCL>::value, bool>;
#endif // KOKKOS_ENABLE_SYCL

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
        case Priority::high: prio = cudaInfo.highPrio_; break;
        case Priority::medium: prio = cudaInfo.mediumPrio_; break;
        case Priority::low: prio = cudaInfo.lowPrio_; break;
        default: throw std::runtime_error("unexpected static Tpetra Space priority");
    }
    CUDA_RUNTIME(cudaStreamCreateWithPriority(&stream, cudaStreamNonBlocking, prio));
    // std::cerr << __FILE__ << ":" << __LINE__ << ": stream " << uintptr_t(stream) << " with prio " << prio << "\n";
    return Kokkos::Cuda(stream, true /*Kokkos will manage this stream*/);
}
#endif // KOKKOS_ENABLE_CUDA

/*! \brief Construct a Kokkos execution space instance with the requested priority
    \tparam ExecSpace the type of Kokkos execution space to produce
    \param prio The Tpetra::Details::Spaces::Priority of the execution space to produce
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
    \tparam ExecSpace the type of Kokkos execution space to manage

    Holds a weak RCP to exec space instances, but provides strong RCPs to callers.
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

        \tparam priority the Spaces::Details::Priority of the provided instance
        \param i Which execution space instance to provide (default = Tpetra::Details::Spaces::Priority::medium)
    */
    template <
      Priority priority = Priority::medium
    >
    rcp_type space_instance(int i = 0) {
        constexpr int p = static_cast<int>(priority);
        static_assert(p < sizeof(instances) / sizeof(instances[0]), "Spaces::Priority enum error");

        if (i < 0) {
            THROW_RUNTIME("requested instance id " << i << " (< 0)");
        }
        if (size_t(i) >= Tpetra::Details::Behavior::spacesIdWarnLimit()) {
            THROW_RUNTIME(
                "requested instance id " << i << " (>= " 
                << Tpetra::Details::Behavior::spacesIdWarnLimit()
                << ") set by TPETRA_SPACES_ID_WARN_LIMT");
        }
        

        // make sure we can store an exec space at index i for priority p
        if (size_t(i) <= instances[p].size()) {
            instances[p].resize(i+1);
        }

        /* no exec space instance i of priority p exists.
           It may have never existed, or all Tpetra objects referencing it have been destructed
           create a new RCP<ExecSpace> and internally store a weak reference, so this space
           will be destructed when all strong references to it are gone, but we can still
           refer to it as long as it lives to prevent recreating
        */
        if (instances[p][i].is_null() || !instances[p][i].is_valid_ptr()) {
            rcp_type r = Teuchos::RCP<const execution_space>(
                new ExecSpace(make_instance<ExecSpace, priority>())
            );
            instances[p][i] = r.create_weak();
            return r; // allow strong rcp to escape so interneral weak one does not immediately go away
        }

        auto r = instances[p][i].create_strong();
        return r;
    }

    /*! \brief Issue a warning if any Tpetra-managed execution space instances survive to the end of static lifetime
    */
    ~InstanceLifetimeManager() {
        for (int i = 0; i < static_cast<int>(Spaces::Priority::NUM_LEVELS); ++i) {
            for (const rcp_type &rcp : instances[i]) {
                if (rcp.is_valid_ptr() && !rcp.is_null()) {
                    // avoid throwing in dtor
                    std::cerr << __FILE__ << ":" << __LINE__ << " execution space instance survived to ~InstanceLifetimeManager. strong_count() = " << rcp.strong_count() << ". Did a Tpetra object live past Kokkos::finalize()?" << std::endl;
                }
            }
        }
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
#ifdef KOKKOS_ENABLE_HIP
extern InstanceLifetimeManager<Kokkos::HIP> HIPSpaces;
#endif
#ifdef KOKKOS_ENABLE_SYCL
extern InstanceLifetimeManager<Kokkos::SYCL> SYCLSpaces;
#endif


#ifdef KOKKOS_ENABLE_CUDA
template <
  typename ExecSpace,
  Priority priority = Priority::medium,
  IsCuda<ExecSpace> = true 
> Teuchos::RCP<const ExecSpace> space_instance(int i = 0) {
    return cudaSpaces.space_instance<priority>(i);
}
#endif
#ifdef KOKKOS_ENABLE_SERIAL
template <
  typename ExecSpace,
  Priority priority = Priority::medium,
  IsSerial<ExecSpace> = true 
> Teuchos::RCP<const ExecSpace> space_instance(int i = 0) {
    return serialSpaces.space_instance<priority>(i);
}
#endif
#ifdef KOKKOS_ENABLE_OPENMP
template <
  typename ExecSpace,
  Priority priority = Priority::medium,
  IsOpenMP<ExecSpace> = true 
> Teuchos::RCP<const ExecSpace> space_instance(int i = 0) {
    return openMPSpaces.space_instance<priority>(i);
}
#endif
#ifdef KOKKOS_ENABLE_HIP
template <
  typename ExecSpace,
  Priority priority = Priority::medium,
  IsOpenMP<ExecSpace> = true 
> Teuchos::RCP<const ExecSpace> space_instance(int i = 0) {
    return HIPSpaces.space_instance<priority>(i);
}
#endif
#ifdef KOKKOS_ENABLE_SYCL
template <
  typename ExecSpace,
  Priority priority = Priority::medium,
  IsOpenMP<ExecSpace> = true 
> Teuchos::RCP<const ExecSpace> space_instance(int i = 0) {
    return SYCLSpaces.space_instance<priority>(i);
}
#endif

/*! \brief
    \param priority the priority of the execution space instance
    \tparam ExecSpace the type of Execution 
    \returns a strong RCP<ExecSpace>
*/
template <typename ExecSpace>
Teuchos::RCP<const ExecSpace> space_instance(const Priority &priority, int i = 0) {
    switch (priority) {
        case Priority::high: return space_instance<ExecSpace, Priority::high>(i);
        case Priority::medium: return space_instance<ExecSpace, Priority::medium>(i);
        case Priority::low: return space_instance<ExecSpace, Priority::low>(i);
        default: throw std::runtime_error("unexpected dynamic Tpetra Space priority in space_instance");
    }
}


/*! \brief cause future work submitted to waiter to wait for the current work in waitee to finish

  \tparam S1 the type of waitee
  \tparam S2 the type of waiter
  \param msg
  \param waitee Future work submitted to this stream will wait for work in \c waiter to finish
  \param waiter Future work submitted to \c waitee will wait for work in this stream to finish

  For Kokkos::Cuda execution spaces, this function may return immediately (i.e., without synchronizing the host).
*/
template <typename S1, typename S2
#ifdef KOKKOS_ENABLE_CUDA
, NotBothCuda<S1, S2> = true
#endif
>
void exec_space_wait(const char *msg, const S1 &waitee, const S2 &/*waiter*/) {
    Tpetra::Details::Range range(msg);
    lazy_init();
    waitee.fence(msg);
}

#ifdef KOKKOS_ENABLE_CUDA
template <typename S1, typename S2,
BothCuda<S1, S2> = true>
void exec_space_wait(const char *msg, const S1 &waitee, const S2 &waiter) {
    Tpetra::Details::Range range(msg);
    lazy_init();

    // TODO: add profiling hooks
    // if they are the same instance, no sync needed
    if (waitee.impl_instance_id() != waiter.impl_instance_id()) { // TODO: use instance operator== once available
        /* cudaStreamWaitEvent is not affected by later calls to cudaEventRecord, even if it overwrites
        the state of a shared event
        this means we only need one event even if many exec_space_waits are in flight at the same time
        */
        CUDA_RUNTIME(cudaEventRecord(cudaInfo.execSpaceWaitEvent_, waitee.cuda_stream()));
        CUDA_RUNTIME(cudaStreamWaitEvent(waiter.cuda_stream(), cudaInfo.execSpaceWaitEvent_, 0 /*flags*/));
    }
}
#endif


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

