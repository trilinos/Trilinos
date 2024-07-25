// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_EXECUTIONSPACES_HPP
#define TPETRA_DETAILS_EXECUTIONSPACES_HPP

#include <iostream>
#include <sstream>
#include <vector>

#include <Kokkos_Core.hpp>

#include <Teuchos_RCP.hpp>

#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_Profiling.hpp"

/*! \file

Interface for Tpetra's managed Kokkos execution spaces
This facility strives to provide:
1. Caching previously-constructed spaces so they can be used on-demand in parts
of Tpetra
2. Spaces of different priorities where supported
3. Fast sync of spaces where supported

For each Kokkos backend, there is a singleton instance manager.
This singleton manages the lifetime of all instances for that Kokkos backend.
These singletons are not meant to be accessed directly, but through
top-level functions in this file.

*/

#define TPETRA_DETAILS_SPACES_THROW(x)                                         \
  {                                                                            \
    std::stringstream ss;                                                      \
    ss << __FILE__ << ":" << __LINE__ << ": " << x;                            \
    throw std::runtime_error(ss.str());                                        \
  }

namespace Tpetra {
namespace Details {
namespace Spaces {

/*! \brief Priority interface for Tpetra's managed execution spaces

    Priority is best-effort. low <= medium <= high, but it may be the case that
   priority levels are equivalent.
*/
enum class Priority {
  low = 0,
  medium = 1,
  high = 2,
  NUM_LEVELS = 3 // not to be used as a priority
};

#if defined(KOKKOS_ENABLE_CUDA)
inline void success_or_throw(cudaError_t err, const char *file,
                             const int line) {
  if (err != cudaSuccess) {
    std::stringstream ss;
    ss << file << ":" << line << ": ";
    ss << cudaGetErrorString(err);
    throw std::runtime_error(ss.str());
  }
}
#define TPETRA_DETAILS_SPACES_CUDA_RUNTIME(x)                                  \
  Tpetra::Details::Spaces::success_or_throw((x), __FILE__, __LINE__)
#endif // KOKKOS_ENABLE_CUDA

/*! \brief Should be called by all functions in the Tpetra::Details::Spaces
   namespace

    * Prepares resources for Kokkos::CUDA exec space instance sync
    * Maps Tpetra::Priority to CUDA stream priorities
*/
void lazy_init();

#if defined(KOKKOS_ENABLE_CUDA)
struct CudaInfo {
  bool initialized_;
  int lowPrio_;
  int mediumPrio_; // same as CUDA default
  int highPrio_;
  cudaEvent_t execSpaceWaitEvent_; // see exec_space_wait

  CudaInfo();
  ~CudaInfo() = default; // execSpaceWaitEvent_ cleaned up by CUDA deinit
  CudaInfo(const CudaInfo &other) = delete;
  CudaInfo(CudaInfo &&other) = delete;
};
extern CudaInfo cudaInfo;
#endif // KOKKOS_ENABLE_CUDA

// Tpetra's managed spaces
#if defined(KOKKOS_ENABLE_CUDA)
template <typename Space>
using IsCuda = std::enable_if_t<std::is_same_v<Space, Kokkos::Cuda>, bool>;
template <typename Space>
using NotCuda = std::enable_if_t<!std::is_same_v<Space, Kokkos::Cuda>, bool>;
template <typename S1, typename S2>
using BothCuda = std::enable_if_t<
    std::is_same_v<S1, Kokkos::Cuda> && std::is_same_v<S2, Kokkos::Cuda>, bool>;
template <typename S1, typename S2>
using NotBothCuda = std::enable_if_t<!std::is_same_v<S1, Kokkos::Cuda> ||
                                         !std::is_same_v<S2, Kokkos::Cuda>,
                                     bool>;
#endif // KOKKOS_ENABLE_CUDA

#if defined(KOKKOS_ENABLE_SERIAL)
///\brief IsSerial<Space> is a type if Space is Kokkos::Serial
template <typename Space>
using IsSerial = std::enable_if_t<std::is_same_v<Space, Kokkos::Serial>, bool>;
#endif // KOKKOS_ENABLE_SERIAL

#if defined(KOKKOS_ENABLE_OPENMP)
///\brief IsOpenMP<Space> is a type if Space is Kokkos::OpenMP
template <typename Space>
using IsOpenMP = std::enable_if_t<std::is_same_v<Space, Kokkos::OpenMP>, bool>;
#endif // KOKKOS_ENABLE_OPENMP

#if defined(KOKKOS_ENABLE_HIP)
///\brief IsHIP<Space> is a type if Space is Kokkos::HIP
template <typename Space>
using IsHIP = std::enable_if_t<std::is_same_v<Space, Kokkos::HIP>, bool>;
#endif // KOKKOS_ENABLE_HIP

#if defined(KOKKOS_ENABLE_SYCL)
///\brief IsSYCL<Space> is a type if Space is Kokkos::Experimental::SYCL
template <typename Space>
using IsSYCL = std::enable_if_t<std::is_same_v<Space, Kokkos::Experimental::SYCL>, bool>;
#endif // KOKKOS_ENABLE_SYCL

/*! \brief Construct a Kokkos execution space instance with the following
   priority

    If priority is not supported, returns the default space instance
*/
template <typename ExecSpace, Priority priority = Priority::medium
#if defined(KOKKOS_ENABLE_CUDA)
          ,
          NotCuda<ExecSpace> = true
#endif // KOKKOS_ENABLE_CUDA
          >
ExecSpace make_instance() {
  return ExecSpace();
}

/*! \brief Construct a Kokkos::Cuda execution space instance with the requested
   priority

    This creates a stream with the requested priority, and attaches it to a new
   space
*/
#if defined(KOKKOS_ENABLE_CUDA)
template <typename ExecSpace, Priority priority = Priority::medium,
          IsCuda<ExecSpace> = true>
Kokkos::Cuda make_instance() {
  lazy_init(); // CUDA priorities
  cudaStream_t stream;
  int prio;
  switch (priority) {
  case Priority::high:
    prio = cudaInfo.highPrio_;
    break;
  case Priority::medium:
    prio = cudaInfo.mediumPrio_;
    break;
  case Priority::low:
    prio = cudaInfo.lowPrio_;
    break;
  default:
    throw std::runtime_error("unexpected static Tpetra Space priority");
  }
  TPETRA_DETAILS_SPACES_CUDA_RUNTIME(
      cudaStreamCreateWithPriority(&stream, cudaStreamNonBlocking, prio));
  return Kokkos::Cuda(stream, true /*Kokkos will manage this stream*/);
}
#endif // KOKKOS_ENABLE_CUDA

/*! \brief Construct a Kokkos execution space instance with the requested
   priority \tparam ExecSpace the type of Kokkos execution space to produce
    \param prio The Tpetra::Details::Spaces::Priority of the execution space to
   produce
*/
template <typename ExecSpace> ExecSpace make_instance(const Priority &prio) {
  switch (prio) {
  case Priority::high:
    return make_instance<ExecSpace, Priority::high>();
  case Priority::medium:
    return make_instance<ExecSpace, Priority::medium>();
  case Priority::low:
    return make_instance<ExecSpace, Priority::low>();
  default:
    throw std::runtime_error("unexpected dynamic Tpetra Space priority");
  }
}

/*! \brief Provides reusable Kokkos execution space instances
    \tparam ExecSpace the type of Kokkos execution space to manage

    Holds a weak RCP to exec space instances, but provides strong RCPs to
   callers. Callers in Tpetra can use this to quickly get an execution space
   instance on-demand.


    When all strong RCP holders go away, the referenced instance will also go
   away. This prevents the spaces from living longer than Kokkos.
*/
template <typename ExecSpace> class InstanceLifetimeManager {
public:
  using execution_space = ExecSpace;
  using rcp_type = Teuchos::RCP<const execution_space>;

  /*! \brief Retrieve a strong `Teuchos::RCP<const ExecSpace>` to instance `i`

      \tparam priority the Spaces::Details::Priority of the provided instance
      \param i Which execution space instance to provide (default =
     Tpetra::Details::Spaces::Priority::medium)
  */
  template <Priority priority = Priority::medium>
  rcp_type space_instance(int i = 0) {
    Tpetra::Details::ProfilingRegion region(
        "Tpetra::Details::Spaces::space_instance");

    constexpr int p = static_cast<int>(priority);
    static_assert(p < sizeof(instances) / sizeof(instances[0]),
                  "Spaces::Priority enum error");

    if (i < 0) {
      TPETRA_DETAILS_SPACES_THROW("requested instance id " << i << " (< 0)");
    }
    if (size_t(i) >= Tpetra::Details::Behavior::spacesIdWarnLimit()) {
      TPETRA_DETAILS_SPACES_THROW(
          "requested instance id "
          << i << " (>= " << Tpetra::Details::Behavior::spacesIdWarnLimit()
          << ") set by TPETRA_SPACES_ID_WARN_LIMIT");
    }

    // make sure we can store an exec space at index i for priority
    // not sure what happens in RCP(), so let's explicitly make it null
    while (size_t(i) >= instances[p].size()) {
      instances[p].push_back(Teuchos::ENull());
    }

    /* no exec space instance i of priority p exists.
       It may have never existed, or all Tpetra objects referencing it have been
       destructed.

       Create a new RCP<ExecSpace> and internally store a weak
       reference, so this space will be destructed when all strong references to
       it are gone, but we can still refer to it as long as it lives to prevent
       recreating
    */
    if (instances[p][i].is_null() || !instances[p][i].is_valid_ptr()) {

      // create a strong RCP to a space
      rcp_type r = Teuchos::RCP<const execution_space>(
          new ExecSpace(make_instance<ExecSpace, priority>()));

      // store a weak RCP to the space
      instances[p][i] = r.create_weak();

      return r; // allow strong rcp to escape so internal weak one does not
                // immediately go away
    }

    auto r = instances[p][i].create_strong();
    return r;
  }

  /*! \brief Issue a warning if any Tpetra-managed execution space instances
   * survive to the end of static lifetime
   */
  ~InstanceLifetimeManager() {
    for (int i = 0; i < static_cast<int>(Spaces::Priority::NUM_LEVELS); ++i) {
      for (const rcp_type &rcp : instances[i]) {
        if (rcp.is_valid_ptr() && !rcp.is_null()) {
          // avoid throwing in dtor
          std::cerr << __FILE__ << ":" << __LINE__
                    << " execution space instance survived to "
                       "~InstanceLifetimeManager. strong_count() = "
                    << rcp.strong_count()
                    << ". Did a Tpetra object live past Kokkos::finalize()?"
                    << std::endl;
        }
      }
    }
  }

private:
  // one vector of instances for each priority level
  std::vector<rcp_type>
      instances[static_cast<int>(Spaces::Priority::NUM_LEVELS)];
};

#if defined(KOKKOS_ENABLE_CUDA)
extern InstanceLifetimeManager<Kokkos::Cuda> cudaSpaces;
#endif
#if defined(KOKKOS_ENABLE_SERIAL)
extern InstanceLifetimeManager<Kokkos::Serial> serialSpaces;
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
extern InstanceLifetimeManager<Kokkos::OpenMP> openMPSpaces;
#endif
#if defined(KOKKOS_ENABLE_HIP)
extern InstanceLifetimeManager<Kokkos::HIP> HIPSpaces;
#endif
#if defined(KOKKOS_ENABLE_SYCL)
extern InstanceLifetimeManager<Kokkos::Experimental::SYCL> SYCLSpaces;
#endif

#if defined(KOKKOS_ENABLE_CUDA)

/*! \brief get a strong Teuchos::RCP to Tpetra-managed Kokkos::Cuda instance \c
 * i
 */
template <typename ExecSpace, Priority priority = Priority::medium,
          IsCuda<ExecSpace> = true>
Teuchos::RCP<const ExecSpace> space_instance(int i = 0) {
  return cudaSpaces.space_instance<priority>(i);
}
#endif

#if defined(KOKKOS_ENABLE_SERIAL)
/*! \brief get a strong Teuchos::RCP to Tpetra-managed Kokkos::Serial instance
 * \c i
 */
template <typename ExecSpace, Priority priority = Priority::medium,
          IsSerial<ExecSpace> = true>
Teuchos::RCP<const ExecSpace> space_instance(int i = 0) {
  return serialSpaces.space_instance<priority>(i);
}
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
/*! \brief get a strong Teuchos::RCP to Tpetra-managed Kokkos::OpenMP instance
 * \c i
 */
template <typename ExecSpace, Priority priority = Priority::medium,
          IsOpenMP<ExecSpace> = true>
Teuchos::RCP<const ExecSpace> space_instance(int i = 0) {
  return openMPSpaces.space_instance<priority>(i);
}
#endif

#if defined(KOKKOS_ENABLE_HIP)
/*! \brief get a strong Teuchos::RCP to Tpetra-managed Kokkos::HIP instance \c i
 */
template <typename ExecSpace, Priority priority = Priority::medium,
          IsHIP<ExecSpace> = true>
Teuchos::RCP<const ExecSpace> space_instance(int i = 0) {
  return HIPSpaces.space_instance<priority>(i);
}
#endif
#if defined(KOKKOS_ENABLE_SYCL)
/*! \brief get a strong Teuchos::RCP to Tpetra-managed Kokkos::Experimental::SYCL instance \c
 * i
 */
template <typename ExecSpace, Priority priority = Priority::medium,
          IsSYCL<ExecSpace> = true>
Teuchos::RCP<const ExecSpace> space_instance(int i = 0) {
  return SYCLSpaces.space_instance<priority>(i);
}
#endif

/*! \brief get a strong Teuchos::RCP to Tpetra-managed Kokkos execution space
   instance \c i \param priority the priority of the execution space instance
    \tparam ExecSpace the type of Execution
    \returns a strong RCP<ExecSpace>
*/
template <typename ExecSpace>
Teuchos::RCP<const ExecSpace> space_instance(const Priority &priority,
                                             int i = 0) {
  switch (priority) {
  case Priority::high:
    return space_instance<ExecSpace, Priority::high>(i);
  case Priority::medium:
    return space_instance<ExecSpace, Priority::medium>(i);
  case Priority::low:
    return space_instance<ExecSpace, Priority::low>(i);
  default:
    throw std::runtime_error(
        "unexpected dynamic Tpetra Space priority in space_instance");
  }
}

/*! \brief cause future work submitted to waiter to wait for the current work in
  waitee to finish

  \tparam S1 the type of waitee
  \tparam S2 the type of waiter
  \param msg
  \param waitee Future work submitted to this stream will wait for work in \c
  waiter to finish \param waiter Future work submitted to \c waitee will wait
  for work in this stream to finish

  For Kokkos::Cuda execution spaces, this function may return immediately (i.e.,
  without synchronizing the host).
*/
template <typename S1, typename S2
#if defined(KOKKOS_ENABLE_CUDA)
          ,
          NotBothCuda<S1, S2> = true
#endif
          >
void exec_space_wait(const char *msg, const S1 &waitee, const S2 & /*waiter*/) {
  Tpetra::Details::ProfilingRegion r(
      "Tpetra::Details::Spaces::exec_space_wait");
  lazy_init();
  waitee.fence(msg);
}

#if defined(KOKKOS_ENABLE_CUDA)
template <typename S1, typename S2, BothCuda<S1, S2> = true>
void exec_space_wait(const char *msg, const S1 &waitee, const S2 &waiter) {
  Tpetra::Details::ProfilingRegion r(
      "Tpetra::Details::Spaces::exec_space_wait");
  lazy_init();

  // if they are the same instance, no sync needed
  if (waitee.impl_instance_id() !=
      waiter
          .impl_instance_id()) { // TODO: use instance operator== once available
    /* cudaStreamWaitEvent is not affected by later calls to cudaEventRecord,
    even if it overwrites the state of a shared event this means we only need
    one event even if many exec_space_waits are in flight at the same time
    */
    TPETRA_DETAILS_SPACES_CUDA_RUNTIME(
        cudaEventRecord(cudaInfo.execSpaceWaitEvent_, waitee.cuda_stream()));
    TPETRA_DETAILS_SPACES_CUDA_RUNTIME(cudaStreamWaitEvent(
        waiter.cuda_stream(), cudaInfo.execSpaceWaitEvent_, 0 /*flags*/));
  }
}
#endif

template <typename S1, typename S2>
void exec_space_wait(const S1 &waitee, const S2 &waiter) {
  Tpetra::Details::ProfilingRegion r(
      "Tpetra::Details::Spaces::exec_space_wait");
  lazy_init();
  exec_space_wait("anonymous", waitee, waiter);
}

template <typename ExecutionSpace>
constexpr KOKKOS_INLINE_FUNCTION bool is_gpu_exec_space() {
  return false;
}

#if defined(KOKKOS_ENABLE_CUDA)
template <>
constexpr KOKKOS_INLINE_FUNCTION bool is_gpu_exec_space<Kokkos::Cuda>() {
  return true;
}
#endif

#if defined(KOKKOS_ENABLE_HIP)
template <>
constexpr KOKKOS_INLINE_FUNCTION bool
is_gpu_exec_space<Kokkos::HIP>() {
  return true;
}
#endif

#if defined(KOKKOS_ENABLE_SYCL)
template <>
constexpr KOKKOS_INLINE_FUNCTION bool
is_gpu_exec_space<Kokkos::Experimental::SYCL>() {
  return true;
}
#endif

} // namespace Spaces
} // namespace Details
} // namespace Tpetra

#undef TPETRA_DETAILS_SPACES_THROW

#endif // TPETRA_DETAILS_EXECUTIONSPACES_HPP
