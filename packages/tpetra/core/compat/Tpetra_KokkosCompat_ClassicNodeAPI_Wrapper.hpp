// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_KOKKOSCOMPAT_CLASSICNODEAPI_WRAPPER_HPP
#define TPETRA_KOKKOSCOMPAT_CLASSICNODEAPI_WRAPPER_HPP

#include "Kokkos_Core.hpp"
#include "TpetraCore_config.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
//
// Dear users: These are just forward declarations.  Please skip
// over them and go down to KokkosDeviceWrapperNode below.  Thank
// you.
//
namespace Teuchos {
  class ParameterList;
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Tpetra {

namespace KokkosCompat {

/// \brief Node that wraps a new Kokkos execution space.
///
/// \tparam ExecutionSpace The type of the Kokkos execution space to wrap.
/// \tparam MemorySpace The Kokkos memory space in which to work.
///   Defaults to the default memory space of ExecutionSpace.
template<class ExecutionSpace,
         class MemorySpace = typename ExecutionSpace::memory_space>
class KokkosDeviceWrapperNode {
public:
  //! The Node's Kokkos execution space.
  typedef ExecutionSpace execution_space;
  //! The Node's Kokkos memory space.
  typedef MemorySpace memory_space;
  /// \brief The Node's Kokkos::Device specialization.
  ///
  /// This is just an (execution space, memory space) pair.
  typedef ::Kokkos::Device<execution_space, memory_space> device_type;

  /// \brief This is NOT a "classic" Node type.
  ///
  /// We will deprecate the "classic" Node types with the 11.14
  /// release of Trilinos, and remove them entirely with the 12.0
  /// release.  This Node type is safe to use.
  static constexpr bool classic = false;

  //! Whether the ExecutionSpace is Kokkos::Serial.
#ifdef KOKKOS_ENABLE_SERIAL
  static constexpr bool is_serial = std::is_same_v<ExecutionSpace, Kokkos::Serial>;
#else
  static constexpr bool is_serial = false;
#endif

  static constexpr bool is_cpu = std::is_same_v<typename ExecutionSpace::memory_space, Kokkos::HostSpace>;

  //! Whether the ExecutionSpace is GPU-like (its default memory space is not HostSpace)
  static constexpr bool is_gpu = !is_cpu;

  KokkosDeviceWrapperNode (Teuchos::ParameterList& /* params */) = delete;
  KokkosDeviceWrapperNode () = delete;

  //! Human-readable name of this Node.
  static std::string name ();
};

#ifdef KOKKOS_ENABLE_SYCL
#ifdef HAVE_TPETRA_SHARED_ALLOCS
  typedef KokkosDeviceWrapperNode<::Kokkos::Experimental::SYCL, ::Kokkos::Experimental::SYCLSharedUSMSpace> KokkosSYCLWrapperNode;
#else
  typedef KokkosDeviceWrapperNode<::Kokkos::Experimental::SYCL, ::Kokkos::Experimental::SYCLDeviceUSMSpace> KokkosSYCLWrapperNode;
#endif
#endif

#ifdef KOKKOS_ENABLE_HIP
#ifdef HAVE_TPETRA_SHARED_ALLOCS
  typedef KokkosDeviceWrapperNode<::Kokkos::HIP, ::Kokkos::HIPManagedSpace> KokkosHIPWrapperNode;
#else
  typedef KokkosDeviceWrapperNode<::Kokkos::HIP, ::Kokkos::HIPSpace> KokkosHIPWrapperNode;
#endif
#endif

#ifdef KOKKOS_ENABLE_CUDA
#ifdef HAVE_TPETRA_SHARED_ALLOCS
  typedef KokkosDeviceWrapperNode<::Kokkos::Cuda,::Kokkos::CudaUVMSpace> KokkosCudaWrapperNode;
#else
  typedef KokkosDeviceWrapperNode<::Kokkos::Cuda> KokkosCudaWrapperNode;
#endif
#endif

#ifdef KOKKOS_ENABLE_OPENMP
  typedef KokkosDeviceWrapperNode<::Kokkos::OpenMP> KokkosOpenMPWrapperNode;
#endif

#ifdef KOKKOS_ENABLE_THREADS
  typedef KokkosDeviceWrapperNode<::Kokkos::Threads> KokkosThreadsWrapperNode;
#endif

#ifdef KOKKOS_ENABLE_SERIAL
  typedef KokkosDeviceWrapperNode<::Kokkos::Serial> KokkosSerialWrapperNode;
#endif // KOKKOS_ENABLE_SERIAL

  // The above definitions / initializations of class (static)
  // variables need to precede the first use of these variables.
  // Otherwise, CUDA 7.5 with GCC 4.8.4 emits a warning ("explicit
  // specialization of member ... must precede its first use").

} // namespace KokkosCompat
} // namespace Tpetra




#endif // TPETRAKOKKOSCOMPAT_CLASSICNODEAPI_WRAPPER_HPP
