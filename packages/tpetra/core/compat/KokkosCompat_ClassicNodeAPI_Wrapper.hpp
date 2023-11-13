#ifndef TPETRAKOKKOSCOMPAT_CLASSICNODEAPI_WRAPPER_HPP
#define TPETRAKOKKOSCOMPAT_CLASSICNODEAPI_WRAPPER_HPP

#include "TpetraCore_config.h"

// the .cpp for this header sets this so we don't get a deprecation warning
// just for that file
#if not defined(TPETRA_DETAILS_INTERNAL_INCLUDE_SILENCE_DEPRECATION)
#if defined(TPETRA_ENABLE_DEPRECATED_CODE)
#warning "The header file Trilinos/packages/tpetra/core/compat/ClassicNodeAPI_Wrapper.hpp is deprecated. Use Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#include "Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#else
#error "The header file Trilinos/packages/tpetra/core/compat/ClassicNodeAPI_Wrapper.hpp is deprecated. Use Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#endif
#endif // not defined(TPETRA_DETAILS_INTERNAL_INCLUDE_SILENCE_DEPRECATION)


namespace Kokkos {

  // can't just do 
  // namespace Compat = Tpetra::KokkosCompat;
  // because Teuchos defines a Kokkos::Compat namespace
  // if no one does this, then we could make the alias

namespace [[deprecated]] Compat {

  template <typename ExecutionSpace,
            typename MemorySpace = 
              typename ExecutionSpace::memory_space>
  using KokkosDeviceWrapperNode [[deprecated]] = Tpetra::KokkosCompat::KokkosDeviceWrapperNode<ExecutionSpace, MemorySpace>;

#ifdef KOKKOS_ENABLE_SYCL
  using KokkosSYCLWrapperNode [[deprecated]] = KokkosDeviceWrapperNode<Kokkos::Experimental::SYCL, Kokkos::Experimental::SYCLDeviceUSMSpace>;
#endif

#ifdef KOKKOS_ENABLE_HIP
  using KokkosHIPWrapperNode [[deprecated]] = Tpetra::KokkosCompat::KokkosHIPWrapperNode;
#endif

#ifdef KOKKOS_ENABLE_CUDA
  using KokkosCudaWrapperNode [[deprecated]] = Tpetra::KokkosCompat::KokkosCudaWrapperNode;
#endif

#ifdef KOKKOS_ENABLE_OPENMP
  using KokkosOpenMPWrapperNode [[deprecated]] = Tpetra::KokkosCompat::KokkosOpenMPWrapperNode;
#endif

#ifdef KOKKOS_ENABLE_THREADS
  using KokkosThreadsWrapperNode [[deprecated]] = Tpetra::KokkosCompat::KokkosThreadsWrapperNode;
#endif

#ifdef KOKKOS_ENABLE_SERIAL
  using KokkosSerialWrapperNode [[deprecated]] = Tpetra::KokkosCompat::KokkosSerialWrapperNode;
#endif // KOKKOS_ENABLE_SERIAL
} // namespace Compat
} // namespace Kokkos

#endif // TPETRAKOKKOSCOMPAT_CLASSICNODEAPI_WRAPPER_HPP
