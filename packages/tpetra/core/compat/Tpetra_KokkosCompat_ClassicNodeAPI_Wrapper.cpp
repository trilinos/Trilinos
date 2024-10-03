// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include "Teuchos_ParameterList.hpp"

namespace Tpetra {
namespace KokkosCompat {

#ifdef KOKKOS_ENABLE_THREADS
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Threads>::name () {
      return "Threads/Wrapper";
    }
#endif

#ifdef KOKKOS_ENABLE_OPENMP
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::OpenMP>::name () {
      return "OpenMP/Wrapper";
    }
#endif

#ifdef KOKKOS_ENABLE_SERIAL
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Serial>::name () {
      return "Serial/Wrapper";
    }
#endif // KOKKOS_ENABLE_SERIAL

#ifdef KOKKOS_ENABLE_CUDA
#ifdef HAVE_TPETRA_SHARED_ALLOCS
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Cuda,Kokkos::CudaUVMSpace>::name() {
      return std::string("Cuda/Wrapper");
    }
#else
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Cuda>::name() {
      return std::string("Cuda/Wrapper");
    }
#endif
#endif // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_HIP
#ifdef HAVE_TPETRA_SHARED_ALLOCS
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::HIP, Kokkos::HIPManagedSpace>::name() {
      return std::string("HIP/Wrapper");
    }
#else
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::HIP, Kokkos::HIPSpace>::name() {
      return std::string("HIP/Wrapper");
    }
#endif
#endif // KOKKOS_ENABLE_HIP

#ifdef KOKKOS_ENABLE_SYCL
#ifdef HAVE_TPETRA_SHARED_ALLOCS
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Experimental::SYCL, Kokkos::Experimental::SYCLSharedUSMSpace>::name() {
      return std::string("SYCL/Wrapper");
    }
#else
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Experimental::SYCL, Kokkos::Experimental::SYCLDeviceUSMSpace>::name() {
      return std::string("SYCL/Wrapper");
    }
#endif
#endif // KOKKOS_ENABLE_SYCL


} // namespace KokkosCompat
} // namespace Tpetra



