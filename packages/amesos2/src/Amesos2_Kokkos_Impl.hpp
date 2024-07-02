// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
  \file   Amesos2_Kokkos_Impl.hpp
  \author
  \date

  \brief  ETI for Solvers using Kokkos adapter
*/

#ifndef AMESOS2_KOKKOS_IMPL_HPP
#define AMESOS2_KOKKOS_IMPL_HPP

#include <type_traits>
#include "Amesos2_KokkosMultiVecAdapter_decl.hpp"
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#define AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(S,LO,NODE_TYPE)                             \
  template class Amesos2::AMESOS2_KOKKOS_IMPL_SOLVER_NAME<KokkosSparse::CrsMatrix<S, LO,         \
    typename NODE_TYPE::device_type>,                                                           \
    Kokkos::View<S**, Kokkos::LayoutLeft, typename NODE_TYPE::device_type> >;

#ifdef KOKKOS_ENABLE_CUDA_UVM
#define AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER_UVM_OFF(S,LO)       \
  template class Amesos2::AMESOS2_KOKKOS_IMPL_SOLVER_NAME<KokkosSparse::CrsMatrix<S, LO,          \
    Kokkos::Device<Kokkos::Cuda,Kokkos::CudaSpace>>,                                              \
    Kokkos::View<S**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Cuda,Kokkos::CudaSpace>> >;
#else
#define AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER_UVM_OFF(S,LO)
#endif

#if defined(KOKKOS_ENABLE_SERIAL)
#ifdef HAVE_TPETRA_INST_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(float, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode)
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(double, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<float>, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<double>, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode)
#endif
#endif

#if defined(KOKKOS_ENABLE_THREADS)
#ifdef HAVE_TPETRA_INST_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(float, int, Tpetra::KokkosCompat::KokkosKokkosThreadsWrapperNode)
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(double, int, Tpetra::KokkosCompat::KokkosKokkosThreadsWrapperNode)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<float>, int, Tpetra::KokkosCompat::KokkosKokkosThreadsWrapperNode)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<double>, int, Tpetra::KokkosCompat::KokkosKokkosThreadsWrapperNode)
#endif
#endif // KOKKOS_ENABLE_THREADS

#if defined(KOKKOS_ENABLE_OPENMP)
#ifdef HAVE_TPETRA_INST_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(float, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode)
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(double, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<float>, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<double>, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode)
#endif
#endif // KOKKOS_ENABLE_OPENMP

#if defined(KOKKOS_ENABLE_CUDA)
#ifdef HAVE_TPETRA_INST_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(float, int, Tpetra::KokkosCompat::KokkosCudaWrapperNode)
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER_UVM_OFF(float, int)
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(double, int, Tpetra::KokkosCompat::KokkosCudaWrapperNode)
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER_UVM_OFF(double, int)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<float>, int, Tpetra::KokkosCompat::KokkosCudaWrapperNode)
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER_UVM_OFF(Kokkos::complex<float>, int)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<double>, int, Tpetra::KokkosCompat::KokkosCudaWrapperNode)
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER_UVM_OFF(Kokkos::complex<double>, int)
#endif
#endif // KOKKOS_ENABLE_CUDA

#if defined(KOKKOS_ENABLE_HIP)
#ifdef HAVE_TPETRA_INST_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(float, int, Tpetra::KokkosCompat::KokkosHIPWrapperNode)
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(double, int, Tpetra::KokkosCompat::KokkosHIPWrapperNode)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<float>, int, Tpetra::KokkosCompat::KokkosHIPWrapperNode)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<double>, int, Tpetra::KokkosCompat::KokkosHIPWrapperNode)
#endif
#endif // KOKKOS_ENABLE_HIP

#if defined(KOKKOS_ENABLE_SYCL)
#ifdef HAVE_TPETRA_INST_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(float, int, Tpetra::KokkosCompat::KokkosSYCLWrapperNode)
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(double, int, Tpetra::KokkosCompat::KokkosSYCLWrapperNode)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<float>, int, Tpetra::KokkosCompat::KokkosSYCLWrapperNode)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<double>, int, Tpetra::KokkosCompat::KokkosSYCLWrapperNode)
#endif
#endif // KOKKOS_ENABLE_SYCL

#endif // AMESOS2_KOKKOS_IMPL_HPP
