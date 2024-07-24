// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_RANDOM_HPP
#define TPETRA_DETAILS_RANDOM_HPP

#include "TpetraCore_config.h"
#include "Kokkos_Random.hpp"

namespace Tpetra {
namespace Details {

template<class ExecutionSpace>
class Static_Random_XorShift64_Pool {
public:
  // The resetPool function will re-initialize the pool based on the system RNG and the MPI rank.
  // On GPU architectures, this will likely involve non-trivial host-to-device transfers.
  static void resetPool(int mpi_rank);

  // The isSet function returns true if resetPool has been callled.
  static bool isSet();
  // The getPool function will return the existing pool.
  static Kokkos::Random_XorShift64_Pool<ExecutionSpace> & getPool();
};


#ifdef KOKKOS_ENABLE_CUDA
template<>
class Static_Random_XorShift64_Pool<typename Kokkos::CudaSpace::execution_space> {
public:
  static void resetPool(int mpi_rank);
  static bool isSet();
  static Kokkos::Random_XorShift64_Pool<typename Kokkos::CudaSpace::execution_space> & getPool();
};
#endif // KOKKOS_ENABLE_CUDA


#ifdef KOKKOS_ENABLE_HIP
template<>
class Static_Random_XorShift64_Pool<typename Kokkos::HIPSpace::execution_space> {
public:
  static void resetPool(int mpi_rank);
  static bool isSet();
  static Kokkos::Random_XorShift64_Pool<typename Kokkos::HIPSpace::execution_space> & getPool();
};
#endif // KOKKOS_ENABLE_HIP


#ifdef KOKKOS_ENABLE_SYCL
template<>
class Static_Random_XorShift64_Pool<typename Kokkos::Experimental::SYCLDeviceUSMSpace::execution_space> {
public:
  static void resetPool(int mpi_rank);
  static bool isSet();
  static Kokkos::Random_XorShift64_Pool<typename Kokkos::Experimental::SYCLDeviceUSMSpace::execution_space> & getPool();
};
#endif // KOKKOS_ENABLE_SYCL

#ifdef KOKKOS_ENABLE_OPENMP
template<>
class Static_Random_XorShift64_Pool<typename Kokkos::OpenMP> {
public:
  static void resetPool(int mpi_rank);
  static bool isSet();
  static Kokkos::Random_XorShift64_Pool<typename Kokkos::OpenMP> & getPool();
};
#endif // KOKKOS_ENABLE_OPENMP

#ifdef KOKKOS_ENABLE_SERIAL
template<>
class Static_Random_XorShift64_Pool<typename Kokkos::Serial> {
public:
  static void resetPool(int mpi_rank);
  static bool isSet();
  static Kokkos::Random_XorShift64_Pool<typename Kokkos::Serial> & getPool();
};
#endif // KOKKOS_ENABLE_SERIAL

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_RANDOM_HPP
