// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
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
