//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_HPP_
#define KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_HPP_

namespace KokkosLapack {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class AMatrix, class SVector, class UMatrix, class VMatrix>
struct svd_tpl_spec_avail {
  enum : bool { value = false };
};

// LAPACK
#if defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK) || defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
#define KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_LAPACK(SCALAR, LAYOUT, EXECSPACE)                                        \
  template <>                                                                                                    \
  struct svd_tpl_spec_avail<                                                                                     \
      EXECSPACE,                                                                                                 \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                               \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                     \
      Kokkos::View<Kokkos::ArithTraits<SCALAR>::mag_type*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                     \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                               \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                     \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                               \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                   \
    enum : bool { value = true };                                                                                \
  };

#if defined(KOKKOS_ENABLE_SERIAL)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Serial)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Serial)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Serial)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Serial)
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_LAPACK(float, Kokkos::LayoutLeft, Kokkos::OpenMP)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_LAPACK(double, Kokkos::LayoutLeft, Kokkos::OpenMP)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::OpenMP)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::OpenMP)
#endif

#if defined(KOKKOS_ENABLE_THREADS)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Threads)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Threads)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Threads)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Threads)
#endif

#endif  // KOKKOSKERNELS_ENABLE_TPL_LAPACK || KOKKOSKERNELS_ENABLE_TPL_MKL

// CUSOLVER
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSOLVER
#define KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_CUSOLVER(SCALAR, LAYOUT, MEMSPACE)                                             \
  template <>                                                                                                          \
  struct svd_tpl_spec_avail<                                                                                           \
      Kokkos::Cuda,                                                                                                    \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<Kokkos::ArithTraits<SCALAR>::mag_type*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,             \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                                           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                         \
    enum : bool { value = true };                                                                                      \
  };

KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)

#if defined(KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif  // CUDAUVMSPACE
#endif  // CUSOLVER

// ROCSOLVER
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER
#define KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_ROCSOLVER(SCALAR, LAYOUT, MEMSPACE)                                           \
  template <>                                                                                                         \
  struct svd_tpl_spec_avail<                                                                                          \
      Kokkos::HIP,                                                                                                    \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<Kokkos::ArithTraits<SCALAR>::mag_type*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEMSPACE>,             \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEMSPACE>,                                           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                        \
    enum : bool { value = true };                                                                                     \
  };

KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_ROCSOLVER(float, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_ROCSOLVER(double, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_ROCSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_ROCSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIPSpace)

#if defined(KOKKOSKERNELS_INST_MEMSPACE_HIPMANAGEDSPACE)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_ROCSOLVER(float, Kokkos::LayoutLeft, Kokkos::HIPManagedSpace)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_ROCSOLVER(double, Kokkos::LayoutLeft, Kokkos::HIPManagedSpace)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_ROCSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIPManagedSpace)
KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_ROCSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIPManagedSpace)
#endif  // HIPMANAGEDSPACE
#endif  // ROCSOLVER

}  // namespace Impl
}  // namespace KokkosLapack

#endif  // KOKKOSLAPACK_SVD_TPL_SPEC_AVAIL_HPP_
