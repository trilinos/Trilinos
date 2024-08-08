/*
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
*/

#ifndef KOKKOSPARSE_SPGEMM_NOREUSE_TPL_SPEC_AVAIL_HPP_
#define KOKKOSPARSE_SPGEMM_NOREUSE_TPL_SPEC_AVAIL_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include "mkl.h"
#endif

namespace KokkosSparse {
namespace Impl {

// Specialization struct which defines whether a specialization exists
template <class CMatrix, class AMatrix, class BMatrix>
struct spgemm_noreuse_tpl_spec_avail {
  enum : bool { value = false };
};

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && (CUDA_VERSION >= 11000)
// For cuSparse 11 and up, use the non-reuse generic interface.
// But for cuSparse 10, there is only one interface
// so just let KokkosSparse::spgemm call the symbolic and numeric wrappers.

#define SPGEMM_NOREUSE_AVAIL_CUSPARSE(SCALAR, MEMSPACE)                                        \
  template <>                                                                                  \
  struct spgemm_noreuse_tpl_spec_avail<                                                        \
      KokkosSparse::CrsMatrix<SCALAR, int, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, void, int>, \
      KokkosSparse::CrsMatrix<const SCALAR, const int, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const int>,             \
      KokkosSparse::CrsMatrix<const SCALAR, const int, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const int>> {           \
    enum : bool { value = true };                                                              \
  };

#define SPGEMM_NOREUSE_AVAIL_CUSPARSE_S(SCALAR)            \
  SPGEMM_NOREUSE_AVAIL_CUSPARSE(SCALAR, Kokkos::CudaSpace) \
  SPGEMM_NOREUSE_AVAIL_CUSPARSE(SCALAR, Kokkos::CudaUVMSpace)

SPGEMM_NOREUSE_AVAIL_CUSPARSE_S(float)
SPGEMM_NOREUSE_AVAIL_CUSPARSE_S(double)
SPGEMM_NOREUSE_AVAIL_CUSPARSE_S(Kokkos::complex<float>)
SPGEMM_NOREUSE_AVAIL_CUSPARSE_S(Kokkos::complex<double>)

#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#define SPGEMM_NOREUSE_AVAIL_MKL(SCALAR, EXEC)                                                          \
  template <>                                                                                           \
  struct spgemm_noreuse_tpl_spec_avail<                                                                 \
      KokkosSparse::CrsMatrix<SCALAR, MKL_INT, Kokkos::Device<EXEC, Kokkos::HostSpace>, void, MKL_INT>, \
      KokkosSparse::CrsMatrix<const SCALAR, const MKL_INT, Kokkos::Device<EXEC, Kokkos::HostSpace>,     \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const MKL_INT>,                  \
      KokkosSparse::CrsMatrix<const SCALAR, const MKL_INT, Kokkos::Device<EXEC, Kokkos::HostSpace>,     \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const MKL_INT>> {                \
    enum : bool { value = true };                                                                       \
  };

#define SPGEMM_NOREUSE_AVAIL_MKL_E(EXEC)                 \
  SPGEMM_NOREUSE_AVAIL_MKL(float, EXEC)                  \
  SPGEMM_NOREUSE_AVAIL_MKL(double, EXEC)                 \
  SPGEMM_NOREUSE_AVAIL_MKL(Kokkos::complex<float>, EXEC) \
  SPGEMM_NOREUSE_AVAIL_MKL(Kokkos::complex<double>, EXEC)

#ifdef KOKKOS_ENABLE_SERIAL
SPGEMM_NOREUSE_AVAIL_MKL_E(Kokkos::Serial)
#endif
#ifdef KOKKOS_ENABLE_OPENMP
SPGEMM_NOREUSE_AVAIL_MKL_E(Kokkos::OpenMP)
#endif
#endif

}  // namespace Impl
}  // namespace KokkosSparse

#endif
