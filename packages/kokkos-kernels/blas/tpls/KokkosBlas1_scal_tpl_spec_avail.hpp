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

#ifndef KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class RV, class AV, class XV, int Xrank = XV::Rank>
struct scal_tpl_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

namespace KokkosBlas {
namespace Impl {

// Generic Host side BLAS (could be MKL or whatever)
#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS)
// double
#define KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_BLAS(SCALAR, LAYOUT, MEMSPACE)         \
  template <class ExecSpace>                                                   \
  struct scal_tpl_spec_avail<                                                  \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      SCALAR,                                                                  \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      1> {                                                                     \
    enum : bool { value = true };                                              \
  };

KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutLeft,
                                     Kokkos::HostSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutLeft,
                                     Kokkos::HostSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>,
                                     Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutLeft,
                                     Kokkos::HostSpace)

#endif

// cuBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
// double
#define KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(SCALAR, LAYOUT, MEMSPACE)       \
  template <class ExecSpace>                                                   \
  struct scal_tpl_spec_avail<                                                  \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      SCALAR,                                                                  \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      1> {                                                                     \
    enum : bool { value = true };                                              \
  };

KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutLeft,
                                       Kokkos::CudaSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutLeft,
                                       Kokkos::CudaSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>,
                                       Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>,
                                       Kokkos::LayoutLeft, Kokkos::CudaSpace)

KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutLeft,
                                       Kokkos::CudaUVMSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutLeft,
                                       Kokkos::CudaUVMSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>,
                                       Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>,
                                       Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)

#endif

// rocBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)

#define KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_ROCBLAS(SCALAR, LAYOUT, MEMSPACE)      \
  template <class ExecSpace>                                                   \
  struct scal_tpl_spec_avail<                                                  \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      SCALAR,                                                                  \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      1> {                                                                     \
    enum : bool { value = true };                                              \
  };

KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_ROCBLAS(double, Kokkos::LayoutLeft,
                                        Kokkos::Experimental::HIPSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_ROCBLAS(float, Kokkos::LayoutLeft,
                                        Kokkos::Experimental::HIPSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<double>,
                                        Kokkos::LayoutLeft,
                                        Kokkos::Experimental::HIPSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<float>,
                                        Kokkos::LayoutLeft,
                                        Kokkos::Experimental::HIPSpace)

#endif

}  // namespace Impl
}  // namespace KokkosBlas
#endif
