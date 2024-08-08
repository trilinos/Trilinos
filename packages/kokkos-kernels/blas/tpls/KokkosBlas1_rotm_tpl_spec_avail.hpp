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

#ifndef KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class VectorView, class ParamView>
struct rotm_tpl_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

namespace KokkosBlas {
namespace Impl {

// Generic Host side BLAS (could be MKL or whatever)
// ARMPL is disabled as it does not detect some corner
// cases correctly which leads to failing unit-tests
#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS)
#define KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_BLAS(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                  \
  template <>                                                                                                        \
  struct rotm_tpl_spec_avail<                                                                                        \
      EXEC_SPACE,                                                                                                    \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<const SCALAR[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                       \
    enum : bool { value = true };                                                                                    \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif
#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#define KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_CUBLAS(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                \
  template <>                                                                                                        \
  struct rotm_tpl_spec_avail<                                                                                        \
      EXEC_SPACE,                                                                                                    \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<const SCALAR[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                       \
    enum : bool { value = true };                                                                                    \
  };

KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#define KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_ROCBLAS(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                               \
  template <>                                                                                                        \
  struct rotm_tpl_spec_avail<                                                                                        \
      EXEC_SPACE,                                                                                                    \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<const SCALAR[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                       \
    enum : bool { value = true };                                                                                    \
  };

KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_ROCBLAS(double, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace)
KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_ROCBLAS(float, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace)
KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_ROCBLAS(double, Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace)
KOKKOSBLAS1_ROTM_TPL_SPEC_AVAIL_ROCBLAS(float, Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace)
#endif

}  // namespace Impl
}  // namespace KokkosBlas
#endif
