// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_HPP_
#define KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_HPP_

namespace KokkosLapack {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class AMatrix, class TauArray, class InfoArray>
struct gegqr_tpl_spec_avail {
  enum : bool { value = false };
};

// Generic Host side LAPACK (could be MKL or whatever)
#if defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK) || defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)

#define KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_LAPACK(SCALAR, LAYOUT, MEMSPACE)                                          \
  template <class ExecSpace>                                                                                        \
  struct gegqr_tpl_spec_avail<                                                                                      \
      ExecSpace,                                                                                                    \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,  \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {   \
    enum : bool { value = true };                                                                                   \
  };

KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_LAPACK(double, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_LAPACK(float, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HostSpace)
#endif
}  // namespace Impl
}  // namespace KokkosLapack

// CUSOLVER
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSOLVER
namespace KokkosLapack {
namespace Impl {

#define KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_CUSOLVER(SCALAR, LAYOUT, MEMSPACE)                                           \
  template <>                                                                                                          \
  struct gegqr_tpl_spec_avail<                                                                                         \
      Kokkos::Cuda,                                                                                                    \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,  \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {   \
    enum : bool { value = true };                                                                                      \
  };

KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaSpace)

#if defined(KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE)
KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // CUSOLVER

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER
#include <rocsolver/rocsolver.h>

namespace KokkosLapack {
namespace Impl {

#define KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_ROCSOLVER(SCALAR, LAYOUT, MEMSPACE)                                         \
  template <>                                                                                                         \
  struct gegqr_tpl_spec_avail<                                                                                        \
      Kokkos::HIP,                                                                                                    \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,  \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {   \
    enum : bool { value = true };                                                                                     \
  };

KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_ROCSOLVER(double, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_ROCSOLVER(float, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_ROCSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_ROCSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIPSpace)

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER

#endif  // KOKKOSLAPACK_GEGQR_TPL_SPEC_AVAIL_HPP_
