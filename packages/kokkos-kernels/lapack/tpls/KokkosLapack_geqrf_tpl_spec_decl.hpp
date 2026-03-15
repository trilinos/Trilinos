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

#ifndef KOKKOSLAPACK_GEQRF_TPL_SPEC_DECL_HPP_
#define KOKKOSLAPACK_GEQRF_TPL_SPEC_DECL_HPP_

#include <iostream>

namespace KokkosLapack {
namespace Impl {
template <class AViewType, class TauViewType, class InfoViewType>
inline void geqrf_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  printf("KokkosLapack::geqrf<> TPL MAGMA specialization for < %s , %s, %s >\n", typeid(AViewType).name(),
         typeid(TauViewType).name(), typeid(InfoViewType).name());
#else
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK
  printf("KokkosLapack::geqrf<> TPL Lapack specialization for < %s , %s, %s >\n", typeid(AViewType).name(),
         typeid(TauViewType).name(), typeid(InfoViewType).name());
#endif
#endif
#endif
}
}  // namespace Impl
}  // namespace KokkosLapack

// Generic Host side LAPACK (could be MKL or whatever)
#if defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK) || defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
#include <KokkosLapack_Host_tpl.hpp>

namespace KokkosLapack {
namespace Impl {

template <class AViewType, class TauViewType, class InfoViewType>
void lapackGeqrfWrapper(const AViewType& A, const TauViewType& Tau, const InfoViewType& Info) {
  using memory_space = typename AViewType::memory_space;
  using Scalar       = typename AViewType::non_const_value_type;
  using ALayout_t    = typename AViewType::array_layout;
  static_assert(std::is_same_v<ALayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - geqrf: A needs to have a Kokkos::LayoutLeft");
  const int m   = A.extent_int(0);
  const int n   = A.extent_int(1);
  const int lda = A.stride(1);

  int lwork = -1;
  // work needs to be at least length 1 to store the returned value for lwork
  Kokkos::View<Scalar*, memory_space> work("work array", 1);

  if constexpr (KokkosKernels::ArithTraits<Scalar>::is_complex) {
    using MagType = typename KokkosKernels::ArithTraits<Scalar>::mag_type;

    HostLapack<std::complex<MagType>>::geqrf(m, n, reinterpret_cast<std::complex<MagType>*>(A.data()), lda,
                                             reinterpret_cast<std::complex<MagType>*>(Tau.data()),
                                             reinterpret_cast<std::complex<MagType>*>(work.data()), lwork, Info.data());

    if (Info[0] < 0) return;

    lwork = static_cast<int>(work(0).real());

    work = Kokkos::View<Scalar*, memory_space>("geqrf work buffer", lwork);

    HostLapack<std::complex<MagType>>::geqrf(m, n, reinterpret_cast<std::complex<MagType>*>(A.data()), lda,
                                             reinterpret_cast<std::complex<MagType>*>(Tau.data()),
                                             reinterpret_cast<std::complex<MagType>*>(work.data()), lwork, Info.data());
  } else {
    HostLapack<Scalar>::geqrf(m, n, A.data(), lda, Tau.data(), work.data(), lwork, Info.data());

    if (Info[0] < 0) return;

    lwork = static_cast<int>(work(0));

    work = Kokkos::View<Scalar*, memory_space>("geqrf work buffer", lwork);

    HostLapack<Scalar>::geqrf(m, n, A.data(), lda, Tau.data(), work.data(), lwork, Info.data());
  }
}

#define KOKKOSLAPACK_GEQRF_LAPACK(SCALAR, LAYOUT, EXECSPACE, MEM_SPACE)                                                \
  template <>                                                                                                          \
  struct GEQRF<                                                                                                        \
      EXECSPACE,                                                                                                       \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, \
      geqrf_eti_spec_avail<EXECSPACE,                                                                                  \
                           Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>,                        \
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                      \
                           Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>,                         \
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                      \
                           Kokkos::View<int*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>,                            \
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                            \
    using AViewType =                                                                                                  \
        Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using TauViewType =                                                                                                \
        Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;  \
    using InfoViewType =                                                                                               \
        Kokkos::View<int*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
                                                                                                                       \
    static void geqrf(const EXECSPACE& /* space */, const AViewType& A, const TauViewType& Tau,                        \
                      const InfoViewType& Info) {                                                                      \
      Kokkos::Profiling::pushRegion("KokkosLapack::geqrf[TPL_LAPACK," #SCALAR "]");                                    \
      geqrf_print_specialization<AViewType, TauViewType, InfoViewType>();                                              \
      lapackGeqrfWrapper(A, Tau, Info);                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#if defined(KOKKOS_ENABLE_SERIAL)
KOKKOSLAPACK_GEQRF_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_GEQRF_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_GEQRF_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_GEQRF_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
KOKKOSLAPACK_GEQRF_LAPACK(float, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_GEQRF_LAPACK(double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_GEQRF_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_GEQRF_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined(KOKKOS_ENABLE_THREADS)
KOKKOSLAPACK_GEQRF_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_GEQRF_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_GEQRF_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_GEQRF_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_LAPACK

// CUSOLVER
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSOLVER
#include "KokkosLapack_cusolver.hpp"

namespace KokkosLapack {
namespace Impl {

template <class ExecutionSpace, class AViewType, class TauViewType, class InfoViewType>
void cusolverGeqrfWrapper(const ExecutionSpace& space, const AViewType& A, const TauViewType& Tau,
                          const InfoViewType& Info) {
  using memory_space = typename AViewType::memory_space;
  using Scalar       = typename AViewType::non_const_value_type;

  using ALayout_t = typename AViewType::array_layout;
  static_assert(std::is_same_v<ALayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - cusolver geqrf: A needs to have a Kokkos::LayoutLeft");
  const int m   = A.extent_int(0);
  const int n   = A.extent_int(1);
  const int lda = A.stride(1);
  int lwork     = 0;

  CudaLapackSingleton& s = CudaLapackSingleton::singleton();
  KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnSetStream(s.handle, space.cuda_stream()));
  if constexpr (std::is_same_v<Scalar, float>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnSgeqrf_bufferSize(s.handle, m, n, A.data(), lda, &lwork));
    Kokkos::View<float*, memory_space> Workspace("cusolver sgeqrf workspace", lwork);

    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnSgeqrf(s.handle, m, n, A.data(), lda, Tau.data(), Workspace.data(), lwork, Info.data()));
  }
  if constexpr (std::is_same_v<Scalar, double>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnDgeqrf_bufferSize(s.handle, m, n, A.data(), lda, &lwork));
    Kokkos::View<double*, memory_space> Workspace("cusolver dgeqrf workspace", lwork);

    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnDgeqrf(s.handle, m, n, A.data(), lda, Tau.data(), Workspace.data(), lwork, Info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnCgeqrf_bufferSize(s.handle, m, n, reinterpret_cast<cuComplex*>(A.data()), lda, &lwork));
    Kokkos::View<cuComplex*, memory_space> Workspace("cusolver cgeqrf workspace", lwork);

    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnCgeqrf(
        s.handle, m, n, reinterpret_cast<cuComplex*>(A.data()), lda, reinterpret_cast<cuComplex*>(Tau.data()),
        reinterpret_cast<cuComplex*>(Workspace.data()), lwork, Info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnZgeqrf_bufferSize(s.handle, m, n, reinterpret_cast<cuDoubleComplex*>(A.data()), lda, &lwork));
    Kokkos::View<cuDoubleComplex*, memory_space> Workspace("cusolver zgeqrf workspace", lwork);

    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnZgeqrf(s.handle, m, n, reinterpret_cast<cuDoubleComplex*>(A.data()),
                                                          lda, reinterpret_cast<cuDoubleComplex*>(Tau.data()),
                                                          reinterpret_cast<cuDoubleComplex*>(Workspace.data()), lwork,
                                                          Info.data()));
  }
  KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnSetStream(s.handle, NULL));
}

#define KOKKOSLAPACK_GEQRF_CUSOLVER(SCALAR, LAYOUT, MEM_SPACE)                                                         \
  template <>                                                                                                          \
  struct GEQRF<                                                                                                        \
      Kokkos::Cuda,                                                                                                    \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                                          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      true,                                                                                                            \
      geqrf_eti_spec_avail<Kokkos::Cuda,                                                                               \
                           Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                     \
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                      \
                           Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                      \
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                      \
                           Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                         \
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                            \
    using AViewType   = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                        \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                         \
    using TauViewType = Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                         \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                         \
    using InfoViewType =                                                                                               \
        Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;  \
                                                                                                                       \
    static void geqrf(const Kokkos::Cuda& space, const AViewType& A, const TauViewType& Tau,                           \
                      const InfoViewType& Info) {                                                                      \
      Kokkos::Profiling::pushRegion("KokkosLapack::geqrf[TPL_CUSOLVER," #SCALAR "]");                                  \
      geqrf_print_specialization<AViewType, TauViewType, InfoViewType>();                                              \
                                                                                                                       \
      cusolverGeqrfWrapper(space, A, Tau, Info);                                                                       \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSLAPACK_GEQRF_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GEQRF_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GEQRF_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GEQRF_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)

#if defined(KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE)
KOKKOSLAPACK_GEQRF_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_GEQRF_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_GEQRF_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_GEQRF_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSOLVER

// ROCSOLVER
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER
#include <KokkosBlas_tpl_spec.hpp>
#include <rocsolver/rocsolver.h>

namespace KokkosLapack {
namespace Impl {

template <class ExecutionSpace, class AViewType, class TauViewType, class InfoViewType>
void rocsolverGeqrfWrapper(const ExecutionSpace& space, const AViewType& A, const TauViewType& Tau,
                           const InfoViewType& Info) {
  using Scalar = typename AViewType::non_const_value_type;

  using ALayout_t = typename AViewType::array_layout;
  static_assert(std::is_same_v<ALayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - rocsolver geqrf: A needs to have a Kokkos::LayoutLeft");
  const rocblas_int m   = static_cast<rocblas_int>(A.extent(0));
  const rocblas_int n   = static_cast<rocblas_int>(A.extent(1));
  const rocblas_int lda = static_cast<rocblas_int>(A.stride(1));

  KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();
  KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));
  if constexpr (std::is_same_v<Scalar, float>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocsolver_sgeqrf(s.handle, m, n, A.data(), lda, Tau.data()));
  }
  if constexpr (std::is_same_v<Scalar, double>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocsolver_dgeqrf(s.handle, m, n, A.data(), lda, Tau.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocsolver_cgeqrf(s.handle, m, n,
                                                       reinterpret_cast<rocblas_float_complex*>(A.data()), lda,
                                                       reinterpret_cast<rocblas_float_complex*>(Tau.data())));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocsolver_zgeqrf(s.handle, m, n,
                                                       reinterpret_cast<rocblas_double_complex*>(A.data()), lda,
                                                       reinterpret_cast<rocblas_double_complex*>(Tau.data())));
  }
  Kokkos::deep_copy(Info, 0);  // Success
  KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));
}

#define KOKKOSLAPACK_GEQRF_ROCSOLVER(SCALAR, LAYOUT, MEM_SPACE)                                                        \
  template <>                                                                                                          \
  struct GEQRF<                                                                                                        \
      Kokkos::HIP,                                                                                                     \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,  \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,     \
      true,                                                                                                            \
      geqrf_eti_spec_avail<Kokkos::HIP,                                                                                \
                           Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                      \
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                      \
                           Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                       \
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                      \
                           Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                          \
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                            \
    using AViewType   = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                         \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                         \
    using TauViewType = Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                          \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                         \
    using InfoViewType =                                                                                               \
        Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;   \
                                                                                                                       \
    static void geqrf(const Kokkos::HIP& space, const AViewType& A, const TauViewType& Tau,                            \
                      const InfoViewType& Info) {                                                                      \
      Kokkos::Profiling::pushRegion("KokkosLapack::geqrf[TPL_ROCSOLVER," #SCALAR "]");                                 \
      geqrf_print_specialization<AViewType, TauViewType, InfoViewType>();                                              \
                                                                                                                       \
      rocsolverGeqrfWrapper(space, A, Tau, Info);                                                                      \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSLAPACK_GEQRF_ROCSOLVER(float, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_GEQRF_ROCSOLVER(double, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_GEQRF_ROCSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_GEQRF_ROCSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIPSpace)

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER

#endif
