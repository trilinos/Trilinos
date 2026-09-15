// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_GETRF_TPL_SPEC_DECL_HPP_
#define KOKKOSLAPACK_GETRF_TPL_SPEC_DECL_HPP_

#include <iostream>

namespace KokkosLapack {
namespace Impl {
template <class AViewType, class TauViewType, class InfoViewType>
inline void getrf_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  printf("KokkosLapack::getrf<> TPL MAGMA specialization for < %s , %s, %s >\n", typeid(AViewType).name(),
         typeid(TauViewType).name(), typeid(InfoViewType).name());
#else
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK
  printf("KokkosLapack::getrf<> TPL Lapack specialization for < %s , %s, %s >\n", typeid(AViewType).name(),
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

template <class AViewType, class IpivView, class InfoView>
void lapackGetrfWrapper(const AViewType& A, const IpivView& Ipiv, const InfoView& Info) {
  using Scalar    = typename AViewType::non_const_value_type;
  using ALayout_t = typename AViewType::array_layout;
  static_assert(std::is_same_v<ALayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - getrf: A needs to have a Kokkos::LayoutLeft");
  const int m   = A.extent_int(0);
  const int n   = A.extent_int(1);
  const int lda = A.stride(1);

  if constexpr (KokkosKernels::ArithTraits<Scalar>::is_complex) {
    using MagType = typename KokkosKernels::ArithTraits<Scalar>::mag_type;

    HostLapack<std::complex<MagType>>::getrf(m, n, reinterpret_cast<std::complex<MagType>*>(A.data()), lda, Ipiv.data(),
                                             Info.data());
  } else {
    HostLapack<Scalar>::getrf(m, n, A.data(), lda, Ipiv.data(), Info.data());
  }
}

#define KOKKOSLAPACK_GETRF_LAPACK(SCALAR, LAYOUT, EXECSPACE, MEM_SPACE)                                                \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct GETRF<                                                                                                        \
      EXECSPACE,                                                                                                       \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,       \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, \
      ETI_SPEC_AVAIL> {                                                                                                \
    using AViewType =                                                                                                  \
        Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using IpivView =                                                                                                   \
        Kokkos::View<int*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    using InfoViewType =                                                                                               \
        Kokkos::View<int*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
                                                                                                                       \
    static void getrf(const EXECSPACE& /* space */, const AViewType& A, const IpivView& Ipiv,                          \
                      const InfoViewType& Info) {                                                                      \
      Kokkos::Profiling::pushRegion("KokkosLapack::getrf[TPL_LAPACK," #SCALAR "]");                                    \
      getrf_print_specialization<AViewType, IpivView, InfoViewType>();                                                 \
      lapackGetrfWrapper(A, Ipiv, Info);                                                                               \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#if defined(KOKKOS_ENABLE_SERIAL)
KOKKOSLAPACK_GETRF_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_GETRF_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_GETRF_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_GETRF_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
KOKKOSLAPACK_GETRF_LAPACK(float, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_GETRF_LAPACK(double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_GETRF_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_GETRF_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined(KOKKOS_ENABLE_THREADS)
KOKKOSLAPACK_GETRF_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_GETRF_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_GETRF_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_GETRF_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_LAPACK

// CUSOLVER
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSOLVER
#include "KokkosLapack_cusolver.hpp"
#include <KokkosKernels_Cuda_Utils.hpp>

namespace KokkosLapack {
namespace Impl {

template <class ExecutionSpace, class AViewType, class IpivViewType, class InfoViewType>
void cusolverGetrfWrapper(const ExecutionSpace& space, const AViewType& A, const IpivViewType& Ipiv,
                          const InfoViewType& Info) {
  using memory_space = typename AViewType::memory_space;
  using Scalar       = typename AViewType::non_const_value_type;

  using ALayout_t = typename AViewType::array_layout;
  static_assert(std::is_same_v<ALayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - cusolver getrf: A needs to have a Kokkos::LayoutLeft");

  const int m   = A.extent_int(0);
  const int n   = A.extent_int(1);
  const int lda = A.stride(1);
  int lwork     = 0;

  CudaLapackSingleton& s = CudaLapackSingleton::singleton();
  KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnSetStream(s.handle, space.cuda_stream()));
  if constexpr (std::is_same_v<Scalar, float>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnSgetrf_bufferSize(s.handle, m, n, A.data(), lda, &lwork));
    Kokkos::View<float*, memory_space> Workspace("cusolver sgetrf workspace", lwork);

    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnSgetrf(s.handle, m, n, A.data(), lda, Workspace.data(), Ipiv.data(), Info.data()));
  }
  if constexpr (std::is_same_v<Scalar, double>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnDgetrf_bufferSize(s.handle, m, n, A.data(), lda, &lwork));
    Kokkos::View<double*, memory_space> Workspace("cusolver dgetrf workspace", lwork);

    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnDgetrf(s.handle, m, n, A.data(), lda, Workspace.data(), Ipiv.data(), Info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnCgetrf_bufferSize(s.handle, m, n, reinterpret_cast<cuComplex*>(A.data()), lda, &lwork));
    Kokkos::View<cuComplex*, memory_space> Workspace("cusolver cgetrf workspace", lwork);

    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnCgetrf(s.handle, m, n, reinterpret_cast<cuComplex*>(A.data()), lda,
                                                          reinterpret_cast<cuComplex*>(Workspace.data()), Ipiv.data(),
                                                          Info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnZgetrf_bufferSize(s.handle, m, n, reinterpret_cast<cuDoubleComplex*>(A.data()), lda, &lwork));
    Kokkos::View<cuDoubleComplex*, memory_space> Workspace("cusolver zgetrf workspace", lwork);

    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnZgetrf(s.handle, m, n, reinterpret_cast<cuDoubleComplex*>(A.data()),
                                                          lda, reinterpret_cast<cuDoubleComplex*>(Workspace.data()),
                                                          Ipiv.data(), Info.data()));
  }
  KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnSetStream(s.handle, NULL));
}

#define KOKKOSLAPACK_GETRF_CUSOLVER(SCALAR, LAYOUT, MEM_SPACE)                                                        \
  template <bool ETI_SPEC_AVAIL>                                                                                      \
  struct GETRF<                                                                                                       \
      Kokkos::Cuda,                                                                                                   \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                         \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                          \
    using IpivViewType =                                                                                              \
        Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using InfoViewType =                                                                                              \
        Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
                                                                                                                      \
    static void getrf(const Kokkos::Cuda& space, const AViewType& A, const IpivViewType& Ipiv,                        \
                      const InfoViewType& Info) {                                                                     \
      Kokkos::Profiling::pushRegion("KokkosLapack::getrf[TPL_CUSOLVER," #SCALAR "]");                                 \
      getrf_print_specialization<AViewType, IpivViewType, InfoViewType>();                                            \
                                                                                                                      \
      cusolverGetrfWrapper(space, A, Ipiv, Info);                                                                     \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

KOKKOSLAPACK_GETRF_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GETRF_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GETRF_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GETRF_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)

#if defined(KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE)
KOKKOSLAPACK_GETRF_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_GETRF_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_GETRF_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_GETRF_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
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

template <class ExecutionSpace, class AViewType, class IpivViewType, class InfoViewType>
void rocsolverGetrfWrapper(const ExecutionSpace& space, const AViewType& A, const IpivViewType& Ipiv,
                           const InfoViewType& Info) {
  using Scalar = typename AViewType::non_const_value_type;

  using ALayout_t = typename AViewType::array_layout;
  static_assert(std::is_same_v<ALayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - rocsolver getrf: A needs to have a Kokkos::LayoutLeft");

  const rocblas_int m   = static_cast<rocblas_int>(A.extent(0));
  const rocblas_int n   = static_cast<rocblas_int>(A.extent(1));
  const rocblas_int lda = static_cast<rocblas_int>(A.stride(1));

  KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();
  KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));
  if constexpr (std::is_same_v<Scalar, float>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocsolver_sgetrf(s.handle, m, n, A.data(), lda, Ipiv.data(), Info.data()));
  }
  if constexpr (std::is_same_v<Scalar, double>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocsolver_dgetrf(s.handle, m, n, A.data(), lda, Ipiv.data(), Info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocsolver_cgetrf(
        s.handle, m, n, reinterpret_cast<rocblas_float_complex*>(A.data()), lda, Ipiv.data(), Info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocsolver_zgetrf(
        s.handle, m, n, reinterpret_cast<rocblas_double_complex*>(A.data()), lda, Ipiv.data(), Info.data()));
  }
  KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));
}

#define KOKKOSLAPACK_GETRF_ROCSOLVER(SCALAR, LAYOUT, MEM_SPACE)                                                        \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct GETRF<                                                                                                        \
      Kokkos::HIP,                                                                                                     \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,     \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,     \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                           \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                           \
    using IpivViewType =                                                                                               \
        Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;   \
    using InfoViewType =                                                                                               \
        Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;   \
                                                                                                                       \
    static void getrf(const Kokkos::HIP& space, const AViewType& A, const IpivViewType& Ipiv,                          \
                      const InfoViewType& Info) {                                                                      \
      Kokkos::Profiling::pushRegion("KokkosLapack::getrf[TPL_ROCSOLVER," #SCALAR "]");                                 \
      getrf_print_specialization<AViewType, IpivViewType, InfoViewType>();                                             \
                                                                                                                       \
      rocsolverGetrfWrapper(space, A, Ipiv, Info);                                                                     \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSLAPACK_GETRF_ROCSOLVER(float, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_GETRF_ROCSOLVER(double, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_GETRF_ROCSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_GETRF_ROCSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIPSpace)

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER

#endif
