// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_GEGQR_TPL_SPEC_DECL_HPP_
#define KOKKOSLAPACK_GEGQR_TPL_SPEC_DECL_HPP_

#include <iostream>
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosLapack {
namespace Impl {
template <class AViewType, class TauViewType, class InfoViewType>
inline void gegqr_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  printf("KokkosLapack::gegqr<> TPL MAGMA specialization for < %s , %s, %s >\n", typeid(AViewType).name(),
         typeid(TauViewType).name(), typeid(InfoViewType).name());
#else
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK
  printf("KokkosLapack::gegqr<> TPL Lapack specialization for < %s , %s, %s >\n", typeid(AViewType).name(),
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
void lapackGegqrWrapper(const int k, const AViewType& A, const TauViewType& Tau, const InfoViewType& Info) {
  using memory_space = typename AViewType::memory_space;
  using Scalar       = typename AViewType::non_const_value_type;
  using ALayout_t    = typename AViewType::array_layout;
  static_assert(std::is_same_v<ALayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - gegqr: A needs to have a Kokkos::LayoutLeft");
  const int m   = A.extent_int(0);
  const int n   = A.extent_int(1);
  const int lda = A.stride(1);

  int lwork = -1;
  // work needs to be at least length 1 to store the returned value for lwork
  Kokkos::View<Scalar*, memory_space> work("work array", 1);

  if constexpr (KokkosKernels::ArithTraits<Scalar>::is_complex) {
    using MagType = typename KokkosKernels::ArithTraits<Scalar>::mag_type;

    HostLapack<std::complex<MagType>>::gegqr(m, n, k, reinterpret_cast<std::complex<MagType>*>(A.data()), lda,
                                             reinterpret_cast<std::complex<MagType>*>(Tau.data()),
                                             reinterpret_cast<std::complex<MagType>*>(work.data()), lwork, Info.data());

    if (Info[0] < 0) return;

    lwork = static_cast<int>(work(0).real());

    work = Kokkos::View<Scalar*, memory_space>("gegqr work buffer", lwork);

    HostLapack<std::complex<MagType>>::gegqr(m, n, k, reinterpret_cast<std::complex<MagType>*>(A.data()), lda,
                                             reinterpret_cast<std::complex<MagType>*>(Tau.data()),
                                             reinterpret_cast<std::complex<MagType>*>(work.data()), lwork, Info.data());
  } else {
    HostLapack<Scalar>::gegqr(m, n, k, A.data(), lda, Tau.data(), work.data(), lwork, Info.data());

    if (Info[0] < 0) return;

    lwork = static_cast<int>(work(0));

    work = Kokkos::View<Scalar*, memory_space>("gegqr work buffer", lwork);

    HostLapack<Scalar>::gegqr(m, n, k, A.data(), lda, Tau.data(), work.data(), lwork, Info.data());
  }
}

#define KOKKOSLAPACK_GEGQR_LAPACK(SCALAR, LAYOUT, EXECSPACE, MEM_SPACE)                                                \
  template <>                                                                                                          \
  struct GEGQR<                                                                                                        \
      EXECSPACE,                                                                                                       \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, \
      gegqr_eti_spec_avail<EXECSPACE,                                                                                  \
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
    static void gegqr(const EXECSPACE& /* space */, const int k, const AViewType& A, const TauViewType& Tau,           \
                      const InfoViewType& Info) {                                                                      \
      Kokkos::Profiling::pushRegion("KokkosLapack::gegqr[TPL_LAPACK," #SCALAR "]");                                    \
      gegqr_print_specialization<AViewType, TauViewType, InfoViewType>();                                              \
      lapackGegqrWrapper(k, A, Tau, Info);                                                                             \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#if defined(KOKKOS_ENABLE_SERIAL)
KOKKOSLAPACK_GEGQR_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_GEGQR_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_GEGQR_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_GEGQR_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
KOKKOSLAPACK_GEGQR_LAPACK(float, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_GEGQR_LAPACK(double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_GEGQR_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_GEGQR_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined(KOKKOS_ENABLE_THREADS)
KOKKOSLAPACK_GEGQR_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_GEGQR_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_GEGQR_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_GEGQR_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
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
void cusolverGegqrWrapper(const ExecutionSpace& space, const int k, const AViewType& A, const TauViewType& Tau,
                          const InfoViewType& Info) {
  using memory_space = typename AViewType::memory_space;
  using Scalar       = typename AViewType::non_const_value_type;

  using ALayout_t = typename AViewType::array_layout;
  static_assert(std::is_same_v<ALayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - cusolver {or,un}mqr: A needs to have a Kokkos::LayoutLeft");
  const int m   = A.extent_int(0);
  const int n   = A.extent_int(1);
  const int lda = A.stride(1);
  int lwork     = 0;

  CudaLapackSingleton& s = CudaLapackSingleton::singleton();
  KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnSetStream(s.handle, space.cuda_stream()));
  if constexpr (std::is_same_v<Scalar, float>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnSorgqr_bufferSize(s.handle, m, n, k, A.data(), lda, Tau.data(), &lwork));
    Kokkos::View<float*, memory_space> Workspace("cusolver sorgqr workspace", lwork);

    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnSorgqr(s.handle, m, n, k, A.data(), lda, Tau.data(), Workspace.data(), lwork, Info.data()));
  }
  if constexpr (std::is_same_v<Scalar, double>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnDorgqr_bufferSize(s.handle, m, n, k, A.data(), lda, Tau.data(), &lwork));
    Kokkos::View<double*, memory_space> Workspace("cusolver dorgqr workspace", lwork);

    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnDorgqr(s.handle, m, n, k, A.data(), lda, Tau.data(), Workspace.data(), lwork, Info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnCungqr_bufferSize(s.handle, m, n, k,
                                                                     reinterpret_cast<cuComplex*>(A.data()), lda,
                                                                     reinterpret_cast<cuComplex*>(Tau.data()), &lwork));
    Kokkos::View<cuComplex*, memory_space> Workspace("cusolver cungqr workspace", lwork);

    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnCungqr(
        s.handle, m, n, k, reinterpret_cast<cuComplex*>(A.data()), lda, reinterpret_cast<cuComplex*>(Tau.data()),
        reinterpret_cast<cuComplex*>(Workspace.data()), lwork, Info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnZungqr_bufferSize(s.handle, m, n, k, reinterpret_cast<cuDoubleComplex*>(A.data()), lda,
                                    reinterpret_cast<cuDoubleComplex*>(Tau.data()), &lwork));
    Kokkos::View<cuDoubleComplex*, memory_space> Workspace("cusolver zungqr workspace", lwork);

    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnZungqr(s.handle, m, n, k, reinterpret_cast<cuDoubleComplex*>(A.data()), lda,
                         reinterpret_cast<cuDoubleComplex*>(Tau.data()),
                         reinterpret_cast<cuDoubleComplex*>(Workspace.data()), lwork, Info.data()));
  }
  KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnSetStream(s.handle, NULL));
}

#define KOKKOSLAPACK_GEGQR_CUSOLVER(SCALAR, LAYOUT, MEM_SPACE)                                                         \
  template <>                                                                                                          \
  struct GEGQR<                                                                                                        \
      Kokkos::Cuda,                                                                                                    \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                                          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      true,                                                                                                            \
      gegqr_eti_spec_avail<Kokkos::Cuda,                                                                               \
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
    static void gegqr(const Kokkos::Cuda& space, const int k, const AViewType& A, const TauViewType& Tau,              \
                      const InfoViewType& Info) {                                                                      \
      Kokkos::Profiling::pushRegion("KokkosLapack::gegqr[TPL_CUSOLVER," #SCALAR "]");                                  \
      gegqr_print_specialization<AViewType, TauViewType, InfoViewType>();                                              \
                                                                                                                       \
      cusolverGegqrWrapper(space, k, A, Tau, Info);                                                                    \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSLAPACK_GEGQR_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GEGQR_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GEGQR_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GEGQR_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)

#if defined(KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE)
KOKKOSLAPACK_GEGQR_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_GEGQR_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_GEGQR_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_GEGQR_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
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
void rocsolverGegqrWrapper(const ExecutionSpace& space, const int k, const AViewType& A, const TauViewType& Tau,
                           const InfoViewType& Info) {
  using Scalar = typename AViewType::non_const_value_type;

  using ALayout_t = typename AViewType::array_layout;
  static_assert(std::is_same_v<ALayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - rocsolver {un,or}gqr: A needs to have a Kokkos::LayoutLeft");
  const rocblas_int m   = static_cast<rocblas_int>(A.extent(0));
  const rocblas_int n   = static_cast<rocblas_int>(A.extent(1));
  const rocblas_int lda = static_cast<rocblas_int>(A.stride(1));

  KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();
  KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));
  if constexpr (std::is_same_v<Scalar, float>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocsolver_sorgqr(s.handle, m, n, k, A.data(), lda, Tau.data()));
  }
  if constexpr (std::is_same_v<Scalar, double>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocsolver_dorgqr(s.handle, m, n, k, A.data(), lda, Tau.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocsolver_cungqr(s.handle, m, n, k,
                                                       reinterpret_cast<rocblas_float_complex*>(A.data()), lda,
                                                       reinterpret_cast<rocblas_float_complex*>(Tau.data())));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocsolver_zungqr(s.handle, m, n, k,
                                                       reinterpret_cast<rocblas_double_complex*>(A.data()), lda,
                                                       reinterpret_cast<rocblas_double_complex*>(Tau.data())));
  }
  Kokkos::deep_copy(Info, 0);  // Success
  KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));
}

#define KOKKOSLAPACK_GEGQR_ROCSOLVER(SCALAR, LAYOUT, MEM_SPACE)                                                        \
  template <>                                                                                                          \
  struct GEGQR<                                                                                                        \
      Kokkos::HIP,                                                                                                     \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,  \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,     \
      true,                                                                                                            \
      gegqr_eti_spec_avail<Kokkos::HIP,                                                                                \
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
    static void gegqr(const Kokkos::HIP& space, const int k, const AViewType& A, const TauViewType& Tau,               \
                      const InfoViewType& Info) {                                                                      \
      Kokkos::Profiling::pushRegion("KokkosLapack::gegqr[TPL_ROCSOLVER," #SCALAR "]");                                 \
      gegqr_print_specialization<AViewType, TauViewType, InfoViewType>();                                              \
                                                                                                                       \
      rocsolverGegqrWrapper(space, k, A, Tau, Info);                                                                   \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSLAPACK_GEGQR_ROCSOLVER(float, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_GEGQR_ROCSOLVER(double, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_GEGQR_ROCSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_GEGQR_ROCSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIPSpace)

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER

#endif
