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

#ifndef KOKKOSLAPACK_GESV_TPL_SPEC_DECL_HPP_
#define KOKKOSLAPACK_GESV_TPL_SPEC_DECL_HPP_

namespace KokkosLapack {
namespace Impl {
template <class AViewType, class BViewType, class PViewType>
inline void gesv_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  printf("KokkosLapack::gesv<> TPL MAGMA specialization for < %s , %s, %s >\n", typeid(AViewType).name(),
         typeid(BViewType).name(), typeid(PViewType).name());
#else
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK
  printf("KokkosLapack::gesv<> TPL Lapack specialization for < %s , %s, %s >\n", typeid(AViewType).name(),
         typeid(BViewType).name(), typeid(PViewType).name());
#endif
#endif
#endif
}
}  // namespace Impl
}  // namespace KokkosLapack

// Generic Host side LAPACK (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK
#include <KokkosLapack_Host_tpl.hpp>

namespace KokkosLapack {
namespace Impl {

template <class AViewType, class BViewType, class IPIVViewType>
void lapackGesvWrapper(const AViewType& A, const BViewType& B, const IPIVViewType& IPIV) {
  using Scalar = typename AViewType::non_const_value_type;

  const bool with_pivot = !((IPIV.extent(0) == 0) && (IPIV.data() == nullptr));

  const int N    = static_cast<int>(A.extent(1));
  const int AST  = static_cast<int>(A.stride(1));
  const int LDA  = (AST == 0) ? 1 : AST;
  const int BST  = static_cast<int>(B.stride(1));
  const int LDB  = (BST == 0) ? 1 : BST;
  const int NRHS = static_cast<int>(B.extent(1));

  int info = 0;

  if (with_pivot) {
    if constexpr (Kokkos::ArithTraits<Scalar>::is_complex) {
      using MagType = typename Kokkos::ArithTraits<Scalar>::mag_type;

      HostLapack<std::complex<MagType>>::gesv(N, NRHS, reinterpret_cast<std::complex<MagType>*>(A.data()), LDA,
                                              IPIV.data(), reinterpret_cast<std::complex<MagType>*>(B.data()), LDB,
                                              info);
    } else {
      HostLapack<Scalar>::gesv(N, NRHS, A.data(), LDA, IPIV.data(), B.data(), LDB, info);
    }
  }
}

#define KOKKOSLAPACK_GESV_LAPACK(SCALAR, LAYOUT, EXECSPACE, MEM_SPACE)                                                 \
  template <>                                                                                                          \
  struct GESV<                                                                                                         \
      EXECSPACE,                                                                                                       \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      true,                                                                                                            \
      gesv_eti_spec_avail<EXECSPACE,                                                                                   \
                          Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>,                         \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                          Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>,                         \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                          Kokkos::View<int*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                     \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                             \
    using AViewType =                                                                                                  \
        Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using BViewType =                                                                                                  \
        Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using PViewType =                                                                                                  \
        Kokkos::View<int*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
                                                                                                                       \
    static void gesv(const EXECSPACE& /* space */, const AViewType& A, const BViewType& B, const PViewType& IPIV) {    \
      Kokkos::Profiling::pushRegion("KokkosLapack::gesv[TPL_LAPACK," #SCALAR "]");                                     \
      gesv_print_specialization<AViewType, BViewType, PViewType>();                                                    \
      lapackGesvWrapper(A, B, IPIV);                                                                                   \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#if defined(KOKKOS_ENABLE_SERIAL)
KOKKOSLAPACK_GESV_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_GESV_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_GESV_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_GESV_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
KOKKOSLAPACK_GESV_LAPACK(float, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_GESV_LAPACK(double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_GESV_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_GESV_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined(KOKKOS_ENABLE_THREADS)
KOKKOSLAPACK_GESV_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_GESV_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_GESV_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_GESV_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_LAPACK

// MAGMA
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
#include <KokkosLapack_magma.hpp>

namespace KokkosLapack {
namespace Impl {

template <class ExecSpace, class AViewType, class BViewType, class IPIVViewType>
void magmaGesvWrapper(const ExecSpace& space, const AViewType& A, const BViewType& B, const IPIVViewType& IPIV) {
  using scalar_type = typename AViewType::non_const_value_type;

  Kokkos::Profiling::pushRegion("KokkosLapack::gesv[TPL_MAGMA," + Kokkos::ArithTraits<scalar_type>::name() + "]");
  gesv_print_specialization<AViewType, BViewType, IPIVViewType>();

  const bool with_pivot = !((IPIV.extent(0) == 0) && (IPIV.data() == nullptr));

  magma_int_t N    = static_cast<magma_int_t>(A.extent(1));
  magma_int_t AST  = static_cast<magma_int_t>(A.stride(1));
  magma_int_t LDA  = (AST == 0) ? 1 : AST;
  magma_int_t BST  = static_cast<magma_int_t>(B.stride(1));
  magma_int_t LDB  = (BST == 0) ? 1 : BST;
  magma_int_t NRHS = static_cast<magma_int_t>(B.extent(1));

  KokkosLapack::Impl::MagmaSingleton& s = KokkosLapack::Impl::MagmaSingleton::singleton();
  magma_int_t info                      = 0;

  space.fence();
  if constexpr (std::is_same_v<scalar_type, float>) {
    if (with_pivot) {
      magma_sgesv_gpu(N, NRHS, reinterpret_cast<magmaFloat_ptr>(A.data()), LDA, IPIV.data(),
                      reinterpret_cast<magmaFloat_ptr>(B.data()), LDB, &info);
    } else {
      magma_sgesv_nopiv_gpu(N, NRHS, reinterpret_cast<magmaFloat_ptr>(A.data()), LDA,
                            reinterpret_cast<magmaFloat_ptr>(B.data()), LDB, &info);
    }
  }

  if constexpr (std::is_same_v<scalar_type, double>) {
    if (with_pivot) {
      magma_dgesv_gpu(N, NRHS, reinterpret_cast<magmaDouble_ptr>(A.data()), LDA, IPIV.data(),
                      reinterpret_cast<magmaDouble_ptr>(B.data()), LDB, &info);
    } else {
      magma_dgesv_nopiv_gpu(N, NRHS, reinterpret_cast<magmaDouble_ptr>(A.data()), LDA,
                            reinterpret_cast<magmaDouble_ptr>(B.data()), LDB, &info);
    }
  }

  if constexpr (std::is_same_v<scalar_type, Kokkos::complex<float>>) {
    if (with_pivot) {
      magma_cgesv_gpu(N, NRHS, reinterpret_cast<magmaFloatComplex_ptr>(A.data()), LDA, IPIV.data(),
                      reinterpret_cast<magmaFloatComplex_ptr>(B.data()), LDB, &info);
    } else {
      magma_cgesv_nopiv_gpu(N, NRHS, reinterpret_cast<magmaFloatComplex_ptr>(A.data()), LDA,
                            reinterpret_cast<magmaFloatComplex_ptr>(B.data()), LDB, &info);
    }
  }

  if constexpr (std::is_same_v<scalar_type, Kokkos::complex<double>>) {
    if (with_pivot) {
      magma_zgesv_gpu(N, NRHS, reinterpret_cast<magmaDoubleComplex_ptr>(A.data()), LDA, IPIV.data(),
                      reinterpret_cast<magmaDoubleComplex_ptr>(B.data()), LDB, &info);
    } else {
      magma_zgesv_nopiv_gpu(N, NRHS, reinterpret_cast<magmaDoubleComplex_ptr>(A.data()), LDA,
                            reinterpret_cast<magmaDoubleComplex_ptr>(B.data()), LDB, &info);
    }
  }
  ExecSpace().fence();
  Kokkos::Profiling::popRegion();
}

#define KOKKOSLAPACK_GESV_MAGMA(SCALAR, LAYOUT, MEM_SPACE)                                                             \
  template <>                                                                                                          \
  struct GESV<Kokkos::Cuda,                                                                                            \
              Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                                  \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                   \
              Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                                  \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                   \
              Kokkos::View<magma_int_t*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                   \
              true,                                                                                                    \
              gesv_eti_spec_avail<Kokkos::Cuda,                                                                        \
                                  Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,              \
                                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                               \
                                  Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,              \
                                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                               \
                                  Kokkos::View<magma_int_t*, LAYOUT,                                                   \
                                               Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>,   \
                                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                     \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                          \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                           \
    using BViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                          \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                           \
    using PViewType =                                                                                                  \
        Kokkos::View<magma_int_t*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>,       \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                         \
                                                                                                                       \
    static void gesv(const Kokkos::Cuda& space, const AViewType& A, const BViewType& B, const PViewType& IPIV) {       \
      magmaGesvWrapper(space, A, B, IPIV);                                                                             \
    }                                                                                                                  \
  };

KOKKOSLAPACK_GESV_MAGMA(float, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GESV_MAGMA(double, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GESV_MAGMA(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GESV_MAGMA(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_MAGMA

// CUSOLVER
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSOLVER
#include "KokkosLapack_cusolver.hpp"

namespace KokkosLapack {
namespace Impl {

template <class ExecutionSpace, class IPIVViewType, class AViewType, class BViewType>
void cusolverGesvWrapper(const ExecutionSpace& space, const IPIVViewType& IPIV, const AViewType& A,
                         const BViewType& B) {
  using memory_space = typename AViewType::memory_space;
  using Scalar       = typename BViewType::non_const_value_type;
  using ALayout_t    = typename AViewType::array_layout;
  using BLayout_t    = typename BViewType::array_layout;

  const int m   = A.extent_int(0);
  const int n   = A.extent_int(1);
  const int lda = std::is_same_v<ALayout_t, Kokkos::LayoutRight> ? A.stride(0) : A.stride(1);

  (void)B;

  const int nrhs = B.extent_int(1);
  const int ldb  = std::is_same_v<BLayout_t, Kokkos::LayoutRight> ? B.stride(0) : B.stride(1);
  int lwork      = 0;
  Kokkos::View<int, memory_space> info("getrf info");

  CudaLapackSingleton& s = CudaLapackSingleton::singleton();
  KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnSetStream(s.handle, space.cuda_stream()));
  if constexpr (std::is_same_v<Scalar, float>) {
    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnSgetrf_bufferSize(s.handle, m, n, A.data(), lda, &lwork));
    Kokkos::View<float*, memory_space> Workspace("getrf workspace", lwork);

    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(
        cusolverDnSgetrf(s.handle, m, n, A.data(), lda, Workspace.data(), IPIV.data(), info.data()));

    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(
        cusolverDnSgetrs(s.handle, CUBLAS_OP_N, m, nrhs, A.data(), lda, IPIV.data(), B.data(), ldb, info.data()));
  }
  if constexpr (std::is_same_v<Scalar, double>) {
    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnDgetrf_bufferSize(s.handle, m, n, A.data(), lda, &lwork));
    Kokkos::View<double*, memory_space> Workspace("getrf workspace", lwork);

    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(
        cusolverDnDgetrf(s.handle, m, n, A.data(), lda, Workspace.data(), IPIV.data(), info.data()));

    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(
        cusolverDnDgetrs(s.handle, CUBLAS_OP_N, m, nrhs, A.data(), lda, IPIV.data(), B.data(), ldb, info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(
        cusolverDnCgetrf_bufferSize(s.handle, m, n, reinterpret_cast<cuComplex*>(A.data()), lda, &lwork));
    Kokkos::View<cuComplex*, memory_space> Workspace("getrf workspace", lwork);

    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnCgetrf(s.handle, m, n, reinterpret_cast<cuComplex*>(A.data()), lda,
                                                    reinterpret_cast<cuComplex*>(Workspace.data()), IPIV.data(),
                                                    info.data()));

    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnCgetrs(s.handle, CUBLAS_OP_N, m, nrhs,
                                                    reinterpret_cast<cuComplex*>(A.data()), lda, IPIV.data(),
                                                    reinterpret_cast<cuComplex*>(B.data()), ldb, info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(
        cusolverDnZgetrf_bufferSize(s.handle, m, n, reinterpret_cast<cuDoubleComplex*>(A.data()), lda, &lwork));
    Kokkos::View<cuDoubleComplex*, memory_space> Workspace("getrf workspace", lwork);

    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnZgetrf(s.handle, m, n, reinterpret_cast<cuDoubleComplex*>(A.data()), lda,
                                                    reinterpret_cast<cuDoubleComplex*>(Workspace.data()), IPIV.data(),
                                                    info.data()));

    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnZgetrs(s.handle, CUBLAS_OP_N, m, nrhs,
                                                    reinterpret_cast<cuDoubleComplex*>(A.data()), lda, IPIV.data(),
                                                    reinterpret_cast<cuDoubleComplex*>(B.data()), ldb, info.data()));
  }
  KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnSetStream(s.handle, NULL));
}

#define KOKKOSLAPACK_GESV_CUSOLVER(SCALAR, LAYOUT, MEM_SPACE)                                                         \
  template <>                                                                                                         \
  struct GESV<                                                                                                        \
      Kokkos::Cuda,                                                                                                   \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      true,                                                                                                           \
      gesv_eti_spec_avail<Kokkos::Cuda,                                                                               \
                          Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                     \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                      \
                          Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                     \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                      \
                          Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                         \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                            \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                         \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                          \
    using BViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                         \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                          \
    using PViewType =                                                                                                 \
        Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
                                                                                                                      \
    static void gesv(const Kokkos::Cuda& space, const AViewType& A, const BViewType& B, const PViewType& IPIV) {      \
      Kokkos::Profiling::pushRegion("KokkosLapack::gesv[TPL_CUSOLVER," #SCALAR "]");                                  \
      gesv_print_specialization<AViewType, BViewType, PViewType>();                                                   \
                                                                                                                      \
      cusolverGesvWrapper(space, IPIV, A, B);                                                                         \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

KOKKOSLAPACK_GESV_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GESV_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GESV_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GESV_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)

#if defined(KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE)
KOKKOSLAPACK_GESV_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_GESV_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_GESV_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_GESV_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
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

template <class ExecutionSpace, class IPIVViewType, class AViewType, class BViewType>
void rocsolverGesvWrapper(const ExecutionSpace& space, const IPIVViewType& IPIV, const AViewType& A,
                          const BViewType& B) {
  using Scalar    = typename BViewType::non_const_value_type;
  using ALayout_t = typename AViewType::array_layout;
  using BLayout_t = typename BViewType::array_layout;

  const rocblas_int N    = static_cast<rocblas_int>(A.extent(0));
  const rocblas_int nrhs = static_cast<rocblas_int>(B.extent(1));
  const rocblas_int lda  = std::is_same_v<ALayout_t, Kokkos::LayoutRight> ? A.stride(0) : A.stride(1);
  const rocblas_int ldb  = std::is_same_v<BLayout_t, Kokkos::LayoutRight> ? B.stride(0) : B.stride(1);
  Kokkos::View<rocblas_int, ExecutionSpace> info("rocsolver info");

  KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();
  KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));
  if constexpr (std::is_same_v<Scalar, float>) {
    KOKKOS_ROCBLAS_SAFE_CALL_IMPL(
        rocsolver_sgesv(s.handle, N, nrhs, A.data(), lda, IPIV.data(), B.data(), ldb, info.data()));
  }
  if constexpr (std::is_same_v<Scalar, double>) {
    KOKKOS_ROCBLAS_SAFE_CALL_IMPL(
        rocsolver_dgesv(s.handle, N, nrhs, A.data(), lda, IPIV.data(), B.data(), ldb, info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
    KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocsolver_cgesv(s.handle, N, nrhs, reinterpret_cast<rocblas_float_complex*>(A.data()),
                                                  lda, IPIV.data(), reinterpret_cast<rocblas_float_complex*>(B.data()),
                                                  ldb, info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
    KOKKOS_ROCBLAS_SAFE_CALL_IMPL(
        rocsolver_zgesv(s.handle, N, nrhs, reinterpret_cast<rocblas_double_complex*>(A.data()), lda, IPIV.data(),
                        reinterpret_cast<rocblas_double_complex*>(B.data()), ldb, info.data()));
  }
  KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));
}

#define KOKKOSLAPACK_GESV_ROCSOLVER(SCALAR, LAYOUT, MEM_SPACE)                                                         \
  template <>                                                                                                          \
  struct GESV<                                                                                                         \
      Kokkos::HIP,                                                                                                     \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<rocblas_int*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      true,                                                                                                            \
      gesv_eti_spec_avail<Kokkos::HIP,                                                                                 \
                          Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                       \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                          Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                       \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                          Kokkos::View<rocblas_int*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                   \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                             \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                           \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                           \
    using BViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                           \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                           \
    using PViewType = Kokkos::View<rocblas_int*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                       \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                           \
                                                                                                                       \
    static void gesv(const Kokkos::HIP& space, const AViewType& A, const BViewType& B, const PViewType& IPIV) {        \
      Kokkos::Profiling::pushRegion("KokkosLapack::gesv[TPL_ROCSOLVER," #SCALAR "]");                                  \
      gesv_print_specialization<AViewType, BViewType, PViewType>();                                                    \
                                                                                                                       \
      rocsolverGesvWrapper(space, IPIV, A, B);                                                                         \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSLAPACK_GESV_ROCSOLVER(float, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_GESV_ROCSOLVER(double, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_GESV_ROCSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_GESV_ROCSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIPSpace)

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER

#endif
