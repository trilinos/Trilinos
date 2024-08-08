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

#ifndef KOKKOSBLAS3_GEMM_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS3_GEMM_TPL_SPEC_DECL_HPP_

#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS)
#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS3_XGEMM_BLAS(SCALAR_TYPE, BASE_SCALAR_TYPE, LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL) \
  template <class ExecSpace>                                                                                        \
  struct GEMM<ExecSpace,                                                                                            \
              Kokkos::View<const SCALAR_TYPE**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                      \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                               \
              Kokkos::View<const SCALAR_TYPE**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                      \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                               \
              Kokkos::View<SCALAR_TYPE**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>,                            \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                               \
              true, ETI_SPEC_AVAIL> {                                                                               \
    typedef SCALAR_TYPE SCALAR;                                                                                     \
    typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                             \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                  \
        AViewType;                                                                                                  \
    typedef Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                             \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                  \
        BViewType;                                                                                                  \
    typedef Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>,                                   \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                  \
        CViewType;                                                                                                  \
                                                                                                                    \
    static void gemm(const ExecSpace& /* space*/, const char transA[], const char transB[],                         \
                     typename AViewType::const_value_type& alpha, const AViewType& A, const BViewType& B,           \
                     typename CViewType::const_value_type& beta, const CViewType& C) {                              \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemm[TPL_BLAS," #SCALAR_TYPE "]");                                 \
      const bool A_t = (transA[0] != 'N') && (transA[0] != 'n');                                                    \
      const KK_INT M = C.extent(0);                                                                                 \
      const KK_INT N = C.extent(1);                                                                                 \
      const KK_INT K = A.extent(A_t ? 0 : 1);                                                                       \
                                                                                                                    \
      bool A_is_lr = std::is_same<Kokkos::LayoutRight, LAYOUTA>::value;                                             \
      bool B_is_lr = std::is_same<Kokkos::LayoutRight, LAYOUTB>::value;                                             \
      bool C_is_lr = std::is_same<Kokkos::LayoutRight, LAYOUTC>::value;                                             \
                                                                                                                    \
      const KK_INT AST = A_is_lr ? A.stride(0) : A.stride(1), LDA = AST == 0 ? 1 : AST;                             \
      const KK_INT BST = B_is_lr ? B.stride(0) : B.stride(1), LDB = BST == 0 ? 1 : BST;                             \
      const KK_INT CST = C_is_lr ? C.stride(0) : C.stride(1), LDC = CST == 0 ? 1 : CST;                             \
                                                                                                                    \
      const BASE_SCALAR_TYPE alpha_val = alpha, beta_val = beta;                                                    \
      if (!A_is_lr && !B_is_lr && !C_is_lr)                                                                         \
        HostBlas<BASE_SCALAR_TYPE>::gemm(transA[0], transB[0], M, N, K, alpha_val,                                  \
                                         reinterpret_cast<const BASE_SCALAR_TYPE*>(A.data()), LDA,                  \
                                         reinterpret_cast<const BASE_SCALAR_TYPE*>(B.data()), LDB, beta_val,        \
                                         reinterpret_cast<BASE_SCALAR_TYPE*>(C.data()), LDC);                       \
      if (A_is_lr && B_is_lr && C_is_lr)                                                                            \
        HostBlas<BASE_SCALAR_TYPE>::gemm(transB[0], transA[0], N, M, K, alpha_val,                                  \
                                         reinterpret_cast<const BASE_SCALAR_TYPE*>(B.data()), LDB,                  \
                                         reinterpret_cast<const BASE_SCALAR_TYPE*>(A.data()), LDA, beta_val,        \
                                         reinterpret_cast<BASE_SCALAR_TYPE*>(C.data()), LDC);                       \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS3_DGEMM_BLAS(LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS3_XGEMM_BLAS(double, double, LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS3_SGEMM_BLAS(LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS3_XGEMM_BLAS(float, float, LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS3_ZGEMM_BLAS(LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL)                          \
  KOKKOSBLAS3_XGEMM_BLAS(Kokkos::complex<double>, std::complex<double>, LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, \
                         ETI_SPEC_AVAIL)

#define KOKKOSBLAS3_CGEMM_BLAS(LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL)                        \
  KOKKOSBLAS3_XGEMM_BLAS(Kokkos::complex<float>, std::complex<float>, LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, \
                         ETI_SPEC_AVAIL)

KOKKOSBLAS3_DGEMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_DGEMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_DGEMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_DGEMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

KOKKOSBLAS3_SGEMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_SGEMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_SGEMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_SGEMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

KOKKOSBLAS3_ZGEMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_ZGEMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_ZGEMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_ZGEMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

KOKKOSBLAS3_CGEMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_CGEMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_CGEMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_CGEMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

// cuBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
#include <KokkosBlas_tpl_spec.hpp>
#include <KokkosBlas3_gemm_dotbased_impl.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS3_XGEMM_CUBLAS(SCALAR_TYPE, CUDA_SCALAR_TYPE, CUBLAS_FN, LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE,      \
                                 ETI_SPEC_AVAIL)                                                                      \
  template <class ExecSpace>                                                                                          \
  struct GEMM<ExecSpace,                                                                                              \
              Kokkos::View<const SCALAR_TYPE**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                        \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                 \
              Kokkos::View<const SCALAR_TYPE**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                        \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                 \
              Kokkos::View<SCALAR_TYPE**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>,                              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                 \
              true, ETI_SPEC_AVAIL> {                                                                                 \
    typedef SCALAR_TYPE SCALAR;                                                                                       \
    typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        AViewType;                                                                                                    \
    typedef Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        BViewType;                                                                                                    \
    typedef Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>,                                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        CViewType;                                                                                                    \
                                                                                                                      \
    static void gemm(const ExecSpace& space, const char transA[], const char transB[],                                \
                     typename AViewType::const_value_type& alpha, const AViewType& A, const BViewType& B,             \
                     typename CViewType::const_value_type& beta, const CViewType& C) {                                \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemm[TPL_CUBLAS," #SCALAR_TYPE "]");                                 \
      const bool A_t = (transA[0] != 'N') && (transA[0] != 'n');                                                      \
      const int M    = static_cast<int>(C.extent(0));                                                                 \
      const int N    = static_cast<int>(C.extent(1));                                                                 \
      const int K    = static_cast<int>(A.extent(A_t ? 0 : 1));                                                       \
                                                                                                                      \
      bool A_is_lr = std::is_same<Kokkos::LayoutRight, LAYOUTA>::value;                                               \
      bool B_is_lr = std::is_same<Kokkos::LayoutRight, LAYOUTB>::value;                                               \
      bool C_is_lr = std::is_same<Kokkos::LayoutRight, LAYOUTC>::value;                                               \
                                                                                                                      \
      const int AST = A_is_lr ? A.stride(0) : A.stride(1), LDA = AST == 0 ? 1 : AST;                                  \
      const int BST = B_is_lr ? B.stride(0) : B.stride(1), LDB = BST == 0 ? 1 : BST;                                  \
      const int CST = C_is_lr ? C.stride(0) : C.stride(1), LDC = CST == 0 ? 1 : CST;                                  \
                                                                                                                      \
      cublasOperation_t transa = trans_mode_kk_to_cublas(transA);                                                     \
      cublasOperation_t transb = trans_mode_kk_to_cublas(transB);                                                     \
                                                                                                                      \
      constexpr int numDotsLayoutLeftThreshold  = 1600;                                                               \
      constexpr int numDotsLayoutRightThreshold = 100;                                                                \
      if ((!A_is_lr && transa != CUBLAS_OP_N && transb == CUBLAS_OP_N && M * N < numDotsLayoutLeftThreshold) ||       \
          (A_is_lr && transa != CUBLAS_OP_N && transb == CUBLAS_OP_N && M * N < numDotsLayoutRightThreshold)) {       \
        DotBasedGEMM<ExecSpace, AViewType, BViewType, CViewType> gemm(alpha, A, B, beta, C);                          \
        bool conjT = (std::is_same<SCALAR, double>::value || std::is_same<SCALAR, float>::value)                      \
                         ? false                                                                                      \
                         : (transa == CUBLAS_OP_C ? true : false);                                                    \
        gemm.run(space, conjT);                                                                                       \
      } else {                                                                                                        \
        KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                    \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                 \
        if (!A_is_lr && !B_is_lr && !C_is_lr)                                                                         \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(CUBLAS_FN(                                                                     \
              s.handle, transa, transb, M, N, K, reinterpret_cast<const CUDA_SCALAR_TYPE*>(&alpha),                   \
              reinterpret_cast<const CUDA_SCALAR_TYPE*>(A.data()), LDA,                                               \
              reinterpret_cast<const CUDA_SCALAR_TYPE*>(B.data()), LDB,                                               \
              reinterpret_cast<const CUDA_SCALAR_TYPE*>(&beta), reinterpret_cast<CUDA_SCALAR_TYPE*>(C.data()), LDC)); \
        if (A_is_lr && B_is_lr && C_is_lr)                                                                            \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(CUBLAS_FN(                                                                     \
              s.handle, transb, transa, N, M, K, reinterpret_cast<const CUDA_SCALAR_TYPE*>(&alpha),                   \
              reinterpret_cast<const CUDA_SCALAR_TYPE*>(B.data()), LDB,                                               \
              reinterpret_cast<const CUDA_SCALAR_TYPE*>(A.data()), LDA,                                               \
              reinterpret_cast<const CUDA_SCALAR_TYPE*>(&beta), reinterpret_cast<CUDA_SCALAR_TYPE*>(C.data()), LDC)); \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                                \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS3_DGEMM_CUBLAS(LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS3_XGEMM_CUBLAS(double, double, cublasDgemm, LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS3_SGEMM_CUBLAS(LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS3_XGEMM_CUBLAS(float, float, cublasSgemm, LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS3_ZGEMM_CUBLAS(LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL)                       \
  KOKKOSBLAS3_XGEMM_CUBLAS(Kokkos::complex<double>, cuDoubleComplex, cublasZgemm, LAYOUTA, LAYOUTB, LAYOUTC, \
                           MEM_SPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS3_CGEMM_CUBLAS(LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL)                           \
  KOKKOSBLAS3_XGEMM_CUBLAS(Kokkos::complex<float>, cuComplex, cublasCgemm, LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, \
                           ETI_SPEC_AVAIL)

KOKKOSBLAS3_DGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_DGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_DGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_DGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_DGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_DGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_DGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_DGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS3_SGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_SGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_SGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_SGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_SGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_SGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_SGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_SGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS3_ZGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_ZGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_ZGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_ZGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_ZGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_ZGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_ZGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_ZGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS3_CGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_CGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_CGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_CGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_CGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_CGEMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_CGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_CGEMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

// rocBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)
#include <KokkosBlas_tpl_spec.hpp>
#include <KokkosBlas3_gemm_dotbased_impl.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS3_XGEMM_ROCBLAS(SCALAR_TYPE, ROCBLAS_SCALAR_TYPE, ROCBLAS_FN, LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL) \
  template <class ExecSpace>                                                                                       \
  struct GEMM<ExecSpace,                                                                                           \
              Kokkos::View<const SCALAR_TYPE**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>,                      \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                              \
              Kokkos::View<const SCALAR_TYPE**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>,                      \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                              \
              Kokkos::View<SCALAR_TYPE**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>,                            \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                              \
              true, ETI_SPEC_AVAIL> {                                                                              \
    typedef SCALAR_TYPE SCALAR;                                                                                    \
    typedef Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>,                             \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                 \
        AViewType;                                                                                                 \
    typedef Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>,                             \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                 \
        BViewType;                                                                                                 \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>,                                   \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                 \
        CViewType;                                                                                                 \
                                                                                                                   \
    static void gemm(const typename CViewType::execution_space& space, const char transA[], const char transB[],   \
                     typename AViewType::const_value_type& alpha, const AViewType& A, const BViewType& B,          \
                     typename CViewType::const_value_type& beta, const CViewType& C) {                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemm[TPL_ROCBLAS," #SCALAR_TYPE "]");                             \
                                                                                                                   \
      const bool A_t = (transA[0] != 'N') && (transA[0] != 'n');                                                   \
      const int M    = static_cast<int>(C.extent(0));                                                              \
      const int N    = static_cast<int>(C.extent(1));                                                              \
      const int K    = static_cast<int>(A.extent(A_t ? 0 : 1));                                                    \
                                                                                                                   \
      bool is_lr = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;                                               \
                                                                                                                   \
      const int AST = is_lr ? A.stride(0) : A.stride(1), LDA = AST == 0 ? 1 : AST;                                 \
      const int BST = is_lr ? B.stride(0) : B.stride(1), LDB = BST == 0 ? 1 : BST;                                 \
      const int CST = is_lr ? C.stride(0) : C.stride(1), LDC = CST == 0 ? 1 : CST;                                 \
                                                                                                                   \
      rocblas_operation transa = trans_mode_kk_to_rocblas(transA);                                                 \
      rocblas_operation transb = trans_mode_kk_to_rocblas(transB);                                                 \
                                                                                                                   \
      constexpr int numDotsLayoutLeftThreshold  = 1600;                                                            \
      constexpr int numDotsLayoutRightThreshold = 100;                                                             \
      if ((!is_lr && transa != rocblas_operation_none && transb == rocblas_operation_none &&                       \
           M * N < numDotsLayoutLeftThreshold) ||                                                                  \
          (is_lr && transa != rocblas_operation_none && transb == rocblas_operation_none &&                        \
           M * N < numDotsLayoutRightThreshold)) {                                                                 \
        DotBasedGEMM<ExecSpace, AViewType, BViewType, CViewType> gemm(alpha, A, B, beta, C);                       \
        bool conjT = (std::is_same<SCALAR, double>::value || std::is_same<SCALAR, float>::value)                   \
                         ? false                                                                                   \
                         : (transa == rocblas_operation_conjugate_transpose ? true : false);                       \
        gemm.run(space, conjT);                                                                                    \
      } else {                                                                                                     \
        KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                   \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));                           \
        if (!is_lr)                                                                                                \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(ROCBLAS_FN(s.handle, transa, transb, M, N, K,                              \
                                                   reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(&alpha),           \
                                                   reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(A.data()), LDA,    \
                                                   reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(B.data()), LDB,    \
                                                   reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(&beta),            \
                                                   reinterpret_cast<ROCBLAS_SCALAR_TYPE*>(C.data()), LDC));        \
        else                                                                                                       \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(ROCBLAS_FN(s.handle, transb, transa, N, M, K,                              \
                                                   reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(&alpha),           \
                                                   reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(B.data()), LDB,    \
                                                   reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(A.data()), LDA,    \
                                                   reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(&beta),            \
                                                   reinterpret_cast<ROCBLAS_SCALAR_TYPE*>(C.data()), LDC));        \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));                                         \
      }                                                                                                            \
      Kokkos::Profiling::popRegion();                                                                              \
    }                                                                                                              \
  };

#define KOKKOSBLAS3_DGEMM_ROCBLAS(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS3_XGEMM_ROCBLAS(double, double, rocblas_dgemm, LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS3_SGEMM_ROCBLAS(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS3_XGEMM_ROCBLAS(float, float, rocblas_sgemm, LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS3_ZGEMM_ROCBLAS(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)                                           \
  KOKKOSBLAS3_XGEMM_ROCBLAS(Kokkos::complex<double>, rocblas_double_complex, rocblas_zgemm, LAYOUT, MEM_SPACE, \
                            ETI_SPEC_AVAIL)

#define KOKKOSBLAS3_CGEMM_ROCBLAS(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)                                         \
  KOKKOSBLAS3_XGEMM_ROCBLAS(Kokkos::complex<float>, rocblas_float_complex, rocblas_cgemm, LAYOUT, MEM_SPACE, \
                            ETI_SPEC_AVAIL)

KOKKOSBLAS3_DGEMM_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIPSpace, true)
KOKKOSBLAS3_DGEMM_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIPSpace, false)
KOKKOSBLAS3_DGEMM_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIPSpace, true)
KOKKOSBLAS3_DGEMM_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIPSpace, false)

KOKKOSBLAS3_SGEMM_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIPSpace, true)
KOKKOSBLAS3_SGEMM_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIPSpace, false)
KOKKOSBLAS3_SGEMM_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIPSpace, true)
KOKKOSBLAS3_SGEMM_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIPSpace, false)

KOKKOSBLAS3_ZGEMM_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIPSpace, true)
KOKKOSBLAS3_ZGEMM_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIPSpace, false)
KOKKOSBLAS3_ZGEMM_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIPSpace, true)
KOKKOSBLAS3_ZGEMM_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIPSpace, false)

KOKKOSBLAS3_CGEMM_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIPSpace, true)
KOKKOSBLAS3_CGEMM_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIPSpace, false)
KOKKOSBLAS3_CGEMM_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIPSpace, true)
KOKKOSBLAS3_CGEMM_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIPSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCBLAS

#endif
