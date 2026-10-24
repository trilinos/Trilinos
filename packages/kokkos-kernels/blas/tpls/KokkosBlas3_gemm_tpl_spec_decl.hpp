// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS3_GEMM_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS3_GEMM_TPL_SPEC_DECL_HPP_

#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS)
#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS3_XGEMM_BLAS(SCALAR_TYPE, BASE_SCALAR_TYPE, LAYOUT)                                            \
  template <typename ExecSpace, bool ETI_SPEC_AVAIL>                                                             \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                \
  struct GEMM<ExecSpace,                                                                                         \
              Kokkos::View<const SCALAR_TYPE**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
              Kokkos::View<const SCALAR_TYPE**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
              Kokkos::View<SCALAR_TYPE**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,    \
              ETI_SPEC_AVAIL> {                                                                                  \
    using SCALAR    = SCALAR_TYPE;                                                                               \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >; \
    using BViewType = Kokkos::View<const SCALAR**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >; \
    using CViewType = Kokkos::View<SCALAR**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;       \
                                                                                                                 \
    static void gemm(const ExecSpace& /* space*/, const char transA[], const char transB[],                      \
                     typename AViewType::const_value_type& alpha, const AViewType& A, const BViewType& B,        \
                     typename CViewType::const_value_type& beta, const CViewType& C) {                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemm[TPL_BLAS," #SCALAR_TYPE "]");                              \
      const bool A_t = (transA[0] != 'N') && (transA[0] != 'n');                                                 \
      const KK_INT M = C.extent(0);                                                                              \
      const KK_INT N = C.extent(1);                                                                              \
      const KK_INT K = A.extent(A_t ? 0 : 1);                                                                    \
                                                                                                                 \
      bool A_is_lr = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;                                           \
      bool B_is_lr = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;                                           \
      bool C_is_lr = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;                                           \
                                                                                                                 \
      const KK_INT AST = A_is_lr ? A.stride(0) : A.stride(1), LDA = AST == 0 ? 1 : AST;                          \
      const KK_INT BST = B_is_lr ? B.stride(0) : B.stride(1), LDB = BST == 0 ? 1 : BST;                          \
      const KK_INT CST = C_is_lr ? C.stride(0) : C.stride(1), LDC = CST == 0 ? 1 : CST;                          \
                                                                                                                 \
      const BASE_SCALAR_TYPE alpha_val = alpha, beta_val = beta;                                                 \
      if (!A_is_lr && !B_is_lr && !C_is_lr)                                                                      \
        HostBlas<BASE_SCALAR_TYPE>::gemm(transA[0], transB[0], M, N, K, alpha_val,                               \
                                         reinterpret_cast<const BASE_SCALAR_TYPE*>(A.data()), LDA,               \
                                         reinterpret_cast<const BASE_SCALAR_TYPE*>(B.data()), LDB, beta_val,     \
                                         reinterpret_cast<BASE_SCALAR_TYPE*>(C.data()), LDC);                    \
      if (A_is_lr && B_is_lr && C_is_lr)                                                                         \
        HostBlas<BASE_SCALAR_TYPE>::gemm(transB[0], transA[0], N, M, K, alpha_val,                               \
                                         reinterpret_cast<const BASE_SCALAR_TYPE*>(B.data()), LDB,               \
                                         reinterpret_cast<const BASE_SCALAR_TYPE*>(A.data()), LDA, beta_val,     \
                                         reinterpret_cast<BASE_SCALAR_TYPE*>(C.data()), LDC);                    \
      Kokkos::Profiling::popRegion();                                                                            \
    }                                                                                                            \
  };

#define KOKKOSBLAS3_DGEMM_BLAS(LAYOUT) KOKKOSBLAS3_XGEMM_BLAS(double, double, LAYOUT)

#define KOKKOSBLAS3_SGEMM_BLAS(LAYOUT) KOKKOSBLAS3_XGEMM_BLAS(float, float, LAYOUT)

#define KOKKOSBLAS3_ZGEMM_BLAS(LAYOUT) KOKKOSBLAS3_XGEMM_BLAS(Kokkos::complex<double>, std::complex<double>, LAYOUT)

#define KOKKOSBLAS3_CGEMM_BLAS(LAYOUT) KOKKOSBLAS3_XGEMM_BLAS(Kokkos::complex<float>, std::complex<float>, LAYOUT)

KOKKOSBLAS3_DGEMM_BLAS(Kokkos::LayoutLeft)
KOKKOSBLAS3_DGEMM_BLAS(Kokkos::LayoutRight)
KOKKOSBLAS3_SGEMM_BLAS(Kokkos::LayoutLeft)
KOKKOSBLAS3_SGEMM_BLAS(Kokkos::LayoutRight)
KOKKOSBLAS3_ZGEMM_BLAS(Kokkos::LayoutLeft)
KOKKOSBLAS3_ZGEMM_BLAS(Kokkos::LayoutRight)
KOKKOSBLAS3_CGEMM_BLAS(Kokkos::LayoutLeft)
KOKKOSBLAS3_CGEMM_BLAS(Kokkos::LayoutRight)

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

// cuBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
#include <KokkosBlas_tpl_spec.hpp>
#include <KokkosBlas3_gemm_dotbased_impl.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS3_XGEMM_CUBLAS(SCALAR_TYPE, CUDA_SCALAR_TYPE, CUBLAS_FN, LAYOUT)                                    \
  template <bool ETI_SPEC_AVAIL>                                                                                      \
  struct GEMM<Kokkos::Cuda,                                                                                           \
              Kokkos::View<const SCALAR_TYPE**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,      \
              Kokkos::View<const SCALAR_TYPE**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,      \
              Kokkos::View<SCALAR_TYPE**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,      \
              ETI_SPEC_AVAIL> {                                                                                       \
    using SCALAR    = SCALAR_TYPE;                                                                                    \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;   \
    using BViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;   \
    using CViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;         \
                                                                                                                      \
    static void gemm(const Kokkos::Cuda& space, const char transA[], const char transB[],                             \
                     typename AViewType::const_value_type& alpha, const AViewType& A, const BViewType& B,             \
                     typename CViewType::const_value_type& beta, const CViewType& C) {                                \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemm[TPL_CUBLAS," #SCALAR_TYPE "]");                                 \
      const bool A_t = (transA[0] != 'N') && (transA[0] != 'n');                                                      \
      const int M    = static_cast<int>(C.extent(0));                                                                 \
      const int N    = static_cast<int>(C.extent(1));                                                                 \
      const int K    = static_cast<int>(A.extent(A_t ? 0 : 1));                                                       \
                                                                                                                      \
      bool A_is_lr = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;                                                \
      bool B_is_lr = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;                                                \
      bool C_is_lr = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;                                                \
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
        DotBasedGEMM<Kokkos::Cuda, AViewType, BViewType, CViewType> gemm(alpha, A, B, beta, C);                       \
        bool conjT = (std::is_same<SCALAR, double>::value || std::is_same<SCALAR, float>::value)                      \
                         ? false                                                                                      \
                         : (transa == CUBLAS_OP_C ? true : false);                                                    \
        gemm.run(space, conjT);                                                                                       \
      } else {                                                                                                        \
        KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                    \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));                             \
        if (!A_is_lr && !B_is_lr && !C_is_lr)                                                                         \
          KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(CUBLAS_FN(                                                                 \
              s.handle, transa, transb, M, N, K, reinterpret_cast<const CUDA_SCALAR_TYPE*>(&alpha),                   \
              reinterpret_cast<const CUDA_SCALAR_TYPE*>(A.data()), LDA,                                               \
              reinterpret_cast<const CUDA_SCALAR_TYPE*>(B.data()), LDB,                                               \
              reinterpret_cast<const CUDA_SCALAR_TYPE*>(&beta), reinterpret_cast<CUDA_SCALAR_TYPE*>(C.data()), LDC)); \
        if (A_is_lr && B_is_lr && C_is_lr)                                                                            \
          KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(CUBLAS_FN(                                                                 \
              s.handle, transb, transa, N, M, K, reinterpret_cast<const CUDA_SCALAR_TYPE*>(&alpha),                   \
              reinterpret_cast<const CUDA_SCALAR_TYPE*>(B.data()), LDB,                                               \
              reinterpret_cast<const CUDA_SCALAR_TYPE*>(A.data()), LDA,                                               \
              reinterpret_cast<const CUDA_SCALAR_TYPE*>(&beta), reinterpret_cast<CUDA_SCALAR_TYPE*>(C.data()), LDC)); \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, NULL));                                            \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS3_DGEMM_CUBLAS(LAYOUT) KOKKOSBLAS3_XGEMM_CUBLAS(double, double, cublasDgemm, LAYOUT)

#define KOKKOSBLAS3_SGEMM_CUBLAS(LAYOUT) KOKKOSBLAS3_XGEMM_CUBLAS(float, float, cublasSgemm, LAYOUT)

#define KOKKOSBLAS3_ZGEMM_CUBLAS(LAYOUT) \
  KOKKOSBLAS3_XGEMM_CUBLAS(Kokkos::complex<double>, cuDoubleComplex, cublasZgemm, LAYOUT)

#define KOKKOSBLAS3_CGEMM_CUBLAS(LAYOUT) \
  KOKKOSBLAS3_XGEMM_CUBLAS(Kokkos::complex<float>, cuComplex, cublasCgemm, LAYOUT)

KOKKOSBLAS3_DGEMM_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS3_DGEMM_CUBLAS(Kokkos::LayoutRight)

KOKKOSBLAS3_SGEMM_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS3_SGEMM_CUBLAS(Kokkos::LayoutRight)

KOKKOSBLAS3_ZGEMM_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS3_ZGEMM_CUBLAS(Kokkos::LayoutRight)

KOKKOSBLAS3_CGEMM_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS3_CGEMM_CUBLAS(Kokkos::LayoutRight)

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

// rocBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)
#include <KokkosBlas_tpl_spec.hpp>
#include <KokkosBlas3_gemm_dotbased_impl.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS3_XGEMM_ROCBLAS(SCALAR_TYPE, ROCBLAS_SCALAR_TYPE, ROCBLAS_FN, LAYOUT)                             \
  template <bool ETI_SPEC_AVAIL>                                                                                    \
  struct GEMM<Kokkos::HIP,                                                                                          \
              Kokkos::View<const SCALAR_TYPE**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,     \
              Kokkos::View<const SCALAR_TYPE**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,     \
              Kokkos::View<SCALAR_TYPE**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,     \
              ETI_SPEC_AVAIL> {                                                                                     \
    using SCALAR    = SCALAR_TYPE;                                                                                  \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;  \
    using BViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;  \
    using CViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;        \
                                                                                                                    \
    static void gemm(const typename CViewType::execution_space& space, const char transA[], const char transB[],    \
                     typename AViewType::const_value_type& alpha, const AViewType& A, const BViewType& B,           \
                     typename CViewType::const_value_type& beta, const CViewType& C) {                              \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemm[TPL_ROCBLAS," #SCALAR_TYPE "]");                              \
                                                                                                                    \
      const bool A_t = (transA[0] != 'N') && (transA[0] != 'n');                                                    \
      const int M    = static_cast<int>(C.extent(0));                                                               \
      const int N    = static_cast<int>(C.extent(1));                                                               \
      const int K    = static_cast<int>(A.extent(A_t ? 0 : 1));                                                     \
                                                                                                                    \
      bool is_lr = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;                                                \
                                                                                                                    \
      const int AST = is_lr ? A.stride(0) : A.stride(1), LDA = AST == 0 ? 1 : AST;                                  \
      const int BST = is_lr ? B.stride(0) : B.stride(1), LDB = BST == 0 ? 1 : BST;                                  \
      const int CST = is_lr ? C.stride(0) : C.stride(1), LDC = CST == 0 ? 1 : CST;                                  \
                                                                                                                    \
      rocblas_operation transa = trans_mode_kk_to_rocblas(transA);                                                  \
      rocblas_operation transb = trans_mode_kk_to_rocblas(transB);                                                  \
                                                                                                                    \
      constexpr int numDotsLayoutLeftThreshold  = 1600;                                                             \
      constexpr int numDotsLayoutRightThreshold = 100;                                                              \
      if ((!is_lr && transa != rocblas_operation_none && transb == rocblas_operation_none &&                        \
           M * N < numDotsLayoutLeftThreshold) ||                                                                   \
          (is_lr && transa != rocblas_operation_none && transb == rocblas_operation_none &&                         \
           M * N < numDotsLayoutRightThreshold)) {                                                                  \
        DotBasedGEMM<Kokkos::HIP, AViewType, BViewType, CViewType> gemm(alpha, A, B, beta, C);                      \
        bool conjT = (std::is_same<SCALAR, double>::value || std::is_same<SCALAR, float>::value)                    \
                         ? false                                                                                    \
                         : (transa == rocblas_operation_conjugate_transpose ? true : false);                        \
        gemm.run(space, conjT);                                                                                     \
      } else {                                                                                                      \
        KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                    \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                        \
        if (!is_lr)                                                                                                 \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(ROCBLAS_FN(s.handle, transa, transb, M, N, K,                           \
                                                       reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(&alpha),        \
                                                       reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(A.data()), LDA, \
                                                       reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(B.data()), LDB, \
                                                       reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(&beta),         \
                                                       reinterpret_cast<ROCBLAS_SCALAR_TYPE*>(C.data()), LDC));     \
        else                                                                                                        \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(ROCBLAS_FN(s.handle, transb, transa, N, M, K,                           \
                                                       reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(&alpha),        \
                                                       reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(B.data()), LDB, \
                                                       reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(A.data()), LDA, \
                                                       reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(&beta),         \
                                                       reinterpret_cast<ROCBLAS_SCALAR_TYPE*>(C.data()), LDC));     \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                      \
      }                                                                                                             \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS3_DGEMM_ROCBLAS(LAYOUT) KOKKOSBLAS3_XGEMM_ROCBLAS(double, double, rocblas_dgemm, LAYOUT)

#define KOKKOSBLAS3_SGEMM_ROCBLAS(LAYOUT) KOKKOSBLAS3_XGEMM_ROCBLAS(float, float, rocblas_sgemm, LAYOUT)

#define KOKKOSBLAS3_ZGEMM_ROCBLAS(LAYOUT) \
  KOKKOSBLAS3_XGEMM_ROCBLAS(Kokkos::complex<double>, rocblas_double_complex, rocblas_zgemm, LAYOUT)

#define KOKKOSBLAS3_CGEMM_ROCBLAS(LAYOUT) \
  KOKKOSBLAS3_XGEMM_ROCBLAS(Kokkos::complex<float>, rocblas_float_complex, rocblas_cgemm, LAYOUT)

KOKKOSBLAS3_DGEMM_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS3_DGEMM_ROCBLAS(Kokkos::LayoutRight)

KOKKOSBLAS3_SGEMM_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS3_SGEMM_ROCBLAS(Kokkos::LayoutRight)

KOKKOSBLAS3_ZGEMM_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS3_ZGEMM_ROCBLAS(Kokkos::LayoutRight)

KOKKOSBLAS3_CGEMM_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS3_CGEMM_ROCBLAS(Kokkos::LayoutRight)

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCBLAS

#endif
