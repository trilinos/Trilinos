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

#ifndef KOKKOSBLAS3_TRMM_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS3_TRMM_TPL_SPEC_DECL_HPP_

// Generic Host side BLAS (could be MKL or anything)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS3_TRMM_BLAS(SCALAR_TYPE, BASE_SCALAR_TYPE, LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)  \
  template <class ExecSpace>                                                                               \
  struct TRMM<ExecSpace,                                                                                   \
              Kokkos::View<const SCALAR_TYPE**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,             \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                      \
              Kokkos::View<SCALAR_TYPE**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                   \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                      \
              true, ETI_SPEC_AVAIL> {                                                                      \
    typedef SCALAR_TYPE SCALAR;                                                                            \
    typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                    \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                         \
        AViewType;                                                                                         \
    typedef Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                          \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                         \
        BViewType;                                                                                         \
                                                                                                           \
    static void trmm(const ExecSpace& /*space*/, const char side[], const char uplo[], const char trans[], \
                     const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A,   \
                     const BViewType& B) {                                                                 \
      Kokkos::Profiling::pushRegion("KokkosBlas::trmm[TPL_BLAS," #SCALAR_TYPE "]");                        \
      const int M = static_cast<int>(B.extent(0));                                                         \
      const int N = static_cast<int>(B.extent(1));                                                         \
                                                                                                           \
      bool A_is_layout_left = std::is_same<Kokkos::LayoutLeft, LAYOUTA>::value;                            \
      bool B_is_layout_left = std::is_same<Kokkos::LayoutLeft, LAYOUTB>::value;                            \
                                                                                                           \
      const int AST = A_is_layout_left ? A.stride(1) : A.stride(0), LDA = (AST == 0) ? 1 : AST;            \
      const int BST = B_is_layout_left ? B.stride(1) : B.stride(0), LDB = (BST == 0) ? 1 : BST;            \
                                                                                                           \
      char side_;                                                                                          \
      char uplo_;                                                                                          \
                                                                                                           \
      if (A_is_layout_left) {                                                                              \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                          \
          side_ = 'L';                                                                                     \
        else                                                                                               \
          side_ = 'R';                                                                                     \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                          \
          uplo_ = 'L';                                                                                     \
        else                                                                                               \
          uplo_ = 'U';                                                                                     \
      } else {                                                                                             \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                          \
          side_ = 'R';                                                                                     \
        else                                                                                               \
          side_ = 'L';                                                                                     \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                          \
          uplo_ = 'U';                                                                                     \
        else                                                                                               \
          uplo_ = 'L';                                                                                     \
      }                                                                                                    \
                                                                                                           \
      if (A_is_layout_left)                                                                                \
        HostBlas<BASE_SCALAR_TYPE>::trmm(side_, uplo_, trans[0], diag[0], M, N, alpha,                     \
                                         reinterpret_cast<const BASE_SCALAR_TYPE*>(A.data()), LDA,         \
                                         reinterpret_cast<BASE_SCALAR_TYPE*>(B.data()), LDB);              \
      else                                                                                                 \
        HostBlas<BASE_SCALAR_TYPE>::trmm(side_, uplo_, trans[0], diag[0], N, M, alpha,                     \
                                         reinterpret_cast<const BASE_SCALAR_TYPE*>(A.data()), LDA,         \
                                         reinterpret_cast<BASE_SCALAR_TYPE*>(B.data()), LDB);              \
      Kokkos::Profiling::popRegion();                                                                      \
    }                                                                                                      \
  };

#define KOKKOSBLAS3_DTRMM_BLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS3_TRMM_BLAS(double, double, LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS3_STRMM_BLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS3_TRMM_BLAS(float, float, LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS3_ZTRMM_BLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS3_TRMM_BLAS(Kokkos::complex<double>, std::complex<double>, LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS3_CTRMM_BLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS3_TRMM_BLAS(Kokkos::complex<float>, std::complex<float>, LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)

// Explicitly define the TRMM class for all permutations listed below

KOKKOSBLAS3_DTRMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_DTRMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_DTRMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_DTRMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

KOKKOSBLAS3_STRMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_STRMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_STRMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_STRMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

KOKKOSBLAS3_ZTRMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_ZTRMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_ZTRMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_ZTRMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

KOKKOSBLAS3_CTRMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_CTRMM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_CTRMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_CTRMM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS3_TRMM_CUBLAS(SCALAR_TYPE, CUDA_SCALAR_TYPE, CUBLAS_FN, LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL) \
  template <class ExecSpace>                                                                                           \
  struct TRMM<ExecSpace,                                                                                               \
              Kokkos::View<const SCALAR_TYPE**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                         \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
              Kokkos::View<SCALAR_TYPE**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                               \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
              true, ETI_SPEC_AVAIL> {                                                                                  \
    typedef SCALAR_TYPE SCALAR;                                                                                        \
    typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        AViewType;                                                                                                     \
    typedef Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                                      \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        BViewType;                                                                                                     \
                                                                                                                       \
    static void trmm(const ExecSpace& space, const char side[], const char uplo[], const char trans[],                 \
                     const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A,               \
                     const BViewType& B) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::trmm[TPL_CUBLAS," #SCALAR_TYPE "]");                                  \
      const int M = static_cast<int>(B.extent(0));                                                                     \
      const int N = static_cast<int>(B.extent(1));                                                                     \
                                                                                                                       \
      bool A_is_layout_left = std::is_same<Kokkos::LayoutLeft, LAYOUTA>::value;                                        \
      bool B_is_layout_left = std::is_same<Kokkos::LayoutLeft, LAYOUTB>::value;                                        \
                                                                                                                       \
      const int AST = A_is_layout_left ? A.stride(1) : A.stride(0), LDA = (AST == 0) ? 1 : AST;                        \
      const int BST = B_is_layout_left ? B.stride(1) : B.stride(0), LDB = (BST == 0) ? 1 : BST;                        \
                                                                                                                       \
      cublasSideMode_t side_;                                                                                          \
      cublasFillMode_t uplo_;                                                                                          \
      cublasOperation_t trans_;                                                                                        \
      cublasDiagType_t diag_;                                                                                          \
                                                                                                                       \
      if (A_is_layout_left) {                                                                                          \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                                      \
          side_ = CUBLAS_SIDE_LEFT;                                                                                    \
        else                                                                                                           \
          side_ = CUBLAS_SIDE_RIGHT;                                                                                   \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                                      \
          uplo_ = CUBLAS_FILL_MODE_LOWER;                                                                              \
        else                                                                                                           \
          uplo_ = CUBLAS_FILL_MODE_UPPER;                                                                              \
      } else {                                                                                                         \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                                      \
          side_ = CUBLAS_SIDE_RIGHT;                                                                                   \
        else                                                                                                           \
          side_ = CUBLAS_SIDE_LEFT;                                                                                    \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                                      \
          uplo_ = CUBLAS_FILL_MODE_UPPER;                                                                              \
        else                                                                                                           \
          uplo_ = CUBLAS_FILL_MODE_LOWER;                                                                              \
      }                                                                                                                \
                                                                                                                       \
      if ((trans[0] == 'N') || (trans[0] == 'n'))                                                                      \
        trans_ = CUBLAS_OP_N;                                                                                          \
      else if ((trans[0] == 'T') || (trans[0] == 't'))                                                                 \
        trans_ = CUBLAS_OP_T;                                                                                          \
      else                                                                                                             \
        trans_ = CUBLAS_OP_C;                                                                                          \
      if ((diag[0] == 'U') || (diag[0] == 'u'))                                                                        \
        diag_ = CUBLAS_DIAG_UNIT;                                                                                      \
      else                                                                                                             \
        diag_ = CUBLAS_DIAG_NON_UNIT;                                                                                  \
                                                                                                                       \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                       \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                    \
      if (A_is_layout_left) {                                                                                          \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(CUBLAS_FN(                                                                        \
            s.handle, side_, uplo_, trans_, diag_, M, N, reinterpret_cast<const CUDA_SCALAR_TYPE*>(&alpha),            \
            reinterpret_cast<const CUDA_SCALAR_TYPE*>(A.data()), LDA, reinterpret_cast<CUDA_SCALAR_TYPE*>(B.data()),   \
            LDB, reinterpret_cast<CUDA_SCALAR_TYPE*>(B.data()), LDB));                                                 \
      } else {                                                                                                         \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(CUBLAS_FN(                                                                        \
            s.handle, side_, uplo_, trans_, diag_, N, M, reinterpret_cast<const CUDA_SCALAR_TYPE*>(&alpha),            \
            reinterpret_cast<const CUDA_SCALAR_TYPE*>(A.data()), LDA, reinterpret_cast<CUDA_SCALAR_TYPE*>(B.data()),   \
            LDB, reinterpret_cast<CUDA_SCALAR_TYPE*>(B.data()), LDB));                                                 \
      }                                                                                                                \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                                   \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS3_DTRMM_CUBLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS3_TRMM_CUBLAS(double, double, cublasDtrmm, LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS3_STRMM_CUBLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS3_TRMM_CUBLAS(float, float, cublasStrmm, LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS3_ZTRMM_CUBLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)                                 \
  KOKKOSBLAS3_TRMM_CUBLAS(Kokkos::complex<double>, cuDoubleComplex, cublasZtrmm, LAYOUTA, LAYOUTB, MEM_SPACE, \
                          ETI_SPEC_AVAIL)

#define KOKKOSBLAS3_CTRMM_CUBLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS3_TRMM_CUBLAS(Kokkos::complex<float>, cuComplex, cublasCtrmm, LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)

// Explicitly define the TRMM class for all permutations listed below

KOKKOSBLAS3_DTRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_DTRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_DTRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_DTRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_DTRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_DTRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_DTRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_DTRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS3_STRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_STRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_STRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_STRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_STRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_STRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_STRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_STRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS3_ZTRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_ZTRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_ZTRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_ZTRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_ZTRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_ZTRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_ZTRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_ZTRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS3_CTRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_CTRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_CTRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_CTRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_CTRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_CTRMM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_CTRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_CTRMM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

#endif  // KOKKOSBLAS3_TRMM_TPL_SPEC_DECL_HPP_
