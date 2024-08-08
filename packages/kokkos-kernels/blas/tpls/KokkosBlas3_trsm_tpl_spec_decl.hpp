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

#ifndef KOKKOSBLAS3_TRSM_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS3_TRSM_TPL_SPEC_DECL_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS3_DTRSM_BLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)                                            \
  template <class ExecSpace>                                                                                           \
  struct TRSM<                                                                                                         \
      ExecSpace,                                                                                                       \
      Kokkos::View<const double**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<double**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    typedef double SCALAR;                                                                                             \
    typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        AViewType;                                                                                                     \
    typedef Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                                      \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        BViewType;                                                                                                     \
                                                                                                                       \
    static void trsm(const ExecSpace& /*space*/, const char side[], const char uplo[], const char trans[],             \
                     const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A,               \
                     const BViewType& B) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::trsm[TPL_BLAS,double]");                                              \
      const int M = static_cast<int>(B.extent(0));                                                                     \
      const int N = static_cast<int>(B.extent(1));                                                                     \
                                                                                                                       \
      bool A_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTA>::value;                                                 \
      bool B_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTB>::value;                                                 \
                                                                                                                       \
      const int AST = A_is_ll ? A.stride(1) : A.stride(0), LDA = (AST == 0) ? 1 : AST;                                 \
      const int BST = B_is_ll ? B.stride(1) : B.stride(0), LDB = (BST == 0) ? 1 : BST;                                 \
                                                                                                                       \
      char side_;                                                                                                      \
      char uplo_;                                                                                                      \
                                                                                                                       \
      if (A_is_ll) {                                                                                                   \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                                      \
          side_ = 'L';                                                                                                 \
        else                                                                                                           \
          side_ = 'R';                                                                                                 \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                                      \
          uplo_ = 'L';                                                                                                 \
        else                                                                                                           \
          uplo_ = 'U';                                                                                                 \
      } else {                                                                                                         \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                                      \
          side_ = 'R';                                                                                                 \
        else                                                                                                           \
          side_ = 'L';                                                                                                 \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                                      \
          uplo_ = 'U';                                                                                                 \
        else                                                                                                           \
          uplo_ = 'L';                                                                                                 \
      }                                                                                                                \
                                                                                                                       \
      if (A_is_ll)                                                                                                     \
        HostBlas<double>::trsm(side_, uplo_, trans[0], diag[0], M, N, alpha, A.data(), LDA, B.data(), LDB);            \
      else                                                                                                             \
        HostBlas<double>::trsm(side_, uplo_, trans[0], diag[0], N, M, alpha, A.data(), LDA, B.data(), LDB);            \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS3_STRSM_BLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)                                           \
  template <class ExecSpace>                                                                                          \
  struct TRSM<                                                                                                        \
      ExecSpace,                                                                                                      \
      Kokkos::View<const float**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<float**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    typedef float SCALAR;                                                                                             \
    typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        AViewType;                                                                                                    \
    typedef Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        BViewType;                                                                                                    \
                                                                                                                      \
    static void trsm(const ExecSpace& /*space*/, const char side[], const char uplo[], const char trans[],            \
                     const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A,              \
                     const BViewType& B) {                                                                            \
      Kokkos::Profiling::pushRegion("KokkosBlas::trsm[TPL_BLAS,float]");                                              \
      const int M = static_cast<int>(B.extent(0));                                                                    \
      const int N = static_cast<int>(B.extent(1));                                                                    \
                                                                                                                      \
      bool A_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTA>::value;                                                \
      bool B_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTB>::value;                                                \
                                                                                                                      \
      const int AST = A_is_ll ? A.stride(1) : A.stride(0), LDA = (AST == 0) ? 1 : AST;                                \
      const int BST = B_is_ll ? B.stride(1) : B.stride(0), LDB = (BST == 0) ? 1 : BST;                                \
                                                                                                                      \
      char side_;                                                                                                     \
      char uplo_;                                                                                                     \
                                                                                                                      \
      if (A_is_ll) {                                                                                                  \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                                     \
          side_ = 'L';                                                                                                \
        else                                                                                                          \
          side_ = 'R';                                                                                                \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                                     \
          uplo_ = 'L';                                                                                                \
        else                                                                                                          \
          uplo_ = 'U';                                                                                                \
      } else {                                                                                                        \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                                     \
          side_ = 'R';                                                                                                \
        else                                                                                                          \
          side_ = 'L';                                                                                                \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                                     \
          uplo_ = 'U';                                                                                                \
        else                                                                                                          \
          uplo_ = 'L';                                                                                                \
      }                                                                                                               \
                                                                                                                      \
      if (A_is_ll)                                                                                                    \
        HostBlas<float>::trsm(side_, uplo_, trans[0], diag[0], M, N, alpha, A.data(), LDA, B.data(), LDB);            \
      else                                                                                                            \
        HostBlas<float>::trsm(side_, uplo_, trans[0], diag[0], N, M, alpha, A.data(), LDA, B.data(), LDB);            \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS3_ZTRSM_BLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)                                 \
  template <class ExecSpace>                                                                                \
  struct TRSM<ExecSpace,                                                                                    \
              Kokkos::View<const Kokkos::complex<double>**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,  \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                       \
              Kokkos::View<Kokkos::complex<double>**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,        \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                       \
              true, ETI_SPEC_AVAIL> {                                                                       \
    typedef Kokkos::complex<double> SCALAR;                                                                 \
    typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                          \
        AViewType;                                                                                          \
    typedef Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                           \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                          \
        BViewType;                                                                                          \
                                                                                                            \
    static void trsm(const ExecSpace& /*space*/, const char side[], const char uplo[], const char trans[],  \
                     const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A,    \
                     const BViewType& B) {                                                                  \
      Kokkos::Profiling::pushRegion("KokkosBlas::trsm[TPL_BLAS,complex<double>]");                          \
      const int M = static_cast<int>(B.extent(0));                                                          \
      const int N = static_cast<int>(B.extent(1));                                                          \
                                                                                                            \
      bool A_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTA>::value;                                      \
      bool B_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTB>::value;                                      \
                                                                                                            \
      const int AST = A_is_ll ? A.stride(1) : A.stride(0), LDA = (AST == 0) ? 1 : AST;                      \
      const int BST = B_is_ll ? B.stride(1) : B.stride(0), LDB = (BST == 0) ? 1 : BST;                      \
                                                                                                            \
      char side_;                                                                                           \
      char uplo_;                                                                                           \
                                                                                                            \
      if (A_is_ll) {                                                                                        \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                           \
          side_ = 'L';                                                                                      \
        else                                                                                                \
          side_ = 'R';                                                                                      \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                           \
          uplo_ = 'L';                                                                                      \
        else                                                                                                \
          uplo_ = 'U';                                                                                      \
      } else {                                                                                              \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                           \
          side_ = 'R';                                                                                      \
        else                                                                                                \
          side_ = 'L';                                                                                      \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                           \
          uplo_ = 'U';                                                                                      \
        else                                                                                                \
          uplo_ = 'L';                                                                                      \
      }                                                                                                     \
                                                                                                            \
      const std::complex<double> alpha_val = alpha;                                                         \
      if (A_is_ll)                                                                                          \
        HostBlas<std::complex<double> >::trsm(side_, uplo_, trans[0], diag[0], M, N, alpha_val,             \
                                              reinterpret_cast<const std::complex<double>*>(A.data()), LDA, \
                                              reinterpret_cast<std::complex<double>*>(B.data()), LDB);      \
      else                                                                                                  \
        HostBlas<std::complex<double> >::trsm(side_, uplo_, trans[0], diag[0], N, M, alpha_val,             \
                                              reinterpret_cast<const std::complex<double>*>(A.data()), LDA, \
                                              reinterpret_cast<std::complex<double>*>(B.data()), LDB);      \
      Kokkos::Profiling::popRegion();                                                                       \
    }                                                                                                       \
  };

#define KOKKOSBLAS3_CTRSM_BLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)                                \
  template <class ExecSpace>                                                                               \
  struct TRSM<ExecSpace,                                                                                   \
              Kokkos::View<const Kokkos::complex<float>**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,  \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                      \
              Kokkos::View<Kokkos::complex<float>**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,        \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                      \
              true, ETI_SPEC_AVAIL> {                                                                      \
    typedef Kokkos::complex<float> SCALAR;                                                                 \
    typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                    \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                         \
        AViewType;                                                                                         \
    typedef Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                          \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                         \
        BViewType;                                                                                         \
                                                                                                           \
    static void trsm(const ExecSpace& /*space*/, const char side[], const char uplo[], const char trans[], \
                     const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A,   \
                     const BViewType& B) {                                                                 \
      Kokkos::Profiling::pushRegion("KokkosBlas::trsm[TPL_BLAS,complex<float>]");                          \
      const int M = static_cast<int>(B.extent(0));                                                         \
      const int N = static_cast<int>(B.extent(1));                                                         \
                                                                                                           \
      bool A_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTA>::value;                                     \
      bool B_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTB>::value;                                     \
                                                                                                           \
      const int AST = A_is_ll ? A.stride(1) : A.stride(0), LDA = (AST == 0) ? 1 : AST;                     \
      const int BST = B_is_ll ? B.stride(1) : B.stride(0), LDB = (BST == 0) ? 1 : BST;                     \
                                                                                                           \
      char side_;                                                                                          \
      char uplo_;                                                                                          \
                                                                                                           \
      if (A_is_ll) {                                                                                       \
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
      const std::complex<float> alpha_val = alpha;                                                         \
      if (A_is_ll)                                                                                         \
        HostBlas<std::complex<float> >::trsm(side_, uplo_, trans[0], diag[0], M, N, alpha_val,             \
                                             reinterpret_cast<const std::complex<float>*>(A.data()), LDA,  \
                                             reinterpret_cast<std::complex<float>*>(B.data()), LDB);       \
      else                                                                                                 \
        HostBlas<std::complex<float> >::trsm(side_, uplo_, trans[0], diag[0], N, M, alpha_val,             \
                                             reinterpret_cast<const std::complex<float>*>(A.data()), LDA,  \
                                             reinterpret_cast<std::complex<float>*>(B.data()), LDB);       \
      Kokkos::Profiling::popRegion();                                                                      \
    }                                                                                                      \
  };

KOKKOSBLAS3_DTRSM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_DTRSM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_DTRSM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_DTRSM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

KOKKOSBLAS3_STRSM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_STRSM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_STRSM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_STRSM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

KOKKOSBLAS3_ZTRSM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_ZTRSM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_ZTRSM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_ZTRSM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

KOKKOSBLAS3_CTRSM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_CTRSM_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_CTRSM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_CTRSM_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS3_DTRSM_CUBLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)                                          \
  template <class ExecSpace>                                                                                           \
  struct TRSM<                                                                                                         \
      ExecSpace,                                                                                                       \
      Kokkos::View<const double**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<double**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    typedef double SCALAR;                                                                                             \
    typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        AViewType;                                                                                                     \
    typedef Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                                      \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        BViewType;                                                                                                     \
                                                                                                                       \
    static void trsm(const ExecSpace& space, const char side[], const char uplo[], const char trans[],                 \
                     const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A,               \
                     const BViewType& B) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::trsm[TPL_CUBLAS,double]");                                            \
      const int M = static_cast<int>(B.extent(0));                                                                     \
      const int N = static_cast<int>(B.extent(1));                                                                     \
                                                                                                                       \
      bool A_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTA>::value;                                                 \
      bool B_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTB>::value;                                                 \
                                                                                                                       \
      const int AST = A_is_ll ? A.stride(1) : A.stride(0), LDA = (AST == 0) ? 1 : AST;                                 \
      const int BST = B_is_ll ? B.stride(1) : B.stride(0), LDB = (BST == 0) ? 1 : BST;                                 \
                                                                                                                       \
      cublasSideMode_t side_;                                                                                          \
      cublasFillMode_t uplo_;                                                                                          \
      cublasOperation_t trans_;                                                                                        \
      cublasDiagType_t diag_;                                                                                          \
                                                                                                                       \
      if (A_is_ll) {                                                                                                   \
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
      if (A_is_ll) {                                                                                                   \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(                                                                                  \
            cublasDtrsm(s.handle, side_, uplo_, trans_, diag_, M, N, &alpha, A.data(), LDA, B.data(), LDB));           \
      } else {                                                                                                         \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(                                                                                  \
            cublasDtrsm(s.handle, side_, uplo_, trans_, diag_, N, M, &alpha, A.data(), LDA, B.data(), LDB));           \
      }                                                                                                                \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                                   \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS3_STRSM_CUBLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)                                         \
  template <class ExecSpace>                                                                                          \
  struct TRSM<                                                                                                        \
      ExecSpace,                                                                                                      \
      Kokkos::View<const float**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<float**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    typedef float SCALAR;                                                                                             \
    typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        AViewType;                                                                                                    \
    typedef Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        BViewType;                                                                                                    \
                                                                                                                      \
    static void trsm(const ExecSpace& space, const char side[], const char uplo[], const char trans[],                \
                     const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A,              \
                     const BViewType& B) {                                                                            \
      Kokkos::Profiling::pushRegion("KokkosBlas::trsm[TPL_CUBLAS,float]");                                            \
      const int M = static_cast<int>(B.extent(0));                                                                    \
      const int N = static_cast<int>(B.extent(1));                                                                    \
                                                                                                                      \
      bool A_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTA>::value;                                                \
      bool B_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTB>::value;                                                \
                                                                                                                      \
      const int AST = A_is_ll ? A.stride(1) : A.stride(0), LDA = (AST == 0) ? 1 : AST;                                \
      const int BST = B_is_ll ? B.stride(1) : B.stride(0), LDB = (BST == 0) ? 1 : BST;                                \
                                                                                                                      \
      cublasSideMode_t side_;                                                                                         \
      cublasFillMode_t uplo_;                                                                                         \
      cublasOperation_t trans_;                                                                                       \
      cublasDiagType_t diag_;                                                                                         \
                                                                                                                      \
      if (A_is_ll) {                                                                                                  \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                                     \
          side_ = CUBLAS_SIDE_LEFT;                                                                                   \
        else                                                                                                          \
          side_ = CUBLAS_SIDE_RIGHT;                                                                                  \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                                     \
          uplo_ = CUBLAS_FILL_MODE_LOWER;                                                                             \
        else                                                                                                          \
          uplo_ = CUBLAS_FILL_MODE_UPPER;                                                                             \
      } else {                                                                                                        \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                                     \
          side_ = CUBLAS_SIDE_RIGHT;                                                                                  \
        else                                                                                                          \
          side_ = CUBLAS_SIDE_LEFT;                                                                                   \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                                     \
          uplo_ = CUBLAS_FILL_MODE_UPPER;                                                                             \
        else                                                                                                          \
          uplo_ = CUBLAS_FILL_MODE_LOWER;                                                                             \
      }                                                                                                               \
                                                                                                                      \
      if ((trans[0] == 'N') || (trans[0] == 'n'))                                                                     \
        trans_ = CUBLAS_OP_N;                                                                                         \
      else if ((trans[0] == 'T') || (trans[0] == 't'))                                                                \
        trans_ = CUBLAS_OP_T;                                                                                         \
      else                                                                                                            \
        trans_ = CUBLAS_OP_C;                                                                                         \
      if ((diag[0] == 'U') || (diag[0] == 'u'))                                                                       \
        diag_ = CUBLAS_DIAG_UNIT;                                                                                     \
      else                                                                                                            \
        diag_ = CUBLAS_DIAG_NON_UNIT;                                                                                 \
                                                                                                                      \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                      \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                   \
      if (A_is_ll) {                                                                                                  \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(                                                                                 \
            cublasStrsm(s.handle, side_, uplo_, trans_, diag_, M, N, &alpha, A.data(), LDA, B.data(), LDB));          \
      } else {                                                                                                        \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(                                                                                 \
            cublasStrsm(s.handle, side_, uplo_, trans_, diag_, N, M, &alpha, A.data(), LDA, B.data(), LDB));          \
      }                                                                                                               \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                                  \
                                                                                                                      \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS3_ZTRSM_CUBLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)                              \
  template <class ExecSpace>                                                                               \
  struct TRSM<ExecSpace,                                                                                   \
              Kokkos::View<const Kokkos::complex<double>**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                      \
              Kokkos::View<Kokkos::complex<double>**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,       \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                      \
              true, ETI_SPEC_AVAIL> {                                                                      \
    typedef Kokkos::complex<double> SCALAR;                                                                \
    typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                    \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                         \
        AViewType;                                                                                         \
    typedef Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                          \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                         \
        BViewType;                                                                                         \
                                                                                                           \
    static void trsm(const ExecSpace& space, const char side[], const char uplo[], const char trans[],     \
                     const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A,   \
                     const BViewType& B) {                                                                 \
      Kokkos::Profiling::pushRegion("KokkosBlas::trsm[TPL_CUBLAS,complex<double>]");                       \
      const int M = static_cast<int>(B.extent(0));                                                         \
      const int N = static_cast<int>(B.extent(1));                                                         \
                                                                                                           \
      bool A_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTA>::value;                                     \
      bool B_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTB>::value;                                     \
                                                                                                           \
      const int AST = A_is_ll ? A.stride(1) : A.stride(0), LDA = (AST == 0) ? 1 : AST;                     \
      const int BST = B_is_ll ? B.stride(1) : B.stride(0), LDB = (BST == 0) ? 1 : BST;                     \
                                                                                                           \
      cublasSideMode_t side_;                                                                              \
      cublasFillMode_t uplo_;                                                                              \
      cublasOperation_t trans_;                                                                            \
      cublasDiagType_t diag_;                                                                              \
                                                                                                           \
      if (A_is_ll) {                                                                                       \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                          \
          side_ = CUBLAS_SIDE_LEFT;                                                                        \
        else                                                                                               \
          side_ = CUBLAS_SIDE_RIGHT;                                                                       \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                          \
          uplo_ = CUBLAS_FILL_MODE_LOWER;                                                                  \
        else                                                                                               \
          uplo_ = CUBLAS_FILL_MODE_UPPER;                                                                  \
      } else {                                                                                             \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                          \
          side_ = CUBLAS_SIDE_RIGHT;                                                                       \
        else                                                                                               \
          side_ = CUBLAS_SIDE_LEFT;                                                                        \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                          \
          uplo_ = CUBLAS_FILL_MODE_UPPER;                                                                  \
        else                                                                                               \
          uplo_ = CUBLAS_FILL_MODE_LOWER;                                                                  \
      }                                                                                                    \
                                                                                                           \
      if ((trans[0] == 'N') || (trans[0] == 'n'))                                                          \
        trans_ = CUBLAS_OP_N;                                                                              \
      else if ((trans[0] == 'T') || (trans[0] == 't'))                                                     \
        trans_ = CUBLAS_OP_T;                                                                              \
      else                                                                                                 \
        trans_ = CUBLAS_OP_C;                                                                              \
      if ((diag[0] == 'U') || (diag[0] == 'u'))                                                            \
        diag_ = CUBLAS_DIAG_UNIT;                                                                          \
      else                                                                                                 \
        diag_ = CUBLAS_DIAG_NON_UNIT;                                                                      \
                                                                                                           \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();           \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                        \
      if (A_is_ll) {                                                                                       \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasZtrsm(s.handle, side_, uplo_, trans_, diag_, M, N,              \
                                                 reinterpret_cast<const cuDoubleComplex*>(&alpha),         \
                                                 reinterpret_cast<const cuDoubleComplex*>(A.data()), LDA,  \
                                                 reinterpret_cast<cuDoubleComplex*>(B.data()), LDB));      \
      } else {                                                                                             \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasZtrsm(s.handle, side_, uplo_, trans_, diag_, N, M,              \
                                                 reinterpret_cast<const cuDoubleComplex*>(&alpha),         \
                                                 reinterpret_cast<const cuDoubleComplex*>(A.data()), LDA,  \
                                                 reinterpret_cast<cuDoubleComplex*>(B.data()), LDB));      \
      }                                                                                                    \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                       \
                                                                                                           \
      Kokkos::Profiling::popRegion();                                                                      \
    }                                                                                                      \
  };

#define KOKKOSBLAS3_CTRSM_CUBLAS(LAYOUTA, LAYOUTB, MEM_SPACE, ETI_SPEC_AVAIL)                                 \
  template <class ExecSpace>                                                                                  \
  struct TRSM<ExecSpace,                                                                                      \
              Kokkos::View<const Kokkos::complex<float>**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,     \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                         \
              Kokkos::View<Kokkos::complex<float>**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,           \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                         \
              true, ETI_SPEC_AVAIL> {                                                                         \
    typedef Kokkos::complex<float> SCALAR;                                                                    \
    typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                       \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                            \
        AViewType;                                                                                            \
    typedef Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>,                             \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                            \
        BViewType;                                                                                            \
                                                                                                              \
    static void trsm(const ExecSpace& space, const char side[], const char uplo[], const char trans[],        \
                     const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A,      \
                     const BViewType& B) {                                                                    \
      Kokkos::Profiling::pushRegion("KokkosBlas::trsm[TPL_CUBLAS,complex<float>]");                           \
      const int M = static_cast<int>(B.extent(0));                                                            \
      const int N = static_cast<int>(B.extent(1));                                                            \
                                                                                                              \
      bool A_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTA>::value;                                        \
      bool B_is_ll = std::is_same<Kokkos::LayoutLeft, LAYOUTB>::value;                                        \
                                                                                                              \
      const int AST = A_is_ll ? A.stride(1) : A.stride(0), LDA = (AST == 0) ? 1 : AST;                        \
      const int BST = B_is_ll ? B.stride(1) : B.stride(0), LDB = (BST == 0) ? 1 : BST;                        \
                                                                                                              \
      cublasSideMode_t side_;                                                                                 \
      cublasFillMode_t uplo_;                                                                                 \
      cublasOperation_t trans_;                                                                               \
      cublasDiagType_t diag_;                                                                                 \
                                                                                                              \
      if (A_is_ll) {                                                                                          \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                             \
          side_ = CUBLAS_SIDE_LEFT;                                                                           \
        else                                                                                                  \
          side_ = CUBLAS_SIDE_RIGHT;                                                                          \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                             \
          uplo_ = CUBLAS_FILL_MODE_LOWER;                                                                     \
        else                                                                                                  \
          uplo_ = CUBLAS_FILL_MODE_UPPER;                                                                     \
      } else {                                                                                                \
        if ((side[0] == 'L') || (side[0] == 'l'))                                                             \
          side_ = CUBLAS_SIDE_RIGHT;                                                                          \
        else                                                                                                  \
          side_ = CUBLAS_SIDE_LEFT;                                                                           \
        if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                             \
          uplo_ = CUBLAS_FILL_MODE_UPPER;                                                                     \
        else                                                                                                  \
          uplo_ = CUBLAS_FILL_MODE_LOWER;                                                                     \
      }                                                                                                       \
                                                                                                              \
      if ((trans[0] == 'N') || (trans[0] == 'n'))                                                             \
        trans_ = CUBLAS_OP_N;                                                                                 \
      else if ((trans[0] == 'T') || (trans[0] == 't'))                                                        \
        trans_ = CUBLAS_OP_T;                                                                                 \
      else                                                                                                    \
        trans_ = CUBLAS_OP_C;                                                                                 \
      if ((diag[0] == 'U') || (diag[0] == 'u'))                                                               \
        diag_ = CUBLAS_DIAG_UNIT;                                                                             \
      else                                                                                                    \
        diag_ = CUBLAS_DIAG_NON_UNIT;                                                                         \
                                                                                                              \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();              \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                           \
      if (A_is_ll) {                                                                                          \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasCtrsm(                                                             \
            s.handle, side_, uplo_, trans_, diag_, M, N, reinterpret_cast<const cuComplex*>(&alpha),          \
            reinterpret_cast<const cuComplex*>(A.data()), LDA, reinterpret_cast<cuComplex*>(B.data()), LDB)); \
      } else {                                                                                                \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasCtrsm(                                                             \
            s.handle, side_, uplo_, trans_, diag_, N, M, reinterpret_cast<const cuComplex*>(&alpha),          \
            reinterpret_cast<const cuComplex*>(A.data()), LDA, reinterpret_cast<cuComplex*>(B.data()), LDB)); \
      }                                                                                                       \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                          \
                                                                                                              \
      Kokkos::Profiling::popRegion();                                                                         \
    }                                                                                                         \
  };

KOKKOSBLAS3_DTRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_DTRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_DTRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_DTRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_DTRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_DTRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_DTRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_DTRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS3_STRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_STRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_STRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_STRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_STRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_STRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_STRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_STRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS3_ZTRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_ZTRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_ZTRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_ZTRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_ZTRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_ZTRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_ZTRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_ZTRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS3_CTRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_CTRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_CTRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_CTRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_CTRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_CTRSM_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_CTRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_CTRSM_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

#endif
