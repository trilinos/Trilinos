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

#ifndef KOKKOSLAPACK_TRTRI_TPL_SPEC_DECL_HPP_
#define KOKKOSLAPACK_TRTRI_TPL_SPEC_DECL_HPP_

#include "KokkosLapack_Host_tpl.hpp"  // trtri prototype

#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
#include "KokkosLapack_magma.hpp"
#endif

namespace KokkosLapack {
namespace Impl {

#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK
#define KOKKOSLAPACK_TRTRI_LAPACK_HOST(SCALAR_TYPE, BASE_SCALAR_TYPE, LAYOUTA, MEM_SPACE, ETI_SPEC_AVAIL)           \
  template <class ExecSpace>                                                                                        \
  struct TRTRI<Kokkos::View<int, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
               Kokkos::View<SCALAR_TYPE**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                           \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                              \
               true, ETI_SPEC_AVAIL> {                                                                              \
    typedef SCALAR_TYPE SCALAR;                                                                                     \
    typedef Kokkos::View<int, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >     \
        RViewType;                                                                                                  \
    typedef Kokkos::View<const SCALAR_TYPE**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                        \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                  \
        AViewType;                                                                                                  \
                                                                                                                    \
    static void trtri(const RViewType& R, const char uplo[], const char diag[], const AViewType& A) {               \
      Kokkos::Profiling::pushRegion("KokkosLapack::trtri[TPL_LAPACK," #SCALAR_TYPE "]");                            \
      const int M = static_cast<int>(A.extent(0));                                                                  \
                                                                                                                    \
      bool A_is_layout_left = std::is_same<Kokkos::LayoutLeft, LAYOUTA>::value;                                     \
                                                                                                                    \
      const int AST = A_is_layout_left ? A.stride(1) : A.stride(0), LDA = (AST == 0) ? 1 : AST;                     \
                                                                                                                    \
      char uplo_;                                                                                                   \
                                                                                                                    \
      if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                                     \
        uplo_ = A_is_layout_left ? 'L' : 'U';                                                                       \
      else                                                                                                          \
        uplo_ = A_is_layout_left ? 'U' : 'L';                                                                       \
                                                                                                                    \
      R() = HostLapack<BASE_SCALAR_TYPE>::trtri(uplo_, diag[0], M,                                                  \
                                                reinterpret_cast<const BASE_SCALAR_TYPE*>(A.data()), LDA);          \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };
#else
#define KOKKOSLAPACK_TRTRI_LAPACK_HOST(SCALAR_TYPE, BASE_SCALAR_TYPE, LAYOUTA, MEM_SPACE, ETI_SPEC_AVAIL)
#endif  // KOKKOSKERNELS_ENABLE_TPL_LAPACK

#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
#define KOKKOSLAPACK_TRTRI_LAPACK_MAGMA(SCALAR_TYPE, BASE_SCALAR_TYPE, MAGMA_FN, LAYOUTA, MEM_SPACE, ETI_SPEC_AVAIL) \
  template <class ExecSpace>                                                                                         \
  struct TRTRI<Kokkos::View<int, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,  \
               Kokkos::View<SCALAR_TYPE**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                            \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                               \
               true, ETI_SPEC_AVAIL> {                                                                               \
    typedef SCALAR_TYPE SCALAR;                                                                                      \
    typedef Kokkos::View<int, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >      \
        RViewType;                                                                                                   \
    typedef Kokkos::View<const SCALAR_TYPE**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>,                         \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                   \
        AViewType;                                                                                                   \
                                                                                                                     \
    static void trtri(const RViewType& R, const char uplo[], const char diag[], const AViewType& A) {                \
      Kokkos::Profiling::pushRegion("KokkosLapack::trtri[TPL_LAPACK," #SCALAR_TYPE "]");                             \
      magma_int_t M = static_cast<magma_int_t>(A.extent(0));                                                         \
                                                                                                                     \
      bool A_is_layout_left = std::is_same<Kokkos::LayoutLeft, LAYOUTA>::value;                                      \
                                                                                                                     \
      magma_int_t AST = A_is_layout_left ? A.stride(1) : A.stride(0), LDA = (AST == 0) ? 1 : AST;                    \
      magma_int_t info = 0;                                                                                          \
      magma_uplo_t uplo_;                                                                                            \
      magma_diag_t diag_;                                                                                            \
                                                                                                                     \
      if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                                                      \
        uplo_ = A_is_layout_left ? MagmaLower : MagmaUpper;                                                          \
      else                                                                                                           \
        uplo_ = A_is_layout_left ? MagmaUpper : MagmaLower;                                                          \
                                                                                                                     \
      if (diag[0] == 'U' || diag[0] == 'u')                                                                          \
        diag_ = MagmaUnit;                                                                                           \
      else                                                                                                           \
        diag_ = MagmaNonUnit;                                                                                        \
                                                                                                                     \
      KokkosLapack::Impl::MagmaSingleton& s = KokkosLapack::Impl::MagmaSingleton::singleton();                       \
      R() = MAGMA_FN(uplo_, diag_, M, reinterpret_cast<BASE_SCALAR_TYPE>(const_cast<SCALAR_TYPE*>(A.data())), LDA,   \
                     &info);                                                                                         \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };
#else
#define KOKKOSLAPACK_TRTRI_LAPACK_MAGMA(SCALAR_TYPE, BASE_SCALAR_TYPE, MAGMA_FN, LAYOUTA, MEM_SPACE, ETI_SPEC_AVAIL)
#endif  // KOKKOSKERNELS_ENABLE_TPL_MAGMA

// Explicitly define the TRTRI class for all permutations listed below

// Handle type and space permutations
#ifdef KOKKOS_ENABLE_CUDA

#define KOKKOSLAPACK_DTRTRI_LAPACK(LAYOUTA, ETI_SPEC_AVAIL)                                                 \
  KOKKOSLAPACK_TRTRI_LAPACK_HOST(double, double, LAYOUTA, Kokkos::HostSpace, ETI_SPEC_AVAIL)                \
  KOKKOSLAPACK_TRTRI_LAPACK_MAGMA(double, magmaDouble_ptr, magma_dtrtri_gpu, LAYOUTA, Kokkos::CudaSpace,    \
                                  ETI_SPEC_AVAIL)                                                           \
  KOKKOSLAPACK_TRTRI_LAPACK_MAGMA(double, magmaDouble_ptr, magma_dtrtri_gpu, LAYOUTA, Kokkos::CudaUVMSpace, \
                                  ETI_SPEC_AVAIL)

#define KOKKOSLAPACK_STRTRI_LAPACK(LAYOUTA, ETI_SPEC_AVAIL)                                                            \
  KOKKOSLAPACK_TRTRI_LAPACK_HOST(float, float, LAYOUTA, Kokkos::HostSpace, ETI_SPEC_AVAIL)                             \
  KOKKOSLAPACK_TRTRI_LAPACK_MAGMA(float, magmaFloat_ptr, magma_strtri_gpu, LAYOUTA, Kokkos::CudaSpace, ETI_SPEC_AVAIL) \
  KOKKOSLAPACK_TRTRI_LAPACK_MAGMA(float, magmaFloat_ptr, magma_strtri_gpu, LAYOUTA, Kokkos::CudaUVMSpace,              \
                                  ETI_SPEC_AVAIL)

#define KOKKOSLAPACK_ZTRTRI_LAPACK(LAYOUTA, ETI_SPEC_AVAIL)                                                   \
  KOKKOSLAPACK_TRTRI_LAPACK_HOST(Kokkos::complex<double>, std::complex<double>, LAYOUTA, Kokkos::HostSpace,   \
                                 ETI_SPEC_AVAIL)                                                              \
  KOKKOSLAPACK_TRTRI_LAPACK_MAGMA(Kokkos::complex<double>, magmaDoubleComplex_ptr, magma_ztrtri_gpu, LAYOUTA, \
                                  Kokkos::CudaSpace, ETI_SPEC_AVAIL)                                          \
  KOKKOSLAPACK_TRTRI_LAPACK_MAGMA(Kokkos::complex<double>, magmaDoubleComplex_ptr, magma_ztrtri_gpu, LAYOUTA, \
                                  Kokkos::CudaUVMSpace, ETI_SPEC_AVAIL)

#define KOKKOSLAPACK_CTRTRI_LAPACK(LAYOUTA, ETI_SPEC_AVAIL)                                                 \
  KOKKOSLAPACK_TRTRI_LAPACK_HOST(Kokkos::complex<float>, std::complex<float>, LAYOUTA, Kokkos::HostSpace,   \
                                 ETI_SPEC_AVAIL)                                                            \
  KOKKOSLAPACK_TRTRI_LAPACK_MAGMA(Kokkos::complex<float>, magmaFloatComplex_ptr, magma_ctrtri_gpu, LAYOUTA, \
                                  Kokkos::CudaSpace, ETI_SPEC_AVAIL)                                        \
  KOKKOSLAPACK_TRTRI_LAPACK_MAGMA(Kokkos::complex<float>, magmaFloatComplex_ptr, magma_ctrtri_gpu, LAYOUTA, \
                                  Kokkos::CudaUVMSpace, ETI_SPEC_AVAIL)

#else

#define KOKKOSLAPACK_DTRTRI_LAPACK(LAYOUTA, ETI_SPEC_AVAIL) \
  KOKKOSLAPACK_TRTRI_LAPACK_HOST(double, double, LAYOUTA, Kokkos::HostSpace, ETI_SPEC_AVAIL)

#define KOKKOSLAPACK_STRTRI_LAPACK(LAYOUTA, ETI_SPEC_AVAIL) \
  KOKKOSLAPACK_TRTRI_LAPACK_HOST(float, float, LAYOUTA, Kokkos::HostSpace, ETI_SPEC_AVAIL)

#define KOKKOSLAPACK_ZTRTRI_LAPACK(LAYOUTA, ETI_SPEC_AVAIL)                                                 \
  KOKKOSLAPACK_TRTRI_LAPACK_HOST(Kokkos::complex<double>, std::complex<double>, LAYOUTA, Kokkos::HostSpace, \
                                 ETI_SPEC_AVAIL)

#define KOKKOSLAPACK_CTRTRI_LAPACK(LAYOUTA, ETI_SPEC_AVAIL)                                               \
  KOKKOSLAPACK_TRTRI_LAPACK_HOST(Kokkos::complex<float>, std::complex<float>, LAYOUTA, Kokkos::HostSpace, \
                                 ETI_SPEC_AVAIL)

#endif

// Handle layout permutations
KOKKOSLAPACK_DTRTRI_LAPACK(Kokkos::LayoutLeft, true)
KOKKOSLAPACK_DTRTRI_LAPACK(Kokkos::LayoutLeft, false)
KOKKOSLAPACK_DTRTRI_LAPACK(Kokkos::LayoutRight, true)
KOKKOSLAPACK_DTRTRI_LAPACK(Kokkos::LayoutRight, false)

KOKKOSLAPACK_STRTRI_LAPACK(Kokkos::LayoutLeft, true)
KOKKOSLAPACK_STRTRI_LAPACK(Kokkos::LayoutLeft, false)
KOKKOSLAPACK_STRTRI_LAPACK(Kokkos::LayoutRight, true)
KOKKOSLAPACK_STRTRI_LAPACK(Kokkos::LayoutRight, false)

KOKKOSLAPACK_ZTRTRI_LAPACK(Kokkos::LayoutLeft, true)
KOKKOSLAPACK_ZTRTRI_LAPACK(Kokkos::LayoutLeft, false)
KOKKOSLAPACK_ZTRTRI_LAPACK(Kokkos::LayoutRight, true)
KOKKOSLAPACK_ZTRTRI_LAPACK(Kokkos::LayoutRight, false)

KOKKOSLAPACK_CTRTRI_LAPACK(Kokkos::LayoutLeft, true)
KOKKOSLAPACK_CTRTRI_LAPACK(Kokkos::LayoutLeft, false)
KOKKOSLAPACK_CTRTRI_LAPACK(Kokkos::LayoutRight, true)
KOKKOSLAPACK_CTRTRI_LAPACK(Kokkos::LayoutRight, false)

}  // namespace Impl
}  // nameSpace KokkosLapack

#endif  // KOKKOSLAPACK_TRTRI_TPL_SPEC_DECL_HPP_
