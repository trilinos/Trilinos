/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSBLAS_TRTRI_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS_TRTRI_TPL_SPEC_DECL_HPP_

#include "KokkosBlas_Host_tpl.hpp"  // trtri prototype
#include "KokkosBlas_tpl_spec.hpp"

namespace KokkosBlas {
namespace Impl {

#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#define KOKKOSBLAS_TRTRI_BLAS_HOST(SCALAR_TYPE, BASE_SCALAR_TYPE, LAYOUTA,     \
                                   MEM_SPACE, ETI_SPEC_AVAIL)                  \
  template <class ExecSpace>                                                   \
  struct TRTRI<Kokkos::View<int, LAYOUTA, Kokkos::HostSpace,                   \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
               Kokkos::View<SCALAR_TYPE**, LAYOUTA,                            \
                            Kokkos::Device<ExecSpace, MEM_SPACE>,              \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
               true, ETI_SPEC_AVAIL> {                                         \
    typedef SCALAR_TYPE SCALAR;                                                \
    typedef Kokkos::View<int, LAYOUTA, Kokkos::HostSpace,                      \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        RViewType;                                                             \
    typedef Kokkos::View<const SCALAR_TYPE**, LAYOUTA,                         \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        AViewType;                                                             \
                                                                               \
    static void trtri(const RViewType& R, const char uplo[],                   \
                      const char diag[], const AViewType& A) {                 \
      Kokkos::Profiling::pushRegion("KokkosBlas::trtri[TPL_BLAS," #SCALAR_TYPE \
                                    "]");                                      \
      const int M = static_cast<int>(A.extent(0));                             \
                                                                               \
      bool A_is_layout_left =                                                  \
          std::is_same<Kokkos::LayoutLeft, LAYOUTA>::value;                    \
                                                                               \
      const int AST = A_is_layout_left ? A.stride(1) : A.stride(0),            \
                LDA = (AST == 0) ? 1 : AST;                                    \
                                                                               \
      char uplo_;                                                              \
                                                                               \
      if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                \
        uplo_ = A_is_layout_left ? 'L' : 'U';                                  \
      else                                                                     \
        uplo_ = A_is_layout_left ? 'U' : 'L';                                  \
                                                                               \
      R() = HostBlas<BASE_SCALAR_TYPE>::trtri(                                 \
          uplo_, diag[0], M,                                                   \
          reinterpret_cast<const BASE_SCALAR_TYPE*>(A.data()), LDA);           \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };
#else
#define KOKKOSBLAS_TRTRI_BLAS_HOST(SCALAR_TYPE, BASE_SCALAR_TYPE, LAYOUTA, \
                                   MEM_SPACE, ETI_SPEC_AVAIL)
#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
#define KOKKOSBLAS_TRTRI_BLAS_MAGMA(SCALAR_TYPE, BASE_SCALAR_TYPE, MAGMA_FN,   \
                                    LAYOUTA, MEM_SPACE, ETI_SPEC_AVAIL)        \
  template <class ExecSpace>                                                   \
  struct TRTRI<Kokkos::View<int, LAYOUTA, Kokkos::HostSpace,                   \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
               Kokkos::View<SCALAR_TYPE**, LAYOUTA,                            \
                            Kokkos::Device<ExecSpace, MEM_SPACE>,              \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
               true, ETI_SPEC_AVAIL> {                                         \
    typedef SCALAR_TYPE SCALAR;                                                \
    typedef Kokkos::View<int, LAYOUTA, Kokkos::HostSpace,                      \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        RViewType;                                                             \
    typedef Kokkos::View<const SCALAR_TYPE**, LAYOUTA,                         \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        AViewType;                                                             \
                                                                               \
    static void trtri(const RViewType& R, const char uplo[],                   \
                      const char diag[], const AViewType& A) {                 \
      Kokkos::Profiling::pushRegion("KokkosBlas::trtri[TPL_BLAS," #SCALAR_TYPE \
                                    "]");                                      \
      magma_int_t M = static_cast<magma_int_t>(A.extent(0));                   \
                                                                               \
      bool A_is_layout_left =                                                  \
          std::is_same<Kokkos::LayoutLeft, LAYOUTA>::value;                    \
                                                                               \
      magma_int_t AST  = A_is_layout_left ? A.stride(1) : A.stride(0),         \
                  LDA  = (AST == 0) ? 1 : AST;                                 \
      magma_int_t info = 0;                                                    \
      magma_uplo_t uplo_;                                                      \
      magma_diag_t diag_;                                                      \
                                                                               \
      if ((uplo[0] == 'L') || (uplo[0] == 'l'))                                \
        uplo_ = A_is_layout_left ? MagmaLower : MagmaUpper;                    \
      else                                                                     \
        uplo_ = A_is_layout_left ? MagmaUpper : MagmaLower;                    \
                                                                               \
      if (diag[0] == 'U' || diag[0] == 'u')                                    \
        diag_ = MagmaUnit;                                                     \
      else                                                                     \
        diag_ = MagmaNonUnit;                                                  \
                                                                               \
      KokkosBlas::Impl::MagmaSingleton& s =                                    \
          KokkosBlas::Impl::MagmaSingleton::singleton();                       \
      R() = MAGMA_FN(uplo_, diag_, M,                                          \
                     reinterpret_cast<BASE_SCALAR_TYPE>(                       \
                         const_cast<SCALAR_TYPE*>(A.data())),                  \
                     LDA, &info);                                              \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };
#else
#define KOKKOSBLAS_TRTRI_BLAS_MAGMA(SCALAR_TYPE, BASE_SCALAR_TYPE, MAGMA_FN, \
                                    LAYOUTA, MEM_SPACE, ETI_SPEC_AVAIL)
#endif  // KOKKOSKERNELS_ENABLE_TPL_MAGMA

// Explicitly define the TRTRI class for all permutations listed below

// Handle type and space permutations
#define KOKKOSBLAS_DTRTRI_BLAS(LAYOUTA, ETI_SPEC_AVAIL)                   \
  KOKKOSBLAS_TRTRI_BLAS_HOST(double, double, LAYOUTA, Kokkos::HostSpace,  \
                             ETI_SPEC_AVAIL)                              \
  KOKKOSBLAS_TRTRI_BLAS_MAGMA(double, magmaDouble_ptr, magma_dtrtri_gpu,  \
                              LAYOUTA, Kokkos::CudaSpace, ETI_SPEC_AVAIL) \
  KOKKOSBLAS_TRTRI_BLAS_MAGMA(double, magmaDouble_ptr, magma_dtrtri_gpu,  \
                              LAYOUTA, Kokkos::CudaUVMSpace, ETI_SPEC_AVAIL)

#define KOKKOSBLAS_STRTRI_BLAS(LAYOUTA, ETI_SPEC_AVAIL)                   \
  KOKKOSBLAS_TRTRI_BLAS_HOST(float, float, LAYOUTA, Kokkos::HostSpace,    \
                             ETI_SPEC_AVAIL)                              \
  KOKKOSBLAS_TRTRI_BLAS_MAGMA(float, magmaFloat_ptr, magma_strtri_gpu,    \
                              LAYOUTA, Kokkos::CudaSpace, ETI_SPEC_AVAIL) \
  KOKKOSBLAS_TRTRI_BLAS_MAGMA(float, magmaFloat_ptr, magma_strtri_gpu,    \
                              LAYOUTA, Kokkos::CudaUVMSpace, ETI_SPEC_AVAIL)

#define KOKKOSBLAS_ZTRTRI_BLAS(LAYOUTA, ETI_SPEC_AVAIL)                        \
  KOKKOSBLAS_TRTRI_BLAS_HOST(Kokkos::complex<double>, std::complex<double>,    \
                             LAYOUTA, Kokkos::HostSpace, ETI_SPEC_AVAIL)       \
  KOKKOSBLAS_TRTRI_BLAS_MAGMA(Kokkos::complex<double>, magmaDoubleComplex_ptr, \
                              magma_ztrtri_gpu, LAYOUTA, Kokkos::CudaSpace,    \
                              ETI_SPEC_AVAIL)                                  \
  KOKKOSBLAS_TRTRI_BLAS_MAGMA(Kokkos::complex<double>, magmaDoubleComplex_ptr, \
                              magma_ztrtri_gpu, LAYOUTA, Kokkos::CudaUVMSpace, \
                              ETI_SPEC_AVAIL)

#define KOKKOSBLAS_CTRTRI_BLAS(LAYOUTA, ETI_SPEC_AVAIL)                        \
  KOKKOSBLAS_TRTRI_BLAS_HOST(Kokkos::complex<float>, std::complex<float>,      \
                             LAYOUTA, Kokkos::HostSpace, ETI_SPEC_AVAIL)       \
  KOKKOSBLAS_TRTRI_BLAS_MAGMA(Kokkos::complex<float>, magmaFloatComplex_ptr,   \
                              magma_ctrtri_gpu, LAYOUTA, Kokkos::CudaSpace,    \
                              ETI_SPEC_AVAIL)                                  \
  KOKKOSBLAS_TRTRI_BLAS_MAGMA(Kokkos::complex<float>, magmaFloatComplex_ptr,   \
                              magma_ctrtri_gpu, LAYOUTA, Kokkos::CudaUVMSpace, \
                              ETI_SPEC_AVAIL)

// Handle layout permutations
KOKKOSBLAS_DTRTRI_BLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS_DTRTRI_BLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS_DTRTRI_BLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS_DTRTRI_BLAS(Kokkos::LayoutRight, false)

KOKKOSBLAS_STRTRI_BLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS_STRTRI_BLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS_STRTRI_BLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS_STRTRI_BLAS(Kokkos::LayoutRight, false)

KOKKOSBLAS_ZTRTRI_BLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS_ZTRTRI_BLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS_ZTRTRI_BLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS_ZTRTRI_BLAS(Kokkos::LayoutRight, false)

KOKKOSBLAS_CTRTRI_BLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS_CTRTRI_BLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS_CTRTRI_BLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS_CTRTRI_BLAS(Kokkos::LayoutRight, false)

}  // namespace Impl
}  // nameSpace KokkosBlas

#endif  // KOKKOSBLAS_TRTRI_TPL_SPEC_DECL_HPP_
