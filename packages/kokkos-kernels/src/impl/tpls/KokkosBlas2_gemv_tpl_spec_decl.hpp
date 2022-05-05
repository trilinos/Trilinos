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

#ifndef KOKKOSBLAS2_GEMV_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS2_GEMV_TPL_SPEC_DECL_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS2_GEMV_DETERMINE_ARGS(LAYOUTA)                             \
  bool A_is_lr      = std::is_same<Kokkos::LayoutRight, LAYOUTA>::value;     \
  const int M       = static_cast<int>(A_is_lr ? A.extent(1) : A.extent(0)); \
  const int N       = static_cast<int>(A_is_lr ? A.extent(0) : A.extent(1)); \
  constexpr int one = 1;                                                     \
  const int LDA     = A_is_lr ? A.stride(0) : A.stride(1);                   \
                                                                             \
  char transa;                                                               \
  if ((trans[0] == 'N') || (trans[0] == 'n'))                                \
    transa = A_is_lr ? 'T' : 'N';                                            \
  else if ((trans[0] == 'T') || (trans[0] == 't'))                           \
    transa = A_is_lr ? 'N' : 'T';                                            \
  else {                                                                     \
    if (A_is_lr)                                                             \
      throw std::runtime_error(                                              \
          "Error: HostBlas::gemv conjugate transpose requires LayoutLeft "   \
          "views.");                                                         \
    transa = 'C';                                                            \
  }

#define KOKKOSBLAS2_DGEMV_BLAS(LAYOUTA, LAYOUTX, LAYOUTY, MEM_SPACE,         \
                               ETI_SPEC_AVAIL)                               \
  template <class ExecSpace>                                                 \
  struct GEMV<                                                               \
      Kokkos::View<const double**, LAYOUTA,                                  \
                   Kokkos::Device<ExecSpace, MEM_SPACE>,                     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                \
      Kokkos::View<const double*, LAYOUTX,                                   \
                   Kokkos::Device<ExecSpace, MEM_SPACE>,                     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                \
      Kokkos::View<double*, LAYOUTY, Kokkos::Device<ExecSpace, MEM_SPACE>,   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                \
      true, ETI_SPEC_AVAIL> {                                                \
    typedef double SCALAR;                                                   \
    typedef Kokkos::View<const SCALAR**, LAYOUTA,                            \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >           \
        AViewType;                                                           \
    typedef Kokkos::View<const SCALAR*, LAYOUTX,                             \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >           \
        XViewType;                                                           \
    typedef Kokkos::View<SCALAR*, LAYOUTY,                                   \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >           \
        YViewType;                                                           \
                                                                             \
    static void gemv(const typename AViewType::execution_space& /* space */, \
                     const char trans[],                                     \
                     typename AViewType::const_value_type& alpha,            \
                     const AViewType& A, const XViewType& X,                 \
                     typename YViewType::const_value_type& beta,             \
                     const YViewType& Y) {                                   \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_BLAS,double]");    \
      KOKKOSBLAS2_GEMV_DETERMINE_ARGS(LAYOUTA);                              \
      HostBlas<double>::gemv(transa, M, N, alpha, A.data(), LDA, X.data(),   \
                             one, beta, Y.data(), one);                      \
      Kokkos::Profiling::popRegion();                                        \
    }                                                                        \
  };

#define KOKKOSBLAS2_SGEMV_BLAS(LAYOUTA, LAYOUTX, LAYOUTY, MEM_SPACE,           \
                               ETI_SPEC_AVAIL)                                 \
  template <class ExecSpace>                                                   \
  struct GEMV<                                                                 \
      Kokkos::View<const float**, LAYOUTA,                                     \
                   Kokkos::Device<ExecSpace, MEM_SPACE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<const float*, LAYOUTX,                                      \
                   Kokkos::Device<ExecSpace, MEM_SPACE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<float*, LAYOUTY, Kokkos::Device<ExecSpace, MEM_SPACE>,      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      true, ETI_SPEC_AVAIL> {                                                  \
    typedef float SCALAR;                                                      \
    typedef Kokkos::View<const SCALAR**, LAYOUTA,                              \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        AViewType;                                                             \
    typedef Kokkos::View<const SCALAR*, LAYOUTX,                               \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XViewType;                                                             \
    typedef Kokkos::View<SCALAR*, LAYOUTY,                                     \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        YViewType;                                                             \
                                                                               \
    static void gemv(const typename AViewType::execution_space& /* space */,   \
                     const char trans[],                                       \
                     typename AViewType::const_value_type& alpha,              \
                     const AViewType& A, const XViewType& X,                   \
                     typename YViewType::const_value_type& beta,               \
                     const YViewType& Y) {                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_BLAS,float]");       \
      KOKKOSBLAS2_GEMV_DETERMINE_ARGS(LAYOUTA);                                \
      HostBlas<float>::gemv(transa, M, N, alpha, A.data(), LDA, X.data(), one, \
                            beta, Y.data(), one);                              \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#define KOKKOSBLAS2_ZGEMV_BLAS(LAYOUTA, LAYOUTX, LAYOUTY, MEM_SPACE,         \
                               ETI_SPEC_AVAIL)                               \
  template <class ExecSpace>                                                 \
  struct GEMV<Kokkos::View<const Kokkos::complex<double>**, LAYOUTA,         \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,             \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,        \
              Kokkos::View<const Kokkos::complex<double>*, LAYOUTX,          \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,             \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,        \
              Kokkos::View<Kokkos::complex<double>*, LAYOUTY,                \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,             \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,        \
              true, ETI_SPEC_AVAIL> {                                        \
    typedef Kokkos::complex<double> SCALAR;                                  \
    typedef Kokkos::View<const SCALAR**, LAYOUTA,                            \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >           \
        AViewType;                                                           \
    typedef Kokkos::View<const SCALAR*, LAYOUTX,                             \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >           \
        XViewType;                                                           \
    typedef Kokkos::View<SCALAR*, LAYOUTY,                                   \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >           \
        YViewType;                                                           \
                                                                             \
    static void gemv(const typename AViewType::execution_space& /* space */, \
                     const char trans[],                                     \
                     typename AViewType::const_value_type& alpha,            \
                     const AViewType& A, const XViewType& X,                 \
                     typename YViewType::const_value_type& beta,             \
                     const YViewType& Y) {                                   \
      Kokkos::Profiling::pushRegion(                                         \
          "KokkosBlas::gemv[TPL_BLAS,complex<double>]");                     \
      KOKKOSBLAS2_GEMV_DETERMINE_ARGS(LAYOUTA);                              \
      const std::complex<double> alpha_val = alpha, beta_val = beta;         \
      HostBlas<std::complex<double> >::gemv(                                 \
          transa, M, N, alpha_val,                                           \
          reinterpret_cast<const std::complex<double>*>(A.data()), LDA,      \
          reinterpret_cast<const std::complex<double>*>(X.data()), one,      \
          beta_val, reinterpret_cast<std::complex<double>*>(Y.data()), one); \
      Kokkos::Profiling::popRegion();                                        \
    }                                                                        \
  };

#define KOKKOSBLAS2_CGEMV_BLAS(LAYOUTA, LAYOUTX, LAYOUTY, MEM_SPACE,         \
                               ETI_SPEC_AVAIL)                               \
  template <class ExecSpace>                                                 \
  struct GEMV<Kokkos::View<const Kokkos::complex<float>**, LAYOUTA,          \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,             \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,        \
              Kokkos::View<const Kokkos::complex<float>*, LAYOUTX,           \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,             \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,        \
              Kokkos::View<Kokkos::complex<float>*, LAYOUTY,                 \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,             \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,        \
              true, ETI_SPEC_AVAIL> {                                        \
    typedef Kokkos::complex<float> SCALAR;                                   \
    typedef Kokkos::View<const SCALAR**, LAYOUTA,                            \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >           \
        AViewType;                                                           \
    typedef Kokkos::View<const SCALAR*, LAYOUTX,                             \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >           \
        XViewType;                                                           \
    typedef Kokkos::View<SCALAR*, LAYOUTY,                                   \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >           \
        YViewType;                                                           \
                                                                             \
    static void gemv(const typename AViewType::execution_space& /* space */, \
                     const char trans[],                                     \
                     typename AViewType::const_value_type& alpha,            \
                     const AViewType& A, const XViewType& X,                 \
                     typename YViewType::const_value_type& beta,             \
                     const YViewType& Y) {                                   \
      Kokkos::Profiling::pushRegion(                                         \
          "KokkosBlas::gemv[TPL_BLAS,complex<float>]");                      \
      KOKKOSBLAS2_GEMV_DETERMINE_ARGS(LAYOUTA);                              \
      const std::complex<float> alpha_val = alpha, beta_val = beta;          \
      HostBlas<std::complex<float> >::gemv(                                  \
          transa, M, N, alpha_val,                                           \
          reinterpret_cast<const std::complex<float>*>(A.data()), LDA,       \
          reinterpret_cast<const std::complex<float>*>(X.data()), one,       \
          beta_val, reinterpret_cast<std::complex<float>*>(Y.data()), one);  \
      Kokkos::Profiling::popRegion();                                        \
    }                                                                        \
  };

KOKKOSBLAS2_DGEMV_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                       Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS2_DGEMV_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                       Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS2_DGEMV_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                       Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS2_DGEMV_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                       Kokkos::LayoutRight, Kokkos::HostSpace, false)

KOKKOSBLAS2_SGEMV_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                       Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS2_SGEMV_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                       Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS2_SGEMV_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                       Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS2_SGEMV_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                       Kokkos::LayoutRight, Kokkos::HostSpace, false)

KOKKOSBLAS2_ZGEMV_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                       Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS2_ZGEMV_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                       Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS2_ZGEMV_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                       Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS2_ZGEMV_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                       Kokkos::LayoutRight, Kokkos::HostSpace, false)

KOKKOSBLAS2_CGEMV_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                       Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS2_CGEMV_BLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                       Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS2_CGEMV_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                       Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS2_CGEMV_BLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                       Kokkos::LayoutRight, Kokkos::HostSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS2_GEMV_CUBLAS_DETERMINE_ARGS(LAYOUTA)                      \
  bool A_is_lr      = std::is_same<Kokkos::LayoutRight, LAYOUTA>::value;     \
  const int M       = static_cast<int>(A_is_lr ? A.extent(1) : A.extent(0)); \
  const int N       = static_cast<int>(A_is_lr ? A.extent(0) : A.extent(1)); \
  constexpr int one = 1;                                                     \
  const int LDA     = A_is_lr ? A.stride(0) : A.stride(1);                   \
                                                                             \
  cublasOperation_t transa;                                                  \
  if ((trans[0] == 'N') || (trans[0] == 'n'))                                \
    transa = A_is_lr ? CUBLAS_OP_T : CUBLAS_OP_N;                            \
  else if ((trans[0] == 'T') || (trans[0] == 't'))                           \
    transa = A_is_lr ? CUBLAS_OP_N : CUBLAS_OP_T;                            \
  else {                                                                     \
    if (A_is_lr)                                                             \
      throw std::runtime_error(                                              \
          "Error: cublas gemv conjugate transpose requires LayoutLeft "      \
          "views.");                                                         \
    transa = CUBLAS_OP_C;                                                    \
  }

#define KOKKOSBLAS2_DGEMV_CUBLAS(LAYOUTA, LAYOUTX, LAYOUTY, MEM_SPACE,         \
                                 ETI_SPEC_AVAIL)                               \
  template <class ExecSpace>                                                   \
  struct GEMV<                                                                 \
      Kokkos::View<const double**, LAYOUTA,                                    \
                   Kokkos::Device<ExecSpace, MEM_SPACE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<const double*, LAYOUTX,                                     \
                   Kokkos::Device<ExecSpace, MEM_SPACE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<double*, LAYOUTY, Kokkos::Device<ExecSpace, MEM_SPACE>,     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      true, ETI_SPEC_AVAIL> {                                                  \
    typedef double SCALAR;                                                     \
    typedef Kokkos::View<const SCALAR**, LAYOUTA,                              \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        AViewType;                                                             \
    typedef Kokkos::View<const SCALAR*, LAYOUTX,                               \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XViewType;                                                             \
    typedef Kokkos::View<SCALAR*, LAYOUTY,                                     \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        YViewType;                                                             \
                                                                               \
    static void gemv(const typename AViewType::execution_space& space,         \
                     const char trans[],                                       \
                     typename AViewType::const_value_type& alpha,              \
                     const AViewType& A, const XViewType& X,                   \
                     typename YViewType::const_value_type& beta,               \
                     const YViewType& Y) {                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_CUBLAS,double]");    \
      KOKKOSBLAS2_GEMV_CUBLAS_DETERMINE_ARGS(LAYOUTA);                         \
      KokkosBlas::Impl::CudaBlasSingleton& s =                                 \
          KokkosBlas::Impl::CudaBlasSingleton::singleton();                    \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(                                            \
          cublasSetStream(s.handle, space.cuda_stream()));                     \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasDgemv(s.handle, transa, M, N, &alpha, \
                                               A.data(), LDA, X.data(), one,   \
                                               &beta, Y.data(), one));         \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));           \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#define KOKKOSBLAS2_SGEMV_CUBLAS(LAYOUTA, LAYOUTX, LAYOUTY, MEM_SPACE,         \
                                 ETI_SPEC_AVAIL)                               \
  template <class ExecSpace>                                                   \
  struct GEMV<                                                                 \
      Kokkos::View<const float**, LAYOUTA,                                     \
                   Kokkos::Device<ExecSpace, MEM_SPACE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<const float*, LAYOUTX,                                      \
                   Kokkos::Device<ExecSpace, MEM_SPACE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<float*, LAYOUTY, Kokkos::Device<ExecSpace, MEM_SPACE>,      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      true, ETI_SPEC_AVAIL> {                                                  \
    typedef float SCALAR;                                                      \
    typedef Kokkos::View<const SCALAR**, LAYOUTA,                              \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        AViewType;                                                             \
    typedef Kokkos::View<const SCALAR*, LAYOUTX,                               \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XViewType;                                                             \
    typedef Kokkos::View<SCALAR*, LAYOUTY,                                     \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        YViewType;                                                             \
                                                                               \
    static void gemv(const typename AViewType::execution_space& space,         \
                     const char trans[],                                       \
                     typename AViewType::const_value_type& alpha,              \
                     const AViewType& A, const XViewType& X,                   \
                     typename YViewType::const_value_type& beta,               \
                     const YViewType& Y) {                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_CUBLAS,float]");     \
      KOKKOSBLAS2_GEMV_CUBLAS_DETERMINE_ARGS(LAYOUTA);                         \
      KokkosBlas::Impl::CudaBlasSingleton& s =                                 \
          KokkosBlas::Impl::CudaBlasSingleton::singleton();                    \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(                                            \
          cublasSetStream(s.handle, space.cuda_stream()));                     \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSgemv(s.handle, transa, M, N, &alpha, \
                                               A.data(), LDA, X.data(), one,   \
                                               &beta, Y.data(), one));         \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));           \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#define KOKKOSBLAS2_ZGEMV_CUBLAS(LAYOUTA, LAYOUTX, LAYOUTY, MEM_SPACE,         \
                                 ETI_SPEC_AVAIL)                               \
  template <class ExecSpace>                                                   \
  struct GEMV<Kokkos::View<const Kokkos::complex<double>**, LAYOUTA,           \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,               \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,          \
              Kokkos::View<const Kokkos::complex<double>*, LAYOUTX,            \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,               \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,          \
              Kokkos::View<Kokkos::complex<double>*, LAYOUTY,                  \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,               \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,          \
              true, ETI_SPEC_AVAIL> {                                          \
    typedef Kokkos::complex<double> SCALAR;                                    \
    typedef Kokkos::View<const SCALAR**, LAYOUTA,                              \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        AViewType;                                                             \
    typedef Kokkos::View<const SCALAR*, LAYOUTX,                               \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XViewType;                                                             \
    typedef Kokkos::View<SCALAR*, LAYOUTY,                                     \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        YViewType;                                                             \
                                                                               \
    static void gemv(const typename AViewType::execution_space& space,         \
                     const char trans[],                                       \
                     typename AViewType::const_value_type& alpha,              \
                     const AViewType& A, const XViewType& X,                   \
                     typename YViewType::const_value_type& beta,               \
                     const YViewType& Y) {                                     \
      Kokkos::Profiling::pushRegion(                                           \
          "KokkosBlas::gemv[TPL_CUBLAS,complex<double>]");                     \
      KOKKOSBLAS2_GEMV_CUBLAS_DETERMINE_ARGS(LAYOUTA);                         \
      KokkosBlas::Impl::CudaBlasSingleton& s =                                 \
          KokkosBlas::Impl::CudaBlasSingleton::singleton();                    \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(                                            \
          cublasSetStream(s.handle, space.cuda_stream()));                     \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(                                            \
          cublasZgemv(s.handle, transa, M, N,                                  \
                      reinterpret_cast<const cuDoubleComplex*>(&alpha),        \
                      reinterpret_cast<const cuDoubleComplex*>(A.data()), LDA, \
                      reinterpret_cast<const cuDoubleComplex*>(X.data()), one, \
                      reinterpret_cast<const cuDoubleComplex*>(&beta),         \
                      reinterpret_cast<cuDoubleComplex*>(Y.data()), one));     \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));           \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#define KOKKOSBLAS2_CGEMV_CUBLAS(LAYOUTA, LAYOUTX, LAYOUTY, MEM_SPACE,        \
                                 ETI_SPEC_AVAIL)                              \
  template <class ExecSpace>                                                  \
  struct GEMV<Kokkos::View<const Kokkos::complex<float>**, LAYOUTA,           \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
              Kokkos::View<const Kokkos::complex<float>*, LAYOUTX,            \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
              Kokkos::View<Kokkos::complex<float>*, LAYOUTY,                  \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
              true, ETI_SPEC_AVAIL> {                                         \
    typedef Kokkos::complex<float> SCALAR;                                    \
    typedef Kokkos::View<const SCALAR**, LAYOUTA,                             \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        AViewType;                                                            \
    typedef Kokkos::View<const SCALAR*, LAYOUTX,                              \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        XViewType;                                                            \
    typedef Kokkos::View<SCALAR*, LAYOUTY,                                    \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        YViewType;                                                            \
                                                                              \
    static void gemv(const typename AViewType::execution_space& space,        \
                     const char trans[],                                      \
                     typename AViewType::const_value_type& alpha,             \
                     const AViewType& A, const XViewType& X,                  \
                     typename YViewType::const_value_type& beta,              \
                     const YViewType& Y) {                                    \
      Kokkos::Profiling::pushRegion(                                          \
          "KokkosBlas::gemv[TPL_CUBLAS,complex<float>]");                     \
      KOKKOSBLAS2_GEMV_CUBLAS_DETERMINE_ARGS(LAYOUTA);                        \
      KokkosBlas::Impl::CudaBlasSingleton& s =                                \
          KokkosBlas::Impl::CudaBlasSingleton::singleton();                   \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(                                           \
          cublasSetStream(s.handle, space.cuda_stream()));                    \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasCgemv(                               \
          s.handle, transa, M, N, reinterpret_cast<const cuComplex*>(&alpha), \
          reinterpret_cast<const cuComplex*>(A.data()), LDA,                  \
          reinterpret_cast<const cuComplex*>(X.data()), one,                  \
          reinterpret_cast<const cuComplex*>(&beta),                          \
          reinterpret_cast<cuComplex*>(Y.data()), one));                      \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));          \
      Kokkos::Profiling::popRegion();                                         \
    }                                                                         \
  };

KOKKOSBLAS2_DGEMV_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                         Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS2_DGEMV_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                         Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS2_DGEMV_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                         Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS2_DGEMV_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                         Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS2_SGEMV_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                         Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS2_SGEMV_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                         Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS2_SGEMV_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                         Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS2_SGEMV_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                         Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS2_ZGEMV_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                         Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS2_ZGEMV_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                         Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS2_ZGEMV_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                         Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS2_ZGEMV_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                         Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS2_CGEMV_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                         Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS2_CGEMV_CUBLAS(Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                         Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS2_CGEMV_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                         Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS2_CGEMV_CUBLAS(Kokkos::LayoutRight, Kokkos::LayoutRight,
                         Kokkos::LayoutRight, Kokkos::CudaSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS2_GEMV_ROCBLAS_DETERMINE_ARGS(LAYOUT)                      \
  bool A_is_lr      = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;      \
  const int M       = static_cast<int>(A_is_lr ? A.extent(1) : A.extent(0)); \
  const int N       = static_cast<int>(A_is_lr ? A.extent(0) : A.extent(1)); \
  constexpr int one = 1;                                                     \
  const int LDA     = A_is_lr ? A.stride(0) : A.stride(1);                   \
                                                                             \
  rocblas_operation transa;                                                  \
  if ((trans[0] == 'N') || (trans[0] == 'n'))                                \
    transa = A_is_lr ? rocblas_operation_transpose : rocblas_operation_none; \
  else if ((trans[0] == 'T') || (trans[0] == 't'))                           \
    transa = A_is_lr ? rocblas_operation_none : rocblas_operation_transpose; \
  else {                                                                     \
    if (A_is_lr)                                                             \
      throw std::runtime_error(                                              \
          "Error: rocBLAS Xgemv conjugate transpose requires LayoutLeft "    \
          "matrix.");                                                        \
    transa = rocblas_operation_conjugate_transpose;                          \
  }

#define KOKKOSBLAS2_DGEMV_ROCBLAS(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)           \
  template <>                                                                  \
  struct GEMV<                                                                 \
      Kokkos::View<const double**, LAYOUT,                                     \
                   Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<const double*, LAYOUT,                                      \
                   Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<double*, LAYOUT,                                            \
                   Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      true, ETI_SPEC_AVAIL> {                                                  \
    typedef double SCALAR;                                                     \
    typedef Kokkos::View<const SCALAR**, LAYOUT,                               \
                         Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        AViewType;                                                             \
    typedef Kokkos::View<const SCALAR*, LAYOUT,                                \
                         Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XViewType;                                                             \
    typedef Kokkos::View<SCALAR*, LAYOUT,                                      \
                         Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        YViewType;                                                             \
                                                                               \
    static void gemv(const typename AViewType::execution_space& space,         \
                     const char trans[],                                       \
                     typename AViewType::const_value_type& alpha,              \
                     const AViewType& A, const XViewType& X,                   \
                     typename YViewType::const_value_type& beta,               \
                     const YViewType& Y) {                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_ROCBLAS,double]");   \
      KOKKOSBLAS2_GEMV_ROCBLAS_DETERMINE_ARGS(LAYOUT);                         \
      KokkosBlas::Impl::RocBlasSingleton& s =                                  \
          KokkosBlas::Impl::RocBlasSingleton::singleton();                     \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(                                           \
          rocblas_dgemv(s.handle, transa, M, N, &alpha, A.data(), LDA,         \
                        X.data(), one, &beta, Y.data(), one));                 \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#define KOKKOSBLAS2_SGEMV_ROCBLAS(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)           \
  template <>                                                                  \
  struct GEMV<                                                                 \
      Kokkos::View<const float**, LAYOUT,                                      \
                   Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<const float*, LAYOUT,                                       \
                   Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<float*, LAYOUT,                                             \
                   Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      true, ETI_SPEC_AVAIL> {                                                  \
    typedef float SCALAR;                                                      \
    typedef Kokkos::View<const SCALAR**, LAYOUT,                               \
                         Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        AViewType;                                                             \
    typedef Kokkos::View<const SCALAR*, LAYOUT,                                \
                         Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XViewType;                                                             \
    typedef Kokkos::View<SCALAR*, LAYOUT,                                      \
                         Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        YViewType;                                                             \
                                                                               \
    static void gemv(const typename AViewType::execution_space& space,         \
                     const char trans[],                                       \
                     typename AViewType::const_value_type& alpha,              \
                     const AViewType& A, const XViewType& X,                   \
                     typename YViewType::const_value_type& beta,               \
                     const YViewType& Y) {                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_ROCBLAS,float]");    \
      KOKKOSBLAS2_GEMV_ROCBLAS_DETERMINE_ARGS(LAYOUT);                         \
      KokkosBlas::Impl::RocBlasSingleton& s =                                  \
          KokkosBlas::Impl::RocBlasSingleton::singleton();                     \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(                                           \
          rocblas_sgemv(s.handle, transa, M, N, &alpha, A.data(), LDA,         \
                        X.data(), one, &beta, Y.data(), one));                 \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#define KOKKOSBLAS2_ZGEMV_ROCBLAS(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)           \
  template <>                                                                  \
  struct GEMV<                                                                 \
      Kokkos::View<const Kokkos::complex<double>**, LAYOUT,                    \
                   Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<const Kokkos::complex<double>*, LAYOUT,                     \
                   Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<Kokkos::complex<double>*, LAYOUT,                           \
                   Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      true, ETI_SPEC_AVAIL> {                                                  \
    typedef Kokkos::complex<double> SCALAR;                                    \
    typedef Kokkos::View<const SCALAR**, LAYOUT,                               \
                         Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        AViewType;                                                             \
    typedef Kokkos::View<const SCALAR*, LAYOUT,                                \
                         Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XViewType;                                                             \
    typedef Kokkos::View<SCALAR*, LAYOUT,                                      \
                         Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        YViewType;                                                             \
                                                                               \
    static void gemv(const typename AViewType::execution_space& space,         \
                     const char trans[],                                       \
                     typename AViewType::const_value_type& alpha,              \
                     const AViewType& A, const XViewType& X,                   \
                     typename YViewType::const_value_type& beta,               \
                     const YViewType& Y) {                                     \
      Kokkos::Profiling::pushRegion(                                           \
          "KokkosBlas::gemv[TPL_ROCBLAS,complex<double>]");                    \
      KOKKOSBLAS2_GEMV_ROCBLAS_DETERMINE_ARGS(LAYOUT);                         \
      KokkosBlas::Impl::RocBlasSingleton& s =                                  \
          KokkosBlas::Impl::RocBlasSingleton::singleton();                     \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_zgemv(                             \
          s.handle, transa, M, N,                                              \
          reinterpret_cast<const rocblas_double_complex*>(&alpha),             \
          reinterpret_cast<const rocblas_double_complex*>(A.data()), LDA,      \
          reinterpret_cast<const rocblas_double_complex*>(X.data()), one,      \
          reinterpret_cast<const rocblas_double_complex*>(&beta),              \
          reinterpret_cast<rocblas_double_complex*>(Y.data()), one));          \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#define KOKKOSBLAS2_CGEMV_ROCBLAS(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)           \
  template <>                                                                  \
  struct GEMV<                                                                 \
      Kokkos::View<const Kokkos::complex<float>**, LAYOUT,                     \
                   Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<const Kokkos::complex<float>*, LAYOUT,                      \
                   Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<Kokkos::complex<float>*, LAYOUT,                            \
                   Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      true, ETI_SPEC_AVAIL> {                                                  \
    typedef Kokkos::complex<float> SCALAR;                                     \
    typedef Kokkos::View<const SCALAR**, LAYOUT,                               \
                         Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        AViewType;                                                             \
    typedef Kokkos::View<const SCALAR*, LAYOUT,                                \
                         Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XViewType;                                                             \
    typedef Kokkos::View<SCALAR*, LAYOUT,                                      \
                         Kokkos::Device<Kokkos::Experimental::HIP, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        YViewType;                                                             \
                                                                               \
    static void gemv(const typename AViewType::execution_space& space,         \
                     const char trans[],                                       \
                     typename AViewType::const_value_type& alpha,              \
                     const AViewType& A, const XViewType& X,                   \
                     typename YViewType::const_value_type& beta,               \
                     const YViewType& Y) {                                     \
      Kokkos::Profiling::pushRegion(                                           \
          "KokkosBlas::gemv[TPL_ROCBLAS,complex<float>]");                     \
      KOKKOSBLAS2_GEMV_ROCBLAS_DETERMINE_ARGS(LAYOUT);                         \
      KokkosBlas::Impl::RocBlasSingleton& s =                                  \
          KokkosBlas::Impl::RocBlasSingleton::singleton();                     \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_cgemv(                             \
          s.handle, transa, M, N,                                              \
          reinterpret_cast<const rocblas_float_complex*>(&alpha),              \
          reinterpret_cast<const rocblas_float_complex*>(A.data()), LDA,       \
          reinterpret_cast<const rocblas_float_complex*>(X.data()), one,       \
          reinterpret_cast<const rocblas_float_complex*>(&beta),               \
          reinterpret_cast<rocblas_float_complex*>(Y.data()), one));           \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

KOKKOSBLAS2_DGEMV_ROCBLAS(Kokkos::LayoutLeft, Kokkos::Experimental::HIPSpace,
                          true)
KOKKOSBLAS2_DGEMV_ROCBLAS(Kokkos::LayoutLeft, Kokkos::Experimental::HIPSpace,
                          false)
KOKKOSBLAS2_DGEMV_ROCBLAS(Kokkos::LayoutRight, Kokkos::Experimental::HIPSpace,
                          true)
KOKKOSBLAS2_DGEMV_ROCBLAS(Kokkos::LayoutRight, Kokkos::Experimental::HIPSpace,
                          false)

KOKKOSBLAS2_SGEMV_ROCBLAS(Kokkos::LayoutLeft, Kokkos::Experimental::HIPSpace,
                          true)
KOKKOSBLAS2_SGEMV_ROCBLAS(Kokkos::LayoutLeft, Kokkos::Experimental::HIPSpace,
                          false)
KOKKOSBLAS2_SGEMV_ROCBLAS(Kokkos::LayoutRight, Kokkos::Experimental::HIPSpace,
                          true)
KOKKOSBLAS2_SGEMV_ROCBLAS(Kokkos::LayoutRight, Kokkos::Experimental::HIPSpace,
                          false)

KOKKOSBLAS2_ZGEMV_ROCBLAS(Kokkos::LayoutLeft, Kokkos::Experimental::HIPSpace,
                          true)
KOKKOSBLAS2_ZGEMV_ROCBLAS(Kokkos::LayoutLeft, Kokkos::Experimental::HIPSpace,
                          false)
KOKKOSBLAS2_ZGEMV_ROCBLAS(Kokkos::LayoutRight, Kokkos::Experimental::HIPSpace,
                          true)
KOKKOSBLAS2_ZGEMV_ROCBLAS(Kokkos::LayoutRight, Kokkos::Experimental::HIPSpace,
                          false)

KOKKOSBLAS2_CGEMV_ROCBLAS(Kokkos::LayoutLeft, Kokkos::Experimental::HIPSpace,
                          true)
KOKKOSBLAS2_CGEMV_ROCBLAS(Kokkos::LayoutLeft, Kokkos::Experimental::HIPSpace,
                          false)
KOKKOSBLAS2_CGEMV_ROCBLAS(Kokkos::LayoutRight, Kokkos::Experimental::HIPSpace,
                          true)
KOKKOSBLAS2_CGEMV_ROCBLAS(Kokkos::LayoutRight, Kokkos::Experimental::HIPSpace,
                          false)

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCBLAS

#endif
