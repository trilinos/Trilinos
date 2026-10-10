// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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

#define KOKKOSBLAS2_DGEMV_BLAS(LAYOUT)                                                                              \
  template <typename ExecSpace, bool ETI_SPEC_AVAIL>                                                                \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                   \
  struct GEMV<ExecSpace, Kokkos::View<const double**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
              Kokkos::View<const double*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
              Kokkos::View<double*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,             \
              ETI_SPEC_AVAIL> {                                                                                     \
    using SCALAR    = double;                                                                                       \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;     \
    using YViewType = Kokkos::View<SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;           \
                                                                                                                    \
    static void gemv(const ExecSpace& /* space */, const char trans[], typename AViewType::const_value_type& alpha, \
                     const AViewType& A, const XViewType& X, typename YViewType::const_value_type& beta,            \
                     const YViewType& Y) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_BLAS,double]");                                           \
      KOKKOSBLAS2_GEMV_DETERMINE_ARGS(LAYOUT);                                                                      \
      HostBlas<double>::gemv(transa, M, N, alpha, A.data(), LDA, X.data(), one, beta, Y.data(), one);               \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS2_SGEMV_BLAS(LAYOUT)                                                                              \
  template <typename ExecSpace, bool ETI_SPEC_AVAIL>                                                                \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                   \
  struct GEMV<ExecSpace, Kokkos::View<const float**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,  \
              Kokkos::View<const float*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,              \
              Kokkos::View<float*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,              \
              ETI_SPEC_AVAIL> {                                                                                     \
    using SCALAR    = float;                                                                                        \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;     \
    using YViewType = Kokkos::View<SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;           \
                                                                                                                    \
    static void gemv(const ExecSpace& /* space */, const char trans[], typename AViewType::const_value_type& alpha, \
                     const AViewType& A, const XViewType& X, typename YViewType::const_value_type& beta,            \
                     const YViewType& Y) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_BLAS,float]");                                            \
      KOKKOSBLAS2_GEMV_DETERMINE_ARGS(LAYOUT);                                                                      \
      HostBlas<float>::gemv(transa, M, N, alpha, A.data(), LDA, X.data(), one, beta, Y.data(), one);                \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS2_ZGEMV_BLAS(LAYOUT)                                                                              \
  template <typename ExecSpace, bool ETI_SPEC_AVAIL>                                                                \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                   \
  struct GEMV<                                                                                                      \
      ExecSpace,                                                                                                    \
      Kokkos::View<const Kokkos::complex<double>**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,   \
      Kokkos::View<const Kokkos::complex<double>*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
      Kokkos::View<Kokkos::complex<double>*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,    \
      ETI_SPEC_AVAIL> {                                                                                             \
    using SCALAR    = Kokkos::complex<double>;                                                                      \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;     \
    using YViewType = Kokkos::View<SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;           \
                                                                                                                    \
    static void gemv(const ExecSpace& /* space */, const char trans[], typename AViewType::const_value_type& alpha, \
                     const AViewType& A, const XViewType& X, typename YViewType::const_value_type& beta,            \
                     const YViewType& Y) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_BLAS,complex<double>]");                                  \
      KOKKOSBLAS2_GEMV_DETERMINE_ARGS(LAYOUT);                                                                      \
      const std::complex<double> alpha_val = alpha, beta_val = beta;                                                \
      HostBlas<std::complex<double> >::gemv(transa, M, N, alpha_val,                                                \
                                            reinterpret_cast<const std::complex<double>*>(A.data()), LDA,           \
                                            reinterpret_cast<const std::complex<double>*>(X.data()), one, beta_val, \
                                            reinterpret_cast<std::complex<double>*>(Y.data()), one);                \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS2_CGEMV_BLAS(LAYOUT)                                                                              \
  template <typename ExecSpace, bool ETI_SPEC_AVAIL>                                                                \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                   \
  struct GEMV<                                                                                                      \
      ExecSpace,                                                                                                    \
      Kokkos::View<const Kokkos::complex<float>**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
      Kokkos::View<const Kokkos::complex<float>*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,     \
      Kokkos::View<Kokkos::complex<float>*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,     \
      ETI_SPEC_AVAIL> {                                                                                             \
    using SCALAR    = Kokkos::complex<float>;                                                                       \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;     \
    using YViewType = Kokkos::View<SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;           \
                                                                                                                    \
    static void gemv(const ExecSpace& /* space */, const char trans[], typename AViewType::const_value_type& alpha, \
                     const AViewType& A, const XViewType& X, typename YViewType::const_value_type& beta,            \
                     const YViewType& Y) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_BLAS,complex<float>]");                                   \
      KOKKOSBLAS2_GEMV_DETERMINE_ARGS(LAYOUT);                                                                      \
      const std::complex<float> alpha_val = alpha, beta_val = beta;                                                 \
      HostBlas<std::complex<float> >::gemv(transa, M, N, alpha_val,                                                 \
                                           reinterpret_cast<const std::complex<float>*>(A.data()), LDA,             \
                                           reinterpret_cast<const std::complex<float>*>(X.data()), one, beta_val,   \
                                           reinterpret_cast<std::complex<float>*>(Y.data()), one);                  \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

KOKKOSBLAS2_DGEMV_BLAS(Kokkos::LayoutLeft)
KOKKOSBLAS2_DGEMV_BLAS(Kokkos::LayoutRight)

KOKKOSBLAS2_SGEMV_BLAS(Kokkos::LayoutLeft)
KOKKOSBLAS2_SGEMV_BLAS(Kokkos::LayoutRight)

KOKKOSBLAS2_ZGEMV_BLAS(Kokkos::LayoutLeft)
KOKKOSBLAS2_ZGEMV_BLAS(Kokkos::LayoutRight)

KOKKOSBLAS2_CGEMV_BLAS(Kokkos::LayoutLeft)
KOKKOSBLAS2_CGEMV_BLAS(Kokkos::LayoutRight)

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

#define KOKKOSBLAS2_DGEMV_CUBLAS(LAYOUT)                                                                             \
  template <bool ETI_SPEC_AVAIL>                                                                                     \
  struct GEMV<                                                                                                       \
      Kokkos::Cuda, Kokkos::View<const double**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
      Kokkos::View<const double*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                   \
      Kokkos::View<double*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true, ETI_SPEC_AVAIL> { \
    using SCALAR    = double;                                                                                        \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;  \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;   \
    using YViewType = Kokkos::View<SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;         \
                                                                                                                     \
    static void gemv(const Kokkos::Cuda& space, const char trans[], typename AViewType::const_value_type& alpha,     \
                     const AViewType& A, const XViewType& X, typename YViewType::const_value_type& beta,             \
                     const YViewType& Y) {                                                                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_CUBLAS,double]");                                          \
      KOKKOSBLAS2_GEMV_CUBLAS_DETERMINE_ARGS(LAYOUT);                                                                \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                     \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));                              \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(                                                                              \
          cublasDgemv(s.handle, transa, M, N, &alpha, A.data(), LDA, X.data(), one, &beta, Y.data(), one));          \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, NULL));                                             \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS2_SGEMV_CUBLAS(LAYOUT)                                                                            \
  template <bool ETI_SPEC_AVAIL>                                                                                    \
  struct GEMV<                                                                                                      \
      Kokkos::Cuda, Kokkos::View<const float**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
      Kokkos::View<const float*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                   \
      Kokkos::View<float*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true, ETI_SPEC_AVAIL> { \
    using SCALAR    = float;                                                                                        \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >; \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;  \
    using YViewType = Kokkos::View<SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;        \
                                                                                                                    \
    static void gemv(const Kokkos::Cuda& space, const char trans[], typename AViewType::const_value_type& alpha,    \
                     const AViewType& A, const XViewType& X, typename YViewType::const_value_type& beta,            \
                     const YViewType& Y) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_CUBLAS,float]");                                          \
      KOKKOSBLAS2_GEMV_CUBLAS_DETERMINE_ARGS(LAYOUT);                                                               \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                    \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));                             \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(                                                                             \
          cublasSgemv(s.handle, transa, M, N, &alpha, A.data(), LDA, X.data(), one, &beta, Y.data(), one));         \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, NULL));                                            \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS2_ZGEMV_CUBLAS(LAYOUT)                                                                               \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct GEMV<                                                                                                         \
      Kokkos::Cuda,                                                                                                    \
      Kokkos::View<const Kokkos::complex<double>**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,   \
      Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
      Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,    \
      ETI_SPEC_AVAIL> {                                                                                                \
    using SCALAR    = Kokkos::complex<double>;                                                                         \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;     \
    using YViewType = Kokkos::View<SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;           \
                                                                                                                       \
    static void gemv(const Kokkos::Cuda& space, const char trans[], typename AViewType::const_value_type& alpha,       \
                     const AViewType& A, const XViewType& X, typename YViewType::const_value_type& beta,               \
                     const YViewType& Y) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_CUBLAS,complex<double>]");                                   \
      KOKKOSBLAS2_GEMV_CUBLAS_DETERMINE_ARGS(LAYOUT);                                                                  \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                       \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));                                \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasZgemv(                                                                    \
          s.handle, transa, M, N, reinterpret_cast<const cuDoubleComplex*>(&alpha),                                    \
          reinterpret_cast<const cuDoubleComplex*>(A.data()), LDA, reinterpret_cast<const cuDoubleComplex*>(X.data()), \
          one, reinterpret_cast<const cuDoubleComplex*>(&beta), reinterpret_cast<cuDoubleComplex*>(Y.data()), one));   \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, NULL));                                               \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS2_CGEMV_CUBLAS(LAYOUT)                                                                               \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct GEMV<                                                                                                         \
      Kokkos::Cuda,                                                                                                    \
      Kokkos::View<const Kokkos::complex<float>**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
      Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,     \
      Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,     \
      ETI_SPEC_AVAIL> {                                                                                                \
    using SCALAR    = Kokkos::complex<float>;                                                                          \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;     \
    using YViewType = Kokkos::View<SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;           \
                                                                                                                       \
    static void gemv(const Kokkos::Cuda& space, const char trans[], typename AViewType::const_value_type& alpha,       \
                     const AViewType& A, const XViewType& X, typename YViewType::const_value_type& beta,               \
                     const YViewType& Y) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_CUBLAS,complex<float>]");                                    \
      KOKKOSBLAS2_GEMV_CUBLAS_DETERMINE_ARGS(LAYOUT);                                                                  \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                       \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));                                \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(                                                                                \
          cublasCgemv(s.handle, transa, M, N, reinterpret_cast<const cuComplex*>(&alpha),                              \
                      reinterpret_cast<const cuComplex*>(A.data()), LDA, reinterpret_cast<const cuComplex*>(X.data()), \
                      one, reinterpret_cast<const cuComplex*>(&beta), reinterpret_cast<cuComplex*>(Y.data()), one));   \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, NULL));                                               \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS2_DGEMV_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS2_DGEMV_CUBLAS(Kokkos::LayoutRight)

KOKKOSBLAS2_SGEMV_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS2_SGEMV_CUBLAS(Kokkos::LayoutRight)

KOKKOSBLAS2_ZGEMV_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS2_ZGEMV_CUBLAS(Kokkos::LayoutRight)

KOKKOSBLAS2_CGEMV_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS2_CGEMV_CUBLAS(Kokkos::LayoutRight)

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

#define KOKKOSBLAS2_DGEMV_ROCBLAS(LAYOUT)                                                                           \
  template <bool ETI_SPEC_AVAIL>                                                                                    \
  struct GEMV<                                                                                                      \
      Kokkos::HIP, Kokkos::View<const double**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,     \
      Kokkos::View<const double*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                   \
      Kokkos::View<double*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true, ETI_SPEC_AVAIL> { \
    using SCALAR    = double;                                                                                       \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;  \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;   \
    using YViewType = Kokkos::View<SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;         \
                                                                                                                    \
    static void gemv(const Kokkos::HIP& space, const char trans[], typename AViewType::const_value_type& alpha,     \
                     const AViewType& A, const XViewType& X, typename YViewType::const_value_type& beta,            \
                     const YViewType& Y) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_ROCBLAS,double]");                                        \
      KOKKOSBLAS2_GEMV_ROCBLAS_DETERMINE_ARGS(LAYOUT);                                                              \
      KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                      \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                          \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                            \
          rocblas_dgemv(s.handle, transa, M, N, &alpha, A.data(), LDA, X.data(), one, &beta, Y.data(), one));       \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                        \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS2_SGEMV_ROCBLAS(LAYOUT)                                                                              \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct GEMV<Kokkos::HIP, Kokkos::View<const float**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
              Kokkos::View<const float*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
              Kokkos::View<float*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,               \
              ETI_SPEC_AVAIL> {                                                                                        \
    using SCALAR    = float;                                                                                           \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;     \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;      \
    using YViewType = Kokkos::View<SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;            \
                                                                                                                       \
    static void gemv(const Kokkos::HIP& space, const char trans[], typename AViewType::const_value_type& alpha,        \
                     const AViewType& A, const XViewType& X, typename YViewType::const_value_type& beta,               \
                     const YViewType& Y) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_ROCBLAS,float]");                                            \
      KOKKOSBLAS2_GEMV_ROCBLAS_DETERMINE_ARGS(LAYOUT);                                                                 \
      KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                         \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                             \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                               \
          rocblas_sgemv(s.handle, transa, M, N, &alpha, A.data(), LDA, X.data(), one, &beta, Y.data(), one));          \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                           \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS2_ZGEMV_ROCBLAS(LAYOUT)                                                                             \
  template <bool ETI_SPEC_AVAIL>                                                                                      \
  struct GEMV<                                                                                                        \
      Kokkos::HIP,                                                                                                    \
      Kokkos::View<const Kokkos::complex<double>**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,   \
      Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
      Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,    \
      ETI_SPEC_AVAIL> {                                                                                               \
    using SCALAR    = Kokkos::complex<double>;                                                                        \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;     \
    using YViewType = Kokkos::View<SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;           \
                                                                                                                      \
    static void gemv(const Kokkos::HIP& space, const char trans[], typename AViewType::const_value_type& alpha,       \
                     const AViewType& A, const XViewType& X, typename YViewType::const_value_type& beta,              \
                     const YViewType& Y) {                                                                            \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_ROCBLAS,complex<double>]");                                 \
      KOKKOSBLAS2_GEMV_ROCBLAS_DETERMINE_ARGS(LAYOUT);                                                                \
      KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                        \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                            \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_zgemv(s.handle, transa, M, N,                                         \
                                                      reinterpret_cast<const rocblas_double_complex*>(&alpha),        \
                                                      reinterpret_cast<const rocblas_double_complex*>(A.data()), LDA, \
                                                      reinterpret_cast<const rocblas_double_complex*>(X.data()), one, \
                                                      reinterpret_cast<const rocblas_double_complex*>(&beta),         \
                                                      reinterpret_cast<rocblas_double_complex*>(Y.data()), one));     \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                          \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS2_CGEMV_ROCBLAS(LAYOUT)                                                                            \
  template <bool ETI_SPEC_AVAIL>                                                                                     \
  struct GEMV<                                                                                                       \
      Kokkos::HIP,                                                                                                   \
      Kokkos::View<const Kokkos::complex<float>**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,   \
      Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
      Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,    \
      ETI_SPEC_AVAIL> {                                                                                              \
    using SCALAR    = Kokkos::complex<float>;                                                                        \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;   \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using YViewType = Kokkos::View<SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;          \
                                                                                                                     \
    static void gemv(const Kokkos::HIP& space, const char trans[], typename AViewType::const_value_type& alpha,      \
                     const AViewType& A, const XViewType& X, typename YViewType::const_value_type& beta,             \
                     const YViewType& Y) {                                                                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::gemv[TPL_ROCBLAS,complex<float>]");                                 \
      KOKKOSBLAS2_GEMV_ROCBLAS_DETERMINE_ARGS(LAYOUT);                                                               \
      KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                       \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                           \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_cgemv(s.handle, transa, M, N,                                        \
                                                      reinterpret_cast<const rocblas_float_complex*>(&alpha),        \
                                                      reinterpret_cast<const rocblas_float_complex*>(A.data()), LDA, \
                                                      reinterpret_cast<const rocblas_float_complex*>(X.data()), one, \
                                                      reinterpret_cast<const rocblas_float_complex*>(&beta),         \
                                                      reinterpret_cast<rocblas_float_complex*>(Y.data()), one));     \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                         \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

KOKKOSBLAS2_DGEMV_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS2_DGEMV_ROCBLAS(Kokkos::LayoutRight)

KOKKOSBLAS2_SGEMV_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS2_SGEMV_ROCBLAS(Kokkos::LayoutRight)

KOKKOSBLAS2_ZGEMV_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS2_ZGEMV_ROCBLAS(Kokkos::LayoutRight)

KOKKOSBLAS2_CGEMV_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS2_CGEMV_ROCBLAS(Kokkos::LayoutRight)

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCBLAS

// ONEMKL
#if defined(KOKKOSKERNELS_ENABLE_TPL_MKL) && defined(KOKKOS_ENABLE_SYCL)
#include <mkl.h>
#include <oneapi/mkl/blas.hpp>
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS2_GEMV_ONEMKL(SCALAR, LAYOUT)                                                                         \
  template <bool ETI_SPEC_AVAIL>                                                                                        \
  struct GEMV<                                                                                                          \
      Kokkos::SYCL, Kokkos::View<const SCALAR**, LAYOUT, Kokkos::SYCL, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,       \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::SYCL, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::SYCL, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true, ETI_SPEC_AVAIL> {    \
    using mem_traits = Kokkos::MemoryTraits<Kokkos::Unmanaged>;                                                         \
    using AViewType  = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::SYCL, mem_traits>;                                  \
    using XViewType  = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::SYCL, mem_traits>;                                   \
    using YViewType  = Kokkos::View<SCALAR*, LAYOUT, Kokkos::SYCL, mem_traits>;                                         \
                                                                                                                        \
    static void gemv(const Kokkos::SYCL& exec, const char kk_trans[], typename AViewType::const_value_type& alpha,      \
                     const AViewType& A, const XViewType& X, typename YViewType::const_value_type& beta,                \
                     const YViewType& Y) {                                                                              \
      if (beta == KokkosKernels::ArithTraits<SCALAR>::zero()) {                                                         \
        Kokkos::deep_copy(Y, KokkosKernels::ArithTraits<SCALAR>::zero());                                               \
      }                                                                                                                 \
                                                                                                                        \
      bool row_major               = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;                                  \
      const std::int64_t M         = A.extent(0);                                                                       \
      const std::int64_t N         = A.extent(1);                                                                       \
      oneapi::mkl::transpose trans = mode_kk_to_onemkl(kk_trans[0]);                                                    \
      const std::int64_t LDA       = row_major ? A.stride(0) : A.stride(1);                                             \
      std::string label            = "KokkosBlas::gemv[TPL_ONEMKL," + KokkosKernels::ArithTraits<SCALAR>::name() + "]"; \
                                                                                                                        \
      Kokkos::Profiling::pushRegion(label);                                                                             \
      using mag_type    = kokkos_to_std_type_map<SCALAR, KokkosKernels::ArithTraits<SCALAR>::is_complex>::type;         \
      const mag_type* a = reinterpret_cast<const mag_type*>(A.data());                                                  \
      const mag_type* x = reinterpret_cast<const mag_type*>(X.data());                                                  \
      mag_type* y       = reinterpret_cast<mag_type*>(Y.data());                                                        \
      if (row_major) {                                                                                                  \
        oneapi::mkl::blas::row_major::gemv(exec.sycl_queue(), trans, M, N, alpha, a, LDA, x, 1, beta, y, 1);            \
      } else {                                                                                                          \
        oneapi::mkl::blas::column_major::gemv(exec.sycl_queue(), trans, M, N, alpha, a, LDA, x, 1, beta, y, 1);         \
      }                                                                                                                 \
      Kokkos::Profiling::popRegion();                                                                                   \
    }                                                                                                                   \
  };

KOKKOSBLAS2_GEMV_ONEMKL(float, Kokkos::LayoutLeft)
KOKKOSBLAS2_GEMV_ONEMKL(float, Kokkos::LayoutRight)
KOKKOSBLAS2_GEMV_ONEMKL(double, Kokkos::LayoutLeft)
KOKKOSBLAS2_GEMV_ONEMKL(double, Kokkos::LayoutRight)
KOKKOSBLAS2_GEMV_ONEMKL(Kokkos::complex<float>, Kokkos::LayoutLeft)
KOKKOSBLAS2_GEMV_ONEMKL(Kokkos::complex<float>, Kokkos::LayoutRight)
KOKKOSBLAS2_GEMV_ONEMKL(Kokkos::complex<double>, Kokkos::LayoutLeft)
KOKKOSBLAS2_GEMV_ONEMKL(Kokkos::complex<double>, Kokkos::LayoutRight)
}  // namespace Impl
}  // namespace KokkosBlas
#endif

#endif
