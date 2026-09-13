// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS2_GER_TPL_SPEC_DECL_CUBLAS_HPP_
#define KOKKOSBLAS2_GER_TPL_SPEC_DECL_CUBLAS_HPP_

#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS2_GER_CUBLAS_DETERMINE_ARGS(LAYOUT)                        \
  bool A_is_ll      = std::is_same<Kokkos::LayoutLeft, LAYOUT>::value;       \
  bool A_is_lr      = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;      \
  const int M       = static_cast<int>(A_is_lr ? A.extent(1) : A.extent(0)); \
  const int N       = static_cast<int>(A_is_lr ? A.extent(0) : A.extent(1)); \
  constexpr int one = 1;                                                     \
  const int LDA     = A_is_lr ? A.stride(0) : A.stride(1);

#define KOKKOSBLAS2_DGER_CUBLAS(LAYOUT, ETI_SPEC_AVAIL)                                                               \
  template <>                                                                                                         \
  struct GER<                                                                                                         \
      Kokkos::Cuda, Kokkos::View<const double*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,      \
      Kokkos::View<const double*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                    \
      Kokkos::View<double**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true, ETI_SPEC_AVAIL> { \
    using SCALAR    = double;                                                                                         \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;         \
                                                                                                                      \
    static void ger(const Kokkos::Cuda& space, const char /*trans*/[], typename AViewType::const_value_type& alpha,   \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_CUBLAS,double]");                                            \
      KOKKOSBLAS2_GER_CUBLAS_DETERMINE_ARGS(LAYOUT);                                                                  \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                      \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));                               \
      if (A_is_ll) {                                                                                                  \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(                                                                             \
            cublasDger(s.handle, M, N, &alpha, X.data(), one, Y.data(), one, A.data(), LDA));                         \
      } else {                                                                                                        \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(                                                                             \
            cublasDger(s.handle, M, N, &alpha, Y.data(), one, X.data(), one, A.data(), LDA));                         \
      }                                                                                                               \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, NULL));                                              \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS2_SGER_CUBLAS(LAYOUT, ETI_SPEC_AVAIL)                                                                \
  template <>                                                                                                          \
  struct GER<Kokkos::Cuda, Kokkos::View<const float*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
             Kokkos::View<const float*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
             Kokkos::View<float**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,              \
             ETI_SPEC_AVAIL> {                                                                                         \
    using SCALAR    = float;                                                                                           \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;     \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;     \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;          \
                                                                                                                       \
    static void ger(const Kokkos::Cuda& space, const char /*trans*/[], typename AViewType::const_value_type& alpha,    \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                      \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_CUBLAS,float]");                                              \
      KOKKOSBLAS2_GER_CUBLAS_DETERMINE_ARGS(LAYOUT);                                                                   \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                       \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));                                \
      if (A_is_ll) {                                                                                                   \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(                                                                              \
            cublasSger(s.handle, M, N, &alpha, X.data(), one, Y.data(), one, A.data(), LDA));                          \
      } else {                                                                                                         \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(                                                                              \
            cublasSger(s.handle, M, N, &alpha, Y.data(), one, X.data(), one, A.data(), LDA));                          \
      }                                                                                                                \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, NULL));                                               \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS2_ZGER_CUBLAS(LAYOUT, ETI_SPEC_AVAIL)                                                              \
  template <>                                                                                                        \
  struct GER<                                                                                                        \
      Kokkos::Cuda,                                                                                                  \
      Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,  \
      Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,  \
      Kokkos::View<Kokkos::complex<double>**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true, \
      ETI_SPEC_AVAIL> {                                                                                              \
    using SCALAR    = Kokkos::complex<double>;                                                                       \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;   \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;   \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;        \
                                                                                                                     \
    static void ger(const Kokkos::Cuda& space, const char trans[], typename AViewType::const_value_type& alpha,      \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                    \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_CUBLAS,complex<double>]");                                  \
      KOKKOSBLAS2_GER_CUBLAS_DETERMINE_ARGS(LAYOUT);                                                                 \
      bool justTranspose                     = (trans[0] == 'T') || (trans[0] == 't');                               \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                     \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));                              \
      if (A_is_ll) {                                                                                                 \
        if (justTranspose) {                                                                                         \
          KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasZgeru(s.handle, M, N,                                               \
                                                       reinterpret_cast<const cuDoubleComplex*>(&alpha),             \
                                                       reinterpret_cast<const cuDoubleComplex*>(X.data()), one,      \
                                                       reinterpret_cast<const cuDoubleComplex*>(Y.data()), one,      \
                                                       reinterpret_cast<cuDoubleComplex*>(A.data()), LDA));          \
        } else {                                                                                                     \
          KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasZgerc(s.handle, M, N,                                               \
                                                       reinterpret_cast<const cuDoubleComplex*>(&alpha),             \
                                                       reinterpret_cast<const cuDoubleComplex*>(X.data()), one,      \
                                                       reinterpret_cast<const cuDoubleComplex*>(Y.data()), one,      \
                                                       reinterpret_cast<cuDoubleComplex*>(A.data()), LDA));          \
        }                                                                                                            \
      } else {                                                                                                       \
        if (justTranspose) {                                                                                         \
          KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasZgeru(s.handle, M, N,                                               \
                                                       reinterpret_cast<const cuDoubleComplex*>(&alpha),             \
                                                       reinterpret_cast<const cuDoubleComplex*>(Y.data()), one,      \
                                                       reinterpret_cast<const cuDoubleComplex*>(X.data()), one,      \
                                                       reinterpret_cast<cuDoubleComplex*>(A.data()), LDA));          \
        } else {                                                                                                     \
          /* cublasZgerc() + ~A_ll => call kokkos-kernels' implementation */                                         \
          GER<Kokkos::Cuda, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::ger(space, trans, alpha, X, Y,  \
                                                                                         A);                         \
        }                                                                                                            \
      }                                                                                                              \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, NULL));                                             \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS2_CGER_CUBLAS(LAYOUT, ETI_SPEC_AVAIL)                                                             \
  template <>                                                                                                       \
  struct GER<                                                                                                       \
      Kokkos::Cuda,                                                                                                 \
      Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,  \
      Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,  \
      Kokkos::View<Kokkos::complex<float>**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true, \
      ETI_SPEC_AVAIL> {                                                                                             \
    using SCALAR    = Kokkos::complex<float>;                                                                       \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;  \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;  \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;       \
                                                                                                                    \
    static void ger(const Kokkos::Cuda& space, const char trans[], typename AViewType::const_value_type& alpha,     \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                   \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_CUBLAS,complex<float>]");                                  \
      KOKKOSBLAS2_GER_CUBLAS_DETERMINE_ARGS(LAYOUT);                                                                \
      bool justTranspose                     = (trans[0] == 'T') || (trans[0] == 't');                              \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                    \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));                             \
      if (A_is_ll) {                                                                                                \
        if (justTranspose) {                                                                                        \
          KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasCgeru(s.handle, M, N, reinterpret_cast<const cuComplex*>(&alpha),  \
                                                       reinterpret_cast<const cuComplex*>(X.data()), one,           \
                                                       reinterpret_cast<const cuComplex*>(Y.data()), one,           \
                                                       reinterpret_cast<cuComplex*>(A.data()), LDA));               \
        } else {                                                                                                    \
          KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasCgerc(s.handle, M, N, reinterpret_cast<const cuComplex*>(&alpha),  \
                                                       reinterpret_cast<const cuComplex*>(X.data()), one,           \
                                                       reinterpret_cast<const cuComplex*>(Y.data()), one,           \
                                                       reinterpret_cast<cuComplex*>(A.data()), LDA));               \
        }                                                                                                           \
      } else {                                                                                                      \
        if (justTranspose) {                                                                                        \
          KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasCgeru(s.handle, M, N, reinterpret_cast<const cuComplex*>(&alpha),  \
                                                       reinterpret_cast<const cuComplex*>(Y.data()), one,           \
                                                       reinterpret_cast<const cuComplex*>(X.data()), one,           \
                                                       reinterpret_cast<cuComplex*>(A.data()), LDA));               \
        } else {                                                                                                    \
          /* cublasCgerc() + ~A_ll => call kokkos-kernels' implementation */                                        \
          GER<Kokkos::Cuda, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::ger(space, trans, alpha, X, Y, \
                                                                                         A);                        \
        }                                                                                                           \
      }                                                                                                             \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, NULL));                                            \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

KOKKOSBLAS2_DGER_CUBLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_DGER_CUBLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_DGER_CUBLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_DGER_CUBLAS(Kokkos::LayoutRight, false)

KOKKOSBLAS2_SGER_CUBLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_SGER_CUBLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_SGER_CUBLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_SGER_CUBLAS(Kokkos::LayoutRight, false)

KOKKOSBLAS2_ZGER_CUBLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_ZGER_CUBLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_ZGER_CUBLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_ZGER_CUBLAS(Kokkos::LayoutRight, false)

KOKKOSBLAS2_CGER_CUBLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_CGER_CUBLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_CGER_CUBLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_CGER_CUBLAS(Kokkos::LayoutRight, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif
