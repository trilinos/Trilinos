// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS2_GER_TPL_SPEC_DECL_ROCBLAS_HPP_
#define KOKKOSBLAS2_GER_TPL_SPEC_DECL_ROCBLAS_HPP_

#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS2_GER_ROCBLAS_DETERMINE_ARGS(LAYOUT)                       \
  bool A_is_ll      = std::is_same<Kokkos::LayoutLeft, LAYOUT>::value;       \
  bool A_is_lr      = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;      \
  const int M       = static_cast<int>(A_is_lr ? A.extent(1) : A.extent(0)); \
  const int N       = static_cast<int>(A_is_lr ? A.extent(0) : A.extent(1)); \
  constexpr int one = 1;                                                     \
  const int LDA     = A_is_lr ? A.stride(0) : A.stride(1);

#define KOKKOSBLAS2_DGER_ROCBLAS(LAYOUT, ETI_SPEC_AVAIL)                                                              \
  template <>                                                                                                         \
  struct GER<Kokkos::HIP, Kokkos::View<const double*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
             Kokkos::View<const double*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,              \
             Kokkos::View<double**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,             \
             ETI_SPEC_AVAIL> {                                                                                        \
    using SCALAR    = double;                                                                                         \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;     \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;     \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;          \
                                                                                                                      \
    static void ger(const Kokkos::HIP& space, const char /*trans*/[], typename AViewType::const_value_type& alpha,    \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_ROCBLAS,double]");                                           \
      KOKKOSBLAS2_GER_ROCBLAS_DETERMINE_ARGS(LAYOUT);                                                                 \
      KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                        \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                            \
      if (A_is_ll) {                                                                                                  \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                            \
            rocblas_dger(s.handle, M, N, &alpha, X.data(), one, Y.data(), one, A.data(), LDA));                       \
      } else {                                                                                                        \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                            \
            rocblas_dger(s.handle, M, N, &alpha, Y.data(), one, X.data(), one, A.data(), LDA));                       \
      }                                                                                                               \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                          \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS2_SGER_ROCBLAS(LAYOUT, ETI_SPEC_AVAIL)                                                             \
  template <>                                                                                                        \
  struct GER<Kokkos::HIP, Kokkos::View<const float*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
             Kokkos::View<const float*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,              \
             Kokkos::View<float**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,             \
             ETI_SPEC_AVAIL> {                                                                                       \
    using SCALAR    = float;                                                                                         \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;         \
                                                                                                                     \
    static void ger(const Kokkos::HIP& space, const char /*trans*/[], typename AViewType::const_value_type& alpha,   \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                    \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_ROCBLAS,float]");                                           \
      KOKKOSBLAS2_GER_ROCBLAS_DETERMINE_ARGS(LAYOUT);                                                                \
      KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                       \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                           \
      if (A_is_ll) {                                                                                                 \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                           \
            rocblas_sger(s.handle, M, N, &alpha, X.data(), one, Y.data(), one, A.data(), LDA));                      \
      } else {                                                                                                       \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                           \
            rocblas_sger(s.handle, M, N, &alpha, Y.data(), one, X.data(), one, A.data(), LDA));                      \
      }                                                                                                              \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                         \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS2_ZGER_ROCBLAS(LAYOUT, ETI_SPEC_AVAIL)                                                               \
  template <>                                                                                                          \
  struct GER<                                                                                                          \
      Kokkos::HIP,                                                                                                     \
      Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,     \
      Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,     \
      Kokkos::View<Kokkos::complex<double>**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,    \
      ETI_SPEC_AVAIL> {                                                                                                \
    using SCALAR    = Kokkos::complex<double>;                                                                         \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;      \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;      \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;           \
                                                                                                                       \
    static void ger(const Kokkos::HIP& space, const char trans[], typename AViewType::const_value_type& alpha,         \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                      \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_ROCBLAS,complex<double>]");                                   \
      KOKKOSBLAS2_GER_ROCBLAS_DETERMINE_ARGS(LAYOUT);                                                                  \
      bool justTranspose                    = (trans[0] == 'T') || (trans[0] == 't');                                  \
      KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                         \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                             \
      if (A_is_ll) {                                                                                                   \
        if (justTranspose) {                                                                                           \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                           \
              rocblas_zgeru(s.handle, M, N, reinterpret_cast<const rocblas_double_complex*>(&alpha),                   \
                            reinterpret_cast<const rocblas_double_complex*>(X.data()), one,                            \
                            reinterpret_cast<const rocblas_double_complex*>(Y.data()), one,                            \
                            reinterpret_cast<rocblas_double_complex*>(A.data()), LDA));                                \
        } else {                                                                                                       \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                           \
              rocblas_zgerc(s.handle, M, N, reinterpret_cast<const rocblas_double_complex*>(&alpha),                   \
                            reinterpret_cast<const rocblas_double_complex*>(X.data()), one,                            \
                            reinterpret_cast<const rocblas_double_complex*>(Y.data()), one,                            \
                            reinterpret_cast<rocblas_double_complex*>(A.data()), LDA));                                \
        }                                                                                                              \
      } else {                                                                                                         \
        if (justTranspose) {                                                                                           \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                           \
              rocblas_zgeru(s.handle, M, N, reinterpret_cast<const rocblas_double_complex*>(&alpha),                   \
                            reinterpret_cast<const rocblas_double_complex*>(Y.data()), one,                            \
                            reinterpret_cast<const rocblas_double_complex*>(X.data()), one,                            \
                            reinterpret_cast<rocblas_double_complex*>(A.data()), LDA));                                \
        } else {                                                                                                       \
          /* rocblas_zgerc() + ~A_ll => call k-kernels' implementation */                                              \
          GER<Kokkos::HIP, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::ger(space, trans, alpha, X, Y, A); \
        }                                                                                                              \
      }                                                                                                                \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                           \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS2_CGER_ROCBLAS(LAYOUT, ETI_SPEC_AVAIL)                                                               \
  template <>                                                                                                          \
  struct GER<                                                                                                          \
      Kokkos::HIP,                                                                                                     \
      Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,      \
      Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,      \
      Kokkos::View<Kokkos::complex<float>**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,     \
      ETI_SPEC_AVAIL> {                                                                                                \
    using SCALAR    = Kokkos::complex<float>;                                                                          \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;      \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;      \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;           \
                                                                                                                       \
    static void ger(const Kokkos::HIP& space, const char trans[], typename AViewType::const_value_type& alpha,         \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                      \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_ROCBLAS,complex<float>]");                                    \
      KOKKOSBLAS2_GER_ROCBLAS_DETERMINE_ARGS(LAYOUT);                                                                  \
      bool justTranspose                    = (trans[0] == 'T') || (trans[0] == 't');                                  \
      KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                         \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                             \
      if (A_is_ll) {                                                                                                   \
        if (justTranspose) {                                                                                           \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                           \
              rocblas_cgeru(s.handle, M, N, reinterpret_cast<const rocblas_float_complex*>(&alpha),                    \
                            reinterpret_cast<const rocblas_float_complex*>(X.data()), one,                             \
                            reinterpret_cast<const rocblas_float_complex*>(Y.data()), one,                             \
                            reinterpret_cast<rocblas_float_complex*>(A.data()), LDA));                                 \
        } else {                                                                                                       \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                           \
              rocblas_cgerc(s.handle, M, N, reinterpret_cast<const rocblas_float_complex*>(&alpha),                    \
                            reinterpret_cast<const rocblas_float_complex*>(X.data()), one,                             \
                            reinterpret_cast<const rocblas_float_complex*>(Y.data()), one,                             \
                            reinterpret_cast<rocblas_float_complex*>(A.data()), LDA));                                 \
        }                                                                                                              \
      } else {                                                                                                         \
        if (justTranspose) {                                                                                           \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                           \
              rocblas_cgeru(s.handle, M, N, reinterpret_cast<const rocblas_float_complex*>(&alpha),                    \
                            reinterpret_cast<const rocblas_float_complex*>(Y.data()), one,                             \
                            reinterpret_cast<const rocblas_float_complex*>(X.data()), one,                             \
                            reinterpret_cast<rocblas_float_complex*>(A.data()), LDA));                                 \
        } else {                                                                                                       \
          /* rocblas_cgerc() + ~A_ll => call k-kernels' implementation */                                              \
          GER<Kokkos::HIP, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::ger(space, trans, alpha, X, Y, A); \
        }                                                                                                              \
      }                                                                                                                \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                           \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS2_DGER_ROCBLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_DGER_ROCBLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_DGER_ROCBLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_DGER_ROCBLAS(Kokkos::LayoutRight, false)

KOKKOSBLAS2_SGER_ROCBLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_SGER_ROCBLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_SGER_ROCBLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_SGER_ROCBLAS(Kokkos::LayoutRight, false)

KOKKOSBLAS2_ZGER_ROCBLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_ZGER_ROCBLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_ZGER_ROCBLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_ZGER_ROCBLAS(Kokkos::LayoutRight, false)

KOKKOSBLAS2_CGER_ROCBLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_CGER_ROCBLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_CGER_ROCBLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_CGER_ROCBLAS(Kokkos::LayoutRight, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif
