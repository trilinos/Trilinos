// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS2_SYR2_TPL_SPEC_DECL_ROCBLAS_HPP_
#define KOKKOSBLAS2_SYR2_TPL_SPEC_DECL_ROCBLAS_HPP_

#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS2_SYR2_ROCBLAS_DETERMINE_ARGS(LAYOUT, uploChar)                \
  bool A_is_ll          = std::is_same<Kokkos::LayoutLeft, LAYOUT>::value;       \
  bool A_is_lr          = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;      \
  const int N           = static_cast<int>(A_is_lr ? A.extent(0) : A.extent(1)); \
  constexpr int one     = 1;                                                     \
  const int LDA         = A_is_lr ? A.stride(0) : A.stride(1);                   \
  rocblas_fill fillMode = (uploChar == 'L' || uploChar == 'l') ? rocblas_fill_lower : rocblas_fill_upper;

#define KOKKOSBLAS2_DSYR2_ROCBLAS(LAYOUT, ETI_SPEC_AVAIL)                                                              \
  template <>                                                                                                          \
  struct SYR2<Kokkos::HIP, Kokkos::View<const double*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
              Kokkos::View<const double*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,              \
              Kokkos::View<double**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,             \
              ETI_SPEC_AVAIL> {                                                                                        \
    using SCALAR    = double;                                                                                          \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;      \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;      \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;           \
                                                                                                                       \
    static void syr2(const typename AViewType::execution_space& space, const char trans[], const char uplo[],          \
                     typename AViewType::const_value_type& alpha, const XViewType& X, const YViewType& Y,              \
                     const AViewType& A) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr2[TPL_ROCBLAS,double]");                                           \
      KOKKOSBLAS2_SYR2_ROCBLAS_DETERMINE_ARGS(LAYOUT, uplo[0]);                                                        \
      if (A_is_ll) {                                                                                                   \
        KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                       \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                           \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                             \
            rocblas_dsyr2(s.handle, fillMode, N, &alpha, X.data(), one, Y.data(), one, A.data(), LDA));                \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                         \
      } else {                                                                                                         \
        /* rocblas_dsyr2() + ~A_ll => call kokkos-kernels' implementation */                                           \
        SYR2<Kokkos::HIP, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, X,  \
                                                                                        Y, A);                         \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS2_SSYR2_ROCBLAS(LAYOUT, ETI_SPEC_AVAIL)                                                             \
  template <>                                                                                                         \
  struct SYR2<Kokkos::HIP, Kokkos::View<const float*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
              Kokkos::View<const float*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,              \
              Kokkos::View<float**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,             \
              ETI_SPEC_AVAIL> {                                                                                       \
    using SCALAR    = float;                                                                                          \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;     \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;     \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;          \
                                                                                                                      \
    static void syr2(const typename AViewType::execution_space& space, const char trans[], const char uplo[],         \
                     typename AViewType::const_value_type& alpha, const XViewType& X, const YViewType& Y,             \
                     const AViewType& A) {                                                                            \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr2[TPL_ROCBLAS,float]");                                           \
      KOKKOSBLAS2_SYR2_ROCBLAS_DETERMINE_ARGS(LAYOUT, uplo[0]);                                                       \
      if (A_is_ll) {                                                                                                  \
        KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                      \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                          \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                            \
            rocblas_ssyr2(s.handle, fillMode, N, &alpha, X.data(), one, Y.data(), one, A.data(), LDA));               \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                        \
      } else {                                                                                                        \
        /* rocblas_ssyr2() + ~A_ll => call kokkos-kernels' implementation */                                          \
        SYR2<Kokkos::HIP, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, X, \
                                                                                        Y, A);                        \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS2_ZSYR2_ROCBLAS(LAYOUT, ETI_SPEC_AVAIL)                                                            \
  template <>                                                                                                        \
  struct SYR2<                                                                                                       \
      Kokkos::HIP,                                                                                                   \
      Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,   \
      Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,   \
      Kokkos::View<Kokkos::complex<double>**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,  \
      ETI_SPEC_AVAIL> {                                                                                              \
    using SCALAR    = Kokkos::complex<double>;                                                                       \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;         \
                                                                                                                     \
    static void syr2(const typename AViewType::execution_space& space, const char trans[], const char uplo[],        \
                     typename AViewType::const_value_type& alpha, const XViewType& X, const YViewType& Y,            \
                     const AViewType& A) {                                                                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr2[TPL_ROCBLAS,complex<double>]");                                \
      KOKKOSBLAS2_SYR2_ROCBLAS_DETERMINE_ARGS(LAYOUT, uplo[0]);                                                      \
      bool justTranspose = (trans[0] == 'T') || (trans[0] == 't');                                                   \
      if (justTranspose) {                                                                                           \
        if (A_is_ll) {                                                                                               \
          KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                   \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                       \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                         \
              rocblas_zsyr2(s.handle, fillMode, N, reinterpret_cast<const rocblas_double_complex*>(&alpha),          \
                            reinterpret_cast<const rocblas_double_complex*>(X.data()), one,                          \
                            reinterpret_cast<const rocblas_double_complex*>(Y.data()), one,                          \
                            reinterpret_cast<rocblas_double_complex*>(A.data()), LDA));                              \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                     \
        } else {                                                                                                     \
          /* rocblas_zsyr2() + ~A_ll => call kokkos-kernels' implementation */                                       \
          SYR2<Kokkos::HIP, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, \
                                                                                          X, Y, A);                  \
        }                                                                                                            \
      } else {                                                                                                       \
        if (A_is_ll && (alpha.imag() == 0.)) {                                                                       \
          KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                   \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                       \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                         \
              rocblas_zher2(s.handle, fillMode, N, reinterpret_cast<const rocblas_double_complex*>(&alpha),          \
                            reinterpret_cast<const rocblas_double_complex*>(X.data()), one,                          \
                            reinterpret_cast<const rocblas_double_complex*>(Y.data()), one,                          \
                            reinterpret_cast<rocblas_double_complex*>(A.data()), LDA));                              \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                     \
        } else {                                                                                                     \
          /* rocblas_zher2() + ~A_ll => call kokkos-kernels' implementation */                                       \
          SYR2<Kokkos::HIP, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, \
                                                                                          X, Y, A);                  \
        }                                                                                                            \
      }                                                                                                              \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS2_CSYR2_ROCBLAS(LAYOUT, ETI_SPEC_AVAIL)                                                            \
  template <>                                                                                                        \
  struct SYR2<                                                                                                       \
      Kokkos::HIP,                                                                                                   \
      Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
      Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
      Kokkos::View<Kokkos::complex<float>**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, true,   \
      ETI_SPEC_AVAIL> {                                                                                              \
    using SCALAR    = Kokkos::complex<float>;                                                                        \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;    \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;         \
                                                                                                                     \
    static void syr2(const typename AViewType::execution_space& space, const char trans[], const char uplo[],        \
                     typename AViewType::const_value_type& alpha, const XViewType& X, const YViewType& Y,            \
                     const AViewType& A) {                                                                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr2[TPL_ROCBLAS,complex<float>]");                                 \
      KOKKOSBLAS2_SYR2_ROCBLAS_DETERMINE_ARGS(LAYOUT, uplo[0]);                                                      \
      bool justTranspose = (trans[0] == 'T') || (trans[0] == 't');                                                   \
      if (justTranspose) {                                                                                           \
        if (A_is_ll) {                                                                                               \
          KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                   \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                       \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                         \
              rocblas_csyr2(s.handle, fillMode, N, reinterpret_cast<const rocblas_float_complex*>(&alpha),           \
                            reinterpret_cast<const rocblas_float_complex*>(X.data()), one,                           \
                            reinterpret_cast<const rocblas_float_complex*>(Y.data()), one,                           \
                            reinterpret_cast<rocblas_float_complex*>(A.data()), LDA));                               \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                     \
        } else {                                                                                                     \
          /* rocblas_csyr2() + ~A_ll => call kokkos-kernels' implementation */                                       \
          SYR2<Kokkos::HIP, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, \
                                                                                          X, Y, A);                  \
        }                                                                                                            \
      } else {                                                                                                       \
        if (A_is_ll && (alpha.imag() == 0.)) {                                                                       \
          KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                   \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                       \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                         \
              rocblas_cher2(s.handle, fillMode, N, reinterpret_cast<const rocblas_float_complex*>(&alpha),           \
                            reinterpret_cast<const rocblas_float_complex*>(X.data()), one,                           \
                            reinterpret_cast<const rocblas_float_complex*>(Y.data()), one,                           \
                            reinterpret_cast<rocblas_float_complex*>(A.data()), LDA));                               \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                     \
        } else {                                                                                                     \
          /* rocblas_cher2() + ~A_ll => call kokkos-kernels' implementation */                                       \
          SYR2<Kokkos::HIP, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, \
                                                                                          X, Y, A);                  \
        }                                                                                                            \
      }                                                                                                              \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

KOKKOSBLAS2_DSYR2_ROCBLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_DSYR2_ROCBLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_DSYR2_ROCBLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_DSYR2_ROCBLAS(Kokkos::LayoutRight, false)

KOKKOSBLAS2_SSYR2_ROCBLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_SSYR2_ROCBLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_SSYR2_ROCBLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_SSYR2_ROCBLAS(Kokkos::LayoutRight, false)

KOKKOSBLAS2_ZSYR2_ROCBLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_ZSYR2_ROCBLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_ZSYR2_ROCBLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_ZSYR2_ROCBLAS(Kokkos::LayoutRight, false)

KOKKOSBLAS2_CSYR2_ROCBLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_CSYR2_ROCBLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_CSYR2_ROCBLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_CSYR2_ROCBLAS(Kokkos::LayoutRight, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif
