// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS2_GER_TPL_SPEC_DECL_BLAS_HPP_
#define KOKKOSBLAS2_GER_TPL_SPEC_DECL_BLAS_HPP_

#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS2_GER_DETERMINE_ARGS(LAYOUT)                               \
  bool A_is_ll      = std::is_same<Kokkos::LayoutLeft, LAYOUT>::value;       \
  bool A_is_lr      = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;      \
  const int M       = static_cast<int>(A_is_lr ? A.extent(1) : A.extent(0)); \
  const int N       = static_cast<int>(A_is_lr ? A.extent(0) : A.extent(1)); \
  constexpr int one = 1;                                                     \
  const int LDA     = A_is_lr ? A.stride(0) : A.stride(1);

#define KOKKOSBLAS2_DGER_BLAS(LAYOUT, ETI_SPEC_AVAIL)                                                                  \
  template <typename ExecSpace>                                                                                        \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                      \
  struct GER<ExecSpace, Kokkos::View<const double*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,       \
             Kokkos::View<const double*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                  \
             Kokkos::View<double**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,                 \
             ETI_SPEC_AVAIL> {                                                                                         \
    using SCALAR    = double;                                                                                          \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;         \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;         \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;              \
                                                                                                                       \
    static void ger(const ExecSpace& /* space */, const char /*trans*/[], typename AViewType::const_value_type& alpha, \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                      \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_BLAS,double]");                                               \
      KOKKOSBLAS2_GER_DETERMINE_ARGS(LAYOUT);                                                                          \
      if (A_is_ll) {                                                                                                   \
        HostBlas<SCALAR>::ger(M, N, alpha, X.data(), one, Y.data(), one, A.data(), LDA);                               \
      } else {                                                                                                         \
        HostBlas<SCALAR>::ger(M, N, alpha, Y.data(), one, X.data(), one, A.data(), LDA);                               \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS2_SGER_BLAS(LAYOUT, ETI_SPEC_AVAIL)                                                           \
  template <typename ExecSpace>                                                                                 \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                               \
  struct GER<ExecSpace, Kokkos::View<const float*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
             Kokkos::View<const float*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,            \
             Kokkos::View<float**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,           \
             ETI_SPEC_AVAIL> {                                                                                  \
    using SCALAR    = float;                                                                                    \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;  \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;  \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
                                                                                                                \
    static void ger(const ExecSpace& /* space */                                                                \
                    ,                                                                                           \
                    const char /*trans*/[], typename AViewType::const_value_type& alpha, const XViewType& X,    \
                    const YViewType& Y, const AViewType& A) {                                                   \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_BLAS,float]");                                         \
      KOKKOSBLAS2_GER_DETERMINE_ARGS(LAYOUT);                                                                   \
      if (A_is_ll) {                                                                                            \
        HostBlas<SCALAR>::ger(M, N, alpha, X.data(), one, Y.data(), one, A.data(), LDA);                        \
      } else {                                                                                                  \
        HostBlas<SCALAR>::ger(M, N, alpha, Y.data(), one, X.data(), one, A.data(), LDA);                        \
      }                                                                                                         \
      Kokkos::Profiling::popRegion();                                                                           \
    }                                                                                                           \
  };

#define KOKKOSBLAS2_ZGER_BLAS(LAYOUT, ETI_SPEC_AVAIL)                                                                  \
  template <typename ExecSpace>                                                                                        \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                      \
  struct GER<ExecSpace,                                                                                                \
             Kokkos::View<const Kokkos::complex<double>*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
             Kokkos::View<const Kokkos::complex<double>*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
             Kokkos::View<Kokkos::complex<double>**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,      \
             true, ETI_SPEC_AVAIL> {                                                                                   \
    using SCALAR    = Kokkos::complex<double>;                                                                         \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;         \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;         \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;              \
                                                                                                                       \
    static void ger(const ExecSpace& space, const char trans[], typename AViewType::const_value_type& alpha,           \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                      \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_BLAS,complex<double>]");                                      \
      KOKKOSBLAS2_GER_DETERMINE_ARGS(LAYOUT);                                                                          \
      const std::complex<double> alpha_val = static_cast<const std::complex<double>>(alpha);                           \
      bool justTranspose                   = (trans[0] == 'T') || (trans[0] == 't');                                   \
      if (A_is_ll) {                                                                                                   \
        if (justTranspose) {                                                                                           \
          HostBlas<std::complex<double>>::geru(M, N, alpha_val,                                                        \
                                               reinterpret_cast<const std::complex<double>*>(X.data()), one,           \
                                               reinterpret_cast<const std::complex<double>*>(Y.data()), one,           \
                                               reinterpret_cast<std::complex<double>*>(A.data()), LDA);                \
        } else {                                                                                                       \
          HostBlas<std::complex<double>>::gerc(M, N, alpha_val,                                                        \
                                               reinterpret_cast<const std::complex<double>*>(X.data()), one,           \
                                               reinterpret_cast<const std::complex<double>*>(Y.data()), one,           \
                                               reinterpret_cast<std::complex<double>*>(A.data()), LDA);                \
        }                                                                                                              \
      } else {                                                                                                         \
        if (justTranspose) {                                                                                           \
          HostBlas<std::complex<double>>::geru(M, N, alpha_val,                                                        \
                                               reinterpret_cast<const std::complex<double>*>(Y.data()), one,           \
                                               reinterpret_cast<const std::complex<double>*>(X.data()), one,           \
                                               reinterpret_cast<std::complex<double>*>(A.data()), LDA);                \
        } else {                                                                                                       \
          /* blasgerc() + ~A_ll => call kokkos-kernels' implementation */                                              \
          GER<ExecSpace, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::ger(space, trans, alpha, X, Y, A);   \
        }                                                                                                              \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS2_CGER_BLAS(LAYOUT, ETI_SPEC_AVAIL)                                                                  \
  template <typename ExecSpace>                                                                                        \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                      \
  struct GER<ExecSpace,                                                                                                \
             Kokkos::View<const Kokkos::complex<float>*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,  \
             Kokkos::View<const Kokkos::complex<float>*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,  \
             Kokkos::View<Kokkos::complex<float>**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, \
             ETI_SPEC_AVAIL> {                                                                                         \
    using SCALAR    = Kokkos::complex<float>;                                                                          \
    using XViewType = Kokkos::View<const SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;         \
    using YViewType = Kokkos::View<const SCALAR*, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;         \
    using AViewType = Kokkos::View<SCALAR**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;              \
                                                                                                                       \
    static void ger(const ExecSpace& space, const char trans[], typename AViewType::const_value_type& alpha,           \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                      \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_BLAS,complex<float>]");                                       \
      KOKKOSBLAS2_GER_DETERMINE_ARGS(LAYOUT);                                                                          \
      const std::complex<float> alpha_val = static_cast<const std::complex<float>>(alpha);                             \
      bool justTranspose                  = (trans[0] == 'T') || (trans[0] == 't');                                    \
      if (A_is_ll) {                                                                                                   \
        if (justTranspose) {                                                                                           \
          HostBlas<std::complex<float>>::geru(M, N, alpha_val, reinterpret_cast<const std::complex<float>*>(X.data()), \
                                              one, reinterpret_cast<const std::complex<float>*>(Y.data()), one,        \
                                              reinterpret_cast<std::complex<float>*>(A.data()), LDA);                  \
        } else {                                                                                                       \
          HostBlas<std::complex<float>>::gerc(M, N, alpha_val, reinterpret_cast<const std::complex<float>*>(X.data()), \
                                              one, reinterpret_cast<const std::complex<float>*>(Y.data()), one,        \
                                              reinterpret_cast<std::complex<float>*>(A.data()), LDA);                  \
        }                                                                                                              \
      } else {                                                                                                         \
        if (justTranspose) {                                                                                           \
          HostBlas<std::complex<float>>::geru(M, N, alpha_val, reinterpret_cast<const std::complex<float>*>(Y.data()), \
                                              one, reinterpret_cast<const std::complex<float>*>(X.data()), one,        \
                                              reinterpret_cast<std::complex<float>*>(A.data()), LDA);                  \
        } else {                                                                                                       \
          /* blasgerc() + ~A_ll => call kokkos-kernels' implementation */                                              \
          GER<ExecSpace, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::ger(space, trans, alpha, X, Y, A);   \
        }                                                                                                              \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS2_DGER_BLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_DGER_BLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_DGER_BLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_DGER_BLAS(Kokkos::LayoutRight, false)

KOKKOSBLAS2_SGER_BLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_SGER_BLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_SGER_BLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_SGER_BLAS(Kokkos::LayoutRight, false)

KOKKOSBLAS2_ZGER_BLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_ZGER_BLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_ZGER_BLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_ZGER_BLAS(Kokkos::LayoutRight, false)

KOKKOSBLAS2_CGER_BLAS(Kokkos::LayoutLeft, true)
KOKKOSBLAS2_CGER_BLAS(Kokkos::LayoutLeft, false)
KOKKOSBLAS2_CGER_BLAS(Kokkos::LayoutRight, true)
KOKKOSBLAS2_CGER_BLAS(Kokkos::LayoutRight, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif
