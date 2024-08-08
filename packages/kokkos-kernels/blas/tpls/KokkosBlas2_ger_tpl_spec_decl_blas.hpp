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

#define KOKKOSBLAS2_DGER_BLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                          \
  template <>                                                                                                         \
  struct GER<                                                                                                         \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<const double*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      Kokkos::View<const double*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      Kokkos::View<double**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    typedef double SCALAR;                                                                                            \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                     \
        XViewType;                                                                                                    \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                     \
        YViewType;                                                                                                    \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                     \
        AViewType;                                                                                                    \
                                                                                                                      \
    static void ger(const EXEC_SPACE& /* space */                                                                     \
                    ,                                                                                                 \
                    const char /*trans*/[], typename AViewType::const_value_type& alpha, const XViewType& X,          \
                    const YViewType& Y, const AViewType& A) {                                                         \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_BLAS,double]");                                              \
      KOKKOSBLAS2_GER_DETERMINE_ARGS(LAYOUT);                                                                         \
      if (A_is_ll) {                                                                                                  \
        HostBlas<SCALAR>::ger(M, N, alpha, X.data(), one, Y.data(), one, A.data(), LDA);                              \
      } else {                                                                                                        \
        HostBlas<SCALAR>::ger(M, N, alpha, Y.data(), one, X.data(), one, A.data(), LDA);                              \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS2_SGER_BLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                         \
  template <>                                                                                                        \
  struct GER<                                                                                                        \
      EXEC_SPACE,                                                                                                    \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                         \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                         \
      Kokkos::View<float**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      true, ETI_SPEC_AVAIL> {                                                                                        \
    typedef float SCALAR;                                                                                            \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                    \
        XViewType;                                                                                                   \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                    \
        YViewType;                                                                                                   \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                    \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                    \
        AViewType;                                                                                                   \
                                                                                                                     \
    static void ger(const EXEC_SPACE& /* space */                                                                    \
                    ,                                                                                                \
                    const char /*trans*/[], typename AViewType::const_value_type& alpha, const XViewType& X,         \
                    const YViewType& Y, const AViewType& A) {                                                        \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_BLAS,float]");                                              \
      KOKKOSBLAS2_GER_DETERMINE_ARGS(LAYOUT);                                                                        \
      if (A_is_ll) {                                                                                                 \
        HostBlas<SCALAR>::ger(M, N, alpha, X.data(), one, Y.data(), one, A.data(), LDA);                             \
      } else {                                                                                                       \
        HostBlas<SCALAR>::ger(M, N, alpha, Y.data(), one, X.data(), one, A.data(), LDA);                             \
      }                                                                                                              \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS2_ZGER_BLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                          \
  template <>                                                                                                         \
  struct GER<EXEC_SPACE,                                                                                              \
             Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,              \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                   \
             Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,              \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                   \
             Kokkos::View<Kokkos::complex<double>**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                   \
             true, ETI_SPEC_AVAIL> {                                                                                  \
    typedef Kokkos::complex<double> SCALAR;                                                                           \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                     \
        XViewType;                                                                                                    \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                     \
        YViewType;                                                                                                    \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                     \
        AViewType;                                                                                                    \
                                                                                                                      \
    static void ger(const EXEC_SPACE& space, const char trans[], typename AViewType::const_value_type& alpha,         \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_BLAS,complex<double>");                                      \
      KOKKOSBLAS2_GER_DETERMINE_ARGS(LAYOUT);                                                                         \
      const std::complex<double> alpha_val = static_cast<const std::complex<double>>(alpha);                          \
      bool justTranspose                   = (trans[0] == 'T') || (trans[0] == 't');                                  \
      if (A_is_ll) {                                                                                                  \
        if (justTranspose) {                                                                                          \
          HostBlas<std::complex<double>>::geru(M, N, alpha_val,                                                       \
                                               reinterpret_cast<const std::complex<double>*>(X.data()), one,          \
                                               reinterpret_cast<const std::complex<double>*>(Y.data()), one,          \
                                               reinterpret_cast<std::complex<double>*>(A.data()), LDA);               \
        } else {                                                                                                      \
          HostBlas<std::complex<double>>::gerc(M, N, alpha_val,                                                       \
                                               reinterpret_cast<const std::complex<double>*>(X.data()), one,          \
                                               reinterpret_cast<const std::complex<double>*>(Y.data()), one,          \
                                               reinterpret_cast<std::complex<double>*>(A.data()), LDA);               \
        }                                                                                                             \
      } else {                                                                                                        \
        if (justTranspose) {                                                                                          \
          HostBlas<std::complex<double>>::geru(M, N, alpha_val,                                                       \
                                               reinterpret_cast<const std::complex<double>*>(Y.data()), one,          \
                                               reinterpret_cast<const std::complex<double>*>(X.data()), one,          \
                                               reinterpret_cast<std::complex<double>*>(A.data()), LDA);               \
        } else {                                                                                                      \
          /* blasgerc() + ~A_ll => call kokkos-kernels' implementation */                                             \
          GER<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::ger(space, trans, alpha, X, Y, A); \
        }                                                                                                             \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS2_CGER_BLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                           \
  template <>                                                                                                          \
  struct GER<EXEC_SPACE,                                                                                               \
             Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             Kokkos::View<Kokkos::complex<float>**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                     \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             true, ETI_SPEC_AVAIL> {                                                                                   \
    typedef Kokkos::complex<float> SCALAR;                                                                             \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                      \
        XViewType;                                                                                                     \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                      \
        YViewType;                                                                                                     \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                      \
        AViewType;                                                                                                     \
                                                                                                                       \
    static void ger(const EXEC_SPACE& space, const char trans[], typename AViewType::const_value_type& alpha,          \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                      \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_BLAS,complex<float>");                                        \
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
          GER<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::ger(space, trans, alpha, X, Y, A);  \
        }                                                                                                              \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSBLAS2_DGER_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_DGER_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS2_DGER_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_DGER_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)

KOKKOSBLAS2_SGER_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_SGER_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS2_SGER_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_SGER_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)

KOKKOSBLAS2_ZGER_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_ZGER_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS2_ZGER_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_ZGER_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)

KOKKOSBLAS2_CGER_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_CGER_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS2_CGER_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_CGER_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSBLAS2_DGER_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_DGER_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS2_DGER_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_DGER_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)

KOKKOSBLAS2_SGER_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_SGER_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS2_SGER_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_SGER_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)

KOKKOSBLAS2_ZGER_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_ZGER_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS2_ZGER_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_ZGER_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)

KOKKOSBLAS2_CGER_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_CGER_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS2_CGER_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_CGER_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)
#endif

}  // namespace Impl
}  // namespace KokkosBlas

#endif
