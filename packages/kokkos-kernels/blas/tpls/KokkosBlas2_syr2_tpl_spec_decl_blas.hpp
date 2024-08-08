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

#ifndef KOKKOSBLAS2_SYR2_TPL_SPEC_DECL_BLAS_HPP_
#define KOKKOSBLAS2_SYR2_TPL_SPEC_DECL_BLAS_HPP_

#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS2_SYR2_DETERMINE_ARGS(LAYOUT)                              \
  bool A_is_ll      = std::is_same<Kokkos::LayoutLeft, LAYOUT>::value;       \
  bool A_is_lr      = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;      \
  const int N       = static_cast<int>(A_is_lr ? A.extent(0) : A.extent(1)); \
  constexpr int one = 1;                                                     \
  const int LDA     = A_is_lr ? A.stride(0) : A.stride(1);

#define KOKKOSBLAS2_DSYR2_BLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                         \
  template <>                                                                                                         \
  struct SYR2<                                                                                                        \
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
    static void syr2(const typename AViewType::execution_space& space, const char trans[], const char uplo[],         \
                     typename AViewType::const_value_type& alpha, const XViewType& X, const YViewType& Y,             \
                     const AViewType& A) {                                                                            \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr2[TPL_BLAS,double]");                                             \
      KOKKOSBLAS2_SYR2_DETERMINE_ARGS(LAYOUT);                                                                        \
      if (A_is_ll) {                                                                                                  \
        HostBlas<SCALAR>::syr2(uplo[0], N, alpha, X.data(), one, Y.data(), one, A.data(), LDA);                       \
      } else {                                                                                                        \
        /* blasDsyr2() + ~A_ll => call kokkos-kernels' implementation */                                              \
        SYR2<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, X,  \
                                                                                       Y, A);                         \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS2_SSYR2_BLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                        \
  template <>                                                                                                        \
  struct SYR2<                                                                                                       \
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
    static void syr2(const typename AViewType::execution_space& space, const char trans[], const char uplo[],        \
                     typename AViewType::const_value_type& alpha, const XViewType& X, const YViewType& Y,            \
                     const AViewType& A) {                                                                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr2[TPL_BLAS,float]");                                             \
      KOKKOSBLAS2_SYR2_DETERMINE_ARGS(LAYOUT);                                                                       \
      if (A_is_ll) {                                                                                                 \
        HostBlas<SCALAR>::syr2(uplo[0], N, alpha, X.data(), one, Y.data(), one, A.data(), LDA);                      \
      } else {                                                                                                       \
        /* blasSsyr2() + ~A_ll => call kokkos-kernels' implementation */                                             \
        SYR2<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, X, \
                                                                                       Y, A);                        \
      }                                                                                                              \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS2_ZSYR2_BLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                          \
  template <>                                                                                                          \
  struct SYR2<EXEC_SPACE,                                                                                              \
              Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                   \
              Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                   \
              Kokkos::View<Kokkos::complex<double>**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                   \
              true, ETI_SPEC_AVAIL> {                                                                                  \
    typedef Kokkos::complex<double> SCALAR;                                                                            \
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
    static void syr2(const typename AViewType::execution_space& space, const char trans[], const char uplo[],          \
                     typename AViewType::const_value_type& alpha, const XViewType& X, const YViewType& Y,              \
                     const AViewType& A) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr2[TPL_BLAS,complex<double>");                                      \
      KOKKOSBLAS2_SYR2_DETERMINE_ARGS(LAYOUT);                                                                         \
      bool justTranspose = (trans[0] == 'T') || (trans[0] == 't');                                                     \
      if (justTranspose) {                                                                                             \
        /* No blasZsyr2() => call kokkos-kernels' implementation */                                                    \
        SYR2<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, X,   \
                                                                                       Y, A);                          \
      } else {                                                                                                         \
        if (A_is_ll) {                                                                                                 \
          HostBlas<std::complex<double>>::her2(uplo[0], N, alpha,                                                      \
                                               reinterpret_cast<const std::complex<double>*>(X.data()), one,           \
                                               reinterpret_cast<const std::complex<double>*>(Y.data()), one,           \
                                               reinterpret_cast<std::complex<double>*>(A.data()), LDA);                \
        } else {                                                                                                       \
          /* blasZher2() + ~A_ll => call kokkos-kernels' implementation */                                             \
          SYR2<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, X, \
                                                                                         Y, A);                        \
        }                                                                                                              \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS2_CSYR2_BLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                          \
  template <>                                                                                                          \
  struct SYR2<EXEC_SPACE,                                                                                              \
              Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,               \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                   \
              Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,               \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                   \
              Kokkos::View<Kokkos::complex<float>**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                    \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                   \
              true, ETI_SPEC_AVAIL> {                                                                                  \
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
    static void syr2(const typename AViewType::execution_space& space, const char trans[], const char uplo[],          \
                     typename AViewType::const_value_type& alpha, const XViewType& X, const YViewType& Y,              \
                     const AViewType& A) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr2[TPL_BLAS,complex<float>");                                       \
      KOKKOSBLAS2_SYR2_DETERMINE_ARGS(LAYOUT);                                                                         \
      bool justTranspose = (trans[0] == 'T') || (trans[0] == 't');                                                     \
      if (justTranspose) {                                                                                             \
        /* No blasCsyr2() => call kokkos-kernels' implementation */                                                    \
        SYR2<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, X,   \
                                                                                       Y, A);                          \
      } else {                                                                                                         \
        if (A_is_ll) {                                                                                                 \
          HostBlas<std::complex<float>>::her2(uplo[0], N, alpha,                                                       \
                                              reinterpret_cast<const std::complex<float>*>(X.data()), one,             \
                                              reinterpret_cast<const std::complex<float>*>(Y.data()), one,             \
                                              reinterpret_cast<std::complex<float>*>(A.data()), LDA);                  \
        } else {                                                                                                       \
          /* blasCher2() + ~A_ll => call kokkos-kernels' implementation */                                             \
          SYR2<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, X, \
                                                                                         Y, A);                        \
        }                                                                                                              \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSBLAS2_DSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_DSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS2_DSYR2_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_DSYR2_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)

KOKKOSBLAS2_SSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_SSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS2_SSYR2_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_SSYR2_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)

KOKKOSBLAS2_ZSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_ZSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS2_ZSYR2_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_ZSYR2_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)

KOKKOSBLAS2_CSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_CSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS2_CSYR2_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_CSYR2_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSBLAS2_DSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_DSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS2_DSYR2_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_DSYR2_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)

KOKKOSBLAS2_SSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_SSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS2_SSYR2_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_SSYR2_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)

KOKKOSBLAS2_ZSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_ZSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS2_ZSYR2_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_ZSYR2_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)

KOKKOSBLAS2_CSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_CSYR2_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS2_CSYR2_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_CSYR2_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)
#endif

}  // namespace Impl
}  // namespace KokkosBlas

#endif
