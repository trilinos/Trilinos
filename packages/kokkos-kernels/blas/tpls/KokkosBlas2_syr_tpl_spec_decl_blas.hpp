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

#ifndef KOKKOSBLAS2_SYR_TPL_SPEC_DECL_BLAS_HPP_
#define KOKKOSBLAS2_SYR_TPL_SPEC_DECL_BLAS_HPP_

#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS2_SYR_DETERMINE_ARGS(LAYOUT)                               \
  bool A_is_ll      = std::is_same<Kokkos::LayoutLeft, LAYOUT>::value;       \
  bool A_is_lr      = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;      \
  const int N       = static_cast<int>(A_is_lr ? A.extent(0) : A.extent(1)); \
  constexpr int one = 1;                                                     \
  const int LDA     = A_is_lr ? A.stride(0) : A.stride(1);

#define KOKKOSBLAS2_DSYR_BLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                          \
  template <>                                                                                                         \
  struct SYR<                                                                                                         \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<const double*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      Kokkos::View<double**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    typedef double SCALAR;                                                                                            \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                     \
        XViewType;                                                                                                    \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                     \
        AViewType;                                                                                                    \
                                                                                                                      \
    static void syr(const typename AViewType::execution_space& space, const char trans[], const char uplo[],          \
                    typename AViewType::const_value_type& alpha, const XViewType& X, const AViewType& A) {            \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr[TPL_BLAS,double]");                                              \
      KOKKOSBLAS2_SYR_DETERMINE_ARGS(LAYOUT);                                                                         \
      if (A_is_ll) {                                                                                                  \
        HostBlas<SCALAR>::syr(uplo[0], N, alpha, X.data(), one, A.data(), LDA);                                       \
      } else {                                                                                                        \
        /* blasDsyr() + ~A_ll => call kokkos-kernels' implementation */                                               \
        SYR<EXEC_SPACE, XViewType, AViewType, false, ETI_SPEC_AVAIL>::syr(space, trans, uplo, alpha, X, A);           \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS2_SSYR_BLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                         \
  template <>                                                                                                        \
  struct SYR<                                                                                                        \
      EXEC_SPACE,                                                                                                    \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                         \
      Kokkos::View<float**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      true, ETI_SPEC_AVAIL> {                                                                                        \
    typedef float SCALAR;                                                                                            \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                    \
        XViewType;                                                                                                   \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                    \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                    \
        AViewType;                                                                                                   \
                                                                                                                     \
    static void syr(const typename AViewType::execution_space& space, const char trans[], const char uplo[],         \
                    typename AViewType::const_value_type& alpha, const XViewType& X, const AViewType& A) {           \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr[TPL_BLAS,float]");                                              \
      KOKKOSBLAS2_SYR_DETERMINE_ARGS(LAYOUT);                                                                        \
      if (A_is_ll) {                                                                                                 \
        HostBlas<SCALAR>::syr(uplo[0], N, alpha, X.data(), one, A.data(), LDA);                                      \
      } else {                                                                                                       \
        /* blasSsyr() + ~A_ll => call kokkos-kernels' implementation */                                              \
        SYR<EXEC_SPACE, XViewType, AViewType, false, ETI_SPEC_AVAIL>::syr(space, trans, uplo, alpha, X, A);          \
      }                                                                                                              \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS2_ZSYR_BLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                        \
  template <>                                                                                                       \
  struct SYR<EXEC_SPACE,                                                                                            \
             Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,            \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                 \
             Kokkos::View<Kokkos::complex<double>**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                 \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                 \
             true, ETI_SPEC_AVAIL> {                                                                                \
    typedef Kokkos::complex<double> SCALAR;                                                                         \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                              \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                   \
        XViewType;                                                                                                  \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                   \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                   \
        AViewType;                                                                                                  \
                                                                                                                    \
    static void syr(const typename AViewType::execution_space& space, const char trans[], const char uplo[],        \
                    typename AViewType::const_value_type& alpha, const XViewType& X, const AViewType& A) {          \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr[TPL_BLAS,complex<double>");                                    \
      KOKKOSBLAS2_SYR_DETERMINE_ARGS(LAYOUT);                                                                       \
      bool justTranspose = (trans[0] == 'T') || (trans[0] == 't');                                                  \
      if (justTranspose) {                                                                                          \
        /* No blasZsyr() => call kokkos-kernels' implementation */                                                  \
        SYR<EXEC_SPACE, XViewType, AViewType, false, ETI_SPEC_AVAIL>::syr(space, trans, uplo, alpha, X, A);         \
      } else {                                                                                                      \
        if (A_is_ll) {                                                                                              \
          HostBlas<std::complex<double>>::her<double>(uplo[0], N, alpha.real(),                                     \
                                                      reinterpret_cast<const std::complex<double>*>(X.data()), one, \
                                                      reinterpret_cast<std::complex<double>*>(A.data()), LDA);      \
        } else {                                                                                                    \
          /* blasZher() + [~A_ll or ~real alpha] => call kokkos-kernels'                                            \
           * implementation */                                                                                      \
          SYR<EXEC_SPACE, XViewType, AViewType, false, ETI_SPEC_AVAIL>::syr(space, trans, uplo, alpha, X, A);       \
        }                                                                                                           \
      }                                                                                                             \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS2_CSYR_BLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                     \
  template <>                                                                                                    \
  struct SYR<EXEC_SPACE,                                                                                         \
             Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,          \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
             Kokkos::View<Kokkos::complex<float>**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,               \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
             true, ETI_SPEC_AVAIL> {                                                                             \
    typedef Kokkos::complex<float> SCALAR;                                                                       \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                           \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                \
        XViewType;                                                                                               \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                \
        AViewType;                                                                                               \
                                                                                                                 \
    static void syr(const typename AViewType::execution_space& space, const char trans[], const char uplo[],     \
                    typename AViewType::const_value_type& alpha, const XViewType& X, const AViewType& A) {       \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr[TPL_BLAS,complex<float>");                                  \
      KOKKOSBLAS2_SYR_DETERMINE_ARGS(LAYOUT);                                                                    \
      bool justTranspose = (trans[0] == 'T') || (trans[0] == 't');                                               \
      if (justTranspose) {                                                                                       \
        /* No blasCsyr() => call kokkos-kernels' implementation */                                               \
        SYR<EXEC_SPACE, XViewType, AViewType, false, ETI_SPEC_AVAIL>::syr(space, trans, uplo, alpha, X, A);      \
      } else {                                                                                                   \
        if (A_is_ll && (alpha.imag() == 0.)) {                                                                   \
          HostBlas<std::complex<float>>::her<float>(uplo[0], N, alpha.real(),                                    \
                                                    reinterpret_cast<const std::complex<float>*>(X.data()), one, \
                                                    reinterpret_cast<std::complex<float>*>(A.data()), LDA);      \
        } else {                                                                                                 \
          /* blasCher() + [~A_ll or ~real alpha] => call kokkos-kernels'                                         \
           * implementation */                                                                                   \
          SYR<EXEC_SPACE, XViewType, AViewType, false, ETI_SPEC_AVAIL>::syr(space, trans, uplo, alpha, X, A);    \
        }                                                                                                        \
      }                                                                                                          \
      Kokkos::Profiling::popRegion();                                                                            \
    }                                                                                                            \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSBLAS2_DSYR_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_DSYR_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS2_DSYR_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_DSYR_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)

KOKKOSBLAS2_SSYR_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_SSYR_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS2_SSYR_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_SSYR_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)

KOKKOSBLAS2_ZSYR_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_ZSYR_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS2_ZSYR_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_ZSYR_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)

KOKKOSBLAS2_CSYR_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_CSYR_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS2_CSYR_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS2_CSYR_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSBLAS2_DSYR_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_DSYR_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS2_DSYR_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_DSYR_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)

KOKKOSBLAS2_SSYR_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_SSYR_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS2_SSYR_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_SSYR_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)

KOKKOSBLAS2_ZSYR_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_ZSYR_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS2_ZSYR_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_ZSYR_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)

KOKKOSBLAS2_CSYR_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_CSYR_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS2_CSYR_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS2_CSYR_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)
#endif

}  // namespace Impl
}  // namespace KokkosBlas

#endif
