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

#define KOKKOSBLAS2_DGER_ROCBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                        \
  template <>                                                                                                          \
  struct GER<                                                                                                          \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<const double*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const double*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<double**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    typedef double SCALAR;                                                                                             \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        XViewType;                                                                                                     \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        YViewType;                                                                                                     \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        AViewType;                                                                                                     \
                                                                                                                       \
    static void ger(const EXEC_SPACE& space, const char /*trans*/[], typename AViewType::const_value_type& alpha,      \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                      \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_ROCBLAS,double]");                                            \
      KOKKOSBLAS2_GER_ROCBLAS_DETERMINE_ARGS(LAYOUT);                                                                  \
      KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                         \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));                                 \
      if (A_is_ll) {                                                                                                   \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(                                                                                 \
            rocblas_dger(s.handle, M, N, &alpha, X.data(), one, Y.data(), one, A.data(), LDA));                        \
      } else {                                                                                                         \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(                                                                                 \
            rocblas_dger(s.handle, M, N, &alpha, Y.data(), one, X.data(), one, A.data(), LDA));                        \
      }                                                                                                                \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));                                               \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS2_SGER_ROCBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                       \
  template <>                                                                                                         \
  struct GER<                                                                                                         \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<float**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    typedef float SCALAR;                                                                                             \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        XViewType;                                                                                                    \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        YViewType;                                                                                                    \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        AViewType;                                                                                                    \
                                                                                                                      \
    static void ger(const EXEC_SPACE& space, const char /*trans*/[], typename AViewType::const_value_type& alpha,     \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_ROCBLAS,float]");                                            \
      KOKKOSBLAS2_GER_ROCBLAS_DETERMINE_ARGS(LAYOUT);                                                                 \
      KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                        \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));                                \
      if (A_is_ll) {                                                                                                  \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(                                                                                \
            rocblas_sger(s.handle, M, N, &alpha, X.data(), one, Y.data(), one, A.data(), LDA));                       \
      } else {                                                                                                        \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(                                                                                \
            rocblas_sger(s.handle, M, N, &alpha, Y.data(), one, X.data(), one, A.data(), LDA));                       \
      }                                                                                                               \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));                                              \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS2_ZGER_ROCBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                       \
  template <>                                                                                                         \
  struct GER<EXEC_SPACE,                                                                                              \
             Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,              \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
             Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,              \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
             Kokkos::View<Kokkos::complex<double>**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
             true, ETI_SPEC_AVAIL> {                                                                                  \
    typedef Kokkos::complex<double> SCALAR;                                                                           \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        XViewType;                                                                                                    \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        YViewType;                                                                                                    \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        AViewType;                                                                                                    \
                                                                                                                      \
    static void ger(const EXEC_SPACE& space, const char trans[], typename AViewType::const_value_type& alpha,         \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_ROCBLAS,complex<double>]");                                  \
      KOKKOSBLAS2_GER_ROCBLAS_DETERMINE_ARGS(LAYOUT);                                                                 \
      bool justTranspose                    = (trans[0] == 'T') || (trans[0] == 't');                                 \
      KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                        \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));                                \
      if (A_is_ll) {                                                                                                  \
        if (justTranspose) {                                                                                          \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_zgeru(s.handle, M, N,                                                 \
                                                      reinterpret_cast<const rocblas_double_complex*>(&alpha),        \
                                                      reinterpret_cast<const rocblas_double_complex*>(X.data()), one, \
                                                      reinterpret_cast<const rocblas_double_complex*>(Y.data()), one, \
                                                      reinterpret_cast<rocblas_double_complex*>(A.data()), LDA));     \
        } else {                                                                                                      \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_zgerc(s.handle, M, N,                                                 \
                                                      reinterpret_cast<const rocblas_double_complex*>(&alpha),        \
                                                      reinterpret_cast<const rocblas_double_complex*>(X.data()), one, \
                                                      reinterpret_cast<const rocblas_double_complex*>(Y.data()), one, \
                                                      reinterpret_cast<rocblas_double_complex*>(A.data()), LDA));     \
        }                                                                                                             \
      } else {                                                                                                        \
        if (justTranspose) {                                                                                          \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_zgeru(s.handle, M, N,                                                 \
                                                      reinterpret_cast<const rocblas_double_complex*>(&alpha),        \
                                                      reinterpret_cast<const rocblas_double_complex*>(Y.data()), one, \
                                                      reinterpret_cast<const rocblas_double_complex*>(X.data()), one, \
                                                      reinterpret_cast<rocblas_double_complex*>(A.data()), LDA));     \
        } else {                                                                                                      \
          /* rocblas_zgerc() + ~A_ll => call k-kernels' implementation */                                             \
          GER<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::ger(space, trans, alpha, X, Y, A); \
        }                                                                                                             \
      }                                                                                                               \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));                                              \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS2_CGER_ROCBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                       \
  template <>                                                                                                         \
  struct GER<EXEC_SPACE,                                                                                              \
             Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,               \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
             Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,               \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
             Kokkos::View<Kokkos::complex<float>**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                    \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
             true, ETI_SPEC_AVAIL> {                                                                                  \
    typedef Kokkos::complex<float> SCALAR;                                                                            \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        XViewType;                                                                                                    \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        YViewType;                                                                                                    \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        AViewType;                                                                                                    \
                                                                                                                      \
    static void ger(const EXEC_SPACE& space, const char trans[], typename AViewType::const_value_type& alpha,         \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_ROCBLAS,complex<float>]");                                   \
      KOKKOSBLAS2_GER_ROCBLAS_DETERMINE_ARGS(LAYOUT);                                                                 \
      bool justTranspose                    = (trans[0] == 'T') || (trans[0] == 't');                                 \
      KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                        \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));                                \
      if (A_is_ll) {                                                                                                  \
        if (justTranspose) {                                                                                          \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_cgeru(s.handle, M, N,                                                 \
                                                      reinterpret_cast<const rocblas_float_complex*>(&alpha),         \
                                                      reinterpret_cast<const rocblas_float_complex*>(X.data()), one,  \
                                                      reinterpret_cast<const rocblas_float_complex*>(Y.data()), one,  \
                                                      reinterpret_cast<rocblas_float_complex*>(A.data()), LDA));      \
        } else {                                                                                                      \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_cgerc(s.handle, M, N,                                                 \
                                                      reinterpret_cast<const rocblas_float_complex*>(&alpha),         \
                                                      reinterpret_cast<const rocblas_float_complex*>(X.data()), one,  \
                                                      reinterpret_cast<const rocblas_float_complex*>(Y.data()), one,  \
                                                      reinterpret_cast<rocblas_float_complex*>(A.data()), LDA));      \
        }                                                                                                             \
      } else {                                                                                                        \
        if (justTranspose) {                                                                                          \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_cgeru(s.handle, M, N,                                                 \
                                                      reinterpret_cast<const rocblas_float_complex*>(&alpha),         \
                                                      reinterpret_cast<const rocblas_float_complex*>(Y.data()), one,  \
                                                      reinterpret_cast<const rocblas_float_complex*>(X.data()), one,  \
                                                      reinterpret_cast<rocblas_float_complex*>(A.data()), LDA));      \
        } else {                                                                                                      \
          /* rocblas_cgerc() + ~A_ll => call k-kernels' implementation */                                             \
          GER<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::ger(space, trans, alpha, X, Y, A); \
        }                                                                                                             \
      }                                                                                                               \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));                                              \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

KOKKOSBLAS2_DGER_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_DGER_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS2_DGER_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_DGER_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

KOKKOSBLAS2_SGER_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_SGER_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS2_SGER_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_SGER_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

KOKKOSBLAS2_ZGER_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_ZGER_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS2_ZGER_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_ZGER_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

KOKKOSBLAS2_CGER_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_CGER_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS2_CGER_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_CGER_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif
