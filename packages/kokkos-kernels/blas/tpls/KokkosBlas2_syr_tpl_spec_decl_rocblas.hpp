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

#ifndef KOKKOSBLAS2_SYR_TPL_SPEC_DECL_ROCBLAS_HPP_
#define KOKKOSBLAS2_SYR_TPL_SPEC_DECL_ROCBLAS_HPP_

#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS2_SYR_ROCBLAS_DETERMINE_ARGS(LAYOUT, uploChar)                 \
  bool A_is_ll          = std::is_same<Kokkos::LayoutLeft, LAYOUT>::value;       \
  bool A_is_lr          = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;      \
  const int N           = static_cast<int>(A_is_lr ? A.extent(0) : A.extent(1)); \
  constexpr int one     = 1;                                                     \
  const int LDA         = A_is_lr ? A.stride(0) : A.stride(1);                   \
  rocblas_fill fillMode = (uploChar == 'L' || uploChar == 'l') ? rocblas_fill_lower : rocblas_fill_upper;

#define KOKKOSBLAS2_DSYR_ROCBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                        \
  template <>                                                                                                          \
  struct SYR<                                                                                                          \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<const double*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<double**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    typedef double SCALAR;                                                                                             \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        XViewType;                                                                                                     \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        AViewType;                                                                                                     \
                                                                                                                       \
    static void syr(const typename AViewType::execution_space& space, const char trans[], const char uplo[],           \
                    typename AViewType::const_value_type& alpha, const XViewType& X, const AViewType& A) {             \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr[TPL_ROCBLAS,double]");                                            \
      KOKKOSBLAS2_SYR_ROCBLAS_DETERMINE_ARGS(LAYOUT, uplo[0]);                                                         \
      if (A_is_ll) {                                                                                                   \
        KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                       \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));                               \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(s.handle, rocblas_pointer_mode_host));                  \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_dsyr(s.handle, fillMode, N, &alpha, X.data(), one, A.data(), LDA));      \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));                                             \
      } else {                                                                                                         \
        /* rocblas_dsyr() + ~A_ll => call kokkos-kernels' implementation */                                            \
        SYR<EXEC_SPACE, XViewType, AViewType, false, ETI_SPEC_AVAIL>::syr(space, trans, uplo, alpha, X, A);            \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS2_SSYR_ROCBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                       \
  template <>                                                                                                         \
  struct SYR<                                                                                                         \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<float**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    typedef float SCALAR;                                                                                             \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        XViewType;                                                                                                    \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        AViewType;                                                                                                    \
                                                                                                                      \
    static void syr(const typename AViewType::execution_space& space, const char trans[], const char uplo[],          \
                    typename AViewType::const_value_type& alpha, const XViewType& X, const AViewType& A) {            \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr[TPL_ROCBLAS,float]");                                            \
      KOKKOSBLAS2_SYR_ROCBLAS_DETERMINE_ARGS(LAYOUT, uplo[0]);                                                        \
      if (A_is_ll) {                                                                                                  \
        KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                      \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));                              \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(s.handle, rocblas_pointer_mode_host));                 \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_ssyr(s.handle, fillMode, N, &alpha, X.data(), one, A.data(), LDA));     \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));                                            \
      } else {                                                                                                        \
        /* rocblas_ssyr() + ~A_ll => call kokkos-kernels' implementation */                                           \
        SYR<EXEC_SPACE, XViewType, AViewType, false, ETI_SPEC_AVAIL>::syr(space, trans, uplo, alpha, X, A);           \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS2_ZSYR_ROCBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                      \
  template <>                                                                                                        \
  struct SYR<EXEC_SPACE,                                                                                             \
             Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,             \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                 \
             Kokkos::View<Kokkos::complex<double>**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                  \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                 \
             true, ETI_SPEC_AVAIL> {                                                                                 \
    typedef Kokkos::complex<double> SCALAR;                                                                          \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                   \
        XViewType;                                                                                                   \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                    \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                   \
        AViewType;                                                                                                   \
                                                                                                                     \
    static void syr(const typename AViewType::execution_space& space, const char trans[], const char uplo[],         \
                    typename AViewType::const_value_type& alpha, const XViewType& X, const AViewType& A) {           \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr[TPL_ROCBLAS,complex<double>]");                                 \
      KOKKOSBLAS2_SYR_ROCBLAS_DETERMINE_ARGS(LAYOUT, uplo[0]);                                                       \
      bool justTranspose = (trans[0] == 'T') || (trans[0] == 't');                                                   \
      if (justTranspose) {                                                                                           \
        if (A_is_ll) {                                                                                               \
          KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                   \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));                           \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(s.handle, rocblas_pointer_mode_host));              \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_zsyr(s.handle, fillMode, N,                                          \
                                                     reinterpret_cast<const rocblas_double_complex*>(&alpha),        \
                                                     reinterpret_cast<const rocblas_double_complex*>(X.data()), one, \
                                                     reinterpret_cast<rocblas_double_complex*>(A.data()), LDA));     \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));                                         \
        } else {                                                                                                     \
          /* rocblas_zsyr() + ~A_ll => call kokkos-kernels' implementation */                                        \
          SYR<EXEC_SPACE, XViewType, AViewType, false, ETI_SPEC_AVAIL>::syr(space, trans, uplo, alpha, X, A);        \
        }                                                                                                            \
      } else {                                                                                                       \
        if (A_is_ll && (alpha.imag() == 0.)) {                                                                       \
          const double alpha_val                = alpha.real();                                                      \
          KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                   \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));                           \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(s.handle, rocblas_pointer_mode_host));              \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_zher(s.handle, fillMode, N, &alpha_val,                              \
                                                     reinterpret_cast<const rocblas_double_complex*>(X.data()), one, \
                                                     reinterpret_cast<rocblas_double_complex*>(A.data()), LDA));     \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));                                         \
        } else {                                                                                                     \
          /* rocblas_zher() + [~A_ll or ~real alpha]=> call kokkos-kernels'                                          \
           * implementation */                                                                                       \
          SYR<EXEC_SPACE, XViewType, AViewType, false, ETI_SPEC_AVAIL>::syr(space, trans, uplo, alpha, X, A);        \
        }                                                                                                            \
      }                                                                                                              \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS2_CSYR_ROCBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                     \
  template <>                                                                                                       \
  struct SYR<EXEC_SPACE,                                                                                            \
             Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,             \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                \
             Kokkos::View<Kokkos::complex<float>**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                  \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                \
             true, ETI_SPEC_AVAIL> {                                                                                \
    typedef Kokkos::complex<float> SCALAR;                                                                          \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                              \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                  \
        XViewType;                                                                                                  \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                   \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                  \
        AViewType;                                                                                                  \
                                                                                                                    \
    static void syr(const typename AViewType::execution_space& space, const char trans[], const char uplo[],        \
                    typename AViewType::const_value_type& alpha, const XViewType& X, const AViewType& A) {          \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr[TPL_ROCBLAS,complex<float>]");                                 \
      KOKKOSBLAS2_SYR_ROCBLAS_DETERMINE_ARGS(LAYOUT, uplo[0]);                                                      \
      bool justTranspose = (trans[0] == 'T') || (trans[0] == 't');                                                  \
      if (justTranspose) {                                                                                          \
        if (A_is_ll) {                                                                                              \
          KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                  \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));                          \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(s.handle, rocblas_pointer_mode_host));             \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_csyr(s.handle, fillMode, N,                                         \
                                                     reinterpret_cast<const rocblas_float_complex*>(&alpha),        \
                                                     reinterpret_cast<const rocblas_float_complex*>(X.data()), one, \
                                                     reinterpret_cast<rocblas_float_complex*>(A.data()), LDA));     \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));                                        \
        } else {                                                                                                    \
          /* rocblas_csyr() + ~A_ll => call kokkos-kernels' implementation */                                       \
          SYR<EXEC_SPACE, XViewType, AViewType, false, ETI_SPEC_AVAIL>::syr(space, trans, uplo, alpha, X, A);       \
        }                                                                                                           \
      } else {                                                                                                      \
        if (A_is_ll && (alpha.imag() == 0.)) {                                                                      \
          const float alpha_val                 = alpha.real();                                                     \
          KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                  \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));                          \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(s.handle, rocblas_pointer_mode_host));             \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_cher(s.handle, fillMode, N, &alpha_val,                             \
                                                     reinterpret_cast<const rocblas_float_complex*>(X.data()), one, \
                                                     reinterpret_cast<rocblas_float_complex*>(A.data()), LDA));     \
          KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));                                        \
        } else {                                                                                                    \
          /* rocblas_cher() + [~A_ll or ~real alpha]=> call kokkos-kernels'                                         \
           * implementation */                                                                                      \
          SYR<EXEC_SPACE, XViewType, AViewType, false, ETI_SPEC_AVAIL>::syr(space, trans, uplo, alpha, X, A);       \
        }                                                                                                           \
      }                                                                                                             \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

KOKKOSBLAS2_DSYR_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_DSYR_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS2_DSYR_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_DSYR_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

KOKKOSBLAS2_SSYR_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_SSYR_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS2_SSYR_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_SSYR_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

KOKKOSBLAS2_ZSYR_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_ZSYR_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS2_ZSYR_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_ZSYR_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

KOKKOSBLAS2_CSYR_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_CSYR_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS2_CSYR_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS2_CSYR_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif
