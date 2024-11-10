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

#ifndef KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {

namespace {
template <class RV, class XV>
inline void nrm2_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
  printf("KokkosBlas1::nrm2<> TPL Blas specialization for < %s , %s >\n", typeid(RV).name(), typeid(XV).name());
#endif
}
}  // namespace
}  // namespace Impl
}  // namespace KokkosBlas

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DNRM2_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                      \
  template <class ExecSpace>                                                                                        \
  struct Nrm2<ExecSpace, Kokkos::View<double, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
              Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                               \
              1, true, ETI_SPEC_AVAIL> {                                                                            \
    typedef Kokkos::View<double, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV;           \
    typedef Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                  \
        XV;                                                                                                         \
    typedef typename XV::size_type size_type;                                                                       \
                                                                                                                    \
    static void nrm2(const ExecSpace& space, RV& R, const XV& X, const bool& take_sqrt) {                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm2[TPL_BLAS,double]");                                           \
      const size_type numElems = X.extent(0);                                                                       \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                             \
        nrm2_print_specialization<RV, XV>();                                                                        \
        int N       = numElems;                                                                                     \
        int int_one = 1;                                                                                            \
        R()         = HostBlas<double>::nrm2(N, X.data(), int_one);                                                 \
        if (!take_sqrt) R() = R() * R();                                                                            \
      } else {                                                                                                      \
        Nrm2<ExecSpace, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm2(space, R, X, take_sqrt);                            \
      }                                                                                                             \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS1_SNRM2_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                     \
  template <class ExecSpace>                                                                                       \
  struct Nrm2<ExecSpace, Kokkos::View<float, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
              Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                              \
              1, true, ETI_SPEC_AVAIL> {                                                                           \
    typedef Kokkos::View<float, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV;           \
    typedef Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                 \
        XV;                                                                                                        \
    typedef typename XV::size_type size_type;                                                                      \
                                                                                                                   \
    static void nrm2(const ExecSpace& space, RV& R, const XV& X, const bool& take_sqrt) {                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm2[TPL_BLAS,float]");                                           \
      const size_type numElems = X.extent(0);                                                                      \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                            \
        nrm2_print_specialization<RV, XV>();                                                                       \
        int N       = numElems;                                                                                    \
        int int_one = 1;                                                                                           \
        R()         = HostBlas<float>::nrm2(N, X.data(), int_one);                                                 \
        if (!take_sqrt) R() = R() * R();                                                                           \
      } else {                                                                                                     \
        Nrm2<ExecSpace, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm2(space, R, X, take_sqrt);                           \
      }                                                                                                            \
      Kokkos::Profiling::popRegion();                                                                              \
    }                                                                                                              \
  };

#define KOKKOSBLAS1_ZNRM2_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                          \
  template <class ExecSpace>                                                                                            \
  struct Nrm2<ExecSpace, Kokkos::View<double, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,     \
              Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                 \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                   \
              1, true, ETI_SPEC_AVAIL> {                                                                                \
    typedef Kokkos::View<double, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV;               \
    typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                   \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                      \
        XV;                                                                                                             \
    typedef typename XV::size_type size_type;                                                                           \
                                                                                                                        \
    static void nrm2(const ExecSpace& space, RV& R, const XV& X, const bool& take_sqrt) {                               \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm2[TPL_BLAS,complex<double>]");                                      \
      const size_type numElems = X.extent(0);                                                                           \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                                 \
        nrm2_print_specialization<RV, XV>();                                                                            \
        int N       = numElems;                                                                                         \
        int int_one = 1;                                                                                                \
        R()         = HostBlas<std::complex<double> >::nrm2(N, reinterpret_cast<const std::complex<double>*>(X.data()), \
                                                            int_one);                                                   \
        if (!take_sqrt) R() = R() * R();                                                                                \
      } else {                                                                                                          \
        Nrm2<ExecSpace, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm2(space, R, X, take_sqrt);                                \
      }                                                                                                                 \
      Kokkos::Profiling::popRegion();                                                                                   \
    }                                                                                                                   \
  };

#define KOKKOSBLAS1_CNRM2_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                        \
  template <class ExecSpace>                                                                                          \
  struct Nrm2<ExecSpace, Kokkos::View<float, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
              Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                 \
              1, true, ETI_SPEC_AVAIL> {                                                                              \
    typedef Kokkos::View<float, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV;              \
    typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                  \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        XV;                                                                                                           \
    typedef typename XV::size_type size_type;                                                                         \
                                                                                                                      \
    static void nrm2(const ExecSpace& space, RV& R, const XV& X, const bool& take_sqrt) {                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm2[TPL_BLAS,complex<float>]");                                     \
      const size_type numElems = X.extent(0);                                                                         \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                               \
        nrm2_print_specialization<RV, XV>();                                                                          \
        int N       = numElems;                                                                                       \
        int int_one = 1;                                                                                              \
        R() =                                                                                                         \
            HostBlas<std::complex<float> >::nrm2(N, reinterpret_cast<const std::complex<float>*>(X.data()), int_one); \
        if (!take_sqrt) R() = R() * R();                                                                              \
      } else {                                                                                                        \
        Nrm2<ExecSpace, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm2(space, R, X, take_sqrt);                              \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

KOKKOSBLAS1_DNRM2_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_DNRM2_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_SNRM2_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_SNRM2_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_ZNRM2_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_ZNRM2_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_CNRM2_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_CNRM2_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_CUBLAS(LAYOUT, KOKKOS_TYPE, TPL_TYPE, EXECSPACE, MEMSPACE, TPL_NRM2,           \
                                              ETI_SPEC_AVAIL)                                                         \
  template <>                                                                                                         \
  struct Nrm2<EXECSPACE,                                                                                              \
              Kokkos::View<Kokkos::ArithTraits<KOKKOS_TYPE>::mag_type, LAYOUT, Kokkos::HostSpace,                     \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                 \
              Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                           \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                 \
              1, true, ETI_SPEC_AVAIL> {                                                                              \
    using RT        = Kokkos::ArithTraits<KOKKOS_TYPE>::mag_type;                                                     \
    using RV        = Kokkos::View<RT, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;          \
    using XV        = Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                   \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >;                                         \
    using size_type = typename XV::size_type;                                                                         \
                                                                                                                      \
    static void nrm2(const EXECSPACE& space, RV& R, const XV& X, const bool& take_sqrt) {                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm2[TPL_CUBLAS," + Kokkos::ArithTraits<KOKKOS_TYPE>::name() + "]"); \
      const size_type numElems = X.extent(0);                                                                         \
      if (numElems <= static_cast<size_type>(std::numeric_limits<int>::max())) {                                      \
        nrm2_print_specialization<RV, XV>();                                                                          \
        const int N                            = static_cast<int>(numElems);                                          \
        KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                    \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                 \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(TPL_NRM2(s.handle, N, reinterpret_cast<const TPL_TYPE*>(X.data()), 1, &R()));    \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                                \
        if (!take_sqrt) R() = R() * R();                                                                              \
      } else {                                                                                                        \
        Nrm2<EXECSPACE, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm2(space, R, X, take_sqrt);                              \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_CUBLAS_EXT(ETI_SPEC_AVAIL)                                                   \
  KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, float, float, Kokkos::Cuda, Kokkos::CudaSpace,          \
                                        cublasSnrm2, ETI_SPEC_AVAIL)                                                \
  KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, double, double, Kokkos::Cuda, Kokkos::CudaSpace,        \
                                        cublasDnrm2, ETI_SPEC_AVAIL)                                                \
  KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::complex<float>, cuComplex, Kokkos::Cuda,        \
                                        Kokkos::CudaSpace, cublasScnrm2, ETI_SPEC_AVAIL)                            \
  KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::complex<double>, cuDoubleComplex, Kokkos::Cuda, \
                                        Kokkos::CudaSpace, cublasDznrm2, ETI_SPEC_AVAIL)

KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_CUBLAS_EXT(true)
KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_CUBLAS_EXT(false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ROCBLAS(LAYOUT, KOKKOS_TYPE, TPL_TYPE, EXECSPACE, MEMSPACE, TPL_NRM2,           \
                                               ETI_SPEC_AVAIL)                                                         \
  template <>                                                                                                          \
  struct Nrm2<EXECSPACE,                                                                                               \
              Kokkos::View<Kokkos::ArithTraits<KOKKOS_TYPE>::mag_type, LAYOUT, Kokkos::HostSpace,                      \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
              Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                            \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
              1, true, ETI_SPEC_AVAIL> {                                                                               \
    using RT        = Kokkos::ArithTraits<KOKKOS_TYPE>::mag_type;                                                      \
    using RV        = Kokkos::View<RT, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;           \
    using XV        = Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                    \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >;                                          \
    using size_type = typename XV::size_type;                                                                          \
                                                                                                                       \
    static void nrm2(const EXECSPACE& space, RV& R, const XV& X, const bool& take_sqrt) {                              \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm2[TPL_ROCBLAS," + Kokkos::ArithTraits<KOKKOS_TYPE>::name() + "]"); \
      const size_type numElems = X.extent(0);                                                                          \
      if (numElems <= static_cast<size_type>(std::numeric_limits<rocblas_int>::max())) {                               \
        nrm2_print_specialization<RV, XV>();                                                                           \
        const rocblas_int N                   = static_cast<rocblas_int>(numElems);                                    \
        KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                       \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));                               \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(TPL_NRM2(s.handle, N, reinterpret_cast<const TPL_TYPE*>(X.data()), 1, &R()));    \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));                                             \
        if (!take_sqrt) R() = R() * R();                                                                               \
      } else {                                                                                                         \
        Nrm2<EXECSPACE, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm2(space, R, X, take_sqrt);                               \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ROCBLAS_EXT(ETI_SPEC_AVAIL)                                            \
  KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, float, float, Kokkos::HIP, Kokkos::HIPSpace,     \
                                         rocblas_snrm2, ETI_SPEC_AVAIL)                                       \
  KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, double, double, Kokkos::HIP, Kokkos::HIPSpace,   \
                                         rocblas_dnrm2, ETI_SPEC_AVAIL)                                       \
  KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::complex<float>, rocblas_float_complex,   \
                                         Kokkos::HIP, Kokkos::HIPSpace, rocblas_scnrm2, ETI_SPEC_AVAIL)       \
  KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::complex<double>, rocblas_double_complex, \
                                         Kokkos::HIP, Kokkos::HIPSpace, rocblas_dznrm2, ETI_SPEC_AVAIL)

KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ROCBLAS_EXT(true)
KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ROCBLAS_EXT(false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_MKL) && defined(KOKKOS_ENABLE_SYCL)
#include <mkl.h>
#include <oneapi/mkl/blas.hpp>
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ONEMKL(LAYOUT, KOKKOS_TYPE, TPL_TYPE, EXECSPACE, MEMSPACE, TPL_NRM2,           \
                                              ETI_SPEC_AVAIL)                                                         \
  template <>                                                                                                         \
  struct Nrm2<EXECSPACE,                                                                                              \
              Kokkos::View<Kokkos::ArithTraits<KOKKOS_TYPE>::mag_type, LAYOUT, Kokkos::HostSpace,                     \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                 \
              Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                           \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                 \
              1, true, ETI_SPEC_AVAIL> {                                                                              \
    using RT        = Kokkos::ArithTraits<KOKKOS_TYPE>::mag_type;                                                     \
    using RV        = Kokkos::View<RT, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;          \
    using XV        = Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                   \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >;                                         \
    using size_type = typename XV::size_type;                                                                         \
                                                                                                                      \
    static void nrm2(const EXECSPACE& space, RV& R, const XV& X, const bool& take_sqrt) {                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm2[TPL_ONEMKL," + Kokkos::ArithTraits<KOKKOS_TYPE>::name() + "]"); \
      const size_type numElems = X.extent(0);                                                                         \
      if (numElems <= static_cast<size_type>(std::numeric_limits<std::int64_t>::max())) {                             \
        nrm2_print_specialization<RV, XV>();                                                                          \
        const std::int64_t N = static_cast<std::int64_t>(numElems);                                                   \
        TPL_NRM2(space.sycl_queue(), N, reinterpret_cast<const TPL_TYPE*>(X.data()), 1, &R());                        \
        if (!take_sqrt) R() = R() * R();                                                                              \
      } else {                                                                                                        \
        Nrm2<EXECSPACE, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm2(space, R, X, take_sqrt);                              \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ONEMKL_EXT(ETI_SPEC_AVAIL)                                                     \
  KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ONEMKL(Kokkos::LayoutLeft, float, float, Kokkos::Experimental::SYCL,                 \
                                        Kokkos::Experimental::SYCLDeviceUSMSpace, oneapi::mkl::blas::row_major::nrm2, \
                                        ETI_SPEC_AVAIL)                                                               \
  KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ONEMKL(Kokkos::LayoutLeft, double, double, Kokkos::Experimental::SYCL,               \
                                        Kokkos::Experimental::SYCLDeviceUSMSpace, oneapi::mkl::blas::row_major::nrm2, \
                                        ETI_SPEC_AVAIL)                                                               \
  KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ONEMKL(Kokkos::LayoutLeft, Kokkos::complex<float>, std::complex<float>,              \
                                        Kokkos::Experimental::SYCL, Kokkos::Experimental::SYCLDeviceUSMSpace,         \
                                        oneapi::mkl::blas::row_major::nrm2, ETI_SPEC_AVAIL)                           \
  KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ONEMKL(Kokkos::LayoutLeft, Kokkos::complex<double>, std::complex<double>,            \
                                        Kokkos::Experimental::SYCL, Kokkos::Experimental::SYCLDeviceUSMSpace,         \
                                        oneapi::mkl::blas::row_major::nrm2, ETI_SPEC_AVAIL)

KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ONEMKL_EXT(true)
KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ONEMKL_EXT(false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSKERNELS_ENABLE_TPL_MKL && KOKKOS_ENABLE_SYCL

#endif
