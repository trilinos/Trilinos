// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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

template <typename ExecSpace, bool ETI_SPEC_AVAIL>
  requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)
struct Nrm2<ExecSpace,
            Kokkos::View<double, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,
            Kokkos::View<const double*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1,
            true, ETI_SPEC_AVAIL> {
  typedef Kokkos::View<double, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV;
  typedef Kokkos::View<const double*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV;
  typedef typename XV::size_type size_type;

  static void nrm2(const ExecSpace& space, RV& R, const XV& X, const bool& take_sqrt) {
    Kokkos::Profiling::pushRegion("KokkosBlas::nrm2[TPL_BLAS,double]");
    const size_type numElems = X.extent(0);
    if (numElems < static_cast<size_type>(INT_MAX)) {
      nrm2_print_specialization<RV, XV>();
      int N       = numElems;
      int int_one = 1;
      R()         = HostBlas<double>::nrm2(N, X.data(), int_one);
      if (!take_sqrt) R() = R() * R();
    } else {
      Nrm2<ExecSpace, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm2(space, R, X, take_sqrt);
    }
    Kokkos::Profiling::popRegion();
  }
};

template <typename ExecSpace, bool ETI_SPEC_AVAIL>
  requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)
struct Nrm2<ExecSpace,
            Kokkos::View<float, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,
            Kokkos::View<const float*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1,
            true, ETI_SPEC_AVAIL> {
  typedef Kokkos::View<float, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV;
  typedef Kokkos::View<const float*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV;
  typedef typename XV::size_type size_type;

  static void nrm2(const ExecSpace& space, RV& R, const XV& X, const bool& take_sqrt) {
    Kokkos::Profiling::pushRegion("KokkosBlas::nrm2[TPL_BLAS,float]");
    const size_type numElems = X.extent(0);
    if (numElems < static_cast<size_type>(INT_MAX)) {
      nrm2_print_specialization<RV, XV>();
      int N       = numElems;
      int int_one = 1;
      R()         = HostBlas<float>::nrm2(N, X.data(), int_one);
      if (!take_sqrt) R() = R() * R();
    } else {
      Nrm2<ExecSpace, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm2(space, R, X, take_sqrt);
    }
    Kokkos::Profiling::popRegion();
  }
};

template <typename ExecSpace, bool ETI_SPEC_AVAIL>
  requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)
struct Nrm2<ExecSpace,
            Kokkos::View<double, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,
            Kokkos::View<const Kokkos::complex<double>*, Kokkos::LayoutLeft, ExecSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >,
            1, true, ETI_SPEC_AVAIL> {
  typedef Kokkos::View<double, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV;
  typedef Kokkos::View<const Kokkos::complex<double>*, Kokkos::LayoutLeft, ExecSpace,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      XV;
  typedef typename XV::size_type size_type;

  static void nrm2(const ExecSpace& space, RV& R, const XV& X, const bool& take_sqrt) {
    Kokkos::Profiling::pushRegion("KokkosBlas::nrm2[TPL_BLAS,complex<double>]");
    const size_type numElems = X.extent(0);
    if (numElems < static_cast<size_type>(INT_MAX)) {
      nrm2_print_specialization<RV, XV>();
      int N       = numElems;
      int int_one = 1;
      R() = HostBlas<std::complex<double> >::nrm2(N, reinterpret_cast<const std::complex<double>*>(X.data()), int_one);
      if (!take_sqrt) R() = R() * R();
    } else {
      Nrm2<ExecSpace, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm2(space, R, X, take_sqrt);
    }
    Kokkos::Profiling::popRegion();
  }
};

template <typename ExecSpace, bool ETI_SPEC_AVAIL>
  requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)
struct Nrm2<ExecSpace,
            Kokkos::View<float, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,
            Kokkos::View<const Kokkos::complex<float>*, Kokkos::LayoutLeft, ExecSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >,
            1, true, ETI_SPEC_AVAIL> {
  typedef Kokkos::View<float, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV;
  typedef Kokkos::View<const Kokkos::complex<float>*, Kokkos::LayoutLeft, ExecSpace,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      XV;
  typedef typename XV::size_type size_type;

  static void nrm2(const ExecSpace& space, RV& R, const XV& X, const bool& take_sqrt) {
    Kokkos::Profiling::pushRegion("KokkosBlas::nrm2[TPL_BLAS,complex<float>]");
    const size_type numElems = X.extent(0);
    if (numElems < static_cast<size_type>(INT_MAX)) {
      nrm2_print_specialization<RV, XV>();
      int N       = numElems;
      int int_one = 1;
      R() = HostBlas<std::complex<float> >::nrm2(N, reinterpret_cast<const std::complex<float>*>(X.data()), int_one);
      if (!take_sqrt) R() = R() * R();
    } else {
      Nrm2<ExecSpace, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm2(space, R, X, take_sqrt);
    }
    Kokkos::Profiling::popRegion();
  }
};

}  // namespace Impl
}  // namespace KokkosBlas

#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_CUBLAS(KOKKOS_TYPE, TPL_TYPE, TPL_NRM2)                                         \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct Nrm2<                                                                                                         \
      Kokkos::Cuda,                                                                                                    \
      Kokkos::View<KokkosKernels::ArithTraits<KOKKOS_TYPE>::mag_type, Kokkos::LayoutLeft, Kokkos::HostSpace,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    using RT = KokkosKernels::ArithTraits<KOKKOS_TYPE>::mag_type;                                                      \
    using RV = Kokkos::View<RT, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;      \
    using XV =                                                                                                         \
        Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;  \
    using size_type = typename XV::size_type;                                                                          \
                                                                                                                       \
    static void nrm2(const Kokkos::Cuda& space, RV& R, const XV& X, const bool& take_sqrt) {                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm2[TPL_CUBLAS," + KokkosKernels::ArithTraits<KOKKOS_TYPE>::name() + \
                                    "]");                                                                              \
      const size_type numElems = X.extent(0);                                                                          \
      if (numElems <= static_cast<size_type>(std::numeric_limits<int>::max())) {                                       \
        nrm2_print_specialization<RV, XV>();                                                                           \
        const int N                            = static_cast<int>(numElems);                                           \
        KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                     \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));                              \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(TPL_NRM2(s.handle, N, reinterpret_cast<const TPL_TYPE*>(X.data()), 1, &R())); \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, NULL));                                             \
        if (!take_sqrt) R() = R() * R();                                                                               \
      } else {                                                                                                         \
        Nrm2<Kokkos::Cuda, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm2(space, R, X, take_sqrt);                            \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_CUBLAS(float, float, cublasSnrm2)
KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_CUBLAS(double, double, cublasDnrm2)
KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<float>, cuComplex, cublasScnrm2)
KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<double>, cuDoubleComplex, cublasDznrm2)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ROCBLAS(KOKKOS_TYPE, TPL_TYPE, TPL_NRM2)                                       \
  template <bool ETI_SPEC_AVAIL>                                                                                      \
  struct Nrm2<                                                                                                        \
      Kokkos::HIP,                                                                                                    \
      Kokkos::View<KokkosKernels::ArithTraits<KOKKOS_TYPE>::mag_type, Kokkos::LayoutLeft, Kokkos::HostSpace,          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    using RT = KokkosKernels::ArithTraits<KOKKOS_TYPE>::mag_type;                                                     \
    using RV = Kokkos::View<RT, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;     \
    using XV =                                                                                                        \
        Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;  \
    using size_type = typename XV::size_type;                                                                         \
                                                                                                                      \
    static void nrm2(const Kokkos::HIP& space, RV& R, const XV& X, const bool& take_sqrt) {                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm2[TPL_ROCBLAS," +                                                 \
                                    KokkosKernels::ArithTraits<KOKKOS_TYPE>::name() + "]");                           \
      const size_type numElems = X.extent(0);                                                                         \
      if (numElems <= static_cast<size_type>(std::numeric_limits<rocblas_int>::max())) {                              \
        nrm2_print_specialization<RV, XV>();                                                                          \
        const rocblas_int N                   = static_cast<rocblas_int>(numElems);                                   \
        KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                      \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                          \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                            \
            TPL_NRM2(s.handle, N, reinterpret_cast<const TPL_TYPE*>(X.data()), 1, &R()));                             \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                        \
        if (!take_sqrt) R() = R() * R();                                                                              \
      } else {                                                                                                        \
        Nrm2<Kokkos::HIP, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm2(space, R, X, take_sqrt);                            \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ROCBLAS(float, float, rocblas_snrm2)
KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ROCBLAS(double, double, rocblas_dnrm2)
KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ROCBLAS(Kokkos::complex<float>, rocblas_float_complex, rocblas_scnrm2)
KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ROCBLAS(Kokkos::complex<double>, rocblas_double_complex, rocblas_dznrm2)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_MKL) && defined(KOKKOS_ENABLE_SYCL)
#include <mkl.h>
#include <oneapi/mkl/blas.hpp>
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ONEMKL(KOKKOS_TYPE, TPL_TYPE, TPL_NRM2)                                         \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct Nrm2<                                                                                                         \
      Kokkos::SYCL,                                                                                                    \
      Kokkos::View<KokkosKernels::ArithTraits<KOKKOS_TYPE>::mag_type, Kokkos::LayoutLeft, Kokkos::HostSpace,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, Kokkos::SYCL, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    using RT = KokkosKernels::ArithTraits<KOKKOS_TYPE>::mag_type;                                                      \
    using RV = Kokkos::View<RT, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;      \
    using XV =                                                                                                         \
        Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, Kokkos::SYCL, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;  \
    using size_type = typename XV::size_type;                                                                          \
                                                                                                                       \
    static void nrm2(const Kokkos::SYCL& space, RV& R, const XV& X, const bool& take_sqrt) {                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm2[TPL_ONEMKL," + KokkosKernels::ArithTraits<KOKKOS_TYPE>::name() + \
                                    "]");                                                                              \
      const size_type numElems = X.extent(0);                                                                          \
      if (numElems <= static_cast<size_type>(std::numeric_limits<std::int64_t>::max())) {                              \
        nrm2_print_specialization<RV, XV>();                                                                           \
        const std::int64_t N = static_cast<std::int64_t>(numElems);                                                    \
        Kokkos::View<typename RV::value_type, Kokkos::SYCL> device_result("device_result");                            \
        TPL_NRM2(space.sycl_queue(), N, reinterpret_cast<const TPL_TYPE*>(X.data()), 1, device_result.data());         \
        Kokkos::deep_copy(R, device_result);                                                                           \
        if (!take_sqrt) R() = R() * R();                                                                               \
      } else {                                                                                                         \
        Nrm2<Kokkos::SYCL, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm2(space, R, X, take_sqrt);                            \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ONEMKL(float, float, oneapi::mkl::blas::row_major::nrm2)
KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ONEMKL(double, double, oneapi::mkl::blas::row_major::nrm2)
KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ONEMKL(Kokkos::complex<float>, std::complex<float>, oneapi::mkl::blas::row_major::nrm2)
KOKKOSBLAS1_NRM2_TPL_SPEC_DECL_ONEMKL(Kokkos::complex<double>, std::complex<double>, oneapi::mkl::blas::row_major::nrm2)

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSKERNELS_ENABLE_TPL_MKL && KOKKOS_ENABLE_SYCL

#endif
