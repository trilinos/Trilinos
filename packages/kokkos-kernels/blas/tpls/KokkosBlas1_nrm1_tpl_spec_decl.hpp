// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {

namespace {
template <class RV, class XV>
inline void nrm1_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
  printf("KokkosBlas1::nrm1<> TPL Blas specialization for < %s , %s >\n", typeid(RV).name(), typeid(XV).name());
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

#define KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(SCALAR, EXECSPACE)                                                         \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct Nrm1<EXECSPACE,                                                                                               \
              Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, Kokkos::LayoutLeft,                  \
                           Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                \
              Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, EXECSPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, 1,  \
              true, ETI_SPEC_AVAIL> {                                                                                  \
    using mag_type = typename KokkosKernels::ArithTraits<SCALAR>::mag_type;                                            \
    using RV = Kokkos::View<mag_type, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using XV = Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, EXECSPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;    \
    using size_type = typename XV::size_type;                                                                          \
                                                                                                                       \
    static void nrm1(const EXECSPACE& space, RV& R, const XV& X) {                                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm1[TPL_BLAS," #SCALAR "]");                                         \
      const size_type numElems = X.extent(0);                                                                          \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                                \
        nrm1_print_specialization<RV, XV>();                                                                           \
        int N   = numElems;                                                                                            \
        int one = 1;                                                                                                   \
        if constexpr (KokkosKernels::ArithTraits<SCALAR>::is_complex) {                                                \
          R() = HostBlas<std::complex<mag_type>>::asum(N, reinterpret_cast<const std::complex<mag_type>*>(X.data()),   \
                                                       one);                                                           \
        } else {                                                                                                       \
          R() = HostBlas<SCALAR>::asum(N, X.data(), one);                                                              \
        }                                                                                                              \
      } else {                                                                                                         \
        Nrm1<EXECSPACE, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm1(space, R, X);                                          \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#if defined(KOKKOS_ENABLE_SERIAL)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(float, Kokkos::Serial)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(double, Kokkos::Serial)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(Kokkos::complex<float>, Kokkos::Serial)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(Kokkos::complex<double>, Kokkos::Serial)
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(float, Kokkos::OpenMP)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(double, Kokkos::OpenMP)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(Kokkos::complex<float>, Kokkos::OpenMP)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(Kokkos::complex<double>, Kokkos::OpenMP)
#endif

#if defined(KOKKOS_ENABLE_THREADS)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(float, Kokkos::Threads)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(double, Kokkos::Threads)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(Kokkos::complex<float>, Kokkos::Threads)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(Kokkos::complex<double>, Kokkos::Threads)
#endif

}  // namespace Impl
}  // namespace KokkosBlas

#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

template <class ExecutionSpace, class XViewType, class RViewType>
void cublasAsumWrapper(const ExecutionSpace& space, RViewType& R, const XViewType& X) {
  using XScalar = typename XViewType::non_const_value_type;

  nrm1_print_specialization<RViewType, XViewType>();
  const int N                            = static_cast<int>(X.extent(0));
  constexpr int one                      = 1;
  KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();

  KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));
  if constexpr (std::is_same_v<XScalar, float>) {
    KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSasum(s.handle, N, X.data(), one, R.data()));
  }
  if constexpr (std::is_same_v<XScalar, double>) {
    KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasDasum(s.handle, N, X.data(), one, R.data()));
  }
  if constexpr (std::is_same_v<XScalar, Kokkos::complex<float>>) {
    KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(
        cublasScasum(s.handle, N, reinterpret_cast<const cuComplex*>(X.data()), one, R.data()));
  }
  if constexpr (std::is_same_v<XScalar, Kokkos::complex<double>>) {
    KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(
        cublasDzasum(s.handle, N, reinterpret_cast<const cuDoubleComplex*>(X.data()), one, R.data()));
  }
  KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, NULL));
}

#define KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_CUBLAS(SCALAR)                                                                  \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct Nrm1<Kokkos::Cuda,                                                                                            \
              Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, Kokkos::LayoutLeft,                  \
                           Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                \
              Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,  \
              1, true, ETI_SPEC_AVAIL> {                                                                               \
    using execution_space = Kokkos::Cuda;                                                                              \
    using RV              = Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, Kokkos::LayoutLeft,    \
                            Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                  \
    using XV = Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using size_type = typename XV::size_type;                                                                          \
                                                                                                                       \
    static void nrm1(const execution_space& space, RV& R, const XV& X) {                                               \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm1[TPL_CUBLAS," #SCALAR "]");                                       \
      const size_type numElems = X.extent(0);                                                                          \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                                \
        cublasAsumWrapper(space, R, X);                                                                                \
      } else {                                                                                                         \
        Nrm1<execution_space, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm1(space, R, X);                                    \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_CUBLAS(float)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_CUBLAS(double)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<float>)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<double>)

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

template <class ExecutionSpace, class XViewType, class RViewType>
void rocblasAsumWrapper(const ExecutionSpace& space, RViewType& R, const XViewType& X) {
  using XScalar = typename XViewType::non_const_value_type;

  nrm1_print_specialization<RViewType, XViewType>();
  const int N                           = static_cast<int>(X.extent(0));
  constexpr int one                     = 1;
  KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();

  KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));
  if constexpr (std::is_same_v<XScalar, float>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_sasum(s.handle, N, X.data(), one, R.data()));
  }
  if constexpr (std::is_same_v<XScalar, double>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_dasum(s.handle, N, X.data(), one, R.data()));
  }
  if constexpr (std::is_same_v<XScalar, Kokkos::complex<float>>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(
        rocblas_scasum(s.handle, N, reinterpret_cast<const rocblas_float_complex*>(X.data()), one, R.data()));
  }
  if constexpr (std::is_same_v<XScalar, Kokkos::complex<double>>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(
        rocblas_dzasum(s.handle, N, reinterpret_cast<const rocblas_double_complex*>(X.data()), one, R.data()));
  }
  KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));
}

#define KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_ROCBLAS(SCALAR)                                                                \
  template <bool ETI_SPEC_AVAIL>                                                                                      \
  struct Nrm1<Kokkos::HIP,                                                                                            \
              Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, Kokkos::LayoutLeft,                 \
                           Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                               \
              Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,  \
              1, true, ETI_SPEC_AVAIL> {                                                                              \
    using RV = Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, Kokkos::LayoutLeft,                \
                            Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                              \
    using XV = Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using size_type = typename XV::size_type;                                                                         \
                                                                                                                      \
    static void nrm1(const Kokkos::HIP& space, RV& R, const XV& X) {                                                  \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm1[TPL_ROCBLAS," #SCALAR "]");                                     \
      const size_type numElems = X.extent(0);                                                                         \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                               \
        rocblasAsumWrapper(space, R, X);                                                                              \
      } else {                                                                                                        \
        Nrm1<Kokkos::HIP, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm1(space, R, X);                                       \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_ROCBLAS(float)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_ROCBLAS(double)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_ROCBLAS(Kokkos::complex<float>)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_ROCBLAS(Kokkos::complex<double>)

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCBLAS

// oneMKL
#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL

#if defined(KOKKOS_ENABLE_SYCL)

#include <KokkosBlas_tpl_spec.hpp>
#include <oneapi/mkl/blas.hpp>

namespace KokkosBlas {
namespace Impl {

template <class ExecutionSpace, class XViewType, class RViewType>
void onemklAsumWrapper(const ExecutionSpace& space, RViewType& R, const XViewType& X) {
  using XScalar  = typename XViewType::non_const_value_type;
  using KAT_X    = KokkosKernels::ArithTraits<XScalar>;
  using layout_t = typename XViewType::array_layout;

  const std::int64_t N = static_cast<std::int64_t>(X.extent(0));

  // Create temp view on device to store the result
  Kokkos::View<typename KokkosKernels::ArithTraits<XScalar>::mag_type, typename XViewType::memory_space> res(
      "sycl asum result");

  // Decide to call row_major or column_major function
  if constexpr (std::is_same_v<Kokkos::LayoutRight, layout_t>) {
    if constexpr (KAT_X::is_complex) {
      oneapi::mkl::blas::row_major::asum(space.sycl_queue(), N,
                                         reinterpret_cast<const std::complex<typename KAT_X::mag_type>*>(X.data()), 1,
                                         res.data());
    } else {
      oneapi::mkl::blas::row_major::asum(space.sycl_queue(), N, X.data(), 1, res.data());
    }
  } else {
    if constexpr (KAT_X::is_complex) {
      oneapi::mkl::blas::column_major::asum(space.sycl_queue(), N,
                                            reinterpret_cast<const std::complex<typename KAT_X::mag_type>*>(X.data()),
                                            1, res.data());
    } else {
      oneapi::mkl::blas::column_major::asum(space.sycl_queue(), X.extent_int(0), X.data(), 1, res.data());
    }
  }
  // Bring result back to host
  Kokkos::deep_copy(space, R, res);
}

#define KOKKOSBLAS1_NRM1_ONEMKL(SCALAR)                                                                                \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct Nrm1<Kokkos::SYCL,                                                                                            \
              Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, Kokkos::LayoutLeft,                  \
                           Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                \
              Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::SYCL, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,  \
              1, true, ETI_SPEC_AVAIL> {                                                                               \
    using execution_space = Kokkos::SYCL;                                                                              \
    using RV              = Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, Kokkos::LayoutLeft,    \
                            Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                  \
    using XV = Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::SYCL, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using size_type = typename XV::size_type;                                                                          \
                                                                                                                       \
    static void nrm1(const execution_space& space, RV& R, const XV& X) {                                               \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm1[TPL_ONEMKL," #SCALAR "]");                                       \
      const size_type numElems = X.extent(0);                                                                          \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                                \
        onemklAsumWrapper(space, R, X);                                                                                \
      } else {                                                                                                         \
        Nrm1<execution_space, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrm1(space, R, X);                                    \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_NRM1_ONEMKL(float)
KOKKOSBLAS1_NRM1_ONEMKL(double)
KOKKOSBLAS1_NRM1_ONEMKL(Kokkos::complex<float>)
KOKKOSBLAS1_NRM1_ONEMKL(Kokkos::complex<double>)

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOS_ENABLE_SYCL
#endif  // KOKKOSKERNELS_ENABLE_TPL_MKL

#endif
