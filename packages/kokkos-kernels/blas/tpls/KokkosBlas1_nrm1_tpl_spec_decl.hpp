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

#define KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(SCALAR, LAYOUT, EXECSPACE, MEMSPACE)                                     \
  template <>                                                                                                        \
  struct Nrm1<EXECSPACE,                                                                                             \
              Kokkos::View<typename Kokkos::ArithTraits<SCALAR>::mag_type, LAYOUT, Kokkos::HostSpace,                \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                 \
              Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                               \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                 \
              1, true,                                                                                               \
              nrm1_eti_spec_avail<EXECSPACE,                                                                         \
                                  Kokkos::View<typename Kokkos::ArithTraits<SCALAR>::mag_type, LAYOUT,               \
                                               Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,          \
                                  Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,           \
                                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                   \
    using mag_type  = typename Kokkos::ArithTraits<SCALAR>::mag_type;                                                \
    using RV        = Kokkos::View<mag_type, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;    \
    using XV        = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                       \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                         \
    using size_type = typename XV::size_type;                                                                        \
                                                                                                                     \
    static void nrm1(const EXECSPACE& space, RV& R, const XV& X) {                                                   \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm1[TPL_BLAS," #SCALAR "]");                                       \
      const size_type numElems = X.extent(0);                                                                        \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                              \
        nrm1_print_specialization<RV, XV>();                                                                         \
        int N   = numElems;                                                                                          \
        int one = 1;                                                                                                 \
        if constexpr (Kokkos::ArithTraits<SCALAR>::is_complex) {                                                     \
          R() = HostBlas<std::complex<mag_type>>::asum(N, reinterpret_cast<const std::complex<mag_type>*>(X.data()), \
                                                       one);                                                         \
        } else {                                                                                                     \
          R() = HostBlas<SCALAR>::asum(N, X.data(), one);                                                            \
        }                                                                                                            \
      } else {                                                                                                       \
        Nrm1<EXECSPACE, RV, XV, 1, false, nrm1_eti_spec_avail<EXECSPACE, RV, XV>::value>::nrm1(space, R, X);         \
      }                                                                                                              \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#if defined(KOKKOS_ENABLE_SERIAL)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(float, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(double, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(float, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined(KOKKOS_ENABLE_THREADS)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(float, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(double, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_BLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
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

  KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));
  if constexpr (std::is_same_v<XScalar, float>) {
    KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSasum(s.handle, N, X.data(), one, R.data()));
  }
  if constexpr (std::is_same_v<XScalar, double>) {
    KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasDasum(s.handle, N, X.data(), one, R.data()));
  }
  if constexpr (std::is_same_v<XScalar, Kokkos::complex<float>>) {
    KOKKOS_CUBLAS_SAFE_CALL_IMPL(
        cublasScasum(s.handle, N, reinterpret_cast<const cuComplex*>(X.data()), one, R.data()));
  }
  if constexpr (std::is_same_v<XScalar, Kokkos::complex<double>>) {
    KOKKOS_CUBLAS_SAFE_CALL_IMPL(
        cublasDzasum(s.handle, N, reinterpret_cast<const cuDoubleComplex*>(X.data()), one, R.data()));
  }
  KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));
}

#define KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_CUBLAS(SCALAR, LAYOUT, MEMSPACE)                                               \
  template <>                                                                                                         \
  struct Nrm1<Kokkos::Cuda,                                                                                           \
              Kokkos::View<typename Kokkos::ArithTraits<SCALAR>::mag_type, LAYOUT, Kokkos::HostSpace,                 \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                  \
              Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                             \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                  \
              1, true,                                                                                                \
              nrm1_eti_spec_avail<Kokkos::Cuda,                                                                       \
                                  Kokkos::View<typename Kokkos::ArithTraits<SCALAR>::mag_type, LAYOUT,                \
                                               Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,           \
                                  Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,         \
                                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                    \
    using execution_space = Kokkos::Cuda;                                                                             \
    using RV              = Kokkos::View<typename Kokkos::ArithTraits<SCALAR>::mag_type, LAYOUT, Kokkos::HostSpace,   \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                    \
    using XV              = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,               \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                    \
    using size_type       = typename XV::size_type;                                                                   \
                                                                                                                      \
    static void nrm1(const execution_space& space, RV& R, const XV& X) {                                              \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm1[TPL_CUBLAS," #SCALAR "]");                                      \
      const size_type numElems = X.extent(0);                                                                         \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                               \
        cublasAsumWrapper(space, R, X);                                                                               \
      } else {                                                                                                        \
        Nrm1<execution_space, RV, XV, 1, false, nrm1_eti_spec_avail<Kokkos::Cuda, RV, XV>::value>::nrm1(space, R, X); \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_CUBLAS(float, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_CUBLAS(double, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)

#if defined(KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_CUBLAS(float, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_CUBLAS(double, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

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

  KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));
  if constexpr (std::is_same_v<XScalar, float>) {
    KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_sasum(s.handle, N, X.data(), one, R.data()));
  }
  if constexpr (std::is_same_v<XScalar, double>) {
    KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_dasum(s.handle, N, X.data(), one, R.data()));
  }
  if constexpr (std::is_same_v<XScalar, Kokkos::complex<float>>) {
    KOKKOS_ROCBLAS_SAFE_CALL_IMPL(
        rocblas_scasum(s.handle, N, reinterpret_cast<const rocblas_float_complex*>(X.data()), one, R.data()));
  }
  if constexpr (std::is_same_v<XScalar, Kokkos::complex<double>>) {
    KOKKOS_ROCBLAS_SAFE_CALL_IMPL(
        rocblas_dzasum(s.handle, N, reinterpret_cast<const rocblas_double_complex*>(X.data()), one, R.data()));
  }
  KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));
}

#define KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_ROCBLAS(SCALAR, LAYOUT, MEMSPACE)                                         \
  template <>                                                                                                    \
  struct Nrm1<Kokkos::HIP,                                                                                       \
              Kokkos::View<typename Kokkos::ArithTraits<SCALAR>::mag_type, LAYOUT, Kokkos::HostSpace,            \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                             \
              Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEMSPACE>,                         \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                             \
              1, true,                                                                                           \
              nrm1_eti_spec_avail<Kokkos::HIP,                                                                   \
                                  Kokkos::View<typename Kokkos::ArithTraits<SCALAR>::mag_type, LAYOUT,           \
                                               Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,      \
                                  Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEMSPACE>,     \
                                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {               \
    using RV        = Kokkos::View<typename Kokkos::ArithTraits<SCALAR>::mag_type, LAYOUT, Kokkos::HostSpace,    \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                     \
    using XV        = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEMSPACE>,                 \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                     \
    using size_type = typename XV::size_type;                                                                    \
                                                                                                                 \
    static void nrm1(const Kokkos::HIP& space, RV& R, const XV& X) {                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm1[TPL_ROCBLAS," #SCALAR "]");                                \
      const size_type numElems = X.extent(0);                                                                    \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                          \
        rocblasAsumWrapper(space, R, X);                                                                         \
      } else {                                                                                                   \
        Nrm1<Kokkos::HIP, RV, XV, 1, false, nrm1_eti_spec_avail<Kokkos::HIP, RV, XV>::value>::nrm1(space, R, X); \
      }                                                                                                          \
      Kokkos::Profiling::popRegion();                                                                            \
    }                                                                                                            \
  };

KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_ROCBLAS(float, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_ROCBLAS(double, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_ROCBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSBLAS1_NRM1_TPL_SPEC_DECL_ROCBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIPSpace)

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
  using KAT_X    = Kokkos::ArithTraits<XScalar>;
  using layout_t = typename XViewType::array_layout;

  const std::int64_t N = static_cast<std::int64_t>(X.extent(0));

  // Create temp view on device to store the result
  Kokkos::View<typename Kokkos::ArithTraits<XScalar>::mag_type, typename XViewType::memory_space> res(
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

#define KOKKOSBLAS1_NRM1_ONEMKL(SCALAR, LAYOUT, MEMSPACE)                                                              \
  template <>                                                                                                          \
  struct Nrm1<                                                                                                         \
      Kokkos::Experimental::SYCL,                                                                                      \
      Kokkos::View<typename Kokkos::ArithTraits<SCALAR>::mag_type, LAYOUT, Kokkos::HostSpace,                          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<Kokkos::Experimental::SYCL, MEMSPACE>,                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      1, true,                                                                                                         \
      nrm1_eti_spec_avail<Kokkos::Experimental::SYCL,                                                                  \
                          Kokkos::View<typename Kokkos::ArithTraits<SCALAR>::mag_type, LAYOUT, Kokkos::HostSpace,      \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                          Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<Kokkos::Experimental::SYCL, MEMSPACE>,    \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                             \
    using execution_space = Kokkos::Experimental::SYCL;                                                                \
    using RV              = Kokkos::View<typename Kokkos::ArithTraits<SCALAR>::mag_type, LAYOUT, Kokkos::HostSpace,    \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                     \
    using XV              = Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<Kokkos::Experimental::SYCL, MEMSPACE>,  \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                     \
    using size_type       = typename XV::size_type;                                                                    \
                                                                                                                       \
    static void nrm1(const execution_space& space, RV& R, const XV& X) {                                               \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrm1[TPL_ONEMKL," #SCALAR "]");                                       \
      const size_type numElems = X.extent(0);                                                                          \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                                \
        onemklAsumWrapper(space, R, X);                                                                                \
      } else {                                                                                                         \
        Nrm1<execution_space, RV, XV, 1, false, nrm1_eti_spec_avail<Kokkos::Experimental::SYCL, RV, XV>::value>::nrm1( \
            space, R, X);                                                                                              \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_NRM1_ONEMKL(float, Kokkos::LayoutLeft, Kokkos::Experimental::SYCLDeviceUSMSpace)
KOKKOSBLAS1_NRM1_ONEMKL(double, Kokkos::LayoutLeft, Kokkos::Experimental::SYCLDeviceUSMSpace)
KOKKOSBLAS1_NRM1_ONEMKL(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Experimental::SYCLDeviceUSMSpace)
KOKKOSBLAS1_NRM1_ONEMKL(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Experimental::SYCLDeviceUSMSpace)

#if defined(KOKKOSKERNELS_INST_MEMSPACE_SYCLSHAREDSPACE)
KOKKOSBLAS1_NRM1_ONEMKL(float, Kokkos::LayoutLeft, Kokkos::Experimental::SYCLSharedUSMSpace)
KOKKOSBLAS1_NRM1_ONEMKL(double, Kokkos::LayoutLeft, Kokkos::Experimental::SYCLSharedUSMSpace)
KOKKOSBLAS1_NRM1_ONEMKL(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Experimental::SYCLSharedUSMSpace)
KOKKOSBLAS1_NRM1_ONEMKL(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Experimental::SYCLSharedUSMSpace)
#endif

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOS_ENABLE_SYCL
#endif  // KOKKOSKERNELS_ENABLE_TPL_MKL

#endif
