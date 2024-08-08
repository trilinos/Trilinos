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

#ifndef KOKKOSLAPACK_SVD_TPL_SPEC_DECL_HPP_
#define KOKKOSLAPACK_SVD_TPL_SPEC_DECL_HPP_

#include "KokkosKernels_Error.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace KokkosLapack {
namespace Impl {
template <class ExecutionSpace, class AMatrix, class SVector, class UMatrix, class VMatrix>
inline void svd_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSOLVER
  if constexpr (std::is_same_v<ExecutionSpace, Kokkos::Cuda>) {
    printf(
        "KokkosLapack::svd<> TPL Cusolver specialization for < %s , %s, %s, %s "
        ">\n",
        typeid(AMatrix).name(), typeid(SVector).name(), typeid(UMatrix).name(), typeid(VMatrix).name());
  }
#endif
#endif
}
}  // namespace Impl
}  // namespace KokkosLapack

// LAPACK
#if defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK) && !defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
#include "KokkosLapack_Host_tpl.hpp"

namespace KokkosLapack {
namespace Impl {

template <class ExecutionSpace, class AMatrix, class SVector, class UMatrix, class VMatrix>
void lapackSvdWrapper(const ExecutionSpace& /* space */, const char jobu[], const char jobvt[], const AMatrix& A,
                      const SVector& S, const UMatrix& U, const VMatrix& Vt) {
  using memory_space = typename AMatrix::memory_space;
  using Scalar       = typename AMatrix::non_const_value_type;
  using Magnitude    = typename SVector::non_const_value_type;
  using ALayout_t    = typename AMatrix::array_layout;
  using ULayout_t    = typename UMatrix::array_layout;
  using VLayout_t    = typename VMatrix::array_layout;

  static_assert(std::is_same_v<ALayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - svd: A needs to have a Kokkos::LayoutLeft");
  static_assert(std::is_same_v<ULayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - svd: U needs to have a Kokkos::LayoutLeft");
  static_assert(std::is_same_v<VLayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - svd: Vt needs to have a Kokkos::LayoutLeft");

  const int m    = A.extent_int(0);
  const int n    = A.extent_int(1);
  const int lda  = A.stride(1);
  const int ldu  = U.stride(1);
  const int ldvt = Vt.stride(1);

  int lwork = -1, info = 0;
  Kokkos::View<Magnitude*, memory_space> rwork("svd rwork buffer", 5 * Kokkos::min(m, n));
  Kokkos::View<Scalar*, memory_space> work("svd work buffer", 1);
  if constexpr (Kokkos::ArithTraits<Scalar>::is_complex) {
    HostLapack<std::complex<Magnitude>>::gesvd(
        jobu[0], jobvt[0], m, n, reinterpret_cast<std::complex<Magnitude>*>(A.data()), lda, S.data(),
        reinterpret_cast<std::complex<Magnitude>*>(U.data()), ldu,
        reinterpret_cast<std::complex<Magnitude>*>(Vt.data()), ldvt,
        reinterpret_cast<std::complex<Magnitude>*>(work.data()), lwork, rwork.data(), info);

    lwork = static_cast<int>(work(0).real());

    work = Kokkos::View<Scalar*, memory_space>("svd work buffer", lwork);
    HostLapack<std::complex<Magnitude>>::gesvd(
        jobu[0], jobvt[0], m, n, reinterpret_cast<std::complex<Magnitude>*>(A.data()), lda, S.data(),
        reinterpret_cast<std::complex<Magnitude>*>(U.data()), ldu,
        reinterpret_cast<std::complex<Magnitude>*>(Vt.data()), ldvt,
        reinterpret_cast<std::complex<Magnitude>*>(work.data()), lwork, rwork.data(), info);
  } else {
    HostLapack<Scalar>::gesvd(jobu[0], jobvt[0], m, n, A.data(), lda, S.data(), U.data(), ldu, Vt.data(), ldvt,
                              work.data(), lwork, rwork.data(), info);

    lwork = static_cast<int>(work(0));

    work = Kokkos::View<Scalar*, memory_space>("svd work buffer", lwork);
    HostLapack<Scalar>::gesvd(jobu[0], jobvt[0], m, n, A.data(), lda, S.data(), U.data(), ldu, Vt.data(), ldvt,
                              work.data(), lwork, rwork.data(), info);
  }
}

#define KOKKOSLAPACK_SVD_LAPACK(SCALAR, LAYOUT, EXEC_SPACE)                                                            \
  template <>                                                                                                          \
  struct SVD<EXEC_SPACE,                                                                                               \
             Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                             \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             Kokkos::View<Kokkos::ArithTraits<SCALAR>::mag_type*, LAYOUT,                                              \
                          Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,     \
             Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                             \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                             \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             true,                                                                                                     \
             svd_eti_spec_avail<                                                                                       \
                 EXEC_SPACE,                                                                                           \
                 Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                         \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                \
                 Kokkos::View<Kokkos::ArithTraits<SCALAR>::mag_type*, LAYOUT,                                          \
                              Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
                 Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                         \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                \
                 Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                         \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                                      \
    using AMatrix = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                      \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                             \
    using SVector =                                                                                                    \
        Kokkos::View<Kokkos::ArithTraits<SCALAR>::mag_type*, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,    \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                         \
    using UMatrix = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                      \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                             \
    using VMatrix = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                      \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                             \
                                                                                                                       \
    static void svd(const EXEC_SPACE& space, const char jobu[], const char jobvt[], const AMatrix& A,                  \
                    const SVector& S, const UMatrix& U, const VMatrix& Vt) {                                           \
      Kokkos::Profiling::pushRegion("KokkosLapack::svd[TPL_LAPACK," #SCALAR "]");                                      \
      svd_print_specialization<EXEC_SPACE, AMatrix, SVector, UMatrix, VMatrix>();                                      \
                                                                                                                       \
      lapackSvdWrapper(space, jobu, jobvt, A, S, U, Vt);                                                               \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#if defined(KOKKOS_ENABLE_SERIAL)
KOKKOSLAPACK_SVD_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Serial)
KOKKOSLAPACK_SVD_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Serial)
KOKKOSLAPACK_SVD_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Serial)
KOKKOSLAPACK_SVD_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Serial)
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
KOKKOSLAPACK_SVD_LAPACK(float, Kokkos::LayoutLeft, Kokkos::OpenMP)
KOKKOSLAPACK_SVD_LAPACK(double, Kokkos::LayoutLeft, Kokkos::OpenMP)
KOKKOSLAPACK_SVD_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::OpenMP)
KOKKOSLAPACK_SVD_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::OpenMP)
#endif

#if defined(KOKKOS_ENABLE_THREADS)
KOKKOSLAPACK_SVD_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Threads)
KOKKOSLAPACK_SVD_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Threads)
KOKKOSLAPACK_SVD_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Threads)
KOKKOSLAPACK_SVD_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Threads)
#endif

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_LAPACK

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include "mkl.h"

namespace KokkosLapack {
namespace Impl {

template <class ExecutionSpace, class AMatrix, class SVector, class UMatrix, class VMatrix>
void mklSvdWrapper(const ExecutionSpace& /* space */, const char jobu[], const char jobvt[], const AMatrix& A,
                   const SVector& S, const UMatrix& U, const VMatrix& Vt) {
  using memory_space = typename AMatrix::memory_space;
  using Scalar       = typename AMatrix::non_const_value_type;
  using Magnitude    = typename SVector::non_const_value_type;
  using ALayout_t    = typename AMatrix::array_layout;
  using ULayout_t    = typename UMatrix::array_layout;
  using VLayout_t    = typename VMatrix::array_layout;

  static_assert(std::is_same_v<ALayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - svd: A needs to have a Kokkos::LayoutLeft");
  static_assert(std::is_same_v<ULayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - svd: U needs to have a Kokkos::LayoutLeft");
  static_assert(std::is_same_v<VLayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - svd: Vt needs to have a Kokkos::LayoutLeft");

  const lapack_int m    = A.extent_int(0);
  const lapack_int n    = A.extent_int(1);
  const lapack_int lda  = A.stride(1);
  const lapack_int ldu  = U.stride(1);
  const lapack_int ldvt = Vt.stride(1);

  Kokkos::View<Magnitude*, memory_space> rwork("svd rwork buffer", Kokkos::min(m, n) - 1);
  lapack_int ret = 0;
  if constexpr (std::is_same_v<Scalar, float>) {
    ret = LAPACKE_sgesvd(LAPACK_COL_MAJOR, jobu[0], jobvt[0], m, n, A.data(), lda, S.data(), U.data(), ldu, Vt.data(),
                         ldvt, rwork.data());
  }
  if constexpr (std::is_same_v<Scalar, double>) {
    ret = LAPACKE_dgesvd(LAPACK_COL_MAJOR, jobu[0], jobvt[0], m, n, A.data(), lda, S.data(), U.data(), ldu, Vt.data(),
                         ldvt, rwork.data());
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
    ret = LAPACKE_cgesvd(LAPACK_COL_MAJOR, jobu[0], jobvt[0], m, n, reinterpret_cast<lapack_complex_float*>(A.data()),
                         lda, S.data(), reinterpret_cast<lapack_complex_float*>(U.data()), ldu,
                         reinterpret_cast<lapack_complex_float*>(Vt.data()), ldvt, rwork.data());
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
    ret = LAPACKE_zgesvd(LAPACK_COL_MAJOR, jobu[0], jobvt[0], m, n, reinterpret_cast<lapack_complex_double*>(A.data()),
                         lda, S.data(), reinterpret_cast<lapack_complex_double*>(U.data()), ldu,
                         reinterpret_cast<lapack_complex_double*>(Vt.data()), ldvt, rwork.data());
  }

  if (ret != 0) {
    std::ostringstream os;
    os << "KokkosLapack::svd: MKL failed with return value: " << ret << "\n";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
}

#define KOKKOSLAPACK_SVD_MKL(SCALAR, LAYOUT, EXEC_SPACE)                                                               \
  template <>                                                                                                          \
  struct SVD<EXEC_SPACE,                                                                                               \
             Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                             \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             Kokkos::View<Kokkos::ArithTraits<SCALAR>::mag_type*, LAYOUT,                                              \
                          Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,     \
             Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                             \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                             \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             true,                                                                                                     \
             svd_eti_spec_avail<                                                                                       \
                 EXEC_SPACE,                                                                                           \
                 Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                         \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                \
                 Kokkos::View<Kokkos::ArithTraits<SCALAR>::mag_type*, LAYOUT,                                          \
                              Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
                 Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                         \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                \
                 Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                         \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                                      \
    using AMatrix = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                      \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                             \
    using SVector =                                                                                                    \
        Kokkos::View<Kokkos::ArithTraits<SCALAR>::mag_type*, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,    \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                         \
    using UMatrix = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                      \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                             \
    using VMatrix = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, Kokkos::HostSpace>,                      \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                             \
                                                                                                                       \
    static void svd(const EXEC_SPACE& space, const char jobu[], const char jobvt[], const AMatrix& A,                  \
                    const SVector& S, const UMatrix& U, const VMatrix& Vt) {                                           \
      Kokkos::Profiling::pushRegion("KokkosLapack::svd[TPL_LAPACK," #SCALAR "]");                                      \
      svd_print_specialization<EXEC_SPACE, AMatrix, SVector, UMatrix, VMatrix>();                                      \
                                                                                                                       \
      mklSvdWrapper(space, jobu, jobvt, A, S, U, Vt);                                                                  \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#if defined(KOKKOS_ENABLE_SERIAL)
KOKKOSLAPACK_SVD_MKL(float, Kokkos::LayoutLeft, Kokkos::Serial)
KOKKOSLAPACK_SVD_MKL(double, Kokkos::LayoutLeft, Kokkos::Serial)
KOKKOSLAPACK_SVD_MKL(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Serial)
KOKKOSLAPACK_SVD_MKL(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Serial)
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
KOKKOSLAPACK_SVD_MKL(float, Kokkos::LayoutLeft, Kokkos::OpenMP)
KOKKOSLAPACK_SVD_MKL(double, Kokkos::LayoutLeft, Kokkos::OpenMP)
KOKKOSLAPACK_SVD_MKL(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::OpenMP)
KOKKOSLAPACK_SVD_MKL(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::OpenMP)
#endif

#if defined(KOKKOS_ENABLE_THREADS)
KOKKOSLAPACK_SVD_MKL(float, Kokkos::LayoutLeft, Kokkos::Threads)
KOKKOSLAPACK_SVD_MKL(double, Kokkos::LayoutLeft, Kokkos::Threads)
KOKKOSLAPACK_SVD_MKL(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Threads)
KOKKOSLAPACK_SVD_MKL(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Threads)
#endif

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_MKL

// CUSOLVER
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSOLVER
#include "KokkosLapack_cusolver.hpp"

namespace KokkosLapack {
namespace Impl {

template <class ExecutionSpace, class AMatrix, class SVector, class UMatrix, class VMatrix>
void cusolverSvdWrapper(const ExecutionSpace& space, const char jobu[], const char jobvt[], const AMatrix& A,
                        const SVector& S, const UMatrix& U, const VMatrix& Vt) {
  using memory_space = typename AMatrix::memory_space;
  using Scalar       = typename AMatrix::non_const_value_type;
  using Magnitude    = typename SVector::non_const_value_type;
  using ALayout_t    = typename AMatrix::array_layout;
  using ULayout_t    = typename UMatrix::array_layout;
  using VLayout_t    = typename VMatrix::array_layout;

  static_assert(std::is_same_v<ALayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - svd: A needs to have a Kokkos::LayoutLeft");
  static_assert(std::is_same_v<ULayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - svd: U needs to have a Kokkos::LayoutLeft");
  static_assert(std::is_same_v<VLayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - svd: Vt needs to have a Kokkos::LayoutLeft");

  const int m    = A.extent_int(0);
  const int n    = A.extent_int(1);
  const int lda  = A.stride(1);
  const int ldu  = U.stride(1);
  const int ldvt = Vt.stride(1);

  int lwork = 0;
  Kokkos::View<int, memory_space> info("svd info");
  Kokkos::View<Magnitude*, memory_space> rwork("svd rwork buffer", Kokkos::min(m, n) - 1);

  CudaLapackSingleton& s = CudaLapackSingleton::singleton();
  KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnSetStream(s.handle, space.cuda_stream()));
  if constexpr (std::is_same_v<Scalar, float>) {
    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnSgesvd_bufferSize(s.handle, m, n, &lwork));
    Kokkos::View<Scalar*, memory_space> work("svd work buffer", lwork);

    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnSgesvd(s.handle, jobu[0], jobvt[0], m, n, A.data(), lda, S.data(),
                                                    U.data(), ldu, Vt.data(), ldvt, work.data(), lwork, rwork.data(),
                                                    info.data()));
  }
  if constexpr (std::is_same_v<Scalar, double>) {
    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnDgesvd_bufferSize(s.handle, m, n, &lwork));
    Kokkos::View<Scalar*, memory_space> work("svd work buffer", lwork);

    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnDgesvd(s.handle, jobu[0], jobvt[0], m, n, A.data(), lda, S.data(),
                                                    U.data(), ldu, Vt.data(), ldvt, work.data(), lwork, rwork.data(),
                                                    info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnCgesvd_bufferSize(s.handle, m, n, &lwork));
    Kokkos::View<Scalar*, memory_space> work("svd work buffer", lwork);

    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(
        cusolverDnCgesvd(s.handle, jobu[0], jobvt[0], m, n, reinterpret_cast<cuComplex*>(A.data()), lda, S.data(),
                         reinterpret_cast<cuComplex*>(U.data()), ldu, reinterpret_cast<cuComplex*>(Vt.data()), ldvt,
                         reinterpret_cast<cuComplex*>(work.data()), lwork, rwork.data(), info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnZgesvd_bufferSize(s.handle, m, n, &lwork));
    Kokkos::View<Scalar*, memory_space> work("svd work buffer", lwork);

    KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnZgesvd(
        s.handle, jobu[0], jobvt[0], m, n, reinterpret_cast<cuDoubleComplex*>(A.data()), lda, S.data(),
        reinterpret_cast<cuDoubleComplex*>(U.data()), ldu, reinterpret_cast<cuDoubleComplex*>(Vt.data()), ldvt,
        reinterpret_cast<cuDoubleComplex*>(work.data()), lwork, rwork.data(), info.data()));
  }
  KOKKOS_CUSOLVER_SAFE_CALL_IMPL(cusolverDnSetStream(s.handle, NULL));
}

#define KOKKOSLAPACK_SVD_CUSOLVER(SCALAR, LAYOUT, MEM_SPACE)                                                           \
  template <>                                                                                                          \
  struct SVD<Kokkos::Cuda,                                                                                             \
             Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                                   \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             Kokkos::View<Kokkos::ArithTraits<SCALAR>::mag_type*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,     \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                                   \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                                   \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             true,                                                                                                     \
             svd_eti_spec_avail<                                                                                       \
                 Kokkos::Cuda,                                                                                         \
                 Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                               \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                \
                 Kokkos::View<Kokkos::ArithTraits<SCALAR>::mag_type*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                \
                 Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                               \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                \
                 Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                               \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                                      \
    using AMatrix = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                            \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                             \
    using SVector = Kokkos::View<Kokkos::ArithTraits<SCALAR>::mag_type*, LAYOUT,                                       \
                                 Kokkos::Device<Kokkos::Cuda, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;    \
    using UMatrix = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                            \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                             \
    using VMatrix = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                            \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                             \
                                                                                                                       \
    static void svd(const Kokkos::Cuda& space, const char jobu[], const char jobvt[], const AMatrix& A,                \
                    const SVector& S, const UMatrix& U, const VMatrix& Vt) {                                           \
      Kokkos::Profiling::pushRegion("KokkosLapack::svd[TPL_CUSOLVER," #SCALAR "]");                                    \
      svd_print_specialization<Kokkos::Cuda, AMatrix, SVector, UMatrix, VMatrix>();                                    \
                                                                                                                       \
      cusolverSvdWrapper(space, jobu, jobvt, A, S, U, Vt);                                                             \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSLAPACK_SVD_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_SVD_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_SVD_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_SVD_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)

#if defined(KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE)
KOKKOSLAPACK_SVD_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_SVD_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_SVD_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_SVD_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSOLVER

// ROCSOLVER
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER
#include <KokkosBlas_tpl_spec.hpp>
#include <rocsolver/rocsolver.h>

namespace KokkosLapack {
namespace Impl {

template <class ExecutionSpace, class AMatrix, class SVector, class UMatrix, class VMatrix>
void rocsolverSvdWrapper(const ExecutionSpace& space, const char jobu[], const char jobvt[], const AMatrix& A,
                         const SVector& S, const UMatrix& U, const VMatrix& Vt) {
  using memory_space = typename AMatrix::memory_space;
  using Scalar       = typename AMatrix::non_const_value_type;
  using Magnitude    = typename SVector::non_const_value_type;
  using ALayout_t    = typename AMatrix::array_layout;
  using ULayout_t    = typename UMatrix::array_layout;
  using VLayout_t    = typename VMatrix::array_layout;

  static_assert(std::is_same_v<ALayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - svd: A needs to have a Kokkos::LayoutLeft");
  static_assert(std::is_same_v<ULayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - svd: U needs to have a Kokkos::LayoutLeft");
  static_assert(std::is_same_v<VLayout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - svd: Vt needs to have a Kokkos::LayoutLeft");

  const rocblas_int m    = A.extent_int(0);
  const rocblas_int n    = A.extent_int(1);
  const rocblas_int lda  = A.stride(1);
  const rocblas_int ldu  = U.stride(1);
  const rocblas_int ldvt = Vt.stride(1);

  rocblas_svect UVecMode = rocblas_svect_all;
  if ((jobu[0] == 'S') || (jobu[0] == 's')) {
    UVecMode = rocblas_svect_singular;
  } else if ((jobu[0] == 'O') || (jobu[0] == 'o')) {
    UVecMode = rocblas_svect_overwrite;
  } else if ((jobu[0] == 'N') || (jobu[0] == 'n')) {
    UVecMode = rocblas_svect_none;
  }
  rocblas_svect VVecMode = rocblas_svect_all;
  if ((jobvt[0] == 'S') || (jobvt[0] == 's')) {
    VVecMode = rocblas_svect_singular;
  } else if ((jobvt[0] == 'O') || (jobvt[0] == 'o')) {
    VVecMode = rocblas_svect_overwrite;
  } else if ((jobvt[0] == 'N') || (jobvt[0] == 'n')) {
    VVecMode = rocblas_svect_none;
  }

  const rocblas_workmode WorkMode = rocblas_outofplace;

  Kokkos::View<rocblas_int, memory_space> info("svd info");
  Kokkos::View<Magnitude*, memory_space> rwork("svd rwork buffer", Kokkos::min(m, n) - 1);

  KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();
  KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));
  if constexpr (std::is_same_v<Scalar, float>) {
    KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocsolver_sgesvd(s.handle, UVecMode, VVecMode, m, n, A.data(), lda, S.data(),
                                                   U.data(), ldu, Vt.data(), ldvt, rwork.data(), WorkMode,
                                                   info.data()));
  }
  if constexpr (std::is_same_v<Scalar, double>) {
    KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocsolver_dgesvd(s.handle, UVecMode, VVecMode, m, n, A.data(), lda, S.data(),
                                                   U.data(), ldu, Vt.data(), ldvt, rwork.data(), WorkMode,
                                                   info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
    KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocsolver_cgesvd(
        s.handle, UVecMode, VVecMode, m, n, reinterpret_cast<rocblas_float_complex*>(A.data()), lda, S.data(),
        reinterpret_cast<rocblas_float_complex*>(U.data()), ldu, reinterpret_cast<rocblas_float_complex*>(Vt.data()),
        ldvt, rwork.data(), WorkMode, info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
    KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocsolver_zgesvd(
        s.handle, UVecMode, VVecMode, m, n, reinterpret_cast<rocblas_double_complex*>(A.data()), lda, S.data(),
        reinterpret_cast<rocblas_double_complex*>(U.data()), ldu, reinterpret_cast<rocblas_double_complex*>(Vt.data()),
        ldvt, rwork.data(), WorkMode, info.data()));
  }
  KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));
}

#define KOKKOSLAPACK_SVD_ROCSOLVER(SCALAR, LAYOUT, MEM_SPACE)                                                          \
  template <>                                                                                                          \
  struct SVD<                                                                                                          \
      Kokkos::HIP,                                                                                                     \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<Kokkos::ArithTraits<SCALAR>::mag_type*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,             \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      true,                                                                                                            \
      svd_eti_spec_avail<                                                                                              \
          Kokkos::HIP,                                                                                                 \
          Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                                       \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                       \
          Kokkos::View<Kokkos::ArithTraits<SCALAR>::mag_type*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,         \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                       \
          Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                                       \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                       \
          Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                                       \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                                             \
    using AMatrix = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                             \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                             \
    using SVector = Kokkos::View<Kokkos::ArithTraits<SCALAR>::mag_type*, LAYOUT,                                       \
                                 Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    using UMatrix = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                             \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                             \
    using VMatrix = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                             \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                             \
                                                                                                                       \
    static void svd(const Kokkos::HIP& space, const char jobu[], const char jobvt[], const AMatrix& A,                 \
                    const SVector& S, const UMatrix& U, const VMatrix& Vt) {                                           \
      Kokkos::Profiling::pushRegion("KokkosLapack::svd[TPL_ROCSOLVER," #SCALAR "]");                                   \
      svd_print_specialization<Kokkos::HIP, AMatrix, SVector, UMatrix, VMatrix>();                                     \
                                                                                                                       \
      rocsolverSvdWrapper(space, jobu, jobvt, A, S, U, Vt);                                                            \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSLAPACK_SVD_ROCSOLVER(float, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_SVD_ROCSOLVER(double, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_SVD_ROCSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_SVD_ROCSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIPSpace)

#if defined(KOKKOSKERNELS_INST_MEMSPACE_HIPMANAGEDSPACE)
KOKKOSLAPACK_SVD_ROCSOLVER(float, Kokkos::LayoutLeft, Kokkos::HIPManagedSpace)
KOKKOSLAPACK_SVD_ROCSOLVER(double, Kokkos::LayoutLeft, Kokkos::HIPManagedSpace)
KOKKOSLAPACK_SVD_ROCSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIPManagedSpace)
KOKKOSLAPACK_SVD_ROCSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIPManagedSpace)
#endif

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER

#endif  // KOKKOSLAPACK_SVD_TPL_SPEC_DECL_HPP_
