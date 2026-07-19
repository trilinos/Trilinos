// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_POTRS_TPL_SPEC_DECL_HPP_
#define KOKKOSLAPACK_POTRS_TPL_SPEC_DECL_HPP_

#include <KokkosKernels_ArithTraits.hpp>
#include <sstream>

namespace KokkosLapack {
namespace Impl {
template <class AViewType, class BViewType>
inline void potrs_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSOLVER
  printf("KokkosLapack::potrs<> TPL cuSolver specialization for < %s, %s >\n", typeid(AViewType).name(),
         typeid(BViewType).name());
#elif defined(KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER)
  printf("KokkosLapack::potrs<> TPL rocSolver specialization for < %s, %s >\n", typeid(AViewType).name(),
         typeid(BViewType).name());
#elif defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK)
  printf("KokkosLapack::potrs<> TPL Lapack specialization for < %s, %s >\n", typeid(AViewType).name(),
         typeid(BViewType).name());
#endif
#endif
}
}  // namespace Impl
}  // namespace KokkosLapack

// Generic Host side LAPACK (could be MKL or whatever)
#if defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK) || defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
#include <KokkosLapack_Host_tpl.hpp>

namespace KokkosLapack {
namespace Impl {

template <class ExecutionSpace, class AViewType, class BViewType>
void lapackPotrsWrapper(const ExecutionSpace& /* space */, const char uplo[], const AViewType& A, BViewType& B) {
  using Scalar   = typename AViewType::non_const_value_type;
  using Layout_t = typename AViewType::array_layout;
  static_assert(std::is_same_v<Layout_t, Kokkos::LayoutLeft>,
                "KokkosLapack - potrs: A needs to have a Kokkos::LayoutLeft");

  const int n    = static_cast<int>(A.extent(0));
  const int nrhs = static_cast<int>(B.extent(1));
  const int lda  = static_cast<int>(A.stride(1));
  const int ldb  = static_cast<int>(B.stride(1));
  int info       = 0;
  if constexpr (KokkosKernels::ArithTraits<Scalar>::is_complex) {
    using MagType = typename KokkosKernels::ArithTraits<Scalar>::mag_type;
    info          = HostLapack<std::complex<MagType>>::potrs(uplo[0], n, nrhs,
                                                             reinterpret_cast<const std::complex<MagType>*>(A.data()), lda,
                                                             reinterpret_cast<std::complex<MagType>*>(B.data()), ldb);
  } else {
    info = HostLapack<Scalar>::potrs(uplo[0], n, nrhs, A.data(), lda, B.data(), ldb);
  }

  if (info != 0) {
    std::ostringstream os;
    if (info < 0)
      os << "KokkosLapack::potrs: LAPACK potrs: argument " << -info << " had an illegal value";
    else
      os << "KokkosLapack::potrs: LAPACK potrs: the leading minor of order " << info << " is not positive definite";
    Kokkos::abort(os.str().c_str());
  }
}

#define KOKKOSLAPACK_POTRS_LAPACK(SCALAR, LAYOUT, EXECSPACE, MEM_SPACE)                                                \
  template <>                                                                                                          \
  struct Potrs<                                                                                                        \
      EXECSPACE,                                                                                                       \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      true,                                                                                                            \
      potrs_eti_spec_avail<EXECSPACE,                                                                                  \
                           Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>,                  \
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                      \
                           Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>,                        \
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                            \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>,                       \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                           \
    using BViewType =                                                                                                  \
        Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    static void potrs(const EXECSPACE& space, const char uplo[], const AViewType& A, BViewType& B) {                   \
      Kokkos::Profiling::pushRegion("KokkosLapack::potrs[TPL_LAPACK," #SCALAR "]");                                    \
      potrs_print_specialization<AViewType, BViewType>();                                                              \
      lapackPotrsWrapper(space, uplo, A, B);                                                                           \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#if defined(KOKKOS_ENABLE_SERIAL)
KOKKOSLAPACK_POTRS_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_POTRS_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_POTRS_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_POTRS_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
KOKKOSLAPACK_POTRS_LAPACK(float, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_POTRS_LAPACK(double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_POTRS_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_POTRS_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined(KOKKOS_ENABLE_THREADS)
KOKKOSLAPACK_POTRS_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_POTRS_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_POTRS_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_POTRS_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_LAPACK || KOKKOSKERNELS_ENABLE_TPL_ACCELERATE

// CUSOLVER
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSOLVER
#include "KokkosLapack_cusolver.hpp"

namespace KokkosLapack {
namespace Impl {

template <class AViewType, class BViewType>
void cusolverPotrsWrapper(const Kokkos::Cuda& space, const char uplo[], const AViewType& A, BViewType& B) {
  using memory_space = typename BViewType::memory_space;
  using Scalar       = typename AViewType::non_const_value_type;

  cublasFillMode_t uplo_t = (uplo[0] == 'U' || uplo[0] == 'u') ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER;
  const int n             = static_cast<int>(A.extent(0));
  const int nrhs          = static_cast<int>(B.extent(1));
  const int lda           = static_cast<int>(A.stride(1));
  const int ldb           = static_cast<int>(B.stride(1));

  Kokkos::View<int, memory_space> info("potrs info");
  CudaLapackSingleton& s = CudaLapackSingleton::singleton();
  KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnSetStream(s.handle, space.cuda_stream()));

  if constexpr (std::is_same_v<Scalar, float>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnSpotrs(s.handle, uplo_t, n, nrhs, A.data(), lda, B.data(), ldb, info.data()));
  }
  if constexpr (std::is_same_v<Scalar, double>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnDpotrs(s.handle, uplo_t, n, nrhs, A.data(), lda, B.data(), ldb, info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnCpotrs(s.handle, uplo_t, n, nrhs,
                                                          reinterpret_cast<const cuComplex*>(A.data()), lda,
                                                          reinterpret_cast<cuComplex*>(B.data()), ldb, info.data()));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
    KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(
        cusolverDnZpotrs(s.handle, uplo_t, n, nrhs, reinterpret_cast<const cuDoubleComplex*>(A.data()), lda,
                         reinterpret_cast<cuDoubleComplex*>(B.data()), ldb, info.data()));
  }
  KOKKOSLAPACK_IMPL_CUSOLVER_SAFE_CALL(cusolverDnSetStream(s.handle, NULL));
  space.fence();
  auto h_info = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, info);
  if (h_info() != 0) {
    std::ostringstream os;
    if (h_info() < 0)
      os << "KokkosLapack::potrs: cuSOLVER potrs: argument " << -h_info() << " had an illegal value";
    else
      os << "KokkosLapack::potrs: cuSOLVER potrs: the leading minor of order " << h_info()
         << " is not positive definite";
    Kokkos::abort(os.str().c_str());
  }
}

#define KOKKOSLAPACK_POTRS_CUSOLVER(SCALAR, LAYOUT, MEM_SPACE)                                                    \
  template <>                                                                                                     \
  struct Potrs<Kokkos::Cuda,                                                                                      \
               Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                      \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                             \
               Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                            \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                             \
               true,                                                                                              \
               potrs_eti_spec_avail<Kokkos::Cuda,                                                                 \
                                    Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>, \
                                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                        \
                                    Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,       \
                                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {              \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,               \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                      \
    using BViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEM_SPACE>,                     \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                      \
    static void potrs(const Kokkos::Cuda& space, const char uplo[], const AViewType& A, BViewType& B) {           \
      Kokkos::Profiling::pushRegion("KokkosLapack::potrs[TPL_CUSOLVER," #SCALAR "]");                             \
      potrs_print_specialization<AViewType, BViewType>();                                                         \
      cusolverPotrsWrapper(space, uplo, A, B);                                                                    \
      Kokkos::Profiling::popRegion();                                                                             \
    }                                                                                                             \
  };

KOKKOSLAPACK_POTRS_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_POTRS_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_POTRS_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_POTRS_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)

#if defined(KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE)
KOKKOSLAPACK_POTRS_CUSOLVER(float, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_POTRS_CUSOLVER(double, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_POTRS_CUSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSLAPACK_POTRS_CUSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
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

template <class AViewType, class BViewType>
void rocsolverPotrsWrapper(const Kokkos::HIP& space, const char uplo[], const AViewType& A, BViewType& B) {
  using Scalar = typename AViewType::non_const_value_type;

  rocblas_fill uplo_t = (uplo[0] == 'U' || uplo[0] == 'u') ? rocblas_fill_upper : rocblas_fill_lower;

  const int n    = static_cast<int>(A.extent(0));
  const int nrhs = static_cast<int>(B.extent(1));
  const int lda  = static_cast<int>(A.stride(1));
  const int ldb  = static_cast<int>(B.stride(1));

  const rocblas_int n_    = static_cast<rocblas_int>(n);
  const rocblas_int nrhs_ = static_cast<rocblas_int>(nrhs);
  const rocblas_int lda_  = static_cast<rocblas_int>(lda);
  const rocblas_int ldb_  = static_cast<rocblas_int>(ldb);

  KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();
  KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));

  if constexpr (std::is_same_v<Scalar, float>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(
        rocsolver_spotrs(s.handle, uplo_t, n_, nrhs_, const_cast<float*>(A.data()), lda_, B.data(), ldb_));
  }
  if constexpr (std::is_same_v<Scalar, double>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(
        rocsolver_dpotrs(s.handle, uplo_t, n_, nrhs_, const_cast<double*>(A.data()), lda_, B.data(), ldb_));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(
        rocsolver_cpotrs(s.handle, uplo_t, n_, nrhs_,
                         const_cast<rocblas_float_complex*>(reinterpret_cast<const rocblas_float_complex*>(A.data())),
                         lda_, reinterpret_cast<rocblas_float_complex*>(B.data()), ldb_));
  }
  if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(
        rocsolver_zpotrs(s.handle, uplo_t, n_, nrhs_,
                         const_cast<rocblas_double_complex*>(reinterpret_cast<const rocblas_double_complex*>(A.data())),
                         lda_, reinterpret_cast<rocblas_double_complex*>(B.data()), ldb_));
  }
  KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));
  space.fence();
}

#define KOKKOSLAPACK_POTRS_ROCSOLVER(SCALAR, LAYOUT, MEM_SPACE)                                                        \
  template <>                                                                                                          \
  struct Potrs<                                                                                                        \
      Kokkos::HIP,                                                                                                     \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                                     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      true,                                                                                                            \
      potrs_eti_spec_avail<Kokkos::HIP,                                                                                \
                           Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                \
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                      \
                           Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                      \
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>>::value> {                            \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                     \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                           \
    using BViewType = Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::HIP, MEM_SPACE>,                           \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                           \
    static void potrs(const Kokkos::HIP& space, const char uplo[], const AViewType& A, BViewType& B) {                 \
      Kokkos::Profiling::pushRegion("KokkosLapack::potrs[TPL_ROCSOLVER," #SCALAR "]");                                 \
      potrs_print_specialization<AViewType, BViewType>();                                                              \
      rocsolverPotrsWrapper(space, uplo, A, B);                                                                        \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSLAPACK_POTRS_ROCSOLVER(float, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_POTRS_ROCSOLVER(double, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_POTRS_ROCSOLVER(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSLAPACK_POTRS_ROCSOLVER(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIPSpace)

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER

#endif  // KOKKOSLAPACK_POTRS_TPL_SPEC_DECL_HPP_
