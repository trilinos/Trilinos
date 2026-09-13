// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class AV, class XMV, int Xrank = XMV::rank>
struct nrm1_tpl_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

namespace KokkosBlas {
namespace Impl {

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
// double
#define KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_BLAS(SCALAR)                                                               \
  template <typename ExecSpace>                                                                                    \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                  \
  struct nrm1_tpl_spec_avail<                                                                                      \
      ExecSpace,                                                                                                   \
      Kokkos::View<typename KokkosKernels::Details::InnerProductSpaceTraits<SCALAR>::mag_type, Kokkos::LayoutLeft, \
                   Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                   \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1> {   \
    enum : bool { value = true };                                                                                  \
  };

KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_BLAS(double)
KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_BLAS(float)
KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>)
KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>)

#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
// double
#define KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_CUBLAS(SCALAR)                                                              \
  template <>                                                                                                       \
  struct nrm1_tpl_spec_avail<                                                                                       \
      Kokkos::Cuda,                                                                                                 \
      Kokkos::View<typename KokkosKernels::Details::InnerProductSpaceTraits<SCALAR>::mag_type, Kokkos::LayoutLeft,  \
                   Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                    \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1> { \
    enum : bool { value = true };                                                                                   \
  };

KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_CUBLAS(double)
KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_CUBLAS(float)
KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>)
KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>)

#endif

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#define KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_ROCBLAS(SCALAR)                                                            \
  template <>                                                                                                      \
  struct nrm1_tpl_spec_avail<                                                                                      \
      Kokkos::HIP,                                                                                                 \
      Kokkos::View<typename KokkosKernels::Details::InnerProductSpaceTraits<SCALAR>::mag_type, Kokkos::LayoutLeft, \
                   Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                   \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1> { \
    enum : bool { value = true };                                                                                  \
  };

KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_ROCBLAS(double)
KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_ROCBLAS(float)
KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<double>)
KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<float>)

#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCBLAS

// oneMKL
#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL

#if defined(KOKKOS_ENABLE_SYCL)

#define KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_MKL_SYCL(SCALAR)                                                            \
  template <>                                                                                                       \
  struct nrm1_tpl_spec_avail<                                                                                       \
      Kokkos::SYCL,                                                                                                 \
      Kokkos::View<typename KokkosKernels::Details::InnerProductSpaceTraits<SCALAR>::mag_type, Kokkos::LayoutLeft,  \
                   Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                    \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::SYCL, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1> { \
    enum : bool { value = true };                                                                                   \
  };

KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_MKL_SYCL(double)
KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_MKL_SYCL(float)
KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_MKL_SYCL(Kokkos::complex<double>)
KOKKOSBLAS1_NRM1_TPL_SPEC_AVAIL_MKL_SYCL(Kokkos::complex<float>)

#endif  // KOKKOS_ENABLE_SYCL
#endif  // KOKKOSKERNELS_ENABLE_TPL_MKL

}  // namespace Impl
}  // namespace KokkosBlas
#endif
