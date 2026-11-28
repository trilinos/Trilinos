// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_NRM2_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS1_NRM2_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class RV, class XMV, int Xrank = XMV::rank>
struct nrm2_tpl_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

namespace KokkosBlas {
namespace Impl {
// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
// double
#define KOKKOSBLAS1_NRM2_TPL_SPEC_AVAIL_BLAS(SCALAR)                                                               \
  template <typename ExecSpace>                                                                                    \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                  \
  struct nrm2_tpl_spec_avail<                                                                                      \
      ExecSpace,                                                                                                   \
      Kokkos::View<typename KokkosKernels::Details::InnerProductSpaceTraits<SCALAR>::mag_type, Kokkos::LayoutLeft, \
                   Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                   \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1> {   \
    enum : bool { value = true };                                                                                  \
  };

KOKKOSBLAS1_NRM2_TPL_SPEC_AVAIL_BLAS(double)
KOKKOSBLAS1_NRM2_TPL_SPEC_AVAIL_BLAS(float)
KOKKOSBLAS1_NRM2_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>)
KOKKOSBLAS1_NRM2_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>)

#endif

#define KOKKOSBLAS1_NRM2_TPL_SPEC(SCALAR, EXEC_SPACE)                                                             \
  template <>                                                                                                     \
  struct nrm2_tpl_spec_avail<                                                                                     \
      EXEC_SPACE,                                                                                                 \
      Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, Kokkos::LayoutLeft, Kokkos::HostSpace,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                     \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, EXEC_SPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1> { \
    enum : bool { value = true };                                                                                 \
  };

#define KOKKOSBLAS1_NRM2_TPL_SPEC_AVAIL(EXEC_SPACE)             \
  KOKKOSBLAS1_NRM2_TPL_SPEC(float, EXEC_SPACE)                  \
  KOKKOSBLAS1_NRM2_TPL_SPEC(double, EXEC_SPACE)                 \
  KOKKOSBLAS1_NRM2_TPL_SPEC(Kokkos::complex<float>, EXEC_SPACE) \
  KOKKOSBLAS1_NRM2_TPL_SPEC(Kokkos::complex<double>, EXEC_SPACE)

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
KOKKOSBLAS1_NRM2_TPL_SPEC_AVAIL(Kokkos::Cuda)
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
KOKKOSBLAS1_NRM2_TPL_SPEC_AVAIL(Kokkos::HIP)
#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_MKL) && defined(KOKKOS_ENABLE_SYCL)
KOKKOSBLAS1_NRM2_TPL_SPEC_AVAIL(Kokkos::SYCL)
#endif

}  // namespace Impl
}  // namespace KokkosBlas
#endif
