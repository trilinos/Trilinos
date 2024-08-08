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

#ifndef KOKKOSBLAS1_NRMINF_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS1_NRMINF_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class AV, class XMV, int Xrank = XMV::rank>
struct nrminf_tpl_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

namespace KokkosBlas {
namespace Impl {

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
// double
#define KOKKOSBLAS1_NRMINF_TPL_SPEC_AVAIL_BLAS(SCALAR, LAYOUT, MEMSPACE)                                          \
  template <class ExecSpace>                                                                                      \
  struct nrminf_tpl_spec_avail<ExecSpace,                                                                         \
                               Kokkos::View<typename Kokkos::Details::InnerProductSpaceTraits<SCALAR>::mag_type,  \
                                            LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                               Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,           \
                                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                            \
                               1> {                                                                               \
    enum : bool { value = true };                                                                                 \
  };

KOKKOSBLAS1_NRMINF_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS1_NRMINF_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS1_NRMINF_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS1_NRMINF_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HostSpace)

#endif

}  // namespace Impl
}  // namespace KokkosBlas
#endif
