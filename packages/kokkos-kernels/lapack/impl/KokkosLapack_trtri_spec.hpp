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
#ifndef KOKKOSLAPACK_TRTRI_SPEC_HPP_
#define KOKKOSLAPACK_TRTRI_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosLapack_trtri_impl.hpp>
#endif

namespace KokkosLapack {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class RVIT, class AVIT>
struct trtri_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosLapack

//
// This Macros provides the ETI specialization of trtri, currently not
// available.
//
#define KOKKOSLAPACK_TRTRI_ETI_SPEC_AVAIL(SCALAR, LAYOUTA, EXEC_SPACE, MEM_SPACE)                          \
  template <>                                                                                              \
  struct trtri_eti_spec_avail<                                                                             \
      Kokkos::View<int, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                               \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                           \
    enum : bool { value = true };                                                                          \
  };

// Include the actual specialization declarations
#include <KokkosLapack_trtri_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosLapack_trtri_eti_spec_avail.hpp>

namespace KokkosLapack {
namespace Impl {

//
// trtri
//

// Unification layer
template <class RVIT, class AVIT, bool tpl_spec_avail = trtri_tpl_spec_avail<RVIT, AVIT>::value,
          bool eti_spec_avail = trtri_eti_spec_avail<RVIT, AVIT>::value>
struct TRTRI {
  static void trtri(const RVIT& R, const char uplo[], const char diag[], const AVIT& A);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
template <class RVIT, class AVIT>
struct TRTRI<RVIT, AVIT, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void trtri(const RVIT& R, const char uplo[], const char diag[], const AVIT& A) {
    static_assert(Kokkos::is_view<AVIT>::value, "AVIT must be a Kokkos::View.");
    static_assert(static_cast<int>(AVIT::rank) == 2, "AVIT must have rank 2.");

    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosLapack::trtri[ETI]"
                                                                     : "KokkosLapack::trtri[noETI]");

    typename AVIT::HostMirror host_A = Kokkos::create_mirror_view(A);
    typename RVIT::HostMirror host_R = Kokkos::create_mirror_view(R);

    Kokkos::deep_copy(host_A, A);

    SerialTrtri_Invoke<typename RVIT::HostMirror, typename AVIT::HostMirror>(R, uplo, diag, host_A);

    Kokkos::deep_copy(A, host_A);

    Kokkos::Profiling::popRegion();
  }
};
#endif  //! defined(KOKKOSKERNELS_ETI_ONLY) ||
        //! KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

}  // namespace Impl
}  // namespace KokkosLapack

//
// These Macros are only included when we are not compiling libkokkoskernels but
// are auto generating files. These macros provide the explicit instantiation
// declaration and definition of TRTRI, potentially reducing user code size. The
// "extern template" skips the implicit instatiation step ensuring that the
// callers code uses this explicit instantiation definition of TRTRI.
//
#define KOKKOSLAPACK_TRTRI_ETI_SPEC_DECL(SCALAR, LAYOUTA, EXEC_SPACE, MEM_SPACE)                           \
  extern template struct TRTRI<                                                                            \
      Kokkos::View<int, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                               \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                              \
      false, true>;

#define KOKKOSLAPACK_TRTRI_ETI_SPEC_INST(SCALAR, LAYOUTA, EXEC_SPACE, MEM_SPACE)                           \
  template struct TRTRI<                                                                                   \
      Kokkos::View<int, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                               \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                              \
      false, true>;

#include <KokkosLapack_trtri_tpl_spec_decl.hpp>

#endif  // KOKKOSLAPACK_TRTRI_SPEC_HPP_
