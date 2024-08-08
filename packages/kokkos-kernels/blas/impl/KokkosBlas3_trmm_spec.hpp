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
#ifndef KOKKOSBLAS3_TRMM_SPEC_HPP_
#define KOKKOSBLAS3_TRMM_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas3_trmm_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class AVIT, class BVIT>
struct trmm_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// This Macro is for readability of the template arguments.
//
#define KOKKOSBLAS3_TRMM_ETI_SPEC_AVAIL_LAYOUT(SCALAR, LAYOUTA, LAYOUTB, EXEC_SPACE, MEM_SPACE)           \
  template <>                                                                                             \
  struct trmm_eti_spec_avail<EXEC_SPACE,                                                                  \
                             Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
                             Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,       \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                   \
    enum : bool { value = true };                                                                         \
  };

//
// This Macros provides the ETI specialization of trmm
//
#define KOKKOSBLAS3_TRMM_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE) \
  KOKKOSBLAS3_TRMM_ETI_SPEC_AVAIL_LAYOUT(SCALAR, LAYOUT, LAYOUT, EXEC_SPACE, MEM_SPACE)

// Include the actual specialization declarations
#include <KokkosBlas3_trmm_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas3_trmm_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

//
// trmm
//

// Unification layer
template <class execution_space, class AVIT, class BVIT,
          bool tpl_spec_avail = trmm_tpl_spec_avail<execution_space, AVIT, BVIT>::value,
          bool eti_spec_avail = trmm_eti_spec_avail<execution_space, AVIT, BVIT>::value>
struct TRMM {
  static void trmm(const execution_space& space, const char side[], const char uplo[], const char trans[],
                   const char diag[], typename BVIT::const_value_type& alpha, const AVIT& A, const BVIT& B);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
template <class execution_space, class AVIT, class BVIT>
struct TRMM<execution_space, AVIT, BVIT, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void trmm(const execution_space& /*space*/, const char side[], const char uplo[], const char trans[],
                   const char diag[], typename BVIT::const_value_type& alpha, const AVIT& A, const BVIT& B) {
    static_assert(Kokkos::is_view<AVIT>::value, "AVIT must be a Kokkos::View.");
    static_assert(Kokkos::is_view<BVIT>::value, "BVIT must be a Kokkos::View.");
    static_assert(static_cast<int>(AVIT::rank) == 2, "AVIT must have rank 2.");
    static_assert(static_cast<int>(BVIT::rank) == 2, "BVIT must have rank 2.");

    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::trmm[ETI]"
                                                                     : "KokkosBlas::trmm[noETI]");

    typename AVIT::HostMirror host_A = Kokkos::create_mirror_view(A);
    typename BVIT::HostMirror host_B = Kokkos::create_mirror_view(B);

    // Copy A to host_A and B to host_B
    // no-op if A and B MemorySpace is HostSpace
    Kokkos::deep_copy(host_A, A);
    Kokkos::deep_copy(host_B, B);

    SerialTrmm_Invoke<typename AVIT::HostMirror, typename BVIT::HostMirror>(side, uplo, trans, diag, alpha, host_A,
                                                                            host_B);

    // Copy host_B to B
    // no-op if B's MemorySpace is HostSpace
    Kokkos::deep_copy(B, host_B);

    Kokkos::Profiling::popRegion();
  }
};
#endif  //! defined(KOKKOSKERNELS_ETI_ONLY) ||
        //! KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

}  // namespace Impl
}  // namespace KokkosBlas

//
// These Macros are for readability.
//
#define KOKKOSBLAS3_TRMM_ETI_SPEC_DECL_LAYOUTS(SCALAR, LAYOUTA, LAYOUTB, EXEC_SPACE, MEM_SPACE)            \
  extern template struct TRMM<EXEC_SPACE,                                                                  \
                              Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
                              Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,       \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
                              false, true>;

#define KOKKOSBLAS3_TRMM_ETI_SPEC_INST_LAYOUTS(SCALAR, LAYOUTA, LAYOUTB, EXEC_SPACE, MEM_SPACE)     \
  template struct TRMM<EXEC_SPACE,                                                                  \
                       Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
                       Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,       \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
                       false, true>;

//
// These Macros are only included when we are not compiling libkokkoskernels but
// are auto generating files. These macros provide the explicit instantiation
// declaration and definition of TRMM, potentially reducing user code size. The
// "extern template" skips the implicit instatiation step ensuring that the
// callers code uses this explicit instantiation definition of TRMM.
//
#define KOKKOSBLAS3_TRMM_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE) \
  KOKKOSBLAS3_TRMM_ETI_SPEC_DECL_LAYOUTS(SCALAR, LAYOUT, LAYOUT, EXEC_SPACE, MEM_SPACE)

#define KOKKOSBLAS3_TRMM_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE) \
  KOKKOSBLAS3_TRMM_ETI_SPEC_INST_LAYOUTS(SCALAR, LAYOUT, LAYOUT, EXEC_SPACE, MEM_SPACE)

#include <KokkosBlas3_trmm_tpl_spec_decl.hpp>

#endif  // KOKKOSBLAS3_TRMM_SPEC_HPP_
