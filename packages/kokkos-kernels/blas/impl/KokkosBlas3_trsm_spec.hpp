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
#ifndef KOKKOSBLAS3_TRSM_SPEC_HPP_
#define KOKKOSBLAS3_TRSM_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"
#include <sstream>

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas3_trsm_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class AVT, class BVT>
struct trsm_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::TRSM.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS3_TRSM_ETI_SPEC_AVAIL_LAYOUT(SCALAR, LAYOUTA, LAYOUTB, EXEC_SPACE, MEM_SPACE)           \
  template <>                                                                                             \
  struct trsm_eti_spec_avail<EXEC_SPACE,                                                                  \
                             Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
                             Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,       \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                   \
    enum : bool { value = true };                                                                         \
  };

#define KOKKOSBLAS3_TRSM_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE) \
  KOKKOSBLAS3_TRSM_ETI_SPEC_AVAIL_LAYOUT(SCALAR, LAYOUT, LAYOUT, EXEC_SPACE, MEM_SPACE)

// Include the actual specialization declarations
#include <KokkosBlas3_trsm_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas3_trsm_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

//
// trsm
//

// Unification layer
template <class execution_space, class AViewType, class BViewType,
          bool tpl_spec_avail = trsm_tpl_spec_avail<execution_space, AViewType, BViewType>::value,
          bool eti_spec_avail = trsm_eti_spec_avail<execution_space, AViewType, BViewType>::value>
struct TRSM {
  static void trsm(const execution_space& space, const char side[], const char uplo[], const char trans[],
                   const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A,
                   const BViewType& B);
};

// Implementation of KokkosBlas::trsm.
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
template <class execution_space, class AViewType, class BViewType>
struct TRSM<execution_space, AViewType, BViewType, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void trsm(const execution_space& /*space*/, const char side[], const char uplo[], const char trans[],
                   const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A,
                   const BViewType& B) {
    static_assert(Kokkos::is_view<AViewType>::value, "AViewType must be a Kokkos::View.");
    static_assert(Kokkos::is_view<BViewType>::value, "BViewType must be a Kokkos::View.");
    static_assert(static_cast<int>(AViewType::rank) == 2, "AViewType must have rank 2.");
    static_assert(static_cast<int>(BViewType::rank) == 2, "BViewType must have rank 2.");

    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::trsm[ETI]"
                                                                     : "KokkosBlas::trsm[noETI]");

    typename AViewType::HostMirror h_A = Kokkos::create_mirror_view(A);
    typename BViewType::HostMirror h_B = Kokkos::create_mirror_view(B);

    Kokkos::deep_copy(h_A, A);
    Kokkos::deep_copy(h_B, B);

    SerialTrsm_Invoke<typename AViewType::HostMirror, typename BViewType::HostMirror>(side, uplo, trans, diag, alpha,
                                                                                      h_A, h_B);

    Kokkos::deep_copy(B, h_B);

    Kokkos::Profiling::popRegion();
  }
};
#endif  //! defined(KOKKOSKERNELS_ETI_ONLY) ||
        //! KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::TRSM.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file. We may spread out definitions (see _DEF macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS3_TRSM_ETI_SPEC_DECL_LAYOUTS(SCALAR, LAYOUTA, LAYOUTB, EXEC_SPACE, MEM_SPACE)            \
  extern template struct TRSM<EXEC_SPACE,                                                                  \
                              Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
                              Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,       \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
                              false, true>;

#define KOKKOSBLAS3_TRSM_ETI_SPEC_INST_LAYOUTS(SCALAR, LAYOUTA, LAYOUTB, EXEC_SPACE, MEM_SPACE)     \
  template struct TRSM<EXEC_SPACE,                                                                  \
                       Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
                       Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,       \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
                       false, true>;

#define KOKKOSBLAS3_TRSM_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE) \
  KOKKOSBLAS3_TRSM_ETI_SPEC_DECL_LAYOUTS(SCALAR, LAYOUT, LAYOUT, EXEC_SPACE, MEM_SPACE)

#define KOKKOSBLAS3_TRSM_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE) \
  KOKKOSBLAS3_TRSM_ETI_SPEC_INST_LAYOUTS(SCALAR, LAYOUT, LAYOUT, EXEC_SPACE, MEM_SPACE)

#include <KokkosBlas3_trsm_tpl_spec_decl.hpp>

#endif  // KOKKOSBLAS3_TRSM_SPEC_HPP_
