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
#ifndef KOKKOS_BLAS1_IMPL_RECIPROCAL_SPEC_HPP_
#define KOKKOS_BLAS1_IMPL_RECIPROCAL_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas1_reciprocal_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class RMV, class XMV, int rank = RMV::rank>
struct reciprocal_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Reciprocal for rank == 1.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_RECIPROCAL_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE,  \
                                              MEM_SPACE)                   \
  template <>                                                              \
  struct reciprocal_eti_spec_avail<                                        \
      EXEC_SPACE,                                                          \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,              \
      Kokkos::View<const SCALAR*, LAYOUT,                                  \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,              \
      1> {                                                                 \
    enum : bool { value = true };                                          \
  };

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Reciprocal for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_RECIPROCAL_MV_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, \
                                                 MEM_SPACE)                  \
  template <>                                                                \
  struct reciprocal_eti_spec_avail<                                          \
      EXEC_SPACE,                                                            \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                \
      Kokkos::View<const SCALAR**, LAYOUT,                                   \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                \
      2> {                                                                   \
    enum : bool { value = true };                                            \
  };

// Include the actual specialization declarations
#include <KokkosBlas1_reciprocal_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_reciprocal_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_reciprocal_mv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

// Unification layer
template <class execution_space, class RMV, class XMV, int rank = RMV::rank,
          bool tpl_spec_avail =
              reciprocal_tpl_spec_avail<execution_space, RMV, XMV>::value,
          bool eti_spec_avail =
              reciprocal_eti_spec_avail<execution_space, RMV, XMV>::value>
struct Reciprocal {
  static void reciprocal(const execution_space& space, const RMV& R,
                         const XMV& X);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of Reciprocal for single vectors (1-D Views).
template <class execution_space, class RMV, class XMV>
struct Reciprocal<execution_space, RMV, XMV, 1, false,
                  KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename XMV::size_type size_type;

  static void reciprocal(const execution_space& space, const RMV& R,
                         const XMV& X) {
    static_assert(Kokkos::is_view<RMV>::value,
                  "KokkosBlas::Impl::"
                  "Reciprocal<1-D>: RMV is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::"
                  "Reciprocal<1-D>: XMV is not a Kokkos::View.");
    static_assert(RMV::rank == 1,
                  "KokkosBlas::Impl::Reciprocal<1-D>: "
                  "RMV is not rank 1.");
    static_assert(XMV::rank == 1,
                  "KokkosBlas::Impl::Reciprocal<1-D>: "
                  "XMV is not rank 1.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
                                      ? "KokkosBlas::reciprocal[ETI]"
                                      : "KokkosBlas::reciprocal[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::reciprocal<> ETI specialization for < %s , %s >\n",
             typeid(RMV).name(), typeid(XMV).name());
    else {
      printf(
          "KokkosBlas1::reciprocal<> non-ETI specialization for < %s , %s >\n",
          typeid(RMV).name(), typeid(XMV).name());
    }
#endif
    const size_type numRows = X.extent(0);

    if (numRows < static_cast<size_type>(INT_MAX)) {
      typedef int index_type;
      V_Reciprocal_Generic<execution_space, RMV, XMV, index_type>(space, R, X);
    } else {
      typedef std::int64_t index_type;
      V_Reciprocal_Generic<execution_space, RMV, XMV, index_type>(space, R, X);
    }
    Kokkos::Profiling::popRegion();
  }
};

template <class execution_space, class RMV, class XMV>
struct Reciprocal<execution_space, RMV, XMV, 2, false,
                  KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename XMV::size_type size_type;

  static void reciprocal(const execution_space& space, const RMV& R,
                         const XMV& X) {
    static_assert(Kokkos::is_view<RMV>::value,
                  "KokkosBlas::Impl::"
                  "Reciprocal<2-D>: RMV is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::"
                  "Reciprocal<2-D>: XMV is not a Kokkos::View.");
    static_assert(RMV::rank == 2,
                  "KokkosBlas::Impl::Reciprocal<2-D>: "
                  "RMV is not rank 2.");
    static_assert(XMV::rank == 2,
                  "KokkosBlas::Impl::Reciprocal<2-D>: "
                  "XMV is not rank 2.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
                                      ? "KokkosBlas::reciprocal[ETI]"
                                      : "KokkosBlas::reciprocal[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::reciprocal<> ETI specialization for < %s , %s >\n",
             typeid(RMV).name(), typeid(XMV).name());
    else {
      printf("KokkosBlas1::asb<> non-ETI specialization for < %s , %s >\n",
             typeid(RMV).name(), typeid(XMV).name());
    }
#endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);
    if (numRows < static_cast<size_type>(INT_MAX) &&
        numRows * numCols < static_cast<size_type>(INT_MAX)) {
      typedef int index_type;
      MV_Reciprocal_Generic<execution_space, RMV, XMV, index_type>(space, R, X);
    } else {
      typedef std::int64_t index_type;
      MV_Reciprocal_Generic<execution_space, RMV, XMV, index_type>(space, R, X);
    }
    Kokkos::Profiling::popRegion();
  }
};
#endif

}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Reciprocal for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_RECIPROCAL_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE,   \
                                             MEM_SPACE)                    \
  extern template struct Reciprocal<                                       \
      EXEC_SPACE,                                                          \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,              \
      Kokkos::View<const SCALAR*, LAYOUT,                                  \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,              \
      1, false, true>;

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Reciprocal for rank == 2.  This is NOT for users!!!  We
// use this macro in one or more .cpp files in this directory.
//
#define KOKKOSBLAS1_RECIPROCAL_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE,   \
                                             MEM_SPACE)                    \
  template struct Reciprocal<                                              \
      EXEC_SPACE,                                                          \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,              \
      Kokkos::View<const SCALAR*, LAYOUT,                                  \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,              \
      1, false, true>;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Reciprocal for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_RECIPROCAL_MV_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, \
                                                MEM_SPACE)                  \
  extern template struct Reciprocal<                                        \
      EXEC_SPACE,                                                           \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      Kokkos::View<const SCALAR**, LAYOUT,                                  \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      2, false, true>;

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Reciprocal for rank == 2.  This is NOT for users!!!  We
// use this macro in one or more .cpp files in this directory.
//
#define KOKKOSBLAS1_RECIPROCAL_MV_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, \
                                                MEM_SPACE)                  \
  template struct Reciprocal<                                               \
      EXEC_SPACE,                                                           \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      Kokkos::View<const SCALAR**, LAYOUT,                                  \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      2, false, true>;

#include <KokkosBlas1_reciprocal_tpl_spec_decl.hpp>

#endif  // KOKKOS_BLAS1_MV_IMPL_RECIPROCAL_HPP_
