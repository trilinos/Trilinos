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
#ifndef KOKKOS_BLAS1_IMPL_ABS_SPEC_HPP_
#define KOKKOS_BLAS1_IMPL_ABS_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas1_abs_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class RMV, class XMV, int rank = RMV::rank>
struct abs_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Abs for rank == 1.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_ABS_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                         \
  template <>                                                                                                         \
  struct abs_eti_spec_avail<                                                                                          \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      1> {                                                                                                            \
    enum : bool { value = true };                                                                                     \
  };

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Abs for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_ABS_MV_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                       \
  template <>                                                                                                          \
  struct abs_eti_spec_avail<                                                                                           \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      2> {                                                                                                             \
    enum : bool { value = true };                                                                                      \
  };

// Include the actual specialization declarations
#include <KokkosBlas1_abs_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_abs_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_abs_mv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

// Unification layer
template <class execution_space, class RMV, class XMV, int rank = RMV::rank,
          bool tpl_spec_avail = abs_tpl_spec_avail<execution_space, RMV, XMV>::value,
          bool eti_spec_avail = abs_eti_spec_avail<execution_space, RMV, XMV>::value>
struct Abs {
  static void abs(const execution_space& space, const RMV& R, const XMV& X);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of Abs for single vectors (1-D Views).
template <class execution_space, class RMV, class XMV>
struct Abs<execution_space, RMV, XMV, 1, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  using size_type = typename XMV::size_type;

  static void abs(const execution_space& space, const RMV& R, const XMV& X) {
    static_assert(Kokkos::is_view<RMV>::value,
                  "KokkosBlas::Impl::"
                  "Abs<1-D>: RMV is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::"
                  "Abs<1-D>: XMV is not a Kokkos::View.");
    static_assert(RMV::rank == 1,
                  "KokkosBlas::Impl::Abs<1-D>: "
                  "RMV is not rank 1.");
    static_assert(XMV::rank == 1,
                  "KokkosBlas::Impl::Abs<1-D>: "
                  "XMV is not rank 1.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::abs[ETI]"
                                                                     : "KokkosBlas::abs[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::abs<> ETI specialization for < %s , %s >\n", typeid(RMV).name(), typeid(XMV).name());
    else {
      printf("KokkosBlas1::abs<> non-ETI specialization for < %s , %s >\n", typeid(RMV).name(), typeid(XMV).name());
    }
#endif
    const size_type numRows = X.extent(0);

    if (numRows < static_cast<size_type>(INT_MAX)) {
      typedef int index_type;
      V_Abs_Generic<execution_space, RMV, XMV, index_type>(space, R, X);
    } else {
      typedef std::int64_t index_type;
      V_Abs_Generic<execution_space, RMV, XMV, index_type>(space, R, X);
    }
    Kokkos::Profiling::popRegion();
  }
};

template <class execution_space, class RMV, class XMV>
struct Abs<execution_space, RMV, XMV, 2, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  using size_type = typename XMV::size_type;

  static void abs(const execution_space& space, const RMV& R, const XMV& X) {
    static_assert(Kokkos::is_view<RMV>::value,
                  "KokkosBlas::Impl::"
                  "Abs<2-D>: RMV is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::"
                  "Abs<2-D>: XMV is not a Kokkos::View.");
    static_assert(RMV::rank == 2,
                  "KokkosBlas::Impl::Abs<2-D>: "
                  "RMV is not rank 2.");
    static_assert(XMV::rank == 2,
                  "KokkosBlas::Impl::Abs<2-D>: "
                  "XMV is not rank 2.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::abs[ETI]"
                                                                     : "KokkosBlas::abs[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::abs<> ETI specialization for < %s , %s >\n", typeid(RMV).name(), typeid(XMV).name());
    else {
      printf("KokkosBlas1::asb<> non-ETI specialization for < %s , %s >\n", typeid(RMV).name(), typeid(XMV).name());
    }
#endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);
    if (numRows < static_cast<size_type>(INT_MAX) && numRows * numCols < static_cast<size_type>(INT_MAX)) {
      typedef int index_type;
      MV_Abs_Generic<execution_space, RMV, XMV, index_type>(space, R, X);
    } else {
      typedef std::int64_t index_type;
      MV_Abs_Generic<execution_space, RMV, XMV, index_type>(space, R, X);
    }
    Kokkos::Profiling::popRegion();
  }
};
#endif

}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Abs for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_ABS_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                          \
  extern template struct Abs<                                                                                         \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      1, false, true>;

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Abs for rank == 2.  This is NOT for users!!!  We
// use this macro in one or more .cpp files in this directory.
//
#define KOKKOSBLAS1_ABS_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                          \
  template struct Abs<                                                                                                \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      1, false, true>;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Abs for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_ABS_MV_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                        \
  extern template struct Abs<                                                                                          \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      2, false, true>;

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Abs for rank == 2.  This is NOT for users!!!  We
// use this macro in one or more .cpp files in this directory.
//
#define KOKKOSBLAS1_ABS_MV_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                        \
  template struct Abs<                                                                                                 \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      2, false, true>;

#include <KokkosBlas1_abs_tpl_spec_decl.hpp>

#endif  // KOKKOS_BLAS1_MV_IMPL_ABS_HPP_
