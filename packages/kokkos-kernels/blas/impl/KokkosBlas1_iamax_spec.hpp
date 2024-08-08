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
#ifndef KOKKOSBLAS1_IAMAX_SPEC_HPP_
#define KOKKOSBLAS1_IAMAX_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas1_iamax_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class RMV, class XMV, int rank = XMV::rank>
struct iamax_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Iamax for rank == 1.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_IAMAX_ETI_SPEC_AVAIL_INDEX(INDEX_TYPE, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                \
  template <>                                                                                                    \
  struct iamax_eti_spec_avail<                                                                                   \
      EXEC_SPACE, Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      1> {                                                                                                       \
    enum : bool { value = true };                                                                                \
  };                                                                                                             \
  template <>                                                                                                    \
  struct iamax_eti_spec_avail<EXEC_SPACE,                                                                        \
                              Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,            \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                            \
                              Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,         \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                            \
                              1> {                                                                               \
    enum : bool { value = true };                                                                                \
  };

#define KOKKOSBLAS1_IAMAX_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                \
  KOKKOSBLAS1_IAMAX_ETI_SPEC_AVAIL_INDEX(unsigned long, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE) \
  KOKKOSBLAS1_IAMAX_ETI_SPEC_AVAIL_INDEX(unsigned int, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)  \
  KOKKOSBLAS1_IAMAX_ETI_SPEC_AVAIL_INDEX(int, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Iamax for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_IAMAX_MV_ETI_SPEC_AVAIL_INDEX_HOST(INDEX_TYPE, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)         \
  template <>                                                                                                     \
  struct iamax_eti_spec_avail<                                                                                    \
      EXEC_SPACE, Kokkos::View<INDEX_TYPE*, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                     \
      2> {                                                                                                        \
    enum : bool { value = true };                                                                                 \
  };                                                                                                              \
  template <>                                                                                                     \
  struct iamax_eti_spec_avail<EXEC_SPACE,                                                                         \
                              Kokkos::View<INDEX_TYPE*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,            \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                             \
                              Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,         \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                             \
                              2> {                                                                                \
    enum : bool { value = true };                                                                                 \
  };

#define KOKKOSBLAS1_IAMAX_MV_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                     \
  KOKKOSBLAS1_IAMAX_MV_ETI_SPEC_AVAIL_INDEX_HOST(unsigned long, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE) \
  KOKKOSBLAS1_IAMAX_MV_ETI_SPEC_AVAIL_INDEX_HOST(unsigned int, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)  \
  KOKKOSBLAS1_IAMAX_MV_ETI_SPEC_AVAIL_INDEX_HOST(int, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)

// Include the actual specialization declarations
#include <KokkosBlas1_iamax_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_iamax_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_iamax_mv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

// Unification layer
template <class execution_space, class RMV, class XMV, int rank = XMV::rank,
          bool tpl_spec_avail = iamax_tpl_spec_avail<execution_space, RMV, XMV>::value,
          bool eti_spec_avail = iamax_eti_spec_avail<execution_space, RMV, XMV>::value>
struct Iamax {
  static void iamax(const execution_space& space, const RMV& R, const XMV& X);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of Iamax for single vectors (1-D Views).
template <class execution_space, class RMV, class XMV>
struct Iamax<execution_space, RMV, XMV, 1, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename XMV::size_type size_type;

  static void iamax(const execution_space& space, const RMV& R, const XMV& X) {
    static_assert(Kokkos::is_view<RMV>::value,
                  "KokkosBlas::Impl::"
                  "Iamax<1-D>: RMV is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::"
                  "Iamax<1-D>: XMV is not a Kokkos::View.");
    static_assert(RMV::rank == 0,
                  "KokkosBlas::Impl::Iamax<1-D>: "
                  "RMV is not rank 0.");
    static_assert(XMV::rank == 1,
                  "KokkosBlas::Impl::Iamax<1-D>: "
                  "XMV is not rank 1.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::iamax[ETI]"
                                                                     : "KokkosBlas::iamax[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::iamax<> ETI specialization for < %s , %s >\n", typeid(RMV).name(), typeid(XMV).name());
    else {
      printf("KokkosBlas1::iamax<> non-ETI specialization for < %s , %s >\n", typeid(RMV).name(), typeid(XMV).name());
    }
#endif
    const size_type numRows = X.extent(0);

    if (numRows < static_cast<size_type>(INT_MAX)) {
      V_Iamax_Invoke<execution_space, RMV, XMV, int>(space, R, X);
    } else {
      typedef std::int64_t index_type;
      V_Iamax_Invoke<execution_space, RMV, XMV, index_type>(space, R, X);
    }
    Kokkos::Profiling::popRegion();
  }
};

template <class execution_space, class RV, class XMV>
struct Iamax<execution_space, RV, XMV, 2, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename XMV::size_type size_type;

  static void iamax(const execution_space& space, const RV& R, const XMV& X) {
    static_assert(Kokkos::is_view<RV>::value,
                  "KokkosBlas::Impl::"
                  "Iamax<2-D>: RV is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::"
                  "Iamax<2-D>: XMV is not a Kokkos::View.");
    static_assert(RV::rank == 1,
                  "KokkosBlas::Impl::Iamax<2-D>: "
                  "RV is not rank 1.");
    static_assert(XMV::rank == 2,
                  "KokkosBlas::Impl::Iamax<2-D>: "
                  "XMV is not rank 2.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::iamax[ETI]"
                                                                     : "KokkosBlas::iamax[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::iamax<> ETI specialization for < %s , %s >\n", typeid(RV).name(), typeid(XMV).name());
    else {
      printf("KokkosBlas1::iamax<> non-ETI specialization for < %s , %s >\n", typeid(RV).name(), typeid(XMV).name());
    }
#endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);
    if (numRows < static_cast<size_type>(INT_MAX) && numRows * numCols < static_cast<size_type>(INT_MAX)) {
      MV_Iamax_Invoke<execution_space, RV, XMV, int>(space, R, X);
    } else {
      typedef std::int64_t index_type;
      MV_Iamax_Invoke<execution_space, RV, XMV, index_type>(space, R, X);
    }
    Kokkos::Profiling::popRegion();
  }
};
#endif

}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Iamax for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_IAMAX_ETI_SPEC_DECL_INDEX(INDEX_TYPE, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                 \
  extern template struct Iamax<                                                                                  \
      EXEC_SPACE, Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      1, false, true>;                                                                                           \
  extern template struct Iamax<EXEC_SPACE,                                                                       \
                               Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,           \
                                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                           \
                               Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,        \
                                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                           \
                               1, false, true>;

#define KOKKOSBLAS1_IAMAX_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                \
  KOKKOSBLAS1_IAMAX_ETI_SPEC_DECL_INDEX(unsigned long, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE) \
  KOKKOSBLAS1_IAMAX_ETI_SPEC_DECL_INDEX(unsigned int, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)  \
  KOKKOSBLAS1_IAMAX_ETI_SPEC_DECL_INDEX(int, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Iamax for rank == 1.  This is NOT for users!!!  We
// use this macro in one or more .cpp files in this directory.
//
#define KOKKOSBLAS1_IAMAX_ETI_SPEC_INST_INDEX(INDEX_TYPE, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                       \
  template struct Iamax<EXEC_SPACE,                                                                                    \
                        Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                     \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                        \
                        1, false, true>;                                                                               \
  template struct Iamax<EXEC_SPACE,                                                                                    \
                        Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                        \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                        \
                        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                     \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                        \
                        1, false, true>;

#define KOKKOSBLAS1_IAMAX_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                \
  KOKKOSBLAS1_IAMAX_ETI_SPEC_INST_INDEX(unsigned long, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE) \
  KOKKOSBLAS1_IAMAX_ETI_SPEC_INST_INDEX(unsigned int, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)  \
  KOKKOSBLAS1_IAMAX_ETI_SPEC_INST_INDEX(int, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Iamax for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_IAMAX_MV_ETI_SPEC_DECL_INDEX(INDEX_TYPE, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)               \
  extern template struct Iamax<                                                                                   \
      EXEC_SPACE, Kokkos::View<INDEX_TYPE*, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                     \
      2, false, true>;                                                                                            \
  extern template struct Iamax<EXEC_SPACE,                                                                        \
                               Kokkos::View<INDEX_TYPE*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,           \
                                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                            \
                               Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,         \
                                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                            \
                               2, false, true>;

#define KOKKOSBLAS1_IAMAX_MV_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                \
  KOKKOSBLAS1_IAMAX_MV_ETI_SPEC_DECL_INDEX(unsigned long, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE) \
  KOKKOSBLAS1_IAMAX_MV_ETI_SPEC_DECL_INDEX(unsigned int, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)  \
  KOKKOSBLAS1_IAMAX_MV_ETI_SPEC_DECL_INDEX(int, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Iamax for rank == 2.  This is NOT for users!!!  We
// use this macro in one or more .cpp files in this directory.
//
#define KOKKOSBLAS1_IAMAX_MV_ETI_SPEC_INST_INDEX(INDEX_TYPE, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)               \
  template struct Iamax<                                                                                          \
      EXEC_SPACE, Kokkos::View<INDEX_TYPE*, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                     \
      2, false, true>;                                                                                            \
  template struct Iamax<EXEC_SPACE,                                                                               \
                        Kokkos::View<INDEX_TYPE*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                  \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                   \
                        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,               \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                   \
                        2, false, true>;

#define KOKKOSBLAS1_IAMAX_MV_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                \
  KOKKOSBLAS1_IAMAX_MV_ETI_SPEC_INST_INDEX(unsigned long, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE) \
  KOKKOSBLAS1_IAMAX_MV_ETI_SPEC_INST_INDEX(unsigned int, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)  \
  KOKKOSBLAS1_IAMAX_MV_ETI_SPEC_INST_INDEX(int, SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)

#include <KokkosBlas1_iamax_tpl_spec_decl.hpp>

#endif  // KOKKOSBLAS1_IAMAX_SPEC_HPP_
