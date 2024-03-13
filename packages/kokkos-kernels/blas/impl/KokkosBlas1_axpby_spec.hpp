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
#ifndef KOKKOS_BLAS1_AXPBY_SPEC_HPP_
#define KOKKOS_BLAS1_AXPBY_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas1_axpby_impl.hpp>
#include <KokkosBlas1_axpby_mv_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class AV, class XMV, class BV, class YMV,
          int rank = YMV::rank>
struct axpby_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Axpby for rank == 1.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_AXPBY_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE,       \
                                         MEM_SPACE)                        \
  template <>                                                              \
  struct axpby_eti_spec_avail<                                             \
      EXEC_SPACE, SCALAR,                                                  \
      Kokkos::View<const SCALAR*, LAYOUT,                                  \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,              \
      SCALAR,                                                              \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,              \
      1> {                                                                 \
    enum : bool { value = true };                                          \
  };

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Axpby for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_AXPBY_MV_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE,     \
                                            MEM_SPACE)                      \
  template <>                                                               \
  struct axpby_eti_spec_avail<                                              \
      EXEC_SPACE, SCALAR,                                                   \
      Kokkos::View<const SCALAR**, LAYOUT,                                  \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      SCALAR,                                                               \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      2> {                                                                  \
    enum : bool { value = true };                                           \
  };                                                                        \
  template <>                                                               \
  struct axpby_eti_spec_avail<                                              \
      EXEC_SPACE,                                                           \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft,                       \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      Kokkos::View<const SCALAR**, LAYOUT,                                  \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft,                       \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      2> {                                                                  \
    enum : bool { value = true };                                           \
  };

// Include the actual specialization declarations
#include <KokkosBlas1_axpby_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_axpby_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_axpby_mv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

//
// axpby
//

/// \brief Implementation of KokkosBlas::axpby for (multi)vectors.
///
/// Compute any of the following, depending on the types of the input
/// arguments of axpxy():
///
/// 1. Y(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j) (if R, X, and Y are 2-D,
///    and av and bv are 1-D)
///
/// 2. Y(i,j) = av*X(i,j) + bv*Y(i,j) (if R, X, and Y are 2-D,
///    and av and bv are scalars)
///
/// 3. Y(i) = av()*X(i) + bv()*Y(i) (if R, X, and Y are 1-D, and av
///    and bv are 0-D Views (not scalars))
///
/// 4. Y(i) = av*X(i) + bv*Y(i) (if R, X, and Y are 1-D, and av and bv
///    are scalars)
///
/// Any <i>scalar</i> coefficient of zero has BLAS semantics of
/// ignoring the corresponding (multi)vector entry.  This does NOT
/// apply to coefficients in av and bv vectors, if they are used.
template <class execution_space, class AV, class XMV, class BV, class YMV,
          int rank = YMV::rank,
          bool tpl_spec_avail =
              axpby_tpl_spec_avail<execution_space, AV, XMV, BV, YMV>::value,
          bool eti_spec_avail =
              axpby_eti_spec_avail<execution_space, AV, XMV, BV, YMV>::value>
struct Axpby {
  static void axpby(const execution_space& space, const AV& av, const XMV& X,
                    const BV& bv, const YMV& Y);
};

template <class execution_space, class AV, class XMV, class BV, class YMV>
struct Axpby<execution_space, AV, XMV, BV, YMV, 0, true, true> {
  static void axpby(const execution_space& /*space*/, const AV& /* av */,
                    const XMV& /* X */, const BV& /* bv */,
                    const YMV& /* Y */) {
    static_assert(YMV::rank == 0, "Oh My God");
  }
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
// Full specialization for XMV and YMV rank-2 Views.
template <class execution_space, class AV, class XMV, class BV, class YMV>
struct Axpby<execution_space, AV, XMV, BV, YMV, 2, false,
             KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename YMV::size_type size_type;

  static void axpby(const execution_space& space, const AV& av, const XMV& X,
                    const BV& bv, const YMV& Y) {
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::"
                  "Axpby<rank-2>::axpby: X is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YMV>::value,
                  "KokkosBlas::Impl::"
                  "Axpby<rank-2>::axpby: Y is not a Kokkos::View.");
    static_assert(std::is_same<typename YMV::value_type,
                               typename YMV::non_const_value_type>::value,
                  "KokkosBlas::Impl::Axpby<rank-2>::axpby: Y is const.  "
                  "It must be nonconst, because it is an output argument "
                  "(we have to be able to write to its entries).");
    static_assert((int)YMV::rank == (int)XMV::rank,
                  "KokkosBlas::Impl::Axpby<rank-2>::axpby (MV): "
                  "X and Y must have the same rank.");
    static_assert(YMV::rank == 2,
                  "KokkosBlas::Impl::Axpby<rank-2>::axpby: "
                  "X and Y must have rank 2.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
                                      ? "KokkosBlas::axpby[ETI]"
                                      : "KokkosBlas::axpby[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf(
          "KokkosBlas1::axpby<> ETI specialization for < %s , %s , %s , %s >\n",
          typeid(AV).name(), typeid(XMV).name(), typeid(BV).name(),
          typeid(YMV).name());
    else {
      printf(
          "KokkosBlas1::axpby<> non-ETI specialization for < %s , %s , %s , %s "
          ">\n",
          typeid(AV).name(), typeid(XMV).name(), typeid(BV).name(),
          typeid(YMV).name());
    }
#endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);
    int a = 2, b = 2;
    if (av.extent(0) == 0) {
      a = 0;
    }
    if (bv.extent(0) == 0) {
      b = 0;
    }

    if (numRows < static_cast<size_type>(INT_MAX) &&
        numRows * numCols < static_cast<size_type>(INT_MAX)) {
      typedef int index_type;
      typedef typename std::conditional<
          std::is_same<typename XMV::array_layout, Kokkos::LayoutLeft>::value,
          Axpby_MV_Invoke_Right<execution_space, AV, XMV, BV, YMV, index_type>,
          Axpby_MV_Invoke_Left<execution_space, AV, XMV, BV, YMV,
                               index_type> >::type Axpby_MV_Invoke_Layout;
      Axpby_MV_Invoke_Layout::run(space, av, X, bv, Y, a, b);
    } else {
      typedef typename XMV::size_type index_type;
      typedef typename std::conditional<
          std::is_same<typename XMV::array_layout, Kokkos::LayoutLeft>::value,
          Axpby_MV_Invoke_Right<execution_space, AV, XMV, BV, YMV, index_type>,
          Axpby_MV_Invoke_Left<execution_space, AV, XMV, BV, YMV,
                               index_type> >::type Axpby_MV_Invoke_Layout;
      Axpby_MV_Invoke_Layout::run(space, av, X, bv, Y, a, b);
    }
    Kokkos::Profiling::popRegion();
  }
};

// Partial specialization for XMV, and YMV rank-2 Views,
// and AV and BV scalars.
template <class execution_space, class XMV, class YMV>
struct Axpby<execution_space, typename XMV::non_const_value_type, XMV,
             typename YMV::non_const_value_type, YMV, 2, false,
             KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename XMV::non_const_value_type AV;
  typedef typename YMV::non_const_value_type BV;
  typedef typename YMV::size_type size_type;
  typedef Kokkos::ArithTraits<typename XMV::non_const_value_type> ATA;
  typedef Kokkos::ArithTraits<typename YMV::non_const_value_type> ATB;

  static void axpby(const execution_space& space, const AV& alpha, const XMV& X,
                    const BV& beta, const YMV& Y) {
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::Axpby::axpby (MV): "
                  "X is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YMV>::value,
                  "KokkosBlas::Impl::Axpby::axpby (MV): "
                  "Y is not a Kokkos::View.");
    static_assert(std::is_same<typename YMV::value_type,
                               typename YMV::non_const_value_type>::value,
                  "KokkosBlas::Impl::Axpby::axpby (MV): Y is const.  "
                  "It must be nonconst, because it is an output argument "
                  "(we have to be able to write to its entries).");
    static_assert((int)YMV::rank == (int)XMV::rank,
                  "KokkosBlas::Impl::Axpby::axpby (MV): "
                  "X and Y must have the same rank.");
    static_assert(YMV::rank == 2,
                  "KokkosBlas::Impl::Axpby::axpby (MV): "
                  "X and Y must have rank 2.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
                                      ? "KokkosBlas::axpby[ETI]"
                                      : "KokkosBlas::axpby[noETI]");

#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf(
          "KokkosBlas1::axpby<> ETI specialization for < %s , %s , %s , %s >\n",
          typeid(AV).name(), typeid(XMV).name(), typeid(BV).name(),
          typeid(YMV).name());
    else {
      printf(
          "KokkosBlas1::axpby<> non-ETI specialization for < %s , %s , %s , %s "
          ">\n",
          typeid(AV).name(), typeid(XMV).name(), typeid(BV).name(),
          typeid(YMV).name());
    }
#endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);
    int a, b;
    if (alpha == ATA::zero()) {
      a = 0;
    }
#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
    else if (alpha == -ATA::one()) {
      a = -1;
    } else if (alpha == ATA::one()) {
      a = 1;
    }
#endif  // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
    else {
      a = 2;
    }
    if (beta == ATB::zero()) {
      b = 0;
    }
#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
    else if (beta == -ATB::one()) {
      b = -1;
    } else if (beta == ATB::one()) {
      b = 1;
    }
#endif  // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
    else {
      b = 2;
    }

    if (numRows < static_cast<size_type>(INT_MAX) &&
        numRows * numCols < static_cast<size_type>(INT_MAX)) {
      typedef int index_type;
      typedef typename std::conditional<
          std::is_same<typename XMV::array_layout, Kokkos::LayoutLeft>::value,
          Axpby_MV_Invoke_Right<execution_space, AV, XMV, BV, YMV, index_type>,
          Axpby_MV_Invoke_Left<execution_space, AV, XMV, BV, YMV,
                               index_type> >::type Axpby_MV_Invoke_Layout;
      Axpby_MV_Invoke_Layout::run(space, alpha, X, beta, Y, a, b);
    } else {
      typedef typename XMV::size_type index_type;
      typedef typename std::conditional<
          std::is_same<typename XMV::array_layout, Kokkos::LayoutLeft>::value,
          Axpby_MV_Invoke_Right<execution_space, AV, XMV, BV, YMV, index_type>,
          Axpby_MV_Invoke_Left<execution_space, AV, XMV, BV, YMV,
                               index_type> >::type Axpby_MV_Invoke_Layout;
      Axpby_MV_Invoke_Layout::run(space, alpha, X, beta, Y, a, b);
    }
    Kokkos::Profiling::popRegion();
  }
};

// Partial specialization for XV and YV rank-1 Views,
// and AV and BV scalars.
template <class execution_space, class XV, class YV>
struct Axpby<execution_space, typename XV::non_const_value_type, XV,
             typename YV::non_const_value_type, YV, 1, false,
             KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename XV::non_const_value_type AV;
  typedef typename YV::non_const_value_type BV;
  typedef typename YV::size_type size_type;
  typedef Kokkos::ArithTraits<typename XV::non_const_value_type> ATA;
  typedef Kokkos::ArithTraits<typename YV::non_const_value_type> ATB;

  static void axpby(const execution_space& space, const AV& alpha, const XV& X,
                    const BV& beta, const YV& Y) {
    static_assert(Kokkos::is_view<XV>::value,
                  "KokkosBlas::Impl::"
                  "Axpby<rank-1>::axpby: X is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YV>::value,
                  "KokkosBlas::Impl::"
                  "Axpby<rank-1>::axpby: Y is not a Kokkos::View.");
    static_assert(std::is_same<typename YV::value_type,
                               typename YV::non_const_value_type>::value,
                  "KokkosBlas::Impl::Axpby<rank-1>::axpby: Y is const.  "
                  "It must be nonconst, because it is an output argument "
                  "(we have to be able to write to its entries).");
    static_assert((int)YV::rank == (int)XV::rank,
                  "KokkosBlas::Impl::"
                  "Axpby<rank-1>::axpby: X and Y must have the same rank.");
    static_assert(YV::rank == 1,
                  "KokkosBlas::Impl::Axpby<rank-1>::axpby: "
                  "X and Y must have rank 1.");

    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
                                      ? "KokkosBlas::axpby[ETI]"
                                      : "KokkosBlas::axpby[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf(
          "KokkosBlas1::axpby<> ETI specialization for < %s , %s , %s , %s >\n",
          typeid(AV).name(), typeid(XV).name(), typeid(BV).name(),
          typeid(YV).name());
    else {
      printf(
          "KokkosBlas1::axpby<> non-ETI specialization for < %s , %s , %s , %s "
          ">\n",
          typeid(AV).name(), typeid(XV).name(), typeid(BV).name(),
          typeid(YV).name());
    }
#endif

    const size_type numRows = X.extent(0);
    int a                   = 2;
    if (alpha == ATA::zero()) {
      a = 0;
    }
#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
    else if (alpha == -ATA::one()) {
      a = -1;
    } else if (alpha == ATA::one()) {
      a = 1;
    }
#endif  // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2

    int b = 2;
    if (beta == ATB::zero()) {
      b = 0;
    }
#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
    else if (beta == -ATB::one()) {
      b = -1;
    } else if (beta == ATB::one()) {
      b = 1;
    }
#endif  // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2

    if (numRows < static_cast<size_type>(INT_MAX)) {
      typedef int index_type;
      Axpby_Generic<execution_space, typename XV::non_const_value_type, XV,
                    typename YV::non_const_value_type, YV, index_type>(
          space, alpha, X, beta, Y, 0, a, b);
    } else {
      typedef typename XV::size_type index_type;
      Axpby_Generic<execution_space, typename XV::non_const_value_type, XV,
                    typename YV::non_const_value_type, YV, index_type>(
          space, alpha, X, beta, Y, 0, a, b);
    }
    Kokkos::Profiling::popRegion();
  }
};
#endif  //! defined(KOKKOSKERNELS_ETI_ONLY) ||
        //! KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Axpby for rank == 1.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _INST macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS1_AXPBY_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE) \
  extern template struct Axpby<                                                \
      EXEC_SPACE, SCALAR,                                                      \
      Kokkos::View<const SCALAR*, LAYOUT,                                      \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      SCALAR,                                                                  \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      1, false, true>;

#define KOKKOSBLAS1_AXPBY_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE) \
  template struct Axpby<                                                       \
      EXEC_SPACE, SCALAR,                                                      \
      Kokkos::View<const SCALAR*, LAYOUT,                                      \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      SCALAR,                                                                  \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      1, false, true>;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Axpby for rank == 2.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _DEF macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS1_AXPBY_MV_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE,      \
                                           MEM_SPACE)                       \
  extern template struct Axpby<                                             \
      EXEC_SPACE, SCALAR,                                                   \
      Kokkos::View<const SCALAR**, LAYOUT,                                  \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      SCALAR,                                                               \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      2, false, true>;                                                      \
  extern template struct Axpby<                                             \
      EXEC_SPACE,                                                           \
      Kokkos::View<const SCALAR*, LAYOUT,                                   \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      Kokkos::View<const SCALAR**, LAYOUT,                                  \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      Kokkos::View<const SCALAR*, LAYOUT,                                   \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      2, false, true>;

#define KOKKOSBLAS1_AXPBY_MV_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE,      \
                                           MEM_SPACE)                       \
  template struct Axpby<                                                    \
      EXEC_SPACE, SCALAR,                                                   \
      Kokkos::View<const SCALAR**, LAYOUT,                                  \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      SCALAR,                                                               \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      2, false, true>;                                                      \
  template struct Axpby<                                                    \
      EXEC_SPACE,                                                           \
      Kokkos::View<const SCALAR*, LAYOUT,                                   \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      Kokkos::View<const SCALAR**, LAYOUT,                                  \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      Kokkos::View<const SCALAR*, LAYOUT,                                   \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      2, false, true>;

#include <KokkosBlas1_axpby_tpl_spec_decl.hpp>

#endif  // KOKKOS_BLAS1_MV_IMPL_AXPBY_HPP_
