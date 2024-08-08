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
#ifndef KOKKOSBLAS1_UPDATE_SPEC_HPP_
#define KOKKOSBLAS1_UPDATE_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas1_update_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class XMV, class YMV, class ZMV, int rank = ZMV::rank>
struct update_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Update for rank == 1.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_UPDATE_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                      \
  template <>                                                                                                         \
  struct update_eti_spec_avail<                                                                                       \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      1> {                                                                                                            \
    enum : bool { value = true };                                                                                     \
  };

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Update for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_UPDATE_MV_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                    \
  template <>                                                                                                          \
  struct update_eti_spec_avail<                                                                                        \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      2> {                                                                                                             \
    enum : bool { value = true };                                                                                      \
  };

// Include the actual specialization declarations
#include <KokkosBlas1_update_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_update_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_update_mv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

//
// update
//

/// \brief Implementation of KokkosBlas::update for single vectors and
///   multivectors.
///
/// Compute
///
/// Z(i,j) = alpha*X(i,j) + beta*Y(i,j) + gamma*Z(i,j),
///
/// with special cases for alpha, beta, or gamma = 0.
template <class execution_space, class XMV, class YMV, class ZMV, int rank = ZMV::rank,
          bool tpl_spec_avail = update_tpl_spec_avail<execution_space, XMV, YMV, ZMV>::value,
          bool eti_spec_avail = update_eti_spec_avail<execution_space, XMV, YMV, ZMV>::value>
struct Update {
  static void update(const execution_space& space, const typename XMV::non_const_value_type& alpha, const XMV& X,
                     const typename YMV::non_const_value_type& beta, const YMV& Y,
                     const typename ZMV::non_const_value_type& gamma, const ZMV& Z);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
// Partial specialization for XMV, YMV, and ZMV rank-2 Views.
template <class execution_space, class XMV, class YMV, class ZMV>
struct Update<execution_space, XMV, YMV, ZMV, 2, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename XMV::size_type size_type;
  typedef Kokkos::ArithTraits<typename XMV::non_const_value_type> ATA;
  typedef Kokkos::ArithTraits<typename YMV::non_const_value_type> ATB;
  typedef Kokkos::ArithTraits<typename ZMV::non_const_value_type> ATC;

  static void update(const execution_space& space, const typename XMV::non_const_value_type& alpha, const XMV& X,
                     const typename YMV::non_const_value_type& beta, const YMV& Y,
                     const typename ZMV::non_const_value_type& gamma, const ZMV& Z) {
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::"
                  "Update<rank 2>::update: X is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YMV>::value,
                  "KokkosBlas::Impl::"
                  "Update<rank 2>::update: Y is not a Kokkos::View.");
    static_assert(Kokkos::is_view<ZMV>::value,
                  "KokkosBlas::Impl::"
                  "Update<rank 2>::update: Z is not a Kokkos::View.");
    static_assert(std::is_same<typename ZMV::value_type, typename ZMV::non_const_value_type>::value,
                  "KokkosBlas::Impl::Update<rank 2>::update: Z is const.  "
                  "It must be nonconst, because it is an output argument "
                  "(we have to be able to write to its entries).");
    // Casting to int avoids compiler warnings about comparing
    // different kinds of enum values.
    static_assert((int)ZMV::rank == (int)XMV::rank && (int)ZMV::rank == (int)YMV::rank,
                  "KokkosBlas::Impl::Update<rank 2>::update: "
                  "X, Y, and Z must have the same rank.");
    static_assert(ZMV::rank == 2,
                  "KokkosBlas::Impl::Update<rank 2>::update: "
                  "XMV, YMV, and ZMV must have rank 2.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::update[ETI]"
                                                                     : "KokkosBlas::update[noETI]");

#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::update<> ETI specialization for < %s , %s , %s >\n", typeid(XMV).name(), typeid(YMV).name(),
             typeid(ZMV).name());
    else {
      printf("KokkosBlas1::update<> non-ETI specialization for < %s , %s , %s >\n", typeid(XMV).name(),
             typeid(YMV).name(), typeid(ZMV).name());
    }
#endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);
    int a = 2, b = 2, c = 2;

    if (alpha == ATA::zero()) {
      a = 0;
    } else {
      a = 2;
    }
    if (beta == ATB::zero()) {
      b = 0;
    } else {
      b = 2;
    }
    if (gamma == ATC::zero()) {
      c = 0;
    } else {
      c = 2;
    }

    if (numCols == static_cast<size_type>(1)) {
      // Special case: ZMV has rank 2, but only 1 column.
      // Dispatch to the rank-1 version for better performance.
      auto X_0 = Kokkos::subview(X, Kokkos::ALL(), 0);
      auto Y_0 = Kokkos::subview(Y, Kokkos::ALL(), 0);
      auto Z_0 = Kokkos::subview(Z, Kokkos::ALL(), 0);

      if (numRows * numCols < static_cast<size_type>(INT_MAX)) {
        typedef int index_type;
        V_Update_Generic<execution_space, decltype(X_0), decltype(Y_0), decltype(Z_0), index_type>(
            space, alpha, X_0, beta, Y_0, gamma, Z_0, a, b, c);
      } else {
        typedef typename XMV::size_type index_type;
        V_Update_Generic<execution_space, decltype(X_0), decltype(Y_0), decltype(Z_0), index_type>(
            space, alpha, X_0, beta, Y_0, gamma, Z_0, a, b, c);
      }
    } else {
      if (numRows * numCols < static_cast<size_type>(INT_MAX)) {
        typedef int index_type;
        MV_Update_Generic<execution_space, XMV, YMV, ZMV, index_type>(space, alpha, X, beta, Y, gamma, Z, a, b, c);
      } else {
        typedef typename XMV::size_type index_type;
        MV_Update_Generic<execution_space, XMV, YMV, ZMV, index_type>(space, alpha, X, beta, Y, gamma, Z, a, b, c);
      }
    }
    Kokkos::Profiling::popRegion();
  }
};

// Partial specialization for XV, YV, and ZV rank-1 Views.
template <class execution_space, class XV, class YV, class ZV>
struct Update<execution_space, XV, YV, ZV, 1, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename XV::size_type size_type;
  typedef Kokkos::ArithTraits<typename XV::non_const_value_type> ATA;
  typedef Kokkos::ArithTraits<typename YV::non_const_value_type> ATB;
  typedef Kokkos::ArithTraits<typename ZV::non_const_value_type> ATC;

  static void update(const execution_space& space, const typename XV::non_const_value_type& alpha, const XV& X,
                     const typename YV::non_const_value_type& beta, const YV& Y,
                     const typename ZV::non_const_value_type& gamma, const ZV& Z) {
    // XV, YV, and ZV must be Kokkos::View specializations.
    static_assert(Kokkos::is_view<XV>::value,
                  "KokkosBlas::Impl::"
                  "Update<rank 1>::update: X is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YV>::value,
                  "KokkosBlas::Impl::"
                  "Update<rank 1>::update: Y is not a Kokkos::View.");
    static_assert(Kokkos::is_view<ZV>::value,
                  "KokkosBlas::Impl::"
                  "Update<rank 1>::update: Z is not a Kokkos::View.");
    // ZV must be nonconst (else it can't be an output argument).
    static_assert(std::is_same<typename ZV::value_type, typename ZV::non_const_value_type>::value,
                  "KokkosBlas::Impl::Update<rank 1>::update: Z is const.  "
                  "It must be nonconst, because it is an output argument "
                  "(we have to be able to write to its entries).");
    static_assert((int)ZV::rank == (int)XV::rank && (int)ZV::rank == (int)YV::rank,
                  "KokkosBlas::Impl::Update<rank 1>::update: "
                  "X, Y, and Z must have the same rank.");
    static_assert(ZV::rank == 1,
                  "KokkosBlas::Impl::Update<rank 1>::update: "
                  "XV, YV, and ZV must have rank 1.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::update[ETI]"
                                                                     : "KokkosBlas::update[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::update<> ETI specialization for < %s , %s , %s >\n", typeid(XV).name(), typeid(YV).name(),
             typeid(ZV).name());
    else {
      printf("KokkosBlas1::update<> non-ETI specialization for < %s , %s , %s >\n", typeid(XV).name(),
             typeid(YV).name(), typeid(ZV).name());
    }
#endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);
    int a = 2, b = 2, c = 2;

    if (alpha == ATA::zero()) {
      a = 0;
    } else {
      a = 2;
    }
    if (beta == ATB::zero()) {
      b = 0;
    } else {
      b = 2;
    }
    if (gamma == ATC::zero()) {
      c = 0;
    } else {
      c = 2;
    }

    if (numRows < static_cast<size_type>(INT_MAX) && numRows * numCols < static_cast<size_type>(INT_MAX)) {
      typedef int index_type;
      V_Update_Generic<execution_space, XV, YV, ZV, index_type>(space, alpha, X, beta, Y, gamma, Z, a, b, c);
    } else {
      typedef typename XV::size_type index_type;
      V_Update_Generic<execution_space, XV, YV, ZV, index_type>(space, alpha, X, beta, Y, gamma, Z, a, b, c);
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
// KokkosBlas::Impl::Update for rank == 1.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _INST macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS1_UPDATE_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                       \
  extern template struct Update<                                                                                      \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      1, false, true>;

#define KOKKOSBLAS1_UPDATE_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                       \
  template struct Update<                                                                                             \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      1, false, true>;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Update for rank == 2.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _DEF macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS1_UPDATE_MV_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                     \
  extern template struct Update<                                                                                       \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      2, false, true>;

#define KOKKOSBLAS1_UPDATE_MV_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                     \
  template struct Update<                                                                                              \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      2, false, true>;

#include <KokkosBlas1_update_tpl_spec_decl.hpp>

#endif  // KOKKOSBLAS1_UPDATE_SPEC_HPP_
