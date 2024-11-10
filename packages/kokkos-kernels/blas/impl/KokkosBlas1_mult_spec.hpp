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
#ifndef KOKKOSBLAS1_MULT_SPEC_HPP_
#define KOKKOSBLAS1_MULT_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas1_mult_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class YMV, class AV, class XMV, int rank = XMV::rank>
struct mult_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Mult for rank == 1.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_MULT_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                        \
  template <>                                                                                                         \
  struct mult_eti_spec_avail<                                                                                         \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      1> {                                                                                                            \
    enum : bool { value = true };                                                                                     \
  };

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Mult for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_MULT_MV_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                      \
  template <>                                                                                                          \
  struct mult_eti_spec_avail<                                                                                          \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      2> {                                                                                                             \
    enum : bool { value = true };                                                                                      \
  };

// Include the actual specialization declarations
#include <KokkosBlas1_mult_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_mult_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_mult_mv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

//
// mult
//
/// \brief Implementation of entry-wise multiply of multivectors or
///   single vectors (depending on the rank template parameter).
///
/// Compute
///
/// Y(i,j) = alpha*A(i,j)*X(i,j) + gamma*Y(i,j)
///
/// with special cases for alpha, or gamma = 0.
template <class execution_space, class YMV, class AV, class XMV, int rank = XMV::rank,
          bool tpl_spec_avail = mult_tpl_spec_avail<execution_space, YMV, AV, XMV>::value,
          bool eti_spec_avail = mult_eti_spec_avail<execution_space, YMV, AV, XMV>::value>
struct Mult {
  static void mult(const execution_space& space, const typename YMV::non_const_value_type& gamma, const YMV& Y,
                   const typename XMV::non_const_value_type& alpha, const AV& A, const XMV& X);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
// Partial specialization for YMV, AV, and XMV rank-2 Views.
template <class execution_space, class YMV, class AV, class XMV>
struct Mult<execution_space, YMV, AV, XMV, 2, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename YMV::size_type size_type;
  typedef typename YMV::non_const_value_type YMV_scalar;
  typedef typename XMV::non_const_value_type XMV_scalar;

  static void mult(const execution_space& space, const YMV_scalar& gamma, const YMV& Y, const XMV_scalar& alpha,
                   const AV& A, const XMV& X) {
    static_assert(Kokkos::is_view<YMV>::value,
                  "KokkosBlas::Impl::"
                  "Mult<rank 2>::mult: Y is not a Kokkos::View.");
    static_assert(Kokkos::is_view<AV>::value,
                  "KokkosBlas::Impl::"
                  "Mult<rank 2>::mult: A is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::"
                  "Mult<rank 2>::mult: X is not a Kokkos::View.");
    static_assert(std::is_same<typename YMV::value_type, typename YMV::non_const_value_type>::value,
                  "KokkosBlas::Impl::Mult<rank 2>::mult: Y is const.  "
                  "It must be nonconst, because it is an output argument "
                  "(we have to be able to write to its entries).");
    // Casting to int avoids compiler warnings about comparing
    // different kinds of enum values.
    static_assert((int)XMV::rank == (int)YMV::rank && (int)XMV::rank == 2,
                  "KokkosBlas::Impl::Mult<rank 2>::mult: "
                  "X, and Y must have the rank 2.");
    static_assert(AV::rank == 1,
                  "KokkosBlas::Impl::Mult<rank 2>::mult: "
                  "AV must have rank 1.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::mult[ETI]"
                                                                     : "KokkosBlas::mult[noETI]");

#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::mult<> ETI specialization for < %s , %s , %s >\n", typeid(YMV).name(), typeid(AV).name(),
             typeid(XMV).name());
    else {
      printf("KokkosBlas1::mult<> non-ETI specialization for < %s , %s , %s >\n", typeid(YMV).name(), typeid(AV).name(),
             typeid(XMV).name());
    }
#endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);

    if (numRows < static_cast<int>(INT_MAX) && numRows * numCols < static_cast<int>(INT_MAX)) {
      MV_Mult_Generic<execution_space, YMV, AV, XMV, int>(space, gamma, Y, alpha, A, X);
    } else {
      MV_Mult_Generic<execution_space, YMV, AV, XMV, int64_t>(space, gamma, Y, alpha, A, X);
    }
    Kokkos::Profiling::popRegion();
  }
};

// Partial specialization for YV, AV, and XV rank-1 Views.
template <class execution_space, class YV, class AV, class XV>
struct Mult<execution_space, YV, AV, XV, 1, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename YV::size_type size_type;
  typedef typename YV::non_const_value_type YV_scalar;
  typedef typename XV::non_const_value_type XV_scalar;

  static void mult(const execution_space& space, const YV_scalar& gamma, const YV& Y, const XV_scalar& alpha,
                   const AV& A, const XV& X) {
    // YV, AV, and XV must be Kokkos::View specializations.
    static_assert(Kokkos::is_view<YV>::value,
                  "KokkosBlas::Impl::"
                  "Mult<rank 1>::mult: Y is not a Kokkos::View.");
    static_assert(Kokkos::is_view<AV>::value,
                  "KokkosBlas::Impl::"
                  "Mult<rank 1>::mult: A is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XV>::value,
                  "KokkosBlas::Impl::"
                  "Mult<rank 1>::mult: X is not a Kokkos::View.");
    // XV must be nonconst (else it can't be an output argument).
    static_assert(std::is_same<typename YV::value_type, typename YV::non_const_value_type>::value,
                  "KokkosBlas::Impl::Mult<rank 1>::mult: Y is const.  "
                  "It must be nonconst, because it is an output argument "
                  "(we have to be able to write to its entries).");
    static_assert((int)XV::rank == (int)YV::rank && (int)AV::rank == 1,
                  "KokkosBlas::Impl::Mult<rank 1>::mult: "
                  "X, Y, and Z must have rank 1.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::mult[ETI]"
                                                                     : "KokkosBlas::mult[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::mult<> ETI specialization for < %s , %s , %s >\n", typeid(YV).name(), typeid(AV).name(),
             typeid(XV).name());
    else {
      printf("KokkosBlas1::mult<> non-ETI specialization for < %s , %s , %s >\n", typeid(YV).name(), typeid(AV).name(),
             typeid(XV).name());
    }
#endif

    const size_type numRows = Y.extent(0);
    if (numRows < static_cast<int>(INT_MAX)) {
      V_Mult_Generic<execution_space, YV, AV, XV, int>(space, gamma, Y, alpha, A, X);
    } else {
      V_Mult_Generic<execution_space, YV, AV, XV, int64_t>(space, gamma, Y, alpha, A, X);
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
// KokkosBlas::Impl::Mult for rank == 1.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _INST macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS1_MULT_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                         \
  extern template struct Mult<                                                                                        \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      1, false, true>;

#define KOKKOSBLAS1_MULT_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                         \
  template struct Mult<                                                                                               \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      1, false, true>;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Mult for rank == 2.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _DEF macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS1_MULT_MV_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                       \
  extern template struct Mult<                                                                                         \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      2, false, true>;

#define KOKKOSBLAS1_MULT_MV_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                       \
  template struct Mult<                                                                                                \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      2, false, true>;

#include <KokkosBlas1_mult_tpl_spec_decl.hpp>

#endif  // KOKKOSBLAS1_MULT_SPEC_HPP_
