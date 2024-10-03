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
#ifndef KOKKOSBLAS1_NRM2W_SPEC_HPP_
#define KOKKOSBLAS1_NRM2W_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas1_nrm2w_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class RMV, class XMV, int rank = XMV::rank>
struct nrm2w_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Nrm2w for rank == 1.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_NRM2W_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                  \
  template <>                                                                                                    \
  struct nrm2w_eti_spec_avail<EXEC_SPACE,                                                                        \
                              Kokkos::View<typename Kokkos::Details::InnerProductSpaceTraits<SCALAR>::mag_type,  \
                                           LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                              Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,         \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                            \
                              1> {                                                                               \
    enum : bool { value = true };                                                                                \
  };

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Nrm2w for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_NRM2W_MV_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                   \
  template <>                                                                                        \
  struct nrm2w_eti_spec_avail<                                                                       \
      EXEC_SPACE,                                                                                    \
      Kokkos::View<typename Kokkos::Details::InnerProductSpaceTraits<SCALAR>::mag_type*, LAYOUT,     \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                        \
      2> {                                                                                           \
    enum : bool { value = true };                                                                    \
  };

// Include the actual specialization declarations
#include <KokkosBlas1_nrm2w_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_nrm2w_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_nrm2w_mv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

// Unification layer
template <class execution_space, class RMV, class XMV, int rank = XMV::rank,
          bool tpl_spec_avail = nrm2w_tpl_spec_avail<execution_space, RMV, XMV>::value,
          bool eti_spec_avail = nrm2w_eti_spec_avail<execution_space, RMV, XMV>::value>
struct Nrm2w {
  static void nrm2w(const execution_space& space, const RMV& R, const XMV& X, const XMV& W, const bool& take_sqrt);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of Nrm2w for single vectors (1-D Views).
template <class execution_space, class RMV, class XMV>
struct Nrm2w<execution_space, RMV, XMV, 1, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  using size_type = typename XMV::size_type;

  static void nrm2w(const execution_space& space, const RMV& R, const XMV& X, const XMV& W, const bool& take_sqrt) {
    static_assert(Kokkos::is_view<RMV>::value,
                  "KokkosBlas::Impl::"
                  "Nrm2w<1-D>: RMV is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::"
                  "Nrm2w<1-D>: XMV is not a Kokkos::View.");
    static_assert(RMV::rank == 0,
                  "KokkosBlas::Impl::Nrm2w<1-D>: "
                  "RMV is not rank 0.");
    static_assert(XMV::rank == 1,
                  "KokkosBlas::Impl::Nrm2w<1-D>: "
                  "XMV is not rank 1.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::nrm2w[ETI]"
                                                                     : "KokkosBlas::nrm2w[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::nrm2w<> ETI specialization for < %s , %s >\n", typeid(RMV).name(), typeid(XMV).name());
    else {
      printf("KokkosBlas1::nrm2w<> non-ETI specialization for < %s , %s >\n", typeid(RMV).name(), typeid(XMV).name());
    }
#endif
    const size_type numRows = X.extent(0);

    if (numRows < static_cast<size_type>(INT_MAX)) {
      V_Nrm2w_Invoke<execution_space, RMV, XMV, int>(space, R, X, W, take_sqrt);
    } else {
      typedef std::int64_t index_type;
      V_Nrm2w_Invoke<execution_space, RMV, XMV, index_type>(space, R, X, W, take_sqrt);
    }
    Kokkos::Profiling::popRegion();
  }
};

template <class execution_space, class RV, class XMV>
struct Nrm2w<execution_space, RV, XMV, 2, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  using size_type = typename XMV::size_type;

  static void nrm2w(const execution_space& space, const RV& R, const XMV& X, const XMV& W, const bool& take_sqrt) {
    static_assert(Kokkos::is_view<RV>::value,
                  "KokkosBlas::Impl::"
                  "Nrm2w<2-D>: RV is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::"
                  "Nrm2w<2-D>: XMV is not a Kokkos::View.");
    static_assert(RV::rank == 1,
                  "KokkosBlas::Impl::Nrm2w<2-D>: "
                  "RV is not rank 1.");
    static_assert(XMV::rank == 2,
                  "KokkosBlas::Impl::Nrm2w<2-D>: "
                  "XMV is not rank 2.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::nrm2w[ETI]"
                                                                     : "KokkosBlas::nrm2w[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::nrm2w<> ETI specialization for < %s , %s >\n", typeid(RV).name(), typeid(XMV).name());
    else {
      printf("KokkosBlas1::nrm2w<> non-ETI specialization for < %s , %s >\n", typeid(RV).name(), typeid(XMV).name());
    }
#endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);
    if (numCols == 1) {
      auto R0 = Kokkos::subview(R, 0);
      auto X0 = Kokkos::subview(X, Kokkos::ALL(), 0);
      auto W0 = Kokkos::subview(W, Kokkos::ALL(), 0);
      if (numRows < static_cast<size_type>(INT_MAX)) {
        V_Nrm2w_Invoke<execution_space, decltype(R0), decltype(X0), int>(space, R0, X0, W0, take_sqrt);
      } else {
        typedef std::int64_t index_type;
        V_Nrm2w_Invoke<execution_space, decltype(R0), decltype(X0), index_type>(space, R0, X0, W0, take_sqrt);
      }
    } else {
      if (numRows < static_cast<size_type>(INT_MAX) && numRows * numCols < static_cast<size_type>(INT_MAX)) {
        MV_Nrm2w_Invoke<execution_space, RV, XMV, int>(space, R, X, W, take_sqrt);
      } else {
        typedef std::int64_t index_type;
        MV_Nrm2w_Invoke<execution_space, RV, XMV, index_type>(space, R, X, W, take_sqrt);
      }
    }
    Kokkos::Profiling::popRegion();
  }
};
#endif

}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Nrm2w for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_NRM2W_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                    \
  extern template struct Nrm2w<Kokkos::View<typename Kokkos::Details::InnerProductSpaceTraits<SCALAR>::mag_type,  \
                                            LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                               Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,         \
                                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                            \
                               1, false, true>;

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Nrm2w for rank == 2.  This is NOT for users!!!  We
// use this macro in one or more .cpp files in this directory.
//
#define KOKKOSBLAS1_NRM2W_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                    \
  template struct Nrm2w<EXEC_SPACE,                                                                               \
                        Kokkos::View<typename Kokkos::Details::InnerProductSpaceTraits<SCALAR>::mag_type, LAYOUT, \
                                     Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                \
                        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                   \
                        1, false, true>;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Nrm2w for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_NRM2W_MV_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                    \
  extern template struct Nrm2w<                                                                      \
      EXEC_SPACE,                                                                                    \
      Kokkos::View<typename Kokkos::Details::InnerProductSpaceTraits<SCALAR>::mag_type*, LAYOUT,     \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                        \
      2, false, true>;

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Nrm2w for rank == 2.  This is NOT for users!!!  We
// use this macro in one or more .cpp files in this directory.
//
#define KOKKOSBLAS1_NRM2W_MV_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                      \
  template struct Nrm2w<EXEC_SPACE,                                                                                    \
                        Kokkos::View<typename Kokkos::Details::InnerProductSpaceTraits<SCALAR>::mag_type*, LAYOUT,     \
                                     Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                    \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                        \
                        2, false, true>;

#include <KokkosBlas1_nrm2w_tpl_spec_decl.hpp>

#endif  // KOKKOSBLAS1_NRM2W_SPEC_HPP_
