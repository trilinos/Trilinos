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
#ifndef KOKKOSBLAS1_ROTMG_SPEC_HPP_
#define KOKKOSBLAS1_ROTMG_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas1_rotmg_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class DXView, class YView, class PView>
struct rotmg_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Rotmg.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_ROTMG_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                     \
  template <>                                                                                                       \
  struct rotmg_eti_spec_avail<                                                                                      \
      EXEC_SPACE,                                                                                                   \
      Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<const SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                        \
      Kokkos::View<SCALAR[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                      \
    enum : bool { value = true };                                                                                   \
  };

// Include the actual specialization declarations
#include <KokkosBlas1_rotmg_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_rotmg_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

// Unification layer
template <class execution_space, class DXView, class YView, class PView,
          bool tpl_spec_avail = rotmg_tpl_spec_avail<execution_space, DXView, YView, PView>::value,
          bool eti_spec_avail = rotmg_eti_spec_avail<execution_space, DXView, YView, PView>::value>
struct Rotmg {
  static void rotmg(execution_space const& space, DXView& d1, DXView& d2, DXView& x1, YView& y1, PView& param);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of Rotmg.
template <class execution_space, class DXView, class YView, class PView>
struct Rotmg<execution_space, DXView, YView, PView, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void rotmg(execution_space const& space, DXView& d1, DXView& d2, DXView& x1, YView& y1, PView& param) {
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::rotmg[ETI]"
                                                                     : "KokkosBlas::rotmg[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::rotmg<> ETI specialization for < %s, %s, %s >\n", typeid(DXView).name(),
             typeid(YView).name(), typeid(PView).name());
    else {
      printf("KokkosBlas1::rotmg<> non-ETI specialization for < %s, %s, %s >\n", typeid(DXView).name(),
             typeid(YView).name(), typeid(PView).name());
    }
#endif
    Rotmg_Invoke<execution_space, DXView, YView, PView>(space, d1, d2, x1, y1, param);
    Kokkos::Profiling::popRegion();
  }
};
#endif

}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Rotmg.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_ROTMG_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                         \
  extern template struct Rotmg<                                                                                        \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      Kokkos::View<const SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<SCALAR[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      false, true>;

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Rotmg.  This is NOT for users!!!  We
// use this macro in one or more .cpp files in this directory.
//
#define KOKKOSBLAS1_ROTMG_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                         \
  template struct Rotmg<                                                                                               \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      Kokkos::View<const SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<SCALAR[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      false, true>;

#include <KokkosBlas1_rotmg_tpl_spec_decl.hpp>

#endif  // KOKKOSBLAS1_ROTMG_SPEC_HPP_
