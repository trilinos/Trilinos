/*
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
*/

#ifndef KOKKOSBLAS1_SWAP_SPEC_HPP_
#define KOKKOSBLAS1_SWAP_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas1_swap_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class XVector, class YVector>
struct swap_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Swap.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_SWAP_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXECSPACE, MEMSPACE)                                         \
  template <>                                                                                                        \
  struct swap_eti_spec_avail<                                                                                        \
      EXECSPACE,                                                                                                     \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> { \
    enum : bool { value = true };                                                                                    \
  };

// Include the actual specialization declarations
#include <KokkosBlas1_swap_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_swap_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

// Unification layer
template <class ExecutionSpace, class XVector, class YVector,
          bool tpl_spec_avail = swap_tpl_spec_avail<ExecutionSpace, XVector, YVector>::value,
          bool eti_spec_avail = swap_eti_spec_avail<ExecutionSpace, XVector, YVector>::value>
struct Swap {
  static void swap(ExecutionSpace const& space, XVector const& X, YVector const& Y);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of Swap.
template <class ExecutionSpace, class XVector, class YVector>
struct Swap<ExecutionSpace, XVector, YVector, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void swap(ExecutionSpace const& space, XVector const& X, YVector const& Y) {
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::swap[ETI]"
                                                                     : "KokkosBlas::swap[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::swap<> ETI specialization for < %s, %s, %s >\n", typeid(ExecutionSpace).name(),
             typeid(XVector).name(), typeid(YVector).name());
    else {
      printf("KokkosBlas1::swap<> non-ETI specialization for < %s, %s, %s >\n", typeid(ExecutionSpace).name(),
             typeid(XVector).name(), typeid(YVector).name());
    }
#endif
    Swap_Invoke<ExecutionSpace, XVector, YVector>(space, X, Y);
    Kokkos::Profiling::popRegion();
  }
};
#endif

}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Swap.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_SWAP_ETI_SPEC_DECL(SCALAR, LAYOUT, EXECSPACE, MEMSPACE)                                        \
  extern template struct Swap<                                                                                     \
      EXECSPACE,                                                                                                   \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      false, true>;

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Swap.  This is NOT for users!!!  We
// use this macro in one or more .cpp files in this directory.
//
#define KOKKOSBLAS1_SWAP_ETI_SPEC_INST(SCALAR, LAYOUT, EXECSPACE, MEMSPACE)                                        \
  template struct Swap<                                                                                            \
      EXECSPACE,                                                                                                   \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      false, true>;

#include <KokkosBlas1_swap_tpl_spec_decl.hpp>

#endif  // KOKKOSBLAS1_SWAP_SPEC_HPP_
