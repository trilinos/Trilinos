// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBLAS1_ROTG_SPEC_HPP_
#define KOKKOSBLAS1_ROTG_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_ArithTraits.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas1_rotg_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class SViewType, class MViewType>
struct rotg_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Rotg.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_ROTG_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXECSPACE, MEMSPACE)                                           \
  template <>                                                                                                          \
  struct rotg_eti_spec_avail<                                                                                          \
      EXECSPACE,                                                                                                       \
      Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,      \
      Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                         \
    enum : bool { value = true };                                                                                      \
  };

// Include the actual specialization declarations
#include <KokkosBlas1_rotg_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_rotg_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

// Unification layer
template <class ExecutionSpace, class SViewType, class MViewType,
          bool tpl_spec_avail = rotg_tpl_spec_avail<ExecutionSpace, SViewType, MViewType>::value,
          bool eti_spec_avail = rotg_eti_spec_avail<ExecutionSpace, SViewType, MViewType>::value>
struct Rotg {
  static void rotg(ExecutionSpace const& space, SViewType const& a, SViewType const& b, MViewType const& c,
                   SViewType const& s);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of Rotg.
template <class ExecutionSpace, class SViewType, class MViewType>
struct Rotg<ExecutionSpace, SViewType, MViewType, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void rotg(ExecutionSpace const& space, SViewType const& a, SViewType const& b, MViewType const& c,
                   SViewType const& s) {
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::rotg[ETI]"
                                                                     : "KokkosBlas::rotg[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::rotg<> ETI specialization for < %s, %s, %s >\n", typeid(ExecutionSpace).name(),
             typeid(SViewType).name(), typeid(MViewType).name());
    else {
      printf("KokkosBlas1::rotg<> non-ETI specialization for < %s, %s, %s >\n", typeid(ExecutionSpace).name(),
             typeid(SViewType).name(), typeid(MViewType).name());
    }
#endif
    Rotg_Invoke<ExecutionSpace, SViewType, MViewType>(space, a, b, c, s);
    Kokkos::Profiling::popRegion();
  }
};
#endif

}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Rotg.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_ROTG_ETI_SPEC_DECL(SCALAR, LAYOUT, EXECSPACE, MEMSPACE)                                            \
  extern template struct Rotg<                                                                                         \
      EXECSPACE,                                                                                                       \
      Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,      \
      Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      false, true>;

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Rotg.  This is NOT for users!!!  We
// use this macro in one or more .cpp files in this directory.
//
#define KOKKOSBLAS1_ROTG_ETI_SPEC_INST(SCALAR, LAYOUT, EXECSPACE, MEMSPACE)                                            \
  template struct Rotg<                                                                                                \
      EXECSPACE,                                                                                                       \
      Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,      \
      Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      false, true>;

#include <KokkosBlas1_rotg_tpl_spec_decl.hpp>
#include <generated_specializations_hpp/KokkosBlas1_rotg_eti_spec_decl.hpp>

#endif  // KOKKOSBLAS1_ROTG_SPEC_HPP_
