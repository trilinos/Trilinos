// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_ROT_SPEC_HPP_
#define KOKKOSBLAS1_ROT_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_ArithTraits.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas1_rot_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class VectorView, class MagnitudeView, class ScalarView>
struct rot_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Rot.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_ROT_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXECSPACE, MEMSPACE)                                            \
  template <>                                                                                                          \
  struct rot_eti_spec_avail<                                                                                           \
      EXECSPACE,                                                                                                       \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,     \
      Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {    \
    enum : bool { value = true };                                                                                      \
  };

// Include the actual specialization declarations
#include <KokkosBlas1_rot_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_rot_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

// Unification layer
template <class ExecutionSpace, class VectorView, class MagnitudeView, class ScalarView,
          bool tpl_spec_avail = rot_tpl_spec_avail<ExecutionSpace, VectorView, MagnitudeView, ScalarView>::value,
          bool eti_spec_avail = rot_eti_spec_avail<ExecutionSpace, VectorView, MagnitudeView, ScalarView>::value>
struct Rot {
  static void rot(ExecutionSpace const& space, VectorView const& X, VectorView const& Y, MagnitudeView const& c,
                  ScalarView const& s);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of Rot.
template <class ExecutionSpace, class VectorView, class MagnitudeView, class ScalarView>
struct Rot<ExecutionSpace, VectorView, MagnitudeView, ScalarView, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void rot(ExecutionSpace const& space, VectorView const& X, VectorView const& Y, MagnitudeView const& c,
                  ScalarView const& s) {
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::rot[ETI]"
                                                                     : "KokkosBlas::rot[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::rot<> ETI specialization for < %s, %s, %s >\n", typeid(ExecutionSpace).name(),
             typeid(VectorView).name(), typeid(ScalarView).name());
    else {
      printf("KokkosBlas1::rot<> non-ETI specialization for < %s, %s, %s >\n", typeid(ExecutionSpace).name(),
             typeid(VectorView).name(), typeid(ScalarView).name());
    }
#endif
    Rot_Invoke<ExecutionSpace, VectorView, MagnitudeView, ScalarView>(space, X, Y, c, s);
    Kokkos::Profiling::popRegion();
  }
};
#endif

}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Rot.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_ROT_ETI_SPEC_DECL(SCALAR, LAYOUT, EXECSPACE, MEMSPACE)                                             \
  extern template struct Rot<                                                                                          \
      EXECSPACE,                                                                                                       \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,     \
      Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,      \
      false, true>;

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Rot.  This is NOT for users!!!  We
// use this macro in one or more .cpp files in this directory.
//
#define KOKKOSBLAS1_ROT_ETI_SPEC_INST(SCALAR, LAYOUT, EXECSPACE, MEMSPACE)                                             \
  template struct Rot<                                                                                                 \
      EXECSPACE,                                                                                                       \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,     \
      Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,      \
      false, true>;

#include <KokkosBlas1_rot_tpl_spec_decl.hpp>
#include <generated_specializations_hpp/KokkosBlas1_rot_eti_spec_decl.hpp>

#endif  // KOKKOSBLAS1_ROT_SPEC_HPP_
