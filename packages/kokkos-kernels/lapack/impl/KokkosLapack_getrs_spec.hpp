// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_IMPL_GETRS_SPEC_HPP_
#define KOKKOSLAPACK_IMPL_GETRS_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_Error.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosLapack_getrs_impl.hpp>
#endif

namespace KokkosLapack {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class AMatrix, class IpivView, class BMatrix, class InfoView>
struct getrs_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosLapack

//
// Macro for declaration of full specialization availability
// KokkosLapack::Impl::GETRS.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSLAPACK_GETRS_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
  template <>                                                                                                      \
  struct getrs_eti_spec_avail<                                                                                     \
      EXEC_SPACE_TYPE,                                                                                             \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,             \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                       \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,             \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                       \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                       \
      Kokkos::View<int *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                     \
    enum : bool { value = true };                                                                                  \
  };

// Include the actual specialization declarations
#include <KokkosLapack_getrs_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosLapack_getrs_eti_spec_avail.hpp>

namespace KokkosLapack {
namespace Impl {

// Unification layer
template <class ExecutionSpace, class AMatrix, class IpivView, class BMatrix, class InfoView,
          bool tpl_spec_avail = getrs_tpl_spec_avail<ExecutionSpace, AMatrix, IpivView, BMatrix, InfoView>::value,
          bool eti_spec_avail = getrs_eti_spec_avail<ExecutionSpace, AMatrix, IpivView, BMatrix, InfoView>::value>
struct GETRS {
  static void getrs(const ExecutionSpace &space, const char trans[], const AMatrix &A, const IpivView &Ipiv,
                    const BMatrix &B, const InfoView &info);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
// Unification layer
template <class ExecutionSpace, class AMatrix, class IpivView, class BMatrix, class InfoView>
struct GETRS<ExecutionSpace, AMatrix, IpivView, BMatrix, InfoView, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void getrs(const ExecutionSpace &space, const char trans[], const AMatrix &A, const IpivView &Ipiv,
                    const BMatrix &B, const InfoView &Info) {
    std::string label = "KokkosLapack::getrs[NATIVE," +
                        KokkosKernels::ArithTraits<typename AMatrix::non_const_value_type>::name() + "]";
    Kokkos::Profiling::pushRegion(label);
    getrs_impl(space, trans, A, Ipiv, B, Info);
    Kokkos::Profiling::popRegion();
  }
};

#endif
}  // namespace Impl
}  // namespace KokkosLapack

//
// Macro for declaration of full specialization of
// KokkosLapack::Impl::GETRS.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSLAPACK_GETRS_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
  extern template struct GETRS<                                                                                   \
      EXEC_SPACE_TYPE,                                                                                            \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                      \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                      \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                      \
      Kokkos::View<int *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                      \
      false, true>;

#define KOKKOSLAPACK_GETRS_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
  template struct GETRS<                                                                                          \
      EXEC_SPACE_TYPE,                                                                                            \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                      \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                      \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                      \
      Kokkos::View<int *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                      \
      false, true>;

#include <KokkosLapack_getrs_tpl_spec_decl.hpp>
#include <generated_specializations_hpp/KokkosLapack_getrs_eti_spec_decl.hpp>

#endif  // KOKKOSLAPACK_IMPL_GETRS_SPEC_HPP_
