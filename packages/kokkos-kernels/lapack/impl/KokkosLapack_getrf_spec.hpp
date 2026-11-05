// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_IMPL_GETRF_SPEC_HPP_
#define KOKKOSLAPACK_IMPL_GETRF_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_ArithTraits.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosLapack_getrf_impl.hpp>
#endif

namespace KokkosLapack {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class AMatrix, class IpivView, class InfoView>
struct getrf_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosLapack

//
// Macro for declaration of full specialization availability
// KokkosLapack::Impl::GETRF.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSLAPACK_GETRF_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
  template <>                                                                                                      \
  struct getrf_eti_spec_avail<                                                                                     \
      EXEC_SPACE_TYPE,                                                                                             \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                       \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                       \
      Kokkos::View<int *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                     \
    enum : bool { value = true };                                                                                  \
  };

// Include the actual specialization declarations
#include <KokkosLapack_getrf_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosLapack_getrf_eti_spec_avail.hpp>

namespace KokkosLapack {
namespace Impl {

// Unification layer
template <class ExecutionSpace, class AMatrix, class IpivView, class InfoView,
          bool tpl_spec_avail = getrf_tpl_spec_avail<ExecutionSpace, AMatrix, IpivView, InfoView>::value,
          bool eti_spec_avail = getrf_eti_spec_avail<ExecutionSpace, AMatrix, IpivView, InfoView>::value>
struct GETRF {
  static void getrf(const ExecutionSpace &space, const AMatrix &A, const IpivView &Ipiv, const InfoView &info);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
// Unification layer
template <class ExecutionSpace, class AMatrix, class IpivView, class InfoView>
struct GETRF<ExecutionSpace, AMatrix, IpivView, InfoView, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void getrf(const ExecutionSpace & /* space */, const AMatrix & /* A */, const IpivView & /* Tau */,
                    const InfoView & /* Info */) {
    // NOTE: Might add the implementation of KokkosLapack::getrf later
    throw std::runtime_error(
        "No fallback implementation of GETRF (general LU factorization) "
        "exists. Enable LAPACK, CUSOLVER, ROCSOLVER or MAGMA TPL.");
  }
};

#endif
}  // namespace Impl
}  // namespace KokkosLapack

//
// Macro for declaration of full specialization of
// KokkosLapack::Impl::GETRF.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSLAPACK_GETRF_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
  extern template struct GETRF<                                                                                   \
      EXEC_SPACE_TYPE,                                                                                            \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                      \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                      \
      Kokkos::View<int *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                      \
      false, true>;

#define KOKKOSLAPACK_GETRF_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE)  \
  template struct GETRF<EXEC_SPACE_TYPE,                                                                           \
                        Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                     \
                        Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                     \
                        Kokkos::View<int *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,          \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                     \
                        false, true>;

#include <KokkosLapack_getrf_tpl_spec_decl.hpp>

#endif  // KOKKOSLAPACK_IMPL_GETRF_SPEC_HPP_
