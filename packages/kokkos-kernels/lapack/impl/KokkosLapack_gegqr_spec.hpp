// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_IMPL_GEGQR_SPEC_HPP_
#define KOKKOSLAPACK_IMPL_GEGQR_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_ArithTraits.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosLapack_gegqr_impl.hpp>
#endif

namespace KokkosLapack {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class AMatrix, class TauArray, class InfoArray>
struct gegqr_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosLapack

//
// Macro for declaration of full specialization availability
// KokkosLapack::Impl::GEGQR.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSLAPACK_GEGQR_ETI_SPEC_AVAIL(SCALAR_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
  template <>                                                                                        \
  struct gegqr_eti_spec_avail<                                                                       \
      EXEC_SPACE_TYPE,                                                                               \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                         \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                         \
      Kokkos::View<int *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                       \
    enum : bool { value = true };                                                                    \
  };

// Include the actual specialization declarations
#include <KokkosLapack_gegqr_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosLapack_gegqr_eti_spec_avail.hpp>

namespace KokkosLapack {
namespace Impl {

// Unification layer
template <class ExecutionSpace, class AMatrix, class TauArray, class InfoArray,
          bool tpl_spec_avail = gegqr_tpl_spec_avail<ExecutionSpace, AMatrix, TauArray, InfoArray>::value,
          bool eti_spec_avail = gegqr_eti_spec_avail<ExecutionSpace, AMatrix, TauArray, InfoArray>::value>
struct GEGQR {
  static void gegqr(const ExecutionSpace &space, const int k, const AMatrix &A, const TauArray &Tau,
                    const InfoArray &info);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
// Unification layer
template <class ExecutionSpace, class AMatrix, class TauArray, class InfoArray>
struct GEGQR<ExecutionSpace, AMatrix, TauArray, InfoArray, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void gegqr(const ExecutionSpace & /* space */, const int /*k*/, const AMatrix & /* A */,
                    const TauArray & /* Tau */, const InfoArray & /* Info */) {
    // NOTE: Might add the implementation of KokkosLapack::gegqr later
    throw std::runtime_error(
        "No fallback implementation of GEGQR (compute Q from QR factorization) "
        "exists. Enable LAPACK, CUSOLVER, ROCSOLVER or MAGMA TPL.");
  }
};

#endif
}  // namespace Impl
}  // namespace KokkosLapack

//
// Macro for declaration of full specialization of
// KokkosLapack::Impl::GEGQR.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSLAPACK_GEGQR_ETI_SPEC_DECL(SCALAR_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
  extern template struct GEGQR<                                                                     \
      EXEC_SPACE_TYPE,                                                                              \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
      Kokkos::View<int *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,             \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
      false, true>;

#define KOKKOSLAPACK_GEGQR_ETI_SPEC_INST(SCALAR_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                \
  template struct GEGQR<EXEC_SPACE_TYPE,                                                                           \
                        Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                     \
                        Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,  \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                     \
                        Kokkos::View<int *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,          \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                     \
                        false, true>;

#include <KokkosLapack_gegqr_tpl_spec_decl.hpp>

#endif  // KOKKOSLAPACK_IMPL_GEGQR_SPEC_HPP_
