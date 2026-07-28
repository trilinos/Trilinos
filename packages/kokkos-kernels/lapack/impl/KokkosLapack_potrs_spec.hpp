// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_POTRS_SPEC_HPP_
#define KOKKOSLAPACK_POTRS_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosLapack_potrs_impl.hpp>
#endif

namespace KokkosLapack {
namespace Impl {
// Specialization struct which defines whether an ETI specialization exists
template <class ExecutionSpace, class AViewType, class BViewType>
struct potrs_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosLapack

//
// Macro for declaration of full specialization availability
// KokkosLapack::Impl::Potrs.  This is NOT for users!!!
//
#define KOKKOSLAPACK_POTRS_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                          \
  template <>                                                                                             \
  struct potrs_eti_spec_avail<EXEC_SPACE,                                                                 \
                              Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                      \
                              Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,       \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                    \
    enum : bool { value = true };                                                                         \
  };

// Include the actual specialization declarations
#include <KokkosLapack_potrs_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosLapack_potrs_eti_spec_avail.hpp>

namespace KokkosLapack {
namespace Impl {

//
// potrs
//

// Unification layer
template <class ExecutionSpace, class AViewType, class BViewType,
          bool tpl_spec_avail = potrs_tpl_spec_avail<ExecutionSpace, AViewType, BViewType>::value,
          bool eti_spec_avail = potrs_eti_spec_avail<ExecutionSpace, AViewType, BViewType>::value>
struct Potrs {
  static void potrs(const ExecutionSpace& space, const char uplo[], const AViewType& A, BViewType& B);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
// Default implementation when no TPL is available
template <class ExecutionSpace, class AViewType, class BViewType>
struct Potrs<ExecutionSpace, AViewType, BViewType, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void potrs(const ExecutionSpace& /* space */, const char uplo[], const AViewType& A, BViewType& B) {
    static_assert(Kokkos::is_view<AViewType>::value, "AViewType must be a Kokkos::View.");
    static_assert(Kokkos::is_view<BViewType>::value, "BViewType must be a Kokkos::View.");

    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosLapack::potrs[ETI]"
                                                                     : "KokkosLapack::potrs[noETI]");

    PotrsImpl<AViewType, BViewType>::potrs(uplo, A, B);

    Kokkos::Profiling::popRegion();
  }
};
#endif  // !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

}  // namespace Impl
}  // namespace KokkosLapack

//
// Macro for declaration of full specialization of
// KokkosLapack::Impl::Potrs.  This is NOT for users!!!
//
#define KOKKOSLAPACK_POTRS_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                       \
  extern template struct Potrs<                                                                                       \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      false, true>;

//
// Macro for definition of full specialization of
// KokkosLapack::Impl::Potrs.  This is NOT for users!!!
//
#define KOKKOSLAPACK_POTRS_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                       \
  template struct Potrs<                                                                                              \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      false, true>;

#include <KokkosLapack_potrs_tpl_spec_decl.hpp>
#include <generated_specializations_hpp/KokkosLapack_potrs_eti_spec_decl.hpp>

#endif  // KOKKOSLAPACK_POTRS_SPEC_HPP_
