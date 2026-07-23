// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_POTRF_SPEC_HPP_
#define KOKKOSLAPACK_POTRF_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosLapack_potrf_impl.hpp>
#endif

namespace KokkosLapack {
namespace Impl {
// Specialization struct which defines whether an ETI specialization exists
template <class ExecutionSpace, class AViewType>
struct potrf_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosLapack

//
// Macro for declaration of full specialization availability
// KokkosLapack::Impl::Potrf.  This is NOT for users!!!
//
#define KOKKOSLAPACK_POTRF_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                \
  template <>                                                                                                   \
  struct potrf_eti_spec_avail<EXEC_SPACE, Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {              \
    enum : bool { value = true };                                                                               \
  };

// Include the actual specialization declarations
#include <KokkosLapack_potrf_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosLapack_potrf_eti_spec_avail.hpp>

namespace KokkosLapack {
namespace Impl {

//
// potrf
//

// Unification layer
template <class ExecutionSpace, class AViewType,
          bool tpl_spec_avail = potrf_tpl_spec_avail<ExecutionSpace, AViewType>::value,
          bool eti_spec_avail = potrf_eti_spec_avail<ExecutionSpace, AViewType>::value>
struct Potrf {
  static void potrf(const ExecutionSpace& space, const char uplo[], AViewType& A);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
// Default implementation when no TPL is available
template <class ExecutionSpace, class AViewType>
struct Potrf<ExecutionSpace, AViewType, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void potrf(const ExecutionSpace& /* space */, const char uplo[], AViewType& A) {
    static_assert(Kokkos::is_view<AViewType>::value, "AViewType must be a Kokkos::View.");

    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosLapack::potrf[ETI]"
                                                                     : "KokkosLapack::potrf[noETI]");

    PotrfImpl<AViewType>::potrf(uplo, A);

    Kokkos::Profiling::popRegion();
  }
};
#endif  // !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

}  // namespace Impl
}  // namespace KokkosLapack

//
// Macro for declaration of full specialization of
// KokkosLapack::Impl::Potrf.  This is NOT for users!!!
//
#define KOKKOSLAPACK_POTRF_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                       \
  extern template struct Potrf<                                                                                       \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      false, true>;

//
// Macro for definition of full specialization of
// KokkosLapack::Impl::Potrf.  This is NOT for users!!!
//
#define KOKKOSLAPACK_POTRF_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                       \
  template struct Potrf<                                                                                              \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      false, true>;

#include <KokkosLapack_potrf_tpl_spec_decl.hpp>
#include <generated_specializations_hpp/KokkosLapack_potrf_eti_spec_decl.hpp>

#endif  // KOKKOSLAPACK_POTRF_SPEC_HPP_
