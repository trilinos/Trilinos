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
#ifndef KOKKOSBLAS_IMPL_GESV_SPEC_HPP_
#define KOKKOSBLAS_IMPL_GESV_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas_gesv_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class AVT, class BVT>
struct gesv_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::GESV.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS_GESV_ETI_SPEC_AVAIL(SCALAR_TYPE, LAYOUT_TYPE,        \
                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
  template <>                                                           \
  struct gesv_eti_spec_avail<                                           \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE,                         \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,           \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE,                         \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {        \
    enum : bool { value = true };                                       \
  };

// Include the actual specialization declarations
#include <KokkosBlas_gesv_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas_gesv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosBlas::gesv.

template <class AMatrix, class BXMV, class IPIVV,
          bool tpl_spec_avail = gesv_tpl_spec_avail<AMatrix, BXMV>::value,
          bool eti_spec_avail = gesv_eti_spec_avail<AMatrix, BXMV>::value>
struct GESV {
  static void gesv(const AMatrix &A, const BXMV &B, const IPIVV &IPIV);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of gesv for multi vectors.
// Unification layer
template <class AMatrix, class BXMV, class IPIVV>
struct GESV<AMatrix, BXMV, IPIVV, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void gesv(const AMatrix & /* A */, const BXMV & /* B */,
                   const IPIVV & /* IPIV */) {
    // NOTE: Might add the implementation of KokkosBlas::gesv later
    throw std::runtime_error(
        "No fallback implementation of GESV (general LU factorization & solve) "
        "exists. Enable BLAS and/or MAGMA TPL.");
  }
};

#endif
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::GESV.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS_GESV_ETI_SPEC_DECL(SCALAR_TYPE, LAYOUT_TYPE,        \
                                      EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
  extern template struct GESV<                                         \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE,                        \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,          \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE,                        \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,          \
      Kokkos::View<int *, LAYOUT_TYPE,                                 \
                   Kokkos::Device<Kokkos::DefaultHostExecutionSpace,   \
                                  Kokkos::HostSpace>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,          \
      false, true>;

#define KOKKOSBLAS_GESV_ETI_SPEC_INST(SCALAR_TYPE, LAYOUT_TYPE,        \
                                      EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
  template struct GESV<                                                \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE,                        \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,          \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE,                        \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,          \
      Kokkos::View<int *, LAYOUT_TYPE,                                 \
                   Kokkos::Device<Kokkos::DefaultHostExecutionSpace,   \
                                  Kokkos::HostSpace>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,          \
      false, true>;

#include <KokkosBlas_gesv_tpl_spec_decl.hpp>
#include <generated_specializations_hpp/KokkosBlas_gesv_eti_spec_decl.hpp>

#endif  // KOKKOSBLAS_IMPL_GESV_SPEC_HPP_
