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
#ifndef KOKKOSLAPACK_IMPL_SVD_SPEC_HPP_
#define KOKKOSLAPACK_IMPL_SVD_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosLapack_svd_impl.hpp>
#endif

namespace KokkosLapack {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class AMatrix, class SVector, class UMatrix, class VMatrix>
struct svd_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosLapack

//
// Macro for declaration of full specialization availability
// KokkosLapack::Impl::SVD.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSLAPACK_SVD_ETI_SPEC_AVAIL(SCALAR_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE)            \
  template <>                                                                                                 \
  struct svd_eti_spec_avail<                                                                                  \
      EXEC_SPACE_TYPE,                                                                                        \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                  \
      Kokkos::View<Kokkos::ArithTraits<SCALAR_TYPE>::mag_type *, LAYOUT_TYPE,                                 \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                  \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                \
    enum : bool { value = true };                                                                             \
  };

// Include the actual specialization declarations
#include <KokkosLapack_svd_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosLapack_svd_eti_spec_avail.hpp>

namespace KokkosLapack {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosLapack::svd.

template <class ExecutionSpace, class AMatrix, class SVector, class UMatrix, class VMatrix,
          bool tpl_spec_avail = svd_tpl_spec_avail<ExecutionSpace, AMatrix, SVector, UMatrix, VMatrix>::value,
          bool eti_spec_avail = svd_eti_spec_avail<ExecutionSpace, AMatrix, SVector, UMatrix, VMatrix>::value>
struct SVD {
  static void svd(const ExecutionSpace &space, const char jobu[], const char jobvt[], const AMatrix &A,
                  const SVector &S, const UMatrix &U, const VMatrix &Vt);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of svd
// Unification layer
template <class ExecutionSpace, class AMatrix, class SVector, class UMatrix, class VMatrix>
struct SVD<ExecutionSpace, AMatrix, SVector, UMatrix, VMatrix, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void svd(const ExecutionSpace & /* space */, const char * /* jobu */, const char * /* jobvt */,
                  const AMatrix & /* A */, const SVector & /* S */, const UMatrix & /* U */, const VMatrix & /* Vt */) {
    // NOTE: Might add the implementation of KokkosLapack::svd later
    throw std::runtime_error(
        "No fallback implementation of SVD (singular value decomposition) "
        "exists. Enable LAPACK, CUSOLVER or ROCSOLVER TPL to use this "
        "function.");
  }
};

#endif
}  // namespace Impl
}  // namespace KokkosLapack

//
// Macro for declaration of full specialization of
// KokkosLapack::Impl::SVD.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSLAPACK_SVD_ETI_SPEC_DECL(SCALAR_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE)             \
  extern template struct SVD<                                                                                 \
      EXEC_SPACE_TYPE,                                                                                        \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                  \
      Kokkos::View<Kokkos::ArithTraits<SCALAR_TYPE>::mag_type *, LAYOUT_TYPE,                                 \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                  \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                  \
      false, true>;

#define KOKKOSLAPACK_SVD_ETI_SPEC_INST(SCALAR_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE)             \
  template struct SVD<                                                                                        \
      EXEC_SPACE_TYPE,                                                                                        \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                  \
      Kokkos::View<Kokkos::ArithTraits<SCALAR_TYPE>::mag_type *, LAYOUT_TYPE,                                 \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                  \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                  \
      false, true>;

#include <KokkosLapack_svd_tpl_spec_decl.hpp>

#endif  // KOKKOSLAPACK_IMPL_SVD_SPEC_HPP_
