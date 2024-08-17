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
#ifndef KOKKOSSPARSE_IMPL_PAR_ILUT_SYMBOLIC_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_PAR_ILUT_SYMBOLIC_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_Handle.hpp"

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosSparse_par_ilut_symbolic_impl.hpp>
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class KernelHandle, class ARowMapType, class AEntriesType, class LRowMapType, class URowMapType>
struct par_ilut_symbolic_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_PAR_ILUT_SYMBOLIC_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,       \
                                                      EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                           \
  template <>                                                                                                    \
  struct par_ilut_symbolic_eti_spec_avail<                                                                       \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> > > {                          \
    enum : bool { value = true };                                                                                \
  };

// Include the actual specialization declarations
#include <KokkosSparse_par_ilut_symbolic_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_par_ilut_symbolic_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosSparse::par_ilut_symbolic

template <
    class KernelHandle, class ARowMapType, class AEntriesType, class LRowMapType, class URowMapType,
    bool tpl_spec_avail =
        par_ilut_symbolic_tpl_spec_avail<KernelHandle, ARowMapType, AEntriesType, LRowMapType, URowMapType>::value,
    bool eti_spec_avail =
        par_ilut_symbolic_eti_spec_avail<KernelHandle, ARowMapType, AEntriesType, LRowMapType, URowMapType>::value>
struct PAR_ILUT_SYMBOLIC {
  static void par_ilut_symbolic(KernelHandle *handle, const ARowMapType &A_row_map, const AEntriesType &A_entries,
                                LRowMapType &L_row_map, URowMapType &U_row_map);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of par_ilut_symbolic
// Unification layer
template <class KernelHandle, class ARowMapType, class AEntriesType, class LRowMapType, class URowMapType>
struct PAR_ILUT_SYMBOLIC<KernelHandle, ARowMapType, AEntriesType, LRowMapType, URowMapType, false,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void par_ilut_symbolic(KernelHandle *handle, const ARowMapType &A_row_map, const AEntriesType &A_entries,
                                LRowMapType &L_row_map, URowMapType &U_row_map) {
    auto par_ilut_handle = handle->get_par_ilut_handle();

    Experimental::ilut_symbolic(*par_ilut_handle, A_row_map, A_entries, L_row_map, U_row_map);
  }
};
#endif
}  // namespace Impl
}  // namespace KokkosSparse

//
// Macro for declaration of full specialization of
// This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSSPARSE_PAR_ILUT_SYMBOLIC_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,        \
                                                     EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                            \
  extern template struct PAR_ILUT_SYMBOLIC<                                                                      \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      false, true>;

#define KOKKOSSPARSE_PAR_ILUT_SYMBOLIC_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,        \
                                                     EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                            \
  template struct PAR_ILUT_SYMBOLIC<                                                                             \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      false, true>;

#include <KokkosSparse_par_ilut_symbolic_tpl_spec_decl.hpp>

#endif
