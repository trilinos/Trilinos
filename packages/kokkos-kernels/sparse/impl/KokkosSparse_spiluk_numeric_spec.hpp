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
#ifndef KOKKOSSPARSE_IMPL_SPILUK_NUMERIC_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_SPILUK_NUMERIC_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_Handle.hpp"

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosSparse_spiluk_numeric_impl.hpp>
#include <KokkosSparse_spiluk_symbolic_impl.hpp>
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class KernelHandle, class ARowMapType, class AEntriesType, class AValuesType,
          class LRowMapType, class LEntriesType, class LValuesType, class URowMapType, class UEntriesType,
          class UValuesType>
struct spiluk_numeric_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPILUK_NUMERIC_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,          \
                                                   EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                              \
  template <>                                                                                                    \
  struct spiluk_numeric_eti_spec_avail<                                                                          \
      EXEC_SPACE_TYPE,                                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> > > {                          \
    enum : bool { value = true };                                                                                \
  };

// Include the actual specialization declarations
#include <KokkosSparse_spiluk_numeric_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_spiluk_numeric_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosSparse::spiluk_numeric

template <
    class ExecutionSpace, class KernelHandle, class ARowMapType, class AEntriesType, class AValuesType,
    class LRowMapType, class LEntriesType, class LValuesType, class URowMapType, class UEntriesType, class UValuesType,
    bool tpl_spec_avail =
        spiluk_numeric_tpl_spec_avail<ExecutionSpace, KernelHandle, ARowMapType, AEntriesType, AValuesType, LRowMapType,
                                      LEntriesType, LValuesType, URowMapType, UEntriesType, UValuesType>::value,
    bool eti_spec_avail =
        spiluk_numeric_eti_spec_avail<ExecutionSpace, KernelHandle, ARowMapType, AEntriesType, AValuesType, LRowMapType,
                                      LEntriesType, LValuesType, URowMapType, UEntriesType, UValuesType>::value>
struct SPILUK_NUMERIC {
  static void spiluk_numeric(KernelHandle *handle, const typename KernelHandle::const_nnz_lno_t &fill_lev,
                             const ARowMapType &A_row_map, const AEntriesType &A_entries, const AValuesType &A_values,
                             LRowMapType &L_row_map, LEntriesType &L_entries, LValuesType &L_values,
                             URowMapType &U_row_map, UEntriesType &U_entries, UValuesType &U_values);
  static void spiluk_numeric_streams(const std::vector<ExecutionSpace> &execspace_v,
                                     std::vector<KernelHandle> &handle_v, const std::vector<ARowMapType> &A_row_map_v,
                                     const std::vector<AEntriesType> &A_entries_v,
                                     const std::vector<AValuesType> &A_values_v,
                                     const std::vector<LRowMapType> &L_row_map_v,
                                     const std::vector<LEntriesType> &L_entries_v, std::vector<LValuesType> &L_values_v,
                                     const std::vector<URowMapType> &U_row_map_v,
                                     const std::vector<UEntriesType> &U_entries_v,
                                     std::vector<UValuesType> &U_values_v);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of spiluk_numeric
// Unification layer
template <class ExecutionSpace, class KernelHandle, class ARowMapType, class AEntriesType, class AValuesType,
          class LRowMapType, class LEntriesType, class LValuesType, class URowMapType, class UEntriesType,
          class UValuesType>
struct SPILUK_NUMERIC<ExecutionSpace, KernelHandle, ARowMapType, AEntriesType, AValuesType, LRowMapType, LEntriesType,
                      LValuesType, URowMapType, UEntriesType, UValuesType, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  using Iluk = Experimental::IlukWrap<typename KernelHandle::SPILUKHandleType>;

  static void spiluk_numeric(KernelHandle *handle, const typename KernelHandle::const_nnz_lno_t & /*fill_lev*/,
                             const ARowMapType &A_row_map, const AEntriesType &A_entries, const AValuesType &A_values,
                             LRowMapType &L_row_map, LEntriesType &L_entries, LValuesType &L_values,
                             URowMapType &U_row_map, UEntriesType &U_entries, UValuesType &U_values) {
    // Call specific algorithm type
    auto spiluk_handle = handle->get_spiluk_handle();

    Iluk::iluk_numeric(*spiluk_handle, A_row_map, A_entries, A_values, L_row_map, L_entries, L_values, U_row_map,
                       U_entries, U_values);
  }

  static void spiluk_numeric_streams(const std::vector<ExecutionSpace> &execspace_v,
                                     std::vector<KernelHandle> &handle_v, const std::vector<ARowMapType> &A_row_map_v,
                                     const std::vector<AEntriesType> &A_entries_v,
                                     const std::vector<AValuesType> &A_values_v,
                                     const std::vector<LRowMapType> &L_row_map_v,
                                     const std::vector<LEntriesType> &L_entries_v, std::vector<LValuesType> &L_values_v,
                                     const std::vector<URowMapType> &U_row_map_v,
                                     const std::vector<UEntriesType> &U_entries_v,
                                     std::vector<UValuesType> &U_values_v) {
    std::vector<typename KernelHandle::SPILUKHandleType *> spiluk_handle_v(execspace_v.size());
    for (int i = 0; i < static_cast<int>(execspace_v.size()); i++) {
      spiluk_handle_v[i] = handle_v[i].get_spiluk_handle();
    }

    Iluk::iluk_numeric_streams(execspace_v, spiluk_handle_v, A_row_map_v, A_entries_v, A_values_v, L_row_map_v,
                               L_entries_v, L_values_v, U_row_map_v, U_entries_v, U_values_v);
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
#define KOKKOSSPARSE_SPILUK_NUMERIC_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,           \
                                                  EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                               \
  extern template struct SPILUK_NUMERIC<                                                                         \
      EXEC_SPACE_TYPE,                                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      false, true>;

#define KOKKOSSPARSE_SPILUK_NUMERIC_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,           \
                                                  EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                               \
  template struct SPILUK_NUMERIC<                                                                                \
      EXEC_SPACE_TYPE,                                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      false, true>;

#include <KokkosSparse_spiluk_numeric_tpl_spec_decl.hpp>

#endif
