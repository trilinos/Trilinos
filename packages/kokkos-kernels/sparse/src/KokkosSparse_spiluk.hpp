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

/// \file KokkosSparse_spiluk.hpp
/// \brief Parallel incomplete LU factorization ILU(k)
///
/// This file provides KokkosSparse::spiluk.  This function performs a
/// local (no MPI) sparse ILU(k) on matrices stored in
/// compressed row sparse ("Crs") format.

#ifndef KOKKOSSPARSE_SPILUK_HPP_
#define KOKKOSSPARSE_SPILUK_HPP_

#include <type_traits>

// #include "KokkosSparse_spiluk_handle.hpp"
#include "KokkosKernels_helpers.hpp"
#include "KokkosKernels_Error.hpp"
#include "KokkosSparse_spiluk_symbolic_spec.hpp"
#include "KokkosSparse_spiluk_numeric_spec.hpp"

namespace KokkosSparse {
namespace Experimental {

#define KOKKOSKERNELS_SPILUK_SAME_TYPE(A, B) \
  std::is_same<typename std::remove_const<A>::type, typename std::remove_const<B>::type>::value

template <typename KernelHandle, typename ARowMapType, typename AEntriesType, typename LRowMapType,
          typename LEntriesType, typename URowMapType, typename UEntriesType>
void spiluk_symbolic(KernelHandle* handle, typename KernelHandle::const_nnz_lno_t fill_lev, ARowMapType& A_rowmap,
                     AEntriesType& A_entries, LRowMapType& L_rowmap, LEntriesType& L_entries, URowMapType& U_rowmap,
                     UEntriesType& U_entries, int nstreams = 1) {
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::nnz_lno_t ordinal_type;

  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename ARowMapType::non_const_value_type, size_type),
                "spiluk_symbolic: A size_type must match KernelHandle "
                "size_type (const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename AEntriesType::non_const_value_type, ordinal_type),
                "spiluk_symbolic: A entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");

  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename LRowMapType::non_const_value_type, size_type),
                "spiluk_symbolic: L size_type must match KernelHandle "
                "size_type (const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename LEntriesType::non_const_value_type, ordinal_type),
                "spiluk_symbolic: L entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");

  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename URowMapType::non_const_value_type, size_type),
                "spiluk_symbolic: U size_type must match KernelHandle "
                "size_type (const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename UEntriesType::non_const_value_type, ordinal_type),
                "spiluk_symbolic: U entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");

  static_assert(Kokkos::is_view<ARowMapType>::value, "spiluk_symbolic: A_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<AEntriesType>::value, "spiluk_symbolic: A_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LRowMapType>::value, "spiluk_symbolic: L_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LEntriesType>::value, "spiluk_symbolic: L_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<URowMapType>::value, "spiluk_symbolic: U_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<UEntriesType>::value, "spiluk_symbolic: U_entries is not a Kokkos::View.");

  static_assert((int)LRowMapType::rank == (int)ARowMapType::rank,
                "spiluk_symbolic: The ranks of L_rowmap and A_rowmap do not match.");
  static_assert((int)LEntriesType::rank == (int)AEntriesType::rank,
                "spiluk_symbolic: The ranks of L_entries and A_entries do not match.");

  static_assert((int)LRowMapType::rank == (int)URowMapType::rank,
                "spiluk_symbolic: The ranks of L_rowmap and U_rowmap do not match.");
  static_assert((int)LEntriesType::rank == (int)UEntriesType::rank,
                "spiluk_symbolic: The ranks of L_entries and U_entries do not match.");

  static_assert(LRowMapType::rank == 1, "spiluk_symbolic: A_rowmap, L_rowmap and U_rowmap must all have rank 1.");
  static_assert(LEntriesType::rank == 1,
                "spiluk_symbolic: A_entries, L_entries and U_entries must all "
                "have rank 1.");

  static_assert(std::is_same<typename LRowMapType::value_type, typename LRowMapType::non_const_value_type>::value,
                "spiluk_symbolic: The output L_rowmap must be nonconst.");
  static_assert(std::is_same<typename LEntriesType::value_type, typename LEntriesType::non_const_value_type>::value,
                "spiluk_symbolic: The output L_entries must be nonconst.");
  static_assert(std::is_same<typename URowMapType::value_type, typename URowMapType::non_const_value_type>::value,
                "spiluk_symbolic: The output U_rowmap must be nonconst.");
  static_assert(std::is_same<typename UEntriesType::value_type, typename UEntriesType::non_const_value_type>::value,
                "spiluk_symbolic: The output U_entries must be nonconst.");

  static_assert(std::is_same<typename LRowMapType::device_type, typename ARowMapType::device_type>::value,
                "spiluk_symbolic: Views LRowMapType and ARowMapType have "
                "different device_types.");
  static_assert(std::is_same<typename LEntriesType::device_type, typename AEntriesType::device_type>::value,
                "spiluk_symbolic: Views LEntriesType and AEntriesType have "
                "different device_types.");

  static_assert(std::is_same<typename LRowMapType::device_type, typename URowMapType::device_type>::value,
                "spiluk_symbolic: Views LRowMapType and URowMapType have "
                "different device_types.");
  static_assert(std::is_same<typename LEntriesType::device_type, typename UEntriesType::device_type>::value,
                "spiluk_symbolic: Views LEntriesType and UEntriesType have "
                "different device_types.");

  static_assert(std::is_same<typename LRowMapType::device_type::execution_space,
                             typename KernelHandle::SPILUKHandleType::execution_space>::value,
                "spiluk_symbolic: KernelHandle and Views have different execution "
                "spaces.");
  static_assert(std::is_same<typename LEntriesType::device_type::execution_space,
                             typename KernelHandle::SPILUKHandleType::execution_space>::value,
                "spiluk_symbolic: KernelHandle and Views have different execution "
                "spaces.");

  static_assert(std::is_same<typename LRowMapType::device_type, typename LEntriesType::device_type>::value,
                "spiluk_symbolic: rowmap and entries have different device types.");

  // Check validity of fill level
  if (fill_lev < 0) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::spiluk_symbolic: fill_lev: " << fill_lev << ". Valid value is >= 0.";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  typedef typename KernelHandle::const_size_type c_size_t;
  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

  typedef typename KernelHandle::HandleExecSpace c_exec_t;
  typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
  typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t,
                                                                    c_persist_t>
      const_handle_type;
  const_handle_type tmp_handle(*handle);

  typedef Kokkos::View<typename ARowMapType::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<ARowMapType>::array_layout,
                       typename ARowMapType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      ARowMap_Internal;

  typedef Kokkos::View<typename AEntriesType::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<AEntriesType>::array_layout,
                       typename AEntriesType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      AEntries_Internal;

  typedef Kokkos::View<typename LRowMapType::non_const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<LRowMapType>::array_layout,
                       typename LRowMapType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      LRowMap_Internal;

  typedef Kokkos::View<typename LEntriesType::non_const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<LEntriesType>::array_layout,
                       typename LEntriesType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      LEntries_Internal;

  typedef Kokkos::View<typename URowMapType::non_const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<URowMapType>::array_layout,
                       typename URowMapType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      URowMap_Internal;

  typedef Kokkos::View<typename UEntriesType::non_const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<UEntriesType>::array_layout,
                       typename UEntriesType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      UEntries_Internal;

  ARowMap_Internal A_rowmap_i   = A_rowmap;
  AEntries_Internal A_entries_i = A_entries;
  LRowMap_Internal L_rowmap_i   = L_rowmap;
  LEntries_Internal L_entries_i = L_entries;
  URowMap_Internal U_rowmap_i   = U_rowmap;
  UEntries_Internal U_entries_i = U_entries;

  KokkosSparse::Impl::SPILUK_SYMBOLIC<const_handle_type, ARowMap_Internal, AEntries_Internal, LRowMap_Internal,
                                      LEntries_Internal, URowMap_Internal,
                                      UEntries_Internal>::spiluk_symbolic(&tmp_handle, fill_lev, A_rowmap_i,
                                                                          A_entries_i, L_rowmap_i, L_entries_i,
                                                                          U_rowmap_i, U_entries_i, nstreams);

}  // spiluk_symbolic

template <typename KernelHandle, typename ARowMapType, typename AEntriesType, typename AValuesType,
          typename LRowMapType, typename LEntriesType, typename LValuesType, typename URowMapType,
          typename UEntriesType, typename UValuesType>
void spiluk_numeric(KernelHandle* handle, typename KernelHandle::const_nnz_lno_t fill_lev, ARowMapType& A_rowmap,
                    AEntriesType& A_entries, AValuesType& A_values, LRowMapType& L_rowmap, LEntriesType& L_entries,
                    LValuesType& L_values, URowMapType& U_rowmap, UEntriesType& U_entries, UValuesType& U_values) {
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::nnz_lno_t ordinal_type;
  typedef typename KernelHandle::nnz_scalar_t scalar_type;

  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename ARowMapType::non_const_value_type, size_type),
                "spiluk_numeric: A size_type must match KernelHandle size_type "
                "(const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename AEntriesType::non_const_value_type, ordinal_type),
                "spiluk_numeric: A entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename AValuesType::value_type, scalar_type),
                "spiluk_numeric: A scalar type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");

  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename LRowMapType::non_const_value_type, size_type),
                "spiluk_numeric: L size_type must match KernelHandle size_type "
                "(const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename LEntriesType::non_const_value_type, ordinal_type),
                "spiluk_numeric: L entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename LValuesType::value_type, scalar_type),
                "spiluk_numeric: L scalar type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");

  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename URowMapType::non_const_value_type, size_type),
                "spiluk_numeric: U size_type must match KernelHandle size_type "
                "(const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename UEntriesType::non_const_value_type, ordinal_type),
                "spiluk_numeric: U entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename UValuesType::value_type, scalar_type),
                "spiluk_numeric: U scalar type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");

  static_assert(Kokkos::is_view<ARowMapType>::value, "spiluk_numeric: A_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<AEntriesType>::value, "spiluk_numeric: A_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<AValuesType>::value, "spiluk_numeric: A_values is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LRowMapType>::value, "spiluk_numeric: L_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LEntriesType>::value, "spiluk_numeric: L_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LValuesType>::value, "spiluk_numeric: L_values is not a Kokkos::View.");
  static_assert(Kokkos::is_view<URowMapType>::value, "spiluk_numeric: U_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<UEntriesType>::value, "spiluk_numeric: U_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<UValuesType>::value, "spiluk_numeric: U_values is not a Kokkos::View.");

  static_assert((int)LRowMapType::rank == (int)ARowMapType::rank,
                "spiluk_numeric: The ranks of L_rowmap and A_rowmap do not match.");
  static_assert((int)LEntriesType::rank == (int)AEntriesType::rank,
                "spiluk_numeric: The ranks of L_entries and A_entries do not match.");
  static_assert((int)LValuesType::rank == (int)AValuesType::rank,
                "spiluk_numeric: The ranks of L_values and A_values do not match.");

  static_assert((int)LRowMapType::rank == (int)URowMapType::rank,
                "spiluk_numeric: The ranks of L_rowmap and U_rowmap do not match.");
  static_assert((int)LEntriesType::rank == (int)UEntriesType::rank,
                "spiluk_numeric: The ranks of L_entries and U_entries do not match.");
  static_assert((int)LValuesType::rank == (int)UValuesType::rank,
                "spiluk_numeric: The ranks of L_values and U_values do not match.");

  static_assert(LRowMapType::rank == 1, "spiluk_numeric: A_rowmap, L_rowmap and U_rowmap must all have rank 1.");
  static_assert(LEntriesType::rank == 1,
                "spiluk_numeric: A_entries, L_entries and U_entries must all "
                "have rank 1.");
  static_assert(LValuesType::rank == 1, "spiluk_numeric: A_values, L_values and U_values must all have rank 1.");

  static_assert(std::is_same<typename LEntriesType::value_type, typename LEntriesType::non_const_value_type>::value,
                "spiluk_numeric: The output L_entries must be nonconst.");
  static_assert(std::is_same<typename LValuesType::value_type, typename LValuesType::non_const_value_type>::value,
                "spiluk_numeric: The output L_values must be nonconst.");
  static_assert(std::is_same<typename UEntriesType::value_type, typename UEntriesType::non_const_value_type>::value,
                "spiluk_numeric: The output U_entries must be nonconst.");
  static_assert(std::is_same<typename UValuesType::value_type, typename UValuesType::non_const_value_type>::value,
                "spiluk_numeric: The output U_values must be nonconst.");

  static_assert(std::is_same<typename LRowMapType::device_type, typename ARowMapType::device_type>::value,
                "spiluk_numeric: Views LRowMapType and ARowMapType have "
                "different device_types.");
  static_assert(std::is_same<typename LEntriesType::device_type, typename AEntriesType::device_type>::value,
                "spiluk_numeric: Views LEntriesType and AEntriesType have "
                "different device_types.");
  static_assert(std::is_same<typename LValuesType::device_type, typename AValuesType::device_type>::value,
                "spiluk_numeric: Views LValuesType and AValuesType have "
                "different device_types.");

  static_assert(std::is_same<typename LRowMapType::device_type, typename URowMapType::device_type>::value,
                "spiluk_numeric: Views LRowMapType and URowMapType have "
                "different device_types.");
  static_assert(std::is_same<typename LEntriesType::device_type, typename UEntriesType::device_type>::value,
                "spiluk_numeric: Views LEntriesType and UEntriesType have "
                "different device_types.");
  static_assert(std::is_same<typename LValuesType::device_type, typename UValuesType::device_type>::value,
                "spiluk_numeric: Views LValuesType and UValuesType have "
                "different device_types.");

  static_assert(std::is_same<typename LRowMapType::device_type::execution_space,
                             typename KernelHandle::SPILUKHandleType::execution_space>::value,
                "spiluk_numeric: KernelHandle and Views have different execution "
                "spaces.");
  static_assert(std::is_same<typename LEntriesType::device_type::execution_space,
                             typename KernelHandle::SPILUKHandleType::execution_space>::value,
                "spiluk_numeric: KernelHandle and Views have different execution "
                "spaces.");
  static_assert(std::is_same<typename LValuesType::device_type::execution_space,
                             typename KernelHandle::SPILUKHandleType::execution_space>::value,
                "spiluk_numeric: KernelHandle and Views have different execution "
                "spaces.");

  static_assert(std::is_same<typename LRowMapType::device_type, typename LEntriesType::device_type>::value,
                "spiluk_numeric: rowmap and entries have different device types.");
  static_assert(std::is_same<typename LRowMapType::device_type, typename LValuesType::device_type>::value,
                "spiluk_numeric: rowmap and values have different device types.");

  // Check validity of fill level
  if (fill_lev < 0) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::spiluk_numeric: fill_lev: " << fill_lev << ". Valid value is >= 0.";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  // Check if symbolic has been called
  if (handle->get_spiluk_handle()->is_symbolic_complete() == false) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::spiluk_numeric: spiluk_symbolic must be "
          "called before spiluk_numeric.";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  typedef typename KernelHandle::const_size_type c_size_t;
  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

  typedef typename KernelHandle::HandleExecSpace c_exec_t;
  typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
  typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t,
                                                                    c_persist_t>
      const_handle_type;
  const_handle_type tmp_handle(*handle);

  typedef Kokkos::View<typename ARowMapType::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<ARowMapType>::array_layout,
                       typename ARowMapType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      ARowMap_Internal;

  typedef Kokkos::View<typename AEntriesType::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<AEntriesType>::array_layout,
                       typename AEntriesType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      AEntries_Internal;

  typedef Kokkos::View<typename AValuesType::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<AValuesType>::array_layout,
                       typename AValuesType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      AValues_Internal;

  typedef Kokkos::View<typename LRowMapType::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<LRowMapType>::array_layout,
                       typename LRowMapType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      LRowMap_Internal;

  typedef Kokkos::View<typename LEntriesType::non_const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<LEntriesType>::array_layout,
                       typename LEntriesType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      LEntries_Internal;

  typedef Kokkos::View<typename LValuesType::non_const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<LValuesType>::array_layout,
                       typename LValuesType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      LValues_Internal;

  typedef Kokkos::View<typename URowMapType::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<URowMapType>::array_layout,
                       typename URowMapType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      URowMap_Internal;

  typedef Kokkos::View<typename UEntriesType::non_const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<UEntriesType>::array_layout,
                       typename UEntriesType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      UEntries_Internal;

  typedef Kokkos::View<typename UValuesType::non_const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<UValuesType>::array_layout,
                       typename UValuesType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      UValues_Internal;

  ARowMap_Internal A_rowmap_i   = A_rowmap;
  AEntries_Internal A_entries_i = A_entries;
  AValues_Internal A_values_i   = A_values;
  LRowMap_Internal L_rowmap_i   = L_rowmap;
  LEntries_Internal L_entries_i = L_entries;
  LValues_Internal L_values_i   = L_values;
  URowMap_Internal U_rowmap_i   = U_rowmap;
  UEntries_Internal U_entries_i = U_entries;
  UValues_Internal U_values_i   = U_values;

  KokkosSparse::Impl::SPILUK_NUMERIC<typename AValuesType::execution_space, const_handle_type, ARowMap_Internal,
                                     AEntries_Internal, AValues_Internal, LRowMap_Internal, LEntries_Internal,
                                     LValues_Internal, URowMap_Internal, UEntries_Internal,
                                     UValues_Internal>::spiluk_numeric(&tmp_handle, fill_lev, A_rowmap_i, A_entries_i,
                                                                       A_values_i, L_rowmap_i, L_entries_i, L_values_i,
                                                                       U_rowmap_i, U_entries_i, U_values_i);
}  // spiluk_numeric

template <class ExecutionSpace, typename KernelHandle, typename ARowMapType, typename AEntriesType,
          typename AValuesType, typename LRowMapType, typename LEntriesType, typename LValuesType, typename URowMapType,
          typename UEntriesType, typename UValuesType>
void spiluk_numeric_streams(const std::vector<ExecutionSpace>& execspace_v, const std::vector<KernelHandle*>& handle_v,
                            typename KernelHandle::const_nnz_lno_t fill_lev, const std::vector<ARowMapType>& A_rowmap_v,
                            const std::vector<AEntriesType>& A_entries_v, const std::vector<AValuesType>& A_values_v,
                            const std::vector<LRowMapType>& L_rowmap_v, const std::vector<LEntriesType>& L_entries_v,
                            std::vector<LValuesType>& L_values_v, const std::vector<URowMapType>& U_rowmap_v,
                            const std::vector<UEntriesType>& U_entries_v, std::vector<UValuesType>& U_values_v) {
  using size_type    = typename KernelHandle::size_type;
  using ordinal_type = typename KernelHandle::nnz_lno_t;
  using scalar_type  = typename KernelHandle::nnz_scalar_t;

  static_assert(Kokkos::is_execution_space<ExecutionSpace>::value, "ExecutionSpace is not valid");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename ARowMapType::memory_space>::accessible,
                "spiluk_numeric_streams: ExecutionSpace cannot access data in "
                "ARowMapType");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename AEntriesType::memory_space>::accessible,
                "spiluk_numeric_streams: ExecutionSpace cannot access data in "
                "AEntriesType");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename AValuesType::memory_space>::accessible,
                "spiluk_numeric_streams: ExecutionSpace cannot access data in "
                "AValuesType");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename LRowMapType::memory_space>::accessible,
                "spiluk_numeric_streams: ExecutionSpace cannot access data in "
                "LRowMapType");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename LEntriesType::memory_space>::accessible,
                "spiluk_numeric_streams: ExecutionSpace cannot access data in "
                "LEntriesType");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename LValuesType::memory_space>::accessible,
                "spiluk_numeric_streams: ExecutionSpace cannot access data in "
                "LValuesType");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename URowMapType::memory_space>::accessible,
                "spiluk_numeric_streams: ExecutionSpace cannot access data in "
                "URowMapType");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename UEntriesType::memory_space>::accessible,
                "spiluk_numeric_streams: ExecutionSpace cannot access data in "
                "UEntriesType");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename UValuesType::memory_space>::accessible,
                "spiluk_numeric_streams: ExecutionSpace cannot access data in "
                "UValuesType");

  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename ARowMapType::non_const_value_type, size_type),
                "spiluk_numeric_streams: A size_type must match KernelHandle size_type "
                "(const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename AEntriesType::non_const_value_type, ordinal_type),
                "spiluk_numeric_streams: A entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename AValuesType::value_type, scalar_type),
                "spiluk_numeric_streams: A scalar type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");

  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename LRowMapType::non_const_value_type, size_type),
                "spiluk_numeric_streams: L size_type must match KernelHandle size_type "
                "(const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename LEntriesType::non_const_value_type, ordinal_type),
                "spiluk_numeric_streams: L entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename LValuesType::value_type, scalar_type),
                "spiluk_numeric_streams: L scalar type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");

  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename URowMapType::non_const_value_type, size_type),
                "spiluk_numeric_streams: U size_type must match KernelHandle size_type "
                "(const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename UEntriesType::non_const_value_type, ordinal_type),
                "spiluk_numeric_streams: U entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPILUK_SAME_TYPE(typename UValuesType::value_type, scalar_type),
                "spiluk_numeric_streams: U scalar type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");

  static_assert(Kokkos::is_view<ARowMapType>::value, "spiluk_numeric_streams: A_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<AEntriesType>::value, "spiluk_numeric_streams: A_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<AValuesType>::value, "spiluk_numeric_streams: A_values is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LRowMapType>::value, "spiluk_numeric_streams: L_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LEntriesType>::value, "spiluk_numeric_streams: L_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LValuesType>::value, "spiluk_numeric_streams: L_values is not a Kokkos::View.");
  static_assert(Kokkos::is_view<URowMapType>::value, "spiluk_numeric_streams: U_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<UEntriesType>::value, "spiluk_numeric_streams: U_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<UValuesType>::value, "spiluk_numeric_streams: U_values is not a Kokkos::View.");

  static_assert((int)LRowMapType::rank == (int)ARowMapType::rank,
                "spiluk_numeric_streams: The ranks of L_rowmap and A_rowmap do "
                "not match.");
  static_assert((int)LEntriesType::rank == (int)AEntriesType::rank,
                "spiluk_numeric_streams: The ranks of L_entries and A_entries "
                "do not match.");
  static_assert((int)LValuesType::rank == (int)AValuesType::rank,
                "spiluk_numeric_streams: The ranks of L_values and A_values do "
                "not match.");

  static_assert((int)LRowMapType::rank == (int)URowMapType::rank,
                "spiluk_numeric_streams: The ranks of L_rowmap and U_rowmap do "
                "not match.");
  static_assert((int)LEntriesType::rank == (int)UEntriesType::rank,
                "spiluk_numeric_streams: The ranks of L_entries and U_entries "
                "do not match.");
  static_assert((int)LValuesType::rank == (int)UValuesType::rank,
                "spiluk_numeric_streams: The ranks of L_values and U_values do "
                "not match.");

  static_assert(LRowMapType::rank == 1,
                "spiluk_numeric_streams: A_rowmap, L_rowmap and U_rowmap must "
                "all have rank 1.");
  static_assert(LEntriesType::rank == 1,
                "spiluk_numeric_streams: A_entries, L_entries and U_entries must all "
                "have rank 1.");
  static_assert(LValuesType::rank == 1,
                "spiluk_numeric_streams: A_values, L_values and U_values must "
                "all have rank 1.");

  static_assert(std::is_same<typename LEntriesType::value_type, typename LEntriesType::non_const_value_type>::value,
                "spiluk_numeric_streams: The output L_entries must be nonconst.");
  static_assert(std::is_same<typename LValuesType::value_type, typename LValuesType::non_const_value_type>::value,
                "spiluk_numeric_streams: The output L_values must be nonconst.");
  static_assert(std::is_same<typename UEntriesType::value_type, typename UEntriesType::non_const_value_type>::value,
                "spiluk_numeric_streams: The output U_entries must be nonconst.");
  static_assert(std::is_same<typename UValuesType::value_type, typename UValuesType::non_const_value_type>::value,
                "spiluk_numeric_streams: The output U_values must be nonconst.");

  static_assert(std::is_same<typename LRowMapType::device_type, typename ARowMapType::device_type>::value,
                "spiluk_numeric_streams: Views LRowMapType and ARowMapType have "
                "different device_types.");
  static_assert(std::is_same<typename LEntriesType::device_type, typename AEntriesType::device_type>::value,
                "spiluk_numeric_streams: Views LEntriesType and AEntriesType have "
                "different device_types.");
  static_assert(std::is_same<typename LValuesType::device_type, typename AValuesType::device_type>::value,
                "spiluk_numeric_streams: Views LValuesType and AValuesType have "
                "different device_types.");

  static_assert(std::is_same<typename LRowMapType::device_type, typename URowMapType::device_type>::value,
                "spiluk_numeric_streams: Views LRowMapType and URowMapType have "
                "different device_types.");
  static_assert(std::is_same<typename LEntriesType::device_type, typename UEntriesType::device_type>::value,
                "spiluk_numeric_streams: Views LEntriesType and UEntriesType have "
                "different device_types.");
  static_assert(std::is_same<typename LValuesType::device_type, typename UValuesType::device_type>::value,
                "spiluk_numeric_streams: Views LValuesType and UValuesType have "
                "different device_types.");

  static_assert(std::is_same<ExecutionSpace, typename KernelHandle::SPILUKHandleType::execution_space>::value,
                "spiluk_numeric_streams: KernelHandle's execution space is different "
                "from "
                "ExecutionSpace.");

  static_assert(std::is_same<typename LRowMapType::device_type::execution_space,
                             typename KernelHandle::SPILUKHandleType::execution_space>::value,
                "spiluk_numeric_streams: KernelHandle and Views have different execution "
                "spaces.");
  static_assert(std::is_same<typename LEntriesType::device_type::execution_space,
                             typename KernelHandle::SPILUKHandleType::execution_space>::value,
                "spiluk_numeric_streams: KernelHandle and Views have different execution "
                "spaces.");
  static_assert(std::is_same<typename LValuesType::device_type::execution_space,
                             typename KernelHandle::SPILUKHandleType::execution_space>::value,
                "spiluk_numeric_streams: KernelHandle and Views have different execution "
                "spaces.");

  static_assert(std::is_same<typename LRowMapType::device_type, typename LEntriesType::device_type>::value,
                "spiluk_numeric_streams: rowmap and entries have different "
                "device types.");
  static_assert(std::is_same<typename LRowMapType::device_type, typename LValuesType::device_type>::value,
                "spiluk_numeric_streams: rowmap and values have different device types.");

  // Check validity of fill level
  if (fill_lev < 0) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::spiluk_numeric_streams: fill_lev: " << fill_lev << ". Valid value is >= 0.";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  // Check sizes of vectors
  if (execspace_v.size() != handle_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::spiluk_numeric_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. handle_v.size() " << handle_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (execspace_v.size() != A_rowmap_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::spiluk_numeric_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. A_rowmap_v.size() " << A_rowmap_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (execspace_v.size() != A_entries_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::spiluk_numeric_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. A_entries_v.size() " << A_entries_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (execspace_v.size() != A_values_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::spiluk_numeric_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. A_values_v.size() " << A_values_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (execspace_v.size() != L_rowmap_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::spiluk_numeric_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. L_rowmap_v.size() " << L_rowmap_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (execspace_v.size() != L_entries_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::spiluk_numeric_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. L_entries_v.size() " << L_entries_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (execspace_v.size() != L_values_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::spiluk_numeric_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. L_values_v.size() " << L_values_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (execspace_v.size() != U_rowmap_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::spiluk_numeric_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. U_rowmap_v.size() " << U_rowmap_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (execspace_v.size() != U_entries_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::spiluk_numeric_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. U_entries_v.size() " << U_entries_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (execspace_v.size() != U_values_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::spiluk_numeric_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. U_values_v.size() " << U_values_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  // Check if symbolic has been called
  for (int i = 0; i < static_cast<int>(execspace_v.size()); i++) {
    if (handle_v[i]->get_spiluk_handle()->is_symbolic_complete() == false) {
      std::ostringstream os;
      os << "KokkosSparse::Experimental::spiluk_numeric_streams: "
            "spiluk_symbolic must be "
            "called before spiluk_numeric_streams -- stream "
         << i;
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }

  using c_size_t    = typename KernelHandle::const_size_type;
  using c_lno_t     = typename KernelHandle::const_nnz_lno_t;
  using c_scalar_t  = typename KernelHandle::const_nnz_scalar_t;
  using c_exec_t    = typename KernelHandle::HandleExecSpace;
  using c_temp_t    = typename KernelHandle::HandleTempMemorySpace;
  using c_persist_t = typename KernelHandle::HandlePersistentMemorySpace;

  using const_handle_type = typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t,
                                                                                      c_exec_t, c_temp_t, c_persist_t>;

  using ARowMap_Internal =
      Kokkos::View<typename ARowMapType::const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<ARowMapType>::array_layout,
                   typename ARowMapType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using AEntries_Internal =
      Kokkos::View<typename AEntriesType::const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<AEntriesType>::array_layout,
                   typename AEntriesType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using AValues_Internal =
      Kokkos::View<typename AValuesType::const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<AValuesType>::array_layout,
                   typename AValuesType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using LRowMap_Internal =
      Kokkos::View<typename LRowMapType::const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<LRowMapType>::array_layout,
                   typename LRowMapType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using LEntries_Internal =
      Kokkos::View<typename LEntriesType::non_const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<LEntriesType>::array_layout,
                   typename LEntriesType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using LValues_Internal =
      Kokkos::View<typename LValuesType::non_const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<LValuesType>::array_layout,
                   typename LValuesType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using URowMap_Internal =
      Kokkos::View<typename URowMapType::const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<URowMapType>::array_layout,
                   typename URowMapType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using UEntries_Internal =
      Kokkos::View<typename UEntriesType::non_const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<UEntriesType>::array_layout,
                   typename UEntriesType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using UValues_Internal =
      Kokkos::View<typename UValuesType::non_const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<UValuesType>::array_layout,
                   typename UValuesType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  std::vector<const_handle_type> handle_i_v(execspace_v.size());
  std::vector<ARowMap_Internal> A_rowmap_i_v(execspace_v.size());
  std::vector<AEntries_Internal> A_entries_i_v(execspace_v.size());
  std::vector<AValues_Internal> A_values_i_v(execspace_v.size());
  std::vector<LRowMap_Internal> L_rowmap_i_v(execspace_v.size());
  std::vector<LEntries_Internal> L_entries_i_v(execspace_v.size());
  std::vector<LValues_Internal> L_values_i_v(execspace_v.size());
  std::vector<URowMap_Internal> U_rowmap_i_v(execspace_v.size());
  std::vector<UEntries_Internal> U_entries_i_v(execspace_v.size());
  std::vector<UValues_Internal> U_values_i_v(execspace_v.size());

  for (int i = 0; i < static_cast<int>(execspace_v.size()); i++) {
    handle_i_v[i]    = const_handle_type(*(handle_v[i]));
    A_rowmap_i_v[i]  = A_rowmap_v[i];
    A_entries_i_v[i] = A_entries_v[i];
    A_values_i_v[i]  = A_values_v[i];
    L_rowmap_i_v[i]  = L_rowmap_v[i];
    L_entries_i_v[i] = L_entries_v[i];
    L_values_i_v[i]  = L_values_v[i];
    U_rowmap_i_v[i]  = U_rowmap_v[i];
    U_entries_i_v[i] = U_entries_v[i];
    U_values_i_v[i]  = U_values_v[i];
  }

  KokkosSparse::Impl::SPILUK_NUMERIC<ExecutionSpace, const_handle_type, ARowMap_Internal, AEntries_Internal,
                                     AValues_Internal, LRowMap_Internal, LEntries_Internal, LValues_Internal,
                                     URowMap_Internal, UEntries_Internal,
                                     UValues_Internal>::spiluk_numeric_streams(execspace_v, handle_i_v, A_rowmap_i_v,
                                                                               A_entries_i_v, A_values_i_v,
                                                                               L_rowmap_i_v, L_entries_i_v,
                                                                               L_values_i_v, U_rowmap_i_v,
                                                                               U_entries_i_v, U_values_i_v);

}  // spiluk_numeric_streams

}  // namespace Experimental
}  // namespace KokkosSparse

#undef KOKKOSKERNELS_SPILUK_SAME_TYPE

#endif  // KOKKOSSPARSE_SPILUK_HPP_
