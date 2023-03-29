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

/// \file KokkosSparse_par_ilut.hpp
/// \brief Parallel threshold incomplete LU factorization ILU(t)
///
/// This file provides KokkosSparse::par_ilut.  This function performs a
/// local (no MPI) sparse ILU(t) on matrices stored in
/// compressed row sparse ("Crs") format. It is expected that symbolic
/// is called before numeric. The numeric function offers a deterministic
/// flag that will force the function to have deterministic results. This
/// is useful for testing but incurs a big performance penalty.
///
/// This algorithm is described in the paper:
/// PARILUT - A New Parallel Threshold ILU Factorization - Anzt, Chow, Dongarra

#ifndef KOKKOSSPARSE_PAR_ILUT_HPP_
#define KOKKOSSPARSE_PAR_ILUT_HPP_

#include <type_traits>

#include "KokkosKernels_helpers.hpp"
#include "KokkosKernels_Error.hpp"
#include "KokkosSparse_par_ilut_symbolic_spec.hpp"
#include "KokkosSparse_par_ilut_numeric_spec.hpp"

namespace KokkosSparse {
namespace Experimental {

#define KOKKOSKERNELS_PAR_ILUT_SAME_TYPE(A, B)      \
  std::is_same<typename std::remove_const<A>::type, \
               typename std::remove_const<B>::type>::value

template <typename KernelHandle, typename ARowMapType, typename AEntriesType,
          typename LRowMapType, typename URowMapType>
void par_ilut_symbolic(KernelHandle* handle, ARowMapType& A_rowmap,
                       AEntriesType& A_entries, LRowMapType& L_rowmap,
                       URowMapType& U_rowmap) {
  using size_type    = typename KernelHandle::size_type;
  using ordinal_type = typename KernelHandle::nnz_lno_t;

  static_assert(KOKKOSKERNELS_PAR_ILUT_SAME_TYPE(
                    typename ARowMapType::non_const_value_type, size_type),
                "par_ilut_symbolic: A size_type must match KernelHandle "
                "size_type (const doesn't matter)");
  static_assert(KOKKOSKERNELS_PAR_ILUT_SAME_TYPE(
                    typename AEntriesType::non_const_value_type, ordinal_type),
                "par_ilut_symbolic: A entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");

  static_assert(KOKKOSKERNELS_PAR_ILUT_SAME_TYPE(
                    typename LRowMapType::non_const_value_type, size_type),
                "par_ilut_symbolic: L size_type must match KernelHandle "
                "size_type (const doesn't matter)");

  static_assert(KOKKOSKERNELS_PAR_ILUT_SAME_TYPE(
                    typename URowMapType::non_const_value_type, size_type),
                "par_ilut_symbolic: U size_type must match KernelHandle "
                "size_type (const doesn't matter)");

  static_assert(Kokkos::is_view<ARowMapType>::value,
                "par_ilut_symbolic: A_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<AEntriesType>::value,
                "par_ilut_symbolic: A_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LRowMapType>::value,
                "par_ilut_symbolic: L_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<URowMapType>::value,
                "par_ilut_symbolic: U_rowmap is not a Kokkos::View.");

  static_assert(
      (int)LRowMapType::rank == (int)ARowMapType::rank,
      "par_ilut_symbolic: The ranks of L_rowmap and A_rowmap do not match.");

  static_assert(
      (int)LRowMapType::rank == (int)URowMapType::rank,
      "par_ilut_symbolic: The ranks of L_rowmap and U_rowmap do not match.");

  static_assert(LRowMapType::rank == 1,
                "par_ilut_symbolic: A_rowmap, L_rowmap and U_rowmap must all "
                "have rank 1.");

  static_assert(std::is_same<typename LRowMapType::value_type,
                             typename LRowMapType::non_const_value_type>::value,
                "par_ilut_symbolic: The output L_rowmap must be nonconst.");
  static_assert(std::is_same<typename URowMapType::value_type,
                             typename URowMapType::non_const_value_type>::value,
                "par_ilut_symbolic: The output U_rowmap must be nonconst.");
  static_assert(std::is_same<typename LRowMapType::device_type,
                             typename ARowMapType::device_type>::value,
                "par_ilut_symbolic: Views LRowMapType and ARowMapType have "
                "different device_types.");
  static_assert(std::is_same<typename LRowMapType::device_type,
                             typename URowMapType::device_type>::value,
                "par_ilut_symbolic: Views LRowMapType and URowMapType have "
                "different device_types.");

  static_assert(
      std::is_same<
          typename LRowMapType::device_type::execution_space,
          typename KernelHandle::PAR_ILUTHandleType::execution_space>::value,
      "par_ilut_symbolic: KernelHandle and Views have different execution "
      "spaces.");

  using c_size_t   = typename KernelHandle::const_size_type;
  using c_lno_t    = typename KernelHandle::const_nnz_lno_t;
  using c_scalar_t = typename KernelHandle::const_nnz_scalar_t;

  using c_exec_t    = typename KernelHandle::HandleExecSpace;
  using c_temp_t    = typename KernelHandle::HandleTempMemorySpace;
  using c_persist_t = typename KernelHandle::HandlePersistentMemorySpace;

  using const_handle_type =
      typename KokkosKernels::Experimental::KokkosKernelsHandle<
          c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t, c_persist_t>;

  const_handle_type tmp_handle(*handle);

  using ARowMap_Internal = Kokkos::View<
      typename ARowMapType::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<ARowMapType>::array_layout,
      typename ARowMapType::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using AEntries_Internal = Kokkos::View<
      typename AEntriesType::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<
          AEntriesType>::array_layout,
      typename AEntriesType::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using LRowMap_Internal = Kokkos::View<
      typename LRowMapType::non_const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<LRowMapType>::array_layout,
      typename LRowMapType::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using URowMap_Internal = Kokkos::View<
      typename URowMapType::non_const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<URowMapType>::array_layout,
      typename URowMapType::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  ARowMap_Internal A_rowmap_i   = A_rowmap;
  AEntries_Internal A_entries_i = A_entries;
  LRowMap_Internal L_rowmap_i   = L_rowmap;
  URowMap_Internal U_rowmap_i   = U_rowmap;

  KokkosSparse::Impl::PAR_ILUT_SYMBOLIC<
      const_handle_type, ARowMap_Internal, AEntries_Internal, LRowMap_Internal,
      URowMap_Internal>::par_ilut_symbolic(&tmp_handle, A_rowmap_i, A_entries_i,
                                           L_rowmap_i, U_rowmap_i);

}  // par_ilut_symbolic

template <typename KernelHandle, typename ARowMapType, typename AEntriesType,
          typename AValuesType, typename LRowMapType, typename LEntriesType,
          typename LValuesType, typename URowMapType, typename UEntriesType,
          typename UValuesType>
void par_ilut_numeric(KernelHandle* handle, ARowMapType& A_rowmap,
                      AEntriesType& A_entries, AValuesType& A_values,
                      LRowMapType& L_rowmap, LEntriesType& L_entries,
                      LValuesType& L_values, URowMapType& U_rowmap,
                      UEntriesType& U_entries, UValuesType& U_values,
                      bool deterministic) {
  using size_type    = typename KernelHandle::size_type;
  using ordinal_type = typename KernelHandle::nnz_lno_t;
  using scalar_type  = typename KernelHandle::nnz_scalar_t;

  static_assert(
      KOKKOSKERNELS_PAR_ILUT_SAME_TYPE(
          typename ARowMapType::non_const_value_type, size_type),
      "par_ilut_numeric: A size_type must match KernelHandle size_type "
      "(const doesn't matter)");
  static_assert(KOKKOSKERNELS_PAR_ILUT_SAME_TYPE(
                    typename AEntriesType::non_const_value_type, ordinal_type),
                "par_ilut_numeric: A entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");
  static_assert(KOKKOSKERNELS_PAR_ILUT_SAME_TYPE(
                    typename AValuesType::value_type, scalar_type),
                "par_ilut_numeric: A scalar type must match KernelHandle entry "
                "type (aka nnz_scalar_t, and const doesn't matter)");

  static_assert(
      KOKKOSKERNELS_PAR_ILUT_SAME_TYPE(
          typename LRowMapType::non_const_value_type, size_type),
      "par_ilut_numeric: L size_type must match KernelHandle size_type "
      "(const doesn't matter)");
  static_assert(KOKKOSKERNELS_PAR_ILUT_SAME_TYPE(
                    typename LEntriesType::non_const_value_type, ordinal_type),
                "par_ilut_numeric: L entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");
  static_assert(KOKKOSKERNELS_PAR_ILUT_SAME_TYPE(
                    typename LValuesType::value_type, scalar_type),
                "par_ilut_numeric: L scalar type must match KernelHandle entry "
                "type (aka nnz_scalar_t, and const doesn't matter)");

  static_assert(
      KOKKOSKERNELS_PAR_ILUT_SAME_TYPE(
          typename URowMapType::non_const_value_type, size_type),
      "par_ilut_numeric: U size_type must match KernelHandle size_type "
      "(const doesn't matter)");
  static_assert(KOKKOSKERNELS_PAR_ILUT_SAME_TYPE(
                    typename UEntriesType::non_const_value_type, ordinal_type),
                "par_ilut_numeric: U entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");
  static_assert(KOKKOSKERNELS_PAR_ILUT_SAME_TYPE(
                    typename UValuesType::value_type, scalar_type),
                "par_ilut_numeric: U scalar type must match KernelHandle entry "
                "type (aka nnz_scalar_t, and const doesn't matter)");

  static_assert(Kokkos::is_view<ARowMapType>::value,
                "par_ilut_numeric: A_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<AEntriesType>::value,
                "par_ilut_numeric: A_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<AValuesType>::value,
                "par_ilut_numeric: A_values is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LRowMapType>::value,
                "par_ilut_numeric: L_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LEntriesType>::value,
                "par_ilut_numeric: L_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LValuesType>::value,
                "par_ilut_numeric: L_values is not a Kokkos::View.");
  static_assert(Kokkos::is_view<URowMapType>::value,
                "par_ilut_numeric: U_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<UEntriesType>::value,
                "par_ilut_numeric: U_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<UValuesType>::value,
                "par_ilut_numeric: U_values is not a Kokkos::View.");

  static_assert(
      (int)LRowMapType::rank == (int)ARowMapType::rank,
      "par_ilut_numeric: The ranks of L_rowmap and A_rowmap do not match.");
  static_assert(
      (int)LEntriesType::rank == (int)AEntriesType::rank,
      "par_ilut_numeric: The ranks of L_entries and A_entries do not match.");
  static_assert(
      (int)LValuesType::rank == (int)AValuesType::rank,
      "par_ilut_numeric: The ranks of L_values and A_values do not match.");

  static_assert(
      (int)LRowMapType::rank == (int)URowMapType::rank,
      "par_ilut_numeric: The ranks of L_rowmap and U_rowmap do not match.");
  static_assert(
      (int)LEntriesType::rank == (int)UEntriesType::rank,
      "par_ilut_numeric: The ranks of L_entries and U_entries do not match.");
  static_assert(
      (int)LValuesType::rank == (int)UValuesType::rank,
      "par_ilut_numeric: The ranks of L_values and U_values do not match.");

  static_assert(LRowMapType::rank == 1,
                "par_ilut_numeric: A_rowmap, L_rowmap and U_rowmap must all "
                "have rank 1.");
  static_assert(LEntriesType::rank == 1,
                "par_ilut_numeric: A_entries, L_entries and U_entries must all "
                "have rank 1.");
  static_assert(LValuesType::rank == 1,
                "par_ilut_numeric: A_values, L_values and U_values must all "
                "have rank 1.");

  static_assert(
      std::is_same<typename LEntriesType::value_type,
                   typename LEntriesType::non_const_value_type>::value,
      "par_ilut_numeric: The output L_entries must be nonconst.");
  static_assert(std::is_same<typename LValuesType::value_type,
                             typename LValuesType::non_const_value_type>::value,
                "par_ilut_numeric: The output L_values must be nonconst.");
  static_assert(
      std::is_same<typename UEntriesType::value_type,
                   typename UEntriesType::non_const_value_type>::value,
      "par_ilut_numeric: The output U_entries must be nonconst.");
  static_assert(std::is_same<typename UValuesType::value_type,
                             typename UValuesType::non_const_value_type>::value,
                "par_ilut_numeric: The output U_values must be nonconst.");

  static_assert(std::is_same<typename LRowMapType::device_type,
                             typename ARowMapType::device_type>::value,
                "par_ilut_numeric: Views LRowMapType and ARowMapType have "
                "different device_types.");
  static_assert(std::is_same<typename LEntriesType::device_type,
                             typename AEntriesType::device_type>::value,
                "par_ilut_numeric: Views LEntriesType and AEntriesType have "
                "different device_types.");
  static_assert(std::is_same<typename LValuesType::device_type,
                             typename AValuesType::device_type>::value,
                "par_ilut_numeric: Views LValuesType and AValuesType have "
                "different device_types.");

  static_assert(std::is_same<typename LRowMapType::device_type,
                             typename URowMapType::device_type>::value,
                "par_ilut_numeric: Views LRowMapType and URowMapType have "
                "different device_types.");
  static_assert(std::is_same<typename LEntriesType::device_type,
                             typename UEntriesType::device_type>::value,
                "par_ilut_numeric: Views LEntriesType and UEntriesType have "
                "different device_types.");
  static_assert(std::is_same<typename LValuesType::device_type,
                             typename UValuesType::device_type>::value,
                "par_ilut_numeric: Views LValuesType and UValuesType have "
                "different device_types.");

  static_assert(
      std::is_same<
          typename LRowMapType::device_type::execution_space,
          typename KernelHandle::PAR_ILUTHandleType::execution_space>::value,
      "par_ilut_numeric: KernelHandle and Views have different execution "
      "spaces.");
  static_assert(
      std::is_same<
          typename LEntriesType::device_type::execution_space,
          typename KernelHandle::PAR_ILUTHandleType::execution_space>::value,
      "par_ilut_numeric: KernelHandle and Views have different execution "
      "spaces.");
  static_assert(
      std::is_same<
          typename LValuesType::device_type::execution_space,
          typename KernelHandle::PAR_ILUTHandleType::execution_space>::value,
      "par_ilut_numeric: KernelHandle and Views have different execution "
      "spaces.");

  static_assert(
      std::is_same<typename LRowMapType::device_type,
                   typename LEntriesType::device_type>::value,
      "par_ilut_numeric: rowmap and entries have different device types.");
  static_assert(
      std::is_same<typename LRowMapType::device_type,
                   typename LValuesType::device_type>::value,
      "par_ilut_numeric: rowmap and values have different device types.");

  // Check if symbolic has been called
  if (handle->get_par_ilut_handle()->is_symbolic_complete() == false) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::par_ilut_numeric: par_ilut_symbolic "
          "must be "
          "called before par_ilut_numeric.";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using c_size_t   = typename KernelHandle::const_size_type;
  using c_lno_t    = typename KernelHandle::const_nnz_lno_t;
  using c_scalar_t = typename KernelHandle::const_nnz_scalar_t;

  using c_exec_t    = typename KernelHandle::HandleExecSpace;
  using c_temp_t    = typename KernelHandle::HandleTempMemorySpace;
  using c_persist_t = typename KernelHandle::HandlePersistentMemorySpace;

  using const_handle_type =
      typename KokkosKernels::Experimental::KokkosKernelsHandle<
          c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t, c_persist_t>;

  const_handle_type tmp_handle(*handle);

  using ARowMap_Internal = Kokkos::View<
      typename ARowMapType::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<ARowMapType>::array_layout,
      typename ARowMapType::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using AEntries_Internal = Kokkos::View<
      typename AEntriesType::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<
          AEntriesType>::array_layout,
      typename AEntriesType::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using AValues_Internal = Kokkos::View<
      typename AValuesType::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<AValuesType>::array_layout,
      typename AValuesType::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using LRowMap_Internal = Kokkos::View<
      typename LRowMapType::non_const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<LRowMapType>::array_layout,
      typename LRowMapType::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using LEntries_Internal =
      Kokkos::View<typename LEntriesType::non_const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<
                       LEntriesType>::array_layout,
                   typename LEntriesType::device_type,
                   Kokkos::MemoryTraits<Kokkos::RandomAccess> >;

  using LValues_Internal = Kokkos::View<
      typename LValuesType::non_const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<LValuesType>::array_layout,
      typename LValuesType::device_type,
      Kokkos::MemoryTraits<Kokkos::RandomAccess> >;

  using URowMap_Internal = Kokkos::View<
      typename URowMapType::non_const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<URowMapType>::array_layout,
      typename URowMapType::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using UEntries_Internal =
      Kokkos::View<typename UEntriesType::non_const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<
                       UEntriesType>::array_layout,
                   typename UEntriesType::device_type,
                   Kokkos::MemoryTraits<Kokkos::RandomAccess> >;

  using UValues_Internal = Kokkos::View<
      typename UValuesType::non_const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<UValuesType>::array_layout,
      typename UValuesType::device_type,
      Kokkos::MemoryTraits<Kokkos::RandomAccess> >;

  ARowMap_Internal A_rowmap_i   = A_rowmap;
  AEntries_Internal A_entries_i = A_entries;
  AValues_Internal A_values_i   = A_values;
  LRowMap_Internal L_rowmap_i   = L_rowmap;
  LEntries_Internal L_entries_i = L_entries;
  LValues_Internal L_values_i   = L_values;
  URowMap_Internal U_rowmap_i   = U_rowmap;
  UEntries_Internal U_entries_i = U_entries;
  UValues_Internal U_values_i   = U_values;

  KokkosSparse::Impl::PAR_ILUT_NUMERIC<
      const_handle_type, ARowMap_Internal, AEntries_Internal, AValues_Internal,
      LRowMap_Internal, LEntries_Internal, LValues_Internal, URowMap_Internal,
      UEntries_Internal,
      UValues_Internal>::par_ilut_numeric(&tmp_handle, A_rowmap_i, A_entries_i,
                                          A_values_i, L_rowmap_i, L_entries_i,
                                          L_values_i, U_rowmap_i, U_entries_i,
                                          U_values_i, deterministic);

  // These may have been resized
  L_entries = L_entries_i;
  L_values  = L_values_i;
  U_entries = U_entries_i;
  U_values  = U_values_i;

}  // par_ilut_numeric

}  // namespace Experimental
}  // namespace KokkosSparse

#undef KOKKOSKERNELS_PAR_ILUT_SAME_TYPE

#endif  // KOKKOSSPARSE_PAR_ILUT_HPP_
