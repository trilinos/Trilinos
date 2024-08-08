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
/// is called before numeric. The handle offers an async_update
/// flag that controls whether asynchronous updates are allowed while computing
/// L U factors. This is useful for testing as it allows for repeatable
/// (deterministic) results but may cause the algorithm to take longer (more
/// iterations) to converge. The par_ilut algorithm will repeat (iterate) until
/// max_iters is hit or the improvement in the residual from iter to iter drops
/// below a certain threshold.
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

// Two types are the same (ignoring const)
template <typename T, typename U>
constexpr bool parilut_same_type = std::is_same_v<typename std::remove_const_t<T>, typename std::remove_const_t<U>>;

/// @brief Performs the symbolic phase of par_ilut.
///        This is a non-blocking function.
///
/// The sparsity pattern of A will be analyzed and L_rowmap and U_rowmap will be
/// populated with the L (lower triangular) and U (upper triagular) non-zero
/// counts respectively. Having a separate symbolic phase allows for reuse when
/// dealing with multiple matrices with the same sparsity pattern. This routine
/// will set some values on handle for symbolic info (row count, nnz counts).
///
/// @tparam KernelHandle Template for the KernelHandle type
/// @tparam ARowMapType Template for A_rowmap type
/// @tparam AEntriesType Template for A_entries type
/// @tparam LRowMapType Template for L_rowmap type
/// @tparam URowMapType Template for U_rowmap type
/// @param handle The kernel handle. It is expected that create_par_ilut_handle
/// has been called on it
/// @param A_rowmap The row map (row nnz offsets) for the A CSR (Input)
/// @param A_entries The entries (column ids) for the A CSR (Input)
/// @param L_rowmap The row map for the L CSR, should already be sized correctly
/// (numRows+1) (Output)
/// @param U_rowmap The row map for the U CSR, should already be sized correctly
/// (numRows+1) (Output)
template <typename KernelHandle, typename ARowMapType, typename AEntriesType, typename LRowMapType,
          typename URowMapType>
void par_ilut_symbolic(KernelHandle* handle, ARowMapType& A_rowmap, AEntriesType& A_entries, LRowMapType& L_rowmap,
                       URowMapType& U_rowmap) {
  using size_type    = typename KernelHandle::size_type;
  using ordinal_type = typename KernelHandle::nnz_lno_t;

  static_assert(parilut_same_type<typename ARowMapType::non_const_value_type, size_type>,
                "par_ilut_symbolic: A size_type must match KernelHandle "
                "size_type (const doesn't matter)");
  static_assert(parilut_same_type<typename AEntriesType::non_const_value_type, ordinal_type>,
                "par_ilut_symbolic: A entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");

  static_assert(parilut_same_type<typename LRowMapType::non_const_value_type, size_type>,
                "par_ilut_symbolic: L size_type must match KernelHandle "
                "size_type (const doesn't matter)");

  static_assert(parilut_same_type<typename URowMapType::non_const_value_type, size_type>,
                "par_ilut_symbolic: U size_type must match KernelHandle "
                "size_type (const doesn't matter)");

  static_assert(Kokkos::is_view<ARowMapType>::value, "par_ilut_symbolic: A_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<AEntriesType>::value, "par_ilut_symbolic: A_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LRowMapType>::value, "par_ilut_symbolic: L_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<URowMapType>::value, "par_ilut_symbolic: U_rowmap is not a Kokkos::View.");

  static_assert((int)LRowMapType::rank == (int)ARowMapType::rank,
                "par_ilut_symbolic: The ranks of L_rowmap and A_rowmap do not match.");

  static_assert((int)LRowMapType::rank == (int)URowMapType::rank,
                "par_ilut_symbolic: The ranks of L_rowmap and U_rowmap do not match.");

  static_assert(LRowMapType::rank == 1,
                "par_ilut_symbolic: A_rowmap, L_rowmap and U_rowmap must all "
                "have rank 1.");

  static_assert(std::is_same_v<typename LRowMapType::value_type, typename LRowMapType::non_const_value_type>,
                "par_ilut_symbolic: The output L_rowmap must be nonconst.");
  static_assert(std::is_same_v<typename URowMapType::value_type, typename URowMapType::non_const_value_type>,
                "par_ilut_symbolic: The output U_rowmap must be nonconst.");
  static_assert(std::is_same_v<typename LRowMapType::device_type, typename ARowMapType::device_type>,
                "par_ilut_symbolic: Views LRowMapType and ARowMapType have "
                "different device_types.");
  static_assert(std::is_same_v<typename LRowMapType::device_type, typename URowMapType::device_type>,
                "par_ilut_symbolic: Views LRowMapType and URowMapType have "
                "different device_types.");

  static_assert(std::is_same_v<typename LRowMapType::device_type::execution_space,
                               typename KernelHandle::PAR_ILUTHandleType::execution_space>,
                "par_ilut_symbolic: KernelHandle and Views have different execution "
                "spaces.");

  if (A_rowmap.extent(0) != 0) {
    KK_REQUIRE_MSG(A_rowmap.extent(0) == L_rowmap.extent(0), "L row map size does not match A row map");
    KK_REQUIRE_MSG(A_rowmap.extent(0) == U_rowmap.extent(0), "U row map size does not match A row map");
  }

  using c_size_t   = typename KernelHandle::const_size_type;
  using c_lno_t    = typename KernelHandle::const_nnz_lno_t;
  using c_scalar_t = typename KernelHandle::const_nnz_scalar_t;

  using c_exec_t    = typename KernelHandle::HandleExecSpace;
  using c_temp_t    = typename KernelHandle::HandleTempMemorySpace;
  using c_persist_t = typename KernelHandle::HandlePersistentMemorySpace;

  using const_handle_type = typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t,
                                                                                      c_exec_t, c_temp_t, c_persist_t>;

  const_handle_type tmp_handle(*handle);

  using ARowMap_Internal =
      Kokkos::View<typename ARowMapType::const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<ARowMapType>::array_layout,
                   typename ARowMapType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

  using AEntries_Internal =
      Kokkos::View<typename AEntriesType::const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<AEntriesType>::array_layout,
                   typename AEntriesType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

  using LRowMap_Internal =
      Kokkos::View<typename LRowMapType::non_const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<LRowMapType>::array_layout,
                   typename LRowMapType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

  using URowMap_Internal =
      Kokkos::View<typename URowMapType::non_const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<URowMapType>::array_layout,
                   typename URowMapType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

  ARowMap_Internal A_rowmap_i   = A_rowmap;
  AEntries_Internal A_entries_i = A_entries;
  LRowMap_Internal L_rowmap_i   = L_rowmap;
  URowMap_Internal U_rowmap_i   = U_rowmap;

  KokkosSparse::Impl::PAR_ILUT_SYMBOLIC<const_handle_type, ARowMap_Internal, AEntries_Internal, LRowMap_Internal,
                                        URowMap_Internal>::par_ilut_symbolic(&tmp_handle, A_rowmap_i, A_entries_i,
                                                                             L_rowmap_i, U_rowmap_i);

}  // par_ilut_symbolic

/// @brief Performs the numeric phase (for specific CSRs, not reusable) of the
/// par_ilut
///        algorithm (described in the header). This is a non-blocking
///        functions. It is expected that par_ilut_symbolic has already been
///        called for the
//         provided KernelHandle.
///
/// @tparam KernelHandle Template for the handle type
/// @tparam ARowMapType Template for the A_rowmap type
/// @tparam AEntriesType Template for the A_entries type
/// @tparam AValuesType Template for the A_values type
/// @tparam LRowMapType Template for the L_rowmap type
/// @tparam LEntriesType Template for the L_entries type
/// @tparam LValuesType Template for the L_values type
/// @tparam URowMapType Template for the U_rowmap type
/// @tparam UEntriesType Template for the U_entries type
/// @tparam UValuesType Template for the U_values type
/// @param handle The kernel handle. It is expected that create_par_ilut_handle
/// has been called on it
/// @param A_rowmap The row map (row nnz offsets) for the A CSR (Input)
/// @param A_entries The entries (column ids) for the A CSR (Input)
/// @param A_values The values (non-zero matrix values) for the A CSR (Input)
/// @param L_rowmap The row map (row nnz offsets) for the L CSR (Input/Output)
/// @param L_entries The entries (column ids) for the L CSR (Output)
/// @param L_values The values (non-zero matrix values) for the L CSR (Output)
/// @param U_rowmap The row map (row nnz offsets) for the U CSR (Input/Output)
/// @param U_entries The entries (column ids) for the U CSR (Output)
/// @param U_values The values (non-zero matrix values) for the U CSR (Output)
template <typename KernelHandle, typename ARowMapType, typename AEntriesType, typename AValuesType,
          typename LRowMapType, typename LEntriesType, typename LValuesType, typename URowMapType,
          typename UEntriesType, typename UValuesType>
void par_ilut_numeric(KernelHandle* handle, ARowMapType& A_rowmap, AEntriesType& A_entries, AValuesType& A_values,
                      LRowMapType& L_rowmap, LEntriesType& L_entries, LValuesType& L_values, URowMapType& U_rowmap,
                      UEntriesType& U_entries, UValuesType& U_values) {
  using size_type    = typename KernelHandle::size_type;
  using ordinal_type = typename KernelHandle::nnz_lno_t;
  using scalar_type  = typename KernelHandle::nnz_scalar_t;

  static_assert(parilut_same_type<typename ARowMapType::non_const_value_type, size_type>,
                "par_ilut_numeric: A size_type must match KernelHandle size_type "
                "(const doesn't matter)");
  static_assert(parilut_same_type<typename AEntriesType::non_const_value_type, ordinal_type>,
                "par_ilut_numeric: A entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");
  static_assert(parilut_same_type<typename AValuesType::value_type, scalar_type>,
                "par_ilut_numeric: A scalar type must match KernelHandle entry "
                "type (aka nnz_scalar_t, and const doesn't matter)");

  static_assert(parilut_same_type<typename LRowMapType::non_const_value_type, size_type>,
                "par_ilut_numeric: L size_type must match KernelHandle size_type "
                "(const doesn't matter)");
  static_assert(parilut_same_type<typename LEntriesType::non_const_value_type, ordinal_type>,
                "par_ilut_numeric: L entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");
  static_assert(parilut_same_type<typename LValuesType::value_type, scalar_type>,
                "par_ilut_numeric: L scalar type must match KernelHandle entry "
                "type (aka nnz_scalar_t, and const doesn't matter)");

  static_assert(parilut_same_type<typename URowMapType::non_const_value_type, size_type>,
                "par_ilut_numeric: U size_type must match KernelHandle size_type "
                "(const doesn't matter)");
  static_assert(parilut_same_type<typename UEntriesType::non_const_value_type, ordinal_type>,
                "par_ilut_numeric: U entry type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");
  static_assert(parilut_same_type<typename UValuesType::value_type, scalar_type>,
                "par_ilut_numeric: U scalar type must match KernelHandle entry "
                "type (aka nnz_scalar_t, and const doesn't matter)");

  static_assert(Kokkos::is_view<ARowMapType>::value, "par_ilut_numeric: A_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<AEntriesType>::value, "par_ilut_numeric: A_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<AValuesType>::value, "par_ilut_numeric: A_values is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LRowMapType>::value, "par_ilut_numeric: L_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LEntriesType>::value, "par_ilut_numeric: L_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<LValuesType>::value, "par_ilut_numeric: L_values is not a Kokkos::View.");
  static_assert(Kokkos::is_view<URowMapType>::value, "par_ilut_numeric: U_rowmap is not a Kokkos::View.");
  static_assert(Kokkos::is_view<UEntriesType>::value, "par_ilut_numeric: U_entries is not a Kokkos::View.");
  static_assert(Kokkos::is_view<UValuesType>::value, "par_ilut_numeric: U_values is not a Kokkos::View.");

  static_assert((int)LRowMapType::rank == (int)ARowMapType::rank,
                "par_ilut_numeric: The ranks of L_rowmap and A_rowmap do not match.");
  static_assert((int)LEntriesType::rank == (int)AEntriesType::rank,
                "par_ilut_numeric: The ranks of L_entries and A_entries do not match.");
  static_assert((int)LValuesType::rank == (int)AValuesType::rank,
                "par_ilut_numeric: The ranks of L_values and A_values do not match.");

  static_assert((int)LRowMapType::rank == (int)URowMapType::rank,
                "par_ilut_numeric: The ranks of L_rowmap and U_rowmap do not match.");
  static_assert((int)LEntriesType::rank == (int)UEntriesType::rank,
                "par_ilut_numeric: The ranks of L_entries and U_entries do not match.");
  static_assert((int)LValuesType::rank == (int)UValuesType::rank,
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

  static_assert(std::is_same_v<typename LEntriesType::value_type, typename LEntriesType::non_const_value_type>,
                "par_ilut_numeric: The output L_entries must be nonconst.");
  static_assert(std::is_same_v<typename LValuesType::value_type, typename LValuesType::non_const_value_type>,
                "par_ilut_numeric: The output L_values must be nonconst.");
  static_assert(std::is_same_v<typename UEntriesType::value_type, typename UEntriesType::non_const_value_type>,
                "par_ilut_numeric: The output U_entries must be nonconst.");
  static_assert(std::is_same_v<typename UValuesType::value_type, typename UValuesType::non_const_value_type>,
                "par_ilut_numeric: The output U_values must be nonconst.");

  static_assert(std::is_same_v<typename LRowMapType::device_type, typename ARowMapType::device_type>,
                "par_ilut_numeric: Views LRowMapType and ARowMapType have "
                "different device_types.");
  static_assert(std::is_same_v<typename LEntriesType::device_type, typename AEntriesType::device_type>,
                "par_ilut_numeric: Views LEntriesType and AEntriesType have "
                "different device_types.");
  static_assert(std::is_same_v<typename LValuesType::device_type, typename AValuesType::device_type>,
                "par_ilut_numeric: Views LValuesType and AValuesType have "
                "different device_types.");

  static_assert(std::is_same_v<typename LRowMapType::device_type, typename URowMapType::device_type>,
                "par_ilut_numeric: Views LRowMapType and URowMapType have "
                "different device_types.");
  static_assert(std::is_same_v<typename LEntriesType::device_type, typename UEntriesType::device_type>,
                "par_ilut_numeric: Views LEntriesType and UEntriesType have "
                "different device_types.");
  static_assert(std::is_same_v<typename LValuesType::device_type, typename UValuesType::device_type>,
                "par_ilut_numeric: Views LValuesType and UValuesType have "
                "different device_types.");

  static_assert(std::is_same_v<typename LRowMapType::device_type::execution_space,
                               typename KernelHandle::PAR_ILUTHandleType::execution_space>,
                "par_ilut_numeric: KernelHandle and Views have different execution "
                "spaces.");
  static_assert(std::is_same_v<typename LEntriesType::device_type::execution_space,
                               typename KernelHandle::PAR_ILUTHandleType::execution_space>,
                "par_ilut_numeric: KernelHandle and Views have different execution "
                "spaces.");
  static_assert(std::is_same_v<typename LValuesType::device_type::execution_space,
                               typename KernelHandle::PAR_ILUTHandleType::execution_space>,
                "par_ilut_numeric: KernelHandle and Views have different execution "
                "spaces.");

  static_assert(std::is_same_v<typename LRowMapType::device_type, typename LEntriesType::device_type>,
                "par_ilut_numeric: rowmap and entries have different device types.");
  static_assert(std::is_same_v<typename LRowMapType::device_type, typename LValuesType::device_type>,
                "par_ilut_numeric: rowmap and values have different device types.");

  // Check if symbolic has been called
  if (handle->get_par_ilut_handle()->is_symbolic_complete() == false) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::par_ilut_numeric: par_ilut_symbolic "
          "must be "
          "called before par_ilut_numeric.";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  KK_REQUIRE_MSG(KokkosSparse::Impl::isCrsGraphSorted(L_rowmap, L_entries), "L is not sorted");
  KK_REQUIRE_MSG(KokkosSparse::Impl::isCrsGraphSorted(U_rowmap, U_entries), "U is not sorted");

  using c_size_t   = typename KernelHandle::const_size_type;
  using c_lno_t    = typename KernelHandle::const_nnz_lno_t;
  using c_scalar_t = typename KernelHandle::const_nnz_scalar_t;

  using c_exec_t    = typename KernelHandle::HandleExecSpace;
  using c_temp_t    = typename KernelHandle::HandleTempMemorySpace;
  using c_persist_t = typename KernelHandle::HandlePersistentMemorySpace;

  using const_handle_type = typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t,
                                                                                      c_exec_t, c_temp_t, c_persist_t>;

  const_handle_type tmp_handle(*handle);

  using ARowMap_Internal =
      Kokkos::View<typename ARowMapType::const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<ARowMapType>::array_layout,
                   typename ARowMapType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

  using AEntries_Internal =
      Kokkos::View<typename AEntriesType::const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<AEntriesType>::array_layout,
                   typename AEntriesType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

  using AValues_Internal =
      Kokkos::View<typename AValuesType::const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<AValuesType>::array_layout,
                   typename AValuesType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

  using LRowMap_Internal =
      Kokkos::View<typename LRowMapType::non_const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<LRowMapType>::array_layout,
                   typename LRowMapType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

  using LEntries_Internal =
      Kokkos::View<typename LEntriesType::non_const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<LEntriesType>::array_layout,
                   typename LEntriesType::device_type, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  using LValues_Internal = Kokkos::View<typename LValuesType::non_const_value_type*,
                                        typename KokkosKernels::Impl::GetUnifiedLayout<LValuesType>::array_layout,
                                        typename LValuesType::device_type, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  using URowMap_Internal =
      Kokkos::View<typename URowMapType::non_const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<URowMapType>::array_layout,
                   typename URowMapType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

  using UEntries_Internal =
      Kokkos::View<typename UEntriesType::non_const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<UEntriesType>::array_layout,
                   typename UEntriesType::device_type, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  using UValues_Internal = Kokkos::View<typename UValuesType::non_const_value_type*,
                                        typename KokkosKernels::Impl::GetUnifiedLayout<UValuesType>::array_layout,
                                        typename UValuesType::device_type, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  ARowMap_Internal A_rowmap_i   = A_rowmap;
  AEntries_Internal A_entries_i = A_entries;
  AValues_Internal A_values_i   = A_values;
  LRowMap_Internal L_rowmap_i   = L_rowmap;
  LEntries_Internal L_entries_i = L_entries;
  LValues_Internal L_values_i   = L_values;
  URowMap_Internal U_rowmap_i   = U_rowmap;
  UEntries_Internal U_entries_i = U_entries;
  UValues_Internal U_values_i   = U_values;

  KokkosSparse::Impl::PAR_ILUT_NUMERIC<const_handle_type, ARowMap_Internal, AEntries_Internal, AValues_Internal,
                                       LRowMap_Internal, LEntries_Internal, LValues_Internal, URowMap_Internal,
                                       UEntries_Internal, UValues_Internal>::par_ilut_numeric(&tmp_handle, A_rowmap_i,
                                                                                              A_entries_i, A_values_i,
                                                                                              L_rowmap_i, L_entries_i,
                                                                                              L_values_i, U_rowmap_i,
                                                                                              U_entries_i, U_values_i);

  // These may have been resized
  L_entries = L_entries_i;
  L_values  = L_values_i;
  U_entries = U_entries_i;
  U_values  = U_values_i;

}  // par_ilut_numeric

}  // namespace Experimental
}  // namespace KokkosSparse

#endif  // KOKKOSSPARSE_PAR_ILUT_HPP_
