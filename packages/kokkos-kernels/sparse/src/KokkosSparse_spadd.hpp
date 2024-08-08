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

#ifndef _KOKKOS_SPADD_HPP
#define _KOKKOS_SPADD_HPP

#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_helpers.hpp"
#include "KokkosBlas1_scal.hpp"
#include "KokkosSparse_spadd_numeric_spec.hpp"
#include "KokkosSparse_spadd_symbolic_spec.hpp"

namespace KokkosSparse {
namespace Experimental {

// Symbolic: count entries in each row in C to produce rowmap
// kernel handle has information about whether it is sorted add or not.
template <typename ExecSpace, typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_,
          typename blno_row_view_t_, typename blno_nnz_view_t_, typename clno_row_view_t_>
void spadd_symbolic(const ExecSpace &exec, KernelHandle *handle,
                    typename KernelHandle::const_nnz_lno_t m,  // same type as column indices
                    typename KernelHandle::const_nnz_lno_t n, const alno_row_view_t_ a_rowmap,
                    const alno_nnz_view_t_ a_entries, const blno_row_view_t_ b_rowmap, const blno_nnz_view_t_ b_entries,
                    clno_row_view_t_ c_rowmap)  // c_rowmap must already be allocated (doesn't
                                                // need to be initialized)
{
  typedef typename KernelHandle::HandleTempMemorySpace MemSpace;
  typedef typename KernelHandle::HandlePersistentMemorySpace PersistentMemSpace;
  typedef typename Kokkos::Device<ExecSpace, MemSpace> DeviceType;

  typedef typename KernelHandle::const_size_type c_size_t;
  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, ExecSpace, MemSpace,
                                                                    PersistentMemSpace>
      ConstKernelHandle;
  ConstKernelHandle tmp_handle(*handle);

  typedef Kokkos::View<typename alno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<alno_row_view_t_>::array_layout, DeviceType,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      Internal_a_rowmap;
  typedef Kokkos::View<typename alno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<alno_nnz_view_t_>::array_layout, DeviceType,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      Internal_a_entries;
  typedef Kokkos::View<typename blno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<blno_row_view_t_>::array_layout, DeviceType,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      Internal_b_rowmap;
  typedef Kokkos::View<typename blno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<blno_nnz_view_t_>::array_layout, DeviceType,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      Internal_b_entries;
  typedef Kokkos::View<typename clno_row_view_t_::value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<clno_row_view_t_>::array_layout, DeviceType,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      Internal_c_rowmap;

  auto addHandle   = handle->get_spadd_handle();
  bool useFallback = !addHandle->is_input_strict_crs();
  if (useFallback) {
    KokkosSparse::Impl::SPADD_SYMBOLIC<ExecSpace, ConstKernelHandle, Internal_a_rowmap, Internal_a_entries,
                                       Internal_b_rowmap, Internal_b_entries, Internal_c_rowmap,
                                       false>::spadd_symbolic(exec, &tmp_handle, m, n,
                                                              Internal_a_rowmap(a_rowmap.data(), a_rowmap.extent(0)),
                                                              Internal_a_entries(a_entries.data(), a_entries.extent(0)),
                                                              Internal_b_rowmap(b_rowmap.data(), b_rowmap.extent(0)),
                                                              Internal_b_entries(b_entries.data(), b_entries.extent(0)),
                                                              Internal_c_rowmap(c_rowmap.data(), c_rowmap.extent(0)));
  } else {
    KokkosSparse::Impl::SPADD_SYMBOLIC<
        ExecSpace, ConstKernelHandle, Internal_a_rowmap, Internal_a_entries, Internal_b_rowmap, Internal_b_entries,
        Internal_c_rowmap>::spadd_symbolic(exec, &tmp_handle, m, n,
                                           Internal_a_rowmap(a_rowmap.data(), a_rowmap.extent(0)),
                                           Internal_a_entries(a_entries.data(), a_entries.extent(0)),
                                           Internal_b_rowmap(b_rowmap.data(), b_rowmap.extent(0)),
                                           Internal_b_entries(b_entries.data(), b_entries.extent(0)),
                                           Internal_c_rowmap(c_rowmap.data(), c_rowmap.extent(0)));
  }
}

// one without an execution space arg
template <typename KernelHandle, typename... Args>
void spadd_symbolic(KernelHandle *handle, Args... args) {
  spadd_symbolic(typename KernelHandle::HandleExecSpace{}, handle, args...);
}

template <typename ExecSpace, typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_,
          typename ascalar_t_, typename ascalar_nnz_view_t_, typename blno_row_view_t_, typename blno_nnz_view_t_,
          typename bscalar_t_, typename bscalar_nnz_view_t_, typename clno_row_view_t_, typename clno_nnz_view_t_,
          typename cscalar_nnz_view_t_>
void spadd_numeric(const ExecSpace &exec, KernelHandle *handle, typename KernelHandle::const_nnz_lno_t m,
                   typename KernelHandle::const_nnz_lno_t n, const alno_row_view_t_ a_rowmap,
                   const alno_nnz_view_t_ a_entries, const ascalar_nnz_view_t_ a_values, const ascalar_t_ alpha,
                   const blno_row_view_t_ b_rowmap, const blno_nnz_view_t_ b_entries,
                   const bscalar_nnz_view_t_ b_values, const bscalar_t_ beta, const clno_row_view_t_ c_rowmap,
                   clno_nnz_view_t_ c_entries, cscalar_nnz_view_t_ c_values) {
  typedef typename KernelHandle::HandleTempMemorySpace MemSpace;
  typedef typename KernelHandle::HandlePersistentMemorySpace PersistentMemSpace;
  typedef typename Kokkos::Device<ExecSpace, MemSpace> DeviceType;

  typedef typename KernelHandle::const_size_type c_size_t;
  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, ExecSpace, MemSpace,
                                                                    PersistentMemSpace>
      ConstKernelHandle;
  ConstKernelHandle tmp_handle(*handle);  // handle->exec_space is also copied

  typedef Kokkos::View<typename alno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<alno_row_view_t_>::array_layout, DeviceType,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      Internal_a_rowmap;
  typedef Kokkos::View<typename alno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<alno_nnz_view_t_>::array_layout, DeviceType,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      Internal_a_entries;
  typedef Kokkos::View<typename ascalar_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<ascalar_nnz_view_t_>::array_layout, DeviceType,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      Internal_a_values;
  typedef Kokkos::View<typename blno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<blno_row_view_t_>::array_layout, DeviceType,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      Internal_b_rowmap;
  typedef Kokkos::View<typename blno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<blno_nnz_view_t_>::array_layout, DeviceType,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      Internal_b_entries;
  typedef Kokkos::View<typename bscalar_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<bscalar_nnz_view_t_>::array_layout, DeviceType,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      Internal_b_values;
  typedef Kokkos::View<typename clno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<clno_row_view_t_>::array_layout, DeviceType,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      Internal_c_rowmap;
  typedef Kokkos::View<typename clno_nnz_view_t_::value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<clno_nnz_view_t_>::array_layout, DeviceType,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      Internal_c_entries;
  typedef Kokkos::View<typename cscalar_nnz_view_t_::value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<cscalar_nnz_view_t_>::array_layout, DeviceType,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      Internal_c_values;

  auto addHandle   = handle->get_spadd_handle();
  bool useFallback = !addHandle->is_input_strict_crs();
  if (useFallback) {
    KokkosSparse::Impl::SPADD_NUMERIC<ExecSpace, ConstKernelHandle, Internal_a_rowmap, Internal_a_entries,
                                      Internal_a_values, Internal_b_rowmap, Internal_b_entries, Internal_b_values,
                                      Internal_c_rowmap, Internal_c_entries, Internal_c_values,
                                      false>::spadd_numeric(exec, &tmp_handle, m, n, alpha,
                                                            Internal_a_rowmap(a_rowmap.data(), a_rowmap.extent(0)),
                                                            Internal_a_entries(a_entries.data(), a_entries.extent(0)),
                                                            Internal_a_values(a_values.data(), a_values.extent(0)),
                                                            beta,
                                                            Internal_b_rowmap(b_rowmap.data(), b_rowmap.extent(0)),
                                                            Internal_b_entries(b_entries.data(), b_entries.extent(0)),
                                                            Internal_b_values(b_values.data(), b_values.extent(0)),
                                                            Internal_c_rowmap(c_rowmap.data(), c_rowmap.extent(0)),
                                                            Internal_c_entries(c_entries.data(), c_entries.extent(0)),
                                                            Internal_c_values(c_values.data(), c_values.extent(0)));
  } else {
    KokkosSparse::Impl::SPADD_NUMERIC<
        ExecSpace, ConstKernelHandle, Internal_a_rowmap, Internal_a_entries, Internal_a_values, Internal_b_rowmap,
        Internal_b_entries, Internal_b_values, Internal_c_rowmap, Internal_c_entries,
        Internal_c_values>::spadd_numeric(exec, &tmp_handle, m, n, alpha,
                                          Internal_a_rowmap(a_rowmap.data(), a_rowmap.extent(0)),
                                          Internal_a_entries(a_entries.data(), a_entries.extent(0)),
                                          Internal_a_values(a_values.data(), a_values.extent(0)), beta,
                                          Internal_b_rowmap(b_rowmap.data(), b_rowmap.extent(0)),
                                          Internal_b_entries(b_entries.data(), b_entries.extent(0)),
                                          Internal_b_values(b_values.data(), b_values.extent(0)),
                                          Internal_c_rowmap(c_rowmap.data(), c_rowmap.extent(0)),
                                          Internal_c_entries(c_entries.data(), c_entries.extent(0)),
                                          Internal_c_values(c_values.data(), c_values.extent(0)));
  }
}

// one without an execution space arg
template <typename KernelHandle, typename... Args>
void spadd_numeric(KernelHandle *handle, Args... args) {
  spadd_numeric(typename KernelHandle::HandleExecSpace{}, handle, args...);
}
}  // namespace Experimental

// Symbolic: count entries in each row in C to produce rowmap
// kernel handle has information about whether it is sorted add or not.
template <typename ExecSpace, typename KernelHandle, typename AMatrix, typename BMatrix, typename CMatrix>
void spadd_symbolic(const ExecSpace &exec, KernelHandle *handle, const AMatrix &A, const BMatrix &B, CMatrix &C) {
  using row_map_type = typename CMatrix::row_map_type::non_const_type;
  using entries_type = typename CMatrix::index_type::non_const_type;
  using values_type  = typename CMatrix::values_type::non_const_type;

  auto addHandle = handle->get_spadd_handle();

  // Create the row_map of C, no need to initialize it
  row_map_type row_mapC(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "row map"), A.numRows() + 1);

  // Shortcuts for special cases as they cause errors in some TPL
  // implementations (e.g., cusparse and hipsparse)
  if (!A.nnz()) {
    Kokkos::deep_copy(exec, row_mapC, B.graph.row_map);
    addHandle->set_c_nnz(B.graph.entries.extent(0));
  } else if (!B.nnz()) {
    Kokkos::deep_copy(exec, row_mapC, A.graph.row_map);
    addHandle->set_c_nnz(A.graph.entries.extent(0));
  } else {
    KokkosSparse::Experimental::spadd_symbolic(exec, handle, A.numRows(), A.numCols(), A.graph.row_map, A.graph.entries,
                                               B.graph.row_map, B.graph.entries, row_mapC);
  }

  // Now create and allocate the entries and values
  // views so we can build a graph and then matrix C
  // and subsequently construct C.
  entries_type entriesC(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "entries"), addHandle->get_c_nnz());
  // Finally since we already have the number of nnz handy
  // we can go ahead and allocate C's values and set them.
  values_type valuesC(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "values"), addHandle->get_c_nnz());

  C = CMatrix("matrix", A.numRows(), A.numCols(), addHandle->get_c_nnz(), valuesC, row_mapC, entriesC);
}

// Numeric: fill the column indices and values
// kernel handle has information about whether it is sorted add or not.
template <typename ExecSpace, typename KernelHandle, typename AScalar, typename AMatrix, typename BScalar,
          typename BMatrix, typename CMatrix>
void spadd_numeric(const ExecSpace &exec, KernelHandle *handle, const AScalar alpha, const AMatrix &A,
                   const BScalar beta, const BMatrix &B, CMatrix &C) {
  if (!A.nnz()) {
    Kokkos::deep_copy(exec, C.graph.entries, B.graph.entries);
    KokkosBlas::scal(exec, C.values, beta, B.values);
  } else if (!B.nnz()) {
    Kokkos::deep_copy(exec, C.graph.entries, A.graph.entries);
    KokkosBlas::scal(exec, C.values, alpha, A.values);
  } else {
    KokkosSparse::Experimental::spadd_numeric(exec, handle, A.numRows(), A.numCols(), A.graph.row_map, A.graph.entries,
                                              A.values, alpha, B.graph.row_map, B.graph.entries, B.values, beta,
                                              C.graph.row_map, C.graph.entries, C.values);
  }
}

// One without an explicit execution space argument
template <typename KernelHandle, typename AMatrix, typename BMatrix, typename CMatrix>
void spadd_symbolic(KernelHandle *handle, const AMatrix &A, const BMatrix &B, CMatrix &C) {
  spadd_symbolic(typename AMatrix::execution_space{}, handle, A, B, C);
}

template <typename KernelHandle, typename AScalar, typename AMatrix, typename BScalar, typename BMatrix,
          typename CMatrix>
void spadd_numeric(KernelHandle *handle, const AScalar alpha, const AMatrix &A, const BScalar beta, const BMatrix &B,
                   CMatrix &C) {
  spadd_numeric(typename AMatrix::execution_space{}, handle, alpha, A, beta, B, C);
}

}  // namespace KokkosSparse

#undef SAME_TYPE

#endif
