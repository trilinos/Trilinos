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
#include "KokkosSparse_spadd_symbolic_spec.hpp"
#include "KokkosSparse_spadd_numeric_spec.hpp"

namespace KokkosSparse {
namespace Experimental {

// Symbolic: count entries in each row in C to produce rowmap
// kernel handle has information about whether it is sorted add or not.
template <typename KernelHandle, typename alno_row_view_t_,
          typename alno_nnz_view_t_, typename blno_row_view_t_,
          typename blno_nnz_view_t_, typename clno_row_view_t_>
void spadd_symbolic(
    KernelHandle* handle, const alno_row_view_t_ a_rowmap,
    const alno_nnz_view_t_ a_entries, const blno_row_view_t_ b_rowmap,
    const blno_nnz_view_t_ b_entries,
    clno_row_view_t_ c_rowmap)  // c_rowmap must already be allocated (doesn't
                                // need to be initialized)
{
  typedef typename KernelHandle::HandleExecSpace ExecSpace;
  typedef typename KernelHandle::HandleTempMemorySpace MemSpace;
  typedef typename KernelHandle::HandlePersistentMemorySpace PersistentMemSpace;
  typedef typename Kokkos::Device<ExecSpace, MemSpace> DeviceType;

  typedef typename KernelHandle::const_size_type c_size_t;
  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<
      c_size_t, c_lno_t, c_scalar_t, ExecSpace, MemSpace, PersistentMemSpace>
      ConstKernelHandle;
  ConstKernelHandle tmp_handle(*handle);

  typedef Kokkos::View<typename alno_row_view_t_::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           alno_row_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_a_rowmap;
  typedef Kokkos::View<typename alno_nnz_view_t_::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           alno_nnz_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_a_entries;
  typedef Kokkos::View<typename blno_row_view_t_::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           blno_row_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_b_rowmap;
  typedef Kokkos::View<typename blno_nnz_view_t_::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           blno_nnz_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_b_entries;
  typedef Kokkos::View<typename clno_row_view_t_::value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           clno_row_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_c_rowmap;
  KokkosSparse::Impl::SPADD_SYMBOLIC<ConstKernelHandle, Internal_a_rowmap,
                                     Internal_a_entries, Internal_b_rowmap,
                                     Internal_b_entries, Internal_c_rowmap>::
      spadd_symbolic(&tmp_handle,
                     Internal_a_rowmap(a_rowmap.data(), a_rowmap.extent(0)),
                     Internal_a_entries(a_entries.data(), a_entries.extent(0)),
                     Internal_b_rowmap(b_rowmap.data(), b_rowmap.extent(0)),
                     Internal_b_entries(b_entries.data(), b_entries.extent(0)),
                     Internal_c_rowmap(c_rowmap.data(), c_rowmap.extent(0)));
}

template <typename KernelHandle, typename alno_row_view_t_,
          typename alno_nnz_view_t_, typename ascalar_t_,
          typename ascalar_nnz_view_t_, typename blno_row_view_t_,
          typename blno_nnz_view_t_, typename bscalar_t_,
          typename bscalar_nnz_view_t_, typename clno_row_view_t_,
          typename clno_nnz_view_t_, typename cscalar_nnz_view_t_>
void spadd_numeric(KernelHandle* handle, const alno_row_view_t_ a_rowmap,
                   const alno_nnz_view_t_ a_entries,
                   const ascalar_nnz_view_t_ a_values, const ascalar_t_ alpha,
                   const blno_row_view_t_ b_rowmap,
                   const blno_nnz_view_t_ b_entries,
                   const bscalar_nnz_view_t_ b_values, const bscalar_t_ beta,
                   const clno_row_view_t_ c_rowmap, clno_nnz_view_t_ c_entries,
                   cscalar_nnz_view_t_ c_values) {
  typedef typename KernelHandle::HandleExecSpace ExecSpace;
  typedef typename KernelHandle::HandleTempMemorySpace MemSpace;
  typedef typename KernelHandle::HandlePersistentMemorySpace PersistentMemSpace;
  typedef typename Kokkos::Device<ExecSpace, MemSpace> DeviceType;

  typedef typename KernelHandle::const_size_type c_size_t;
  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<
      c_size_t, c_lno_t, c_scalar_t, ExecSpace, MemSpace, PersistentMemSpace>
      ConstKernelHandle;
  ConstKernelHandle tmp_handle(*handle);

  typedef Kokkos::View<typename alno_row_view_t_::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           alno_row_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_a_rowmap;
  typedef Kokkos::View<typename alno_nnz_view_t_::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           alno_nnz_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_a_entries;
  typedef Kokkos::View<typename ascalar_nnz_view_t_::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           ascalar_nnz_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_a_values;
  typedef Kokkos::View<typename blno_row_view_t_::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           blno_row_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_b_rowmap;
  typedef Kokkos::View<typename blno_nnz_view_t_::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           blno_nnz_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_b_entries;
  typedef Kokkos::View<typename bscalar_nnz_view_t_::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           bscalar_nnz_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_b_values;
  typedef Kokkos::View<typename clno_row_view_t_::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           clno_row_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_c_rowmap;
  typedef Kokkos::View<typename clno_nnz_view_t_::value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           clno_nnz_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_c_entries;
  typedef Kokkos::View<typename cscalar_nnz_view_t_::value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           cscalar_nnz_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_c_values;
  KokkosSparse::Impl::SPADD_NUMERIC<ConstKernelHandle, Internal_a_rowmap,
                                    Internal_a_entries, Internal_a_values,
                                    Internal_b_rowmap, Internal_b_entries,
                                    Internal_b_values, Internal_c_rowmap,
                                    Internal_c_entries, Internal_c_values>::
      spadd_numeric(&tmp_handle, alpha,
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
}
}  // namespace Experimental

// Symbolic: count entries in each row in C to produce rowmap
// kernel handle has information about whether it is sorted add or not.
template <typename KernelHandle, typename AMatrix, typename BMatrix,
          typename CMatrix>
void spadd_symbolic(KernelHandle* handle, const AMatrix& A, const BMatrix& B,
                    CMatrix& C) {
  using row_map_type = typename CMatrix::row_map_type::non_const_type;
  using entries_type = typename CMatrix::index_type::non_const_type;
  using values_type  = typename CMatrix::values_type::non_const_type;

  // Create the row_map of C, no need to initialize it
  row_map_type row_mapC(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "row map"),
      A.numRows() + 1);
  KokkosSparse::Experimental::spadd_symbolic(handle, A.graph.row_map,
                                             A.graph.entries, B.graph.row_map,
                                             B.graph.entries, row_mapC);

  // Now create and allocate the entries and values
  // views so we can build a graph and then matrix C
  // and subsequently construct C.
  auto addHandle = handle->get_spadd_handle();
  entries_type entriesC(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "entries"),
      addHandle->get_c_nnz());
  // Finally since we already have the number of nnz handy
  // we can go ahead and allocate C's values and set them.
  values_type valuesC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "values"),
                      addHandle->get_c_nnz());

  C = CMatrix("matrix", A.numRows(), A.numCols(), addHandle->get_c_nnz(),
              valuesC, row_mapC, entriesC);
}

// Symbolic: count entries in each row in C to produce rowmap
// kernel handle has information about whether it is sorted add or not.
template <typename KernelHandle, typename AScalar, typename AMatrix,
          typename BScalar, typename BMatrix, typename CMatrix>
void spadd_numeric(KernelHandle* handle, const AScalar alpha, const AMatrix& A,
                   const BScalar beta, const BMatrix& B, CMatrix& C) {
  KokkosSparse::Experimental::spadd_numeric(
      handle, A.graph.row_map, A.graph.entries, A.values, alpha,
      B.graph.row_map, B.graph.entries, B.values, beta, C.graph.row_map,
      C.graph.entries, C.values);
}

}  // namespace KokkosSparse

#undef SAME_TYPE

#endif
