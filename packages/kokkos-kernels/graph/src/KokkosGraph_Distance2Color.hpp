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

#ifndef _KOKKOS_GRAPH_COLORDISTANCE2_HPP
#define _KOKKOS_GRAPH_COLORDISTANCE2_HPP

#include "KokkosGraph_Distance2ColorHandle.hpp"
#include "KokkosGraph_Distance2Color_impl.hpp"
#include "KokkosKernels_Utils.hpp"
#include "KokkosSparse_Utils.hpp"

namespace KokkosGraph {

namespace Experimental {

/**
 * Compute the distance-2 coloring of an undirected graph.
 *
 * The graph must be symmetric, but it is not required to have
 * diagonal entries. The coloring will not have distance-1 or distance-2
 * conflicts.
 *
 * @param[in]  handle         The Kernel Handle
 * @param[in]  num_vertices   Number of vertices in the graph
 * @param[in]  row_map        Row map
 * @param[in]  row_entries    Row entries
 *
 * \post
 * <code>handle->get_distance2_graph_coloring_handle()->get_vertex_colors()</code>
 *    will return a view of length num_vertices, containing the colors.
 */

template <class KernelHandle, typename InRowmap, typename InEntries>
void graph_color_distance2(KernelHandle *handle, typename KernelHandle::nnz_lno_t num_verts, InRowmap row_map,
                           InEntries row_entries) {
  using size_type       = typename KernelHandle::size_type;
  using lno_t           = typename KernelHandle::nnz_lno_t;
  using InternalRowmap  = Kokkos::View<const size_type *, Kokkos::LayoutLeft, typename InRowmap::device_type,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using InternalEntries = Kokkos::View<const lno_t *, Kokkos::LayoutLeft, typename InEntries::device_type,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  Kokkos::Timer timer;
  size_type nnz = row_entries.extent(0);
  InternalRowmap rowmap_internal(row_map.data(), row_map.extent(0));
  InternalEntries rowentries_internal(row_entries.data(), nnz);
  auto gch_d2 = handle->get_distance2_graph_coloring_handle();
  // note: last template argument 'false' means do distance-2, not bipartite
  KokkosGraph::Impl::GraphColorDistance2<typename KernelHandle::GraphColorDistance2HandleType, InternalRowmap,
                                         InternalEntries, false>
      gc(num_verts, num_verts, rowmap_internal, rowentries_internal, rowmap_internal, rowentries_internal, gch_d2);
  gc.compute_distance2_color();
  gch_d2->add_to_overall_coloring_time(timer.seconds());
  gch_d2->set_coloring_time(timer.seconds());
}

/**
 * Color the left part (rows) of a bipartite graph: rows r1 and r2 can have the
 * same color if there is no column c such that edges (r1, c) and (r2, c) exist.
 * This means only conflicts over a path exactly two edges long are avoided.
 *
 * This problem is equivalent to grouping the matrix rows into a minimal number
 * of sets, so that within each set, the intersection of any two rows' entries
 * is empty.
 *
 * Distance-1 conflicts (where r1 and c are neighbors) are not avoided,
 * since columns are not colored. In general, there is no particular
 * relationship between a row and column that happen to have the same index.
 *
 * However, if the input graph is symmetric and has diagonal entries in every
 * row, then rows and columns are equivalent and distance-1 conflicts are
 * present through edges (r1, r1) and (r1, r2).
 *
 * @param[in]  handle         The Kernel Handle
 * @param[in]  num_rows       Number of "rows" (vertices in the left part of the
 * graph)
 * @param[in]  num_columns    Number of "columns" (vertices in the right part of
 * the graph)
 * @param[in]  row_map        Row map (CRS format)
 * @param[in]  row_entries    Row entries (CRS format)
 * @param[in]  is_symmetric   Whether (rowmap, row_entries) is known to be
 * symmetric. If it is, this saves computing the transpose internally.
 *
 * \post
 * <code>handle->get_distance2_graph_coloring_handle()->get_vertex_colors()</code>
 *    will return a view of length num_rows, containing the colors.
 */

template <class KernelHandle, typename InRowmap, typename InEntries>
void bipartite_color_rows(KernelHandle *handle, typename KernelHandle::nnz_lno_t num_rows,
                          typename KernelHandle::nnz_lno_t num_columns, InRowmap row_map, InEntries row_entries,
                          bool is_symmetric = false) {
  using execution_space = typename KernelHandle::HandleExecSpace;
  using size_type       = typename KernelHandle::size_type;
  using lno_t           = typename KernelHandle::nnz_lno_t;
  using InternalRowmap  = Kokkos::View<const size_type *, Kokkos::LayoutLeft, typename InRowmap::device_type,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using InternalEntries = Kokkos::View<const lno_t *, Kokkos::LayoutLeft, typename InEntries::device_type,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using TRowmap         = Kokkos::View<size_type *, Kokkos::LayoutLeft, typename InRowmap::device_type>;
  using TEntries        = Kokkos::View<lno_t *, Kokkos::LayoutLeft, typename InEntries::device_type>;
  Kokkos::Timer timer;
  size_type nnz = row_entries.extent(0);
  TRowmap col_map;
  TEntries col_entries;
  if (!is_symmetric) {
    // Compute the transpose
    col_map     = TRowmap("Col map", num_columns + 1);
    col_entries = TEntries("Col entries", nnz);
    KokkosSparse::Impl::transpose_graph<InRowmap, InEntries, TRowmap, TEntries, TRowmap, execution_space>(
        num_rows, num_columns, row_map, row_entries, col_map, col_entries);
  }
  InternalRowmap rowmap_internal(row_map.data(), row_map.extent(0));
  InternalEntries rowentries_internal(row_entries.data(), nnz);
  InternalRowmap colmap_internal;
  InternalEntries colentries_internal;
  if (is_symmetric) {
    colmap_internal     = InternalRowmap(row_map.data(), row_map.extent(0));
    colentries_internal = InternalEntries(row_entries.data(), nnz);
  } else {
    colmap_internal     = InternalRowmap(col_map.data(), col_map.extent(0));
    colentries_internal = InternalEntries(col_entries.data(), nnz);
  }
  auto gch_d2 = handle->get_distance2_graph_coloring_handle();
  // note: last template argument 'true' means do bipartite one-sided
  KokkosGraph::Impl::GraphColorDistance2<typename KernelHandle::GraphColorDistance2HandleType, InternalRowmap,
                                         InternalEntries, true>
      gc(num_rows, num_columns, rowmap_internal, rowentries_internal, colmap_internal, colentries_internal, gch_d2);
  gc.compute_distance2_color();
  gch_d2->add_to_overall_coloring_time(timer.seconds());
  gch_d2->set_coloring_time(timer.seconds());
}

/**
 * Color the right part (columns) of a bipartite graph: columns c1 and c2 can
 * have the same color if there is no row r such that edges (r, c1) and (r, c2)
 * exist.
 *
 * This problem is equivalent to grouping the matrix columns into a minimal
 * number of sets, so that within each set, no two columns appear together in
 * any row's entries. This can be used for computing a compressed Jacobian
 * matrix.
 *
 * Note that the input to this function is still a CRS (row-wise) graph. If you
 * have a CCS (column-wise) or a symmetric graph, use bipartite_color_rows()
 * instead. Calling that with the column-wise graph is equivalent to calling
 * this with the row-wise graph, and that way the transpose will be computed
 * automatically as needed.
 *
 * @param[in]  handle         The Kernel Handle
 * @param[in]  num_rows       Number of "rows" (vertices in the left part of the
 * graph)
 * @param[in]  num_columns    Number of "columns" (vertices in the right part of
 * the graph)
 * @param[in]  row_map        Row map
 * @param[in]  row_entries    Row entries
 *
 * \post handle->get_distance2_graph_coloring_handle()->get_vertex_colors() will
 * return a view of length num_columns, containing the colors.
 */
template <class KernelHandle, typename InRowmap, typename InEntries>
void bipartite_color_columns(KernelHandle *handle, typename KernelHandle::nnz_lno_t num_rows,
                             typename KernelHandle::nnz_lno_t num_columns, InRowmap row_map, InEntries row_entries) {
  using execution_space = typename KernelHandle::HandleExecSpace;
  using size_type       = typename KernelHandle::size_type;
  using lno_t           = typename KernelHandle::nnz_lno_t;
  using InternalRowmap  = Kokkos::View<const size_type *, Kokkos::LayoutLeft, typename InRowmap::device_type,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using InternalEntries = Kokkos::View<const lno_t *, Kokkos::LayoutLeft, typename InEntries::device_type,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using TRowmap         = Kokkos::View<size_type *, Kokkos::LayoutLeft, typename InRowmap::device_type>;
  using TEntries        = Kokkos::View<lno_t *, Kokkos::LayoutLeft, typename InEntries::device_type>;
  Kokkos::Timer timer;
  size_type nnz = row_entries.extent(0);
  // Compute the transpose
  TRowmap col_map("Col map", num_columns + 1);
  TEntries col_entries(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Col entries"), nnz);
  KokkosSparse::Impl::transpose_graph<InRowmap, InEntries, TRowmap, TEntries, TRowmap, execution_space>(
      num_rows, num_columns, row_map, row_entries, col_map, col_entries);
  // Get unmanaged views for both graph and its transpose
  InternalRowmap colmap_internal(col_map.data(), col_map.extent(0));
  InternalEntries colentries_internal(col_entries.data(), nnz);
  InternalRowmap rowmap_internal(row_map.data(), row_map.extent(0));
  InternalEntries rowentries_internal(row_entries.data(), nnz);
  auto gch_d2 = handle->get_distance2_graph_coloring_handle();
  // note: last template argument 'true' means do bipartite one-sided
  KokkosGraph::Impl::GraphColorDistance2<typename KernelHandle::GraphColorDistance2HandleType, InternalRowmap,
                                         InternalEntries, true>
      gc(num_columns, num_rows, colmap_internal, colentries_internal, rowmap_internal, rowentries_internal, gch_d2);
  gc.compute_distance2_color();
  gch_d2->add_to_overall_coloring_time(timer.seconds());
  gch_d2->set_coloring_time(timer.seconds());
}

}  // end namespace Experimental
}  // end namespace KokkosGraph

#endif  //_KOKKOS_GRAPH_COLORDISTANCE2_HPP
