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
#ifndef _KOKKOSKERNELS_SPARSEUTILS_HPP
#define _KOKKOSKERNELS_SPARSEUTILS_HPP
#include <vector>

#include "Kokkos_Core.hpp"
#include "KokkosKernels_SimpleUtils.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"
#include "KokkosKernels_PrintUtils.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_BsrMatrix.hpp"
#include "Kokkos_Bitset.hpp"
#include "KokkosGraph_RCM.hpp"

#ifdef KOKKOSKERNELS_HAVE_PARALLEL_GNUSORT
#include <parallel/algorithm>
#endif

namespace KokkosSparse {

enum class SparseMatrixFormat {
  BSR,
  CRS,
};

namespace Impl {

/* create a block-sparse version of a CrsMatrix
 */
template <typename in_row_view_t, typename in_nnz_view_t, typename in_val_view_t, typename out_row_view_t,
          typename out_nnz_view_t, typename out_val_view_t>
void kk_create_bsr_formated_point_crsmatrix(int block_size, size_t num_rows, size_t num_cols, in_row_view_t in_xadj,
                                            in_nnz_view_t in_adj, in_val_view_t in_vals, size_t &out_num_rows,
                                            size_t &out_num_cols, out_row_view_t &out_xadj, out_nnz_view_t &out_adj,
                                            out_val_view_t &out_vals) {
  typedef typename in_nnz_view_t::non_const_value_type lno_t;
  typedef typename in_row_view_t::non_const_value_type size_type;
  typedef typename in_val_view_t::non_const_value_type scalar_t;

  typename in_row_view_t::HostMirror hr = Kokkos::create_mirror_view(in_xadj);
  Kokkos::deep_copy(hr, in_xadj);
  typename in_nnz_view_t::HostMirror he = Kokkos::create_mirror_view(in_adj);
  Kokkos::deep_copy(he, in_adj);
  typename in_val_view_t::HostMirror hv = Kokkos::create_mirror_view(in_vals);
  Kokkos::deep_copy(hv, in_vals);

  out_num_rows = (num_rows / block_size) * block_size;
  if (num_rows % block_size) out_num_rows += block_size;

  out_num_cols = (num_cols / block_size) * block_size;
  if (num_cols % block_size) out_num_cols += block_size;

  std::vector<size_type> block_rows_xadj(out_num_rows + 1, 0);
  std::vector<lno_t> block_adj;      //(in_adj.extent(0), 0);
  std::vector<scalar_t> block_vals;  // (in_adj.extent(0), 0);

  std::vector<lno_t> block_columns(out_num_cols, 0);
  std::vector<scalar_t> block_accumulators(out_num_cols, 0);
  std::vector<bool> block_flags(out_num_cols, false);

  // loop over first rows of each block-row
  for (lno_t i = 0; i < lno_t(num_rows); i += block_size) {
    lno_t outputrowsize = 0;

    // loop over rows in block
    for (lno_t block_ind = 0; block_ind < block_size; ++block_ind) {
      const lno_t row_ind = block_ind + i;
      if (row_ind < lno_t(num_rows)) {
        size_type adj_begin = hr(row_ind);
        size_type adj_end   = hr(row_ind + 1);

        lno_t row_size = adj_end - adj_begin;

        for (lno_t col_ind = 0; col_ind < row_size; ++col_ind) {
          lno_t colid = he(col_ind + adj_begin);
          // scalar_t colval = hv(col_ind);

          lno_t block_column_start = (colid / block_size) * block_size;

          for (lno_t kk = 0; kk < block_size; ++kk) {
            lno_t col_id_to_insert = block_column_start + kk;
            // std::cout << colid << " " << block_column_start << " " <<
            // col_id_to_insert << " " << block_flags[col_id_to_insert] << " ##
            // ";

            if (block_flags[col_id_to_insert] == false) {
              block_flags[col_id_to_insert] = true;
              // block_adj[output_index + outputrowsize++] = col_id_to_insert;
              // block_adj.push_back(col_id_to_insert);
              block_columns[outputrowsize++] = col_id_to_insert;
            }
          }
        }
      } else {
        lno_t colid = row_ind;
        // scalar_t colval = hv(col_ind);

        lno_t block_column_start = (colid / block_size) * block_size;

        for (lno_t kk = 0; kk < block_size; ++kk) {
          lno_t col_id_to_insert = block_column_start + kk;
          if (block_flags[col_id_to_insert] == false) {
            block_flags[col_id_to_insert] = true;
            // block_adj[output_index + outputrowsize++] = col_id_to_insert;
            // block_adj.push_back(col_id_to_insert);
            block_columns[outputrowsize++] = col_id_to_insert;
          }
        }
      }
    }
    std::sort(block_columns.begin(), block_columns.begin() + outputrowsize);
    // std::cout << "\nrow:" << i << " outputrowsize:" << outputrowsize <<
    // std::endl;
    for (lno_t kk = 0; kk < outputrowsize; ++kk) {
      block_flags[block_columns[kk]] = false;
      // std::cout << block_columns[kk] << " ";
    }
    // std::cout << std::endl;

    for (lno_t block_ind = 0; block_ind < block_size; ++block_ind) {
      lno_t row_ind = block_ind + i;
      if (row_ind < lno_t(num_rows)) {
        size_type adj_begin = hr(row_ind);
        size_type adj_end   = hr(row_ind + 1);

        lno_t row_size = adj_end - adj_begin;

        for (lno_t col_ind = 0; col_ind < row_size; ++col_ind) {
          lno_t colid               = he(col_ind + adj_begin);
          scalar_t colval           = hv(col_ind + adj_begin);
          block_accumulators[colid] = colval;
        }
      } else {
        block_accumulators[row_ind] = 1;
      }

      for (lno_t kk = 0; kk < outputrowsize; ++kk) {
        lno_t outcol = block_columns[kk];
        block_adj.push_back(outcol);
        block_vals.push_back(block_accumulators[outcol]);
        block_accumulators[outcol] = 0;
      }
      block_rows_xadj[row_ind + 1] = block_rows_xadj[row_ind] + outputrowsize;
    }
  }

  out_xadj = out_row_view_t("BlockedPointCRS XADJ", out_num_rows + 1);
  out_adj  = out_nnz_view_t("BlockedPointCRS ADJ", block_adj.size());
  out_vals = out_val_view_t("BlockedPointCRS VALS", block_vals.size());

  typename out_row_view_t::HostMirror hor = Kokkos::create_mirror_view(out_xadj);
  typename out_nnz_view_t::HostMirror hoe = Kokkos::create_mirror_view(out_adj);
  typename out_val_view_t::HostMirror hov = Kokkos::create_mirror_view(out_vals);

  for (lno_t i = 0; i < lno_t(out_num_rows) + 1; ++i) {
    hor(i) = block_rows_xadj[i];
  }

  size_type ne = block_adj.size();
  for (size_type i = 0; i < ne; ++i) {
    hoe(i) = block_adj[i];
  }
  for (size_type i = 0; i < ne; ++i) {
    hov(i) = block_vals[i];
  }

  Kokkos::deep_copy(out_xadj, hor);
  Kokkos::deep_copy(out_adj, hoe);
  Kokkos::deep_copy(out_vals, hov);
}

/* can be used in a parallel for to copy between compatible views
 */
template <typename T, typename U>
struct ViewConverter {
  ViewConverter(const T &_dst, const U &_src) : dst(_dst), src(_src) {}
  KOKKOS_INLINE_FUNCTION void operator()(size_t i) const { dst[i] = src[i]; }
  T dst;
  U src;
};

/* Create output row pointer, col index, and value arrays
   for BSR-format data from CRS data consistent with BSR format

*/
template <typename in_row_view_t, typename in_nnz_view_t, typename in_val_view_t,
          typename out_row_view_t /*output row_map view type*/, typename out_nnz_view_t, typename out_val_view_t>
void kk_create_bsr_from_bsr_formatted_point_crs(int block_size, size_t num_rows, size_t num_cols,
                                                in_row_view_t in_xadj,  // row pointer (CrsMatrix::graph.row_map)
                                                in_nnz_view_t in_adj,   // col index (CrsMatrix::graph.entries)
                                                in_val_view_t in_vals,  // values CrsMatrix::values
                                                size_t &out_num_rows,   // rows of blocks in output
                                                size_t &out_num_cols,   // cols of blocks in output
                                                out_row_view_t &out_xadj, out_nnz_view_t &out_adj,
                                                out_val_view_t &out_vals) {
  typedef typename in_nnz_view_t::non_const_value_type in_ordinal_type;
  typedef typename in_val_view_t::non_const_value_type in_scalar_type;
  typedef typename in_nnz_view_t::device_type in_device_type;
  typedef typename out_nnz_view_t::non_const_value_type out_ordinal_type;
  typedef typename out_val_view_t::non_const_value_type out_scalar_type;
  typedef typename out_nnz_view_t::device_type out_device_type;

  // in_row_view_t and out_row_view_t may not be the same, so use ViewConverter
  // to do the conversion
  typedef KokkosSparse::CrsMatrix<in_scalar_type, in_ordinal_type, in_device_type> InMatrix;
  typedef KokkosSparse::Experimental::BsrMatrix<out_scalar_type, out_ordinal_type, out_device_type> OutMatrix;

  // in_rowmap <- in_xadj
  Kokkos::View<typename InMatrix::non_const_size_type *, in_device_type> in_rowmap("", in_xadj.size());
  Kokkos::parallel_for("", in_xadj.size(), ViewConverter<decltype(in_rowmap), in_row_view_t>(in_rowmap, in_xadj));

  // reconstruct original CrsMatrix
  InMatrix in("", num_rows, num_cols, in_vals.size(), in_vals, in_rowmap, in_adj);

  // convert to BsrMatrix
  OutMatrix out(in, block_size);

  // out_xadj <- out.graph.row_map
  Kokkos::resize(out_xadj, out.graph.row_map.size());
  Kokkos::parallel_for("", out_xadj.size(),
                       ViewConverter<out_row_view_t, decltype(out.graph.row_map)>(out_xadj, out.graph.row_map));

  out_adj      = out.graph.entries;
  out_vals     = out.values;
  out_num_rows = out.numRows();
  out_num_cols = out.numCols();
}

template <typename in_row_view_t, typename in_nnz_view_t, typename in_scalar_view_t, typename out_row_view_t,
          typename out_nnz_view_t, typename out_scalar_view_t, typename tempwork_row_view_t, typename MyExecSpace>
struct TransposeMatrix {
  struct CountTag {};
  struct FillTag {};

  using team_count_policy_t = Kokkos::TeamPolicy<CountTag, MyExecSpace>;
  using team_fill_policy_t  = Kokkos::TeamPolicy<FillTag, MyExecSpace>;

  using team_count_member_t = typename team_count_policy_t::member_type;
  using team_fill_member_t  = typename team_fill_policy_t::member_type;

  using nnz_lno_t = typename in_nnz_view_t::non_const_value_type;
  using size_type = typename in_row_view_t::non_const_value_type;

  nnz_lno_t num_rows;
  nnz_lno_t num_cols;
  in_row_view_t xadj;
  in_nnz_view_t adj;
  in_scalar_view_t vals;
  out_row_view_t t_xadj;     // allocated
  out_nnz_view_t t_adj;      // allocated
  out_scalar_view_t t_vals;  // allocated
  tempwork_row_view_t tmp_txadj;
  bool transpose_values;
  nnz_lno_t team_work_size;

  TransposeMatrix(nnz_lno_t num_rows_, nnz_lno_t num_cols_, in_row_view_t xadj_, in_nnz_view_t adj_,
                  in_scalar_view_t vals_, out_row_view_t t_xadj_, out_nnz_view_t t_adj_, out_scalar_view_t t_vals_,
                  tempwork_row_view_t tmp_txadj_, bool transpose_values_, nnz_lno_t team_row_work_size_)
      : num_rows(num_rows_),
        num_cols(num_cols_),
        xadj(xadj_),
        adj(adj_),
        vals(vals_),
        t_xadj(t_xadj_),
        t_adj(t_adj_),
        t_vals(t_vals_),
        tmp_txadj(tmp_txadj_),
        transpose_values(transpose_values_),
        team_work_size(team_row_work_size_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const CountTag &, const team_count_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, num_rows);
    // TODO we dont need to go over rows
    // just go over nonzeroes.
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_index) {
                           const size_type col_begin = xadj[row_index];
                           const size_type col_end   = xadj[row_index + 1];
                           const nnz_lno_t left_work = col_end - col_begin;
                           Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, left_work), [&](nnz_lno_t i) {
                             const size_type adjind   = i + col_begin;
                             const nnz_lno_t colIndex = adj[adjind];
                             typedef typename std::remove_reference<decltype(t_xadj(0))>::type atomic_incr_type;
                             Kokkos::atomic_fetch_add(&(t_xadj(colIndex)), atomic_incr_type(1));
                           });
                         });
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const FillTag &, const team_fill_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, num_rows);

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&](const nnz_lno_t &row_index) {
          // const nnz_lno_t teamsize = teamMember.team_size();
          // for (nnz_lno_t row_index = team_row_begin + teamMember.team_rank();
          // row_index < team_row_end; row_index += teamsize){
          const size_type col_begin = xadj[row_index];
          const size_type col_end   = xadj[row_index + 1];
          const nnz_lno_t left_work = col_end - col_begin;
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, left_work), [&](nnz_lno_t i) {
            const size_type adjind   = i + col_begin;
            const nnz_lno_t colIndex = adj[adjind];
            typedef typename std::remove_reference<decltype(tmp_txadj(0))>::type atomic_incr_type;
            const size_type pos = Kokkos::atomic_fetch_add(&(tmp_txadj(colIndex)), atomic_incr_type(1));

            t_adj(pos) = row_index;
            if (transpose_values) {
              t_vals(pos) = vals[adjind];
            }
          });
          //}
        });
  }
};

template <typename in_row_view_t, typename in_nnz_view_t, typename in_scalar_view_t, typename out_row_view_t,
          typename out_nnz_view_t, typename out_scalar_view_t, typename tempwork_row_view_t, typename MyExecSpace>
void transpose_matrix(typename in_nnz_view_t::non_const_value_type num_rows,
                      typename in_nnz_view_t::non_const_value_type num_cols, in_row_view_t xadj, in_nnz_view_t adj,
                      in_scalar_view_t vals,
                      out_row_view_t t_xadj,    // pre-allocated -- initialized with 0
                      out_nnz_view_t t_adj,     // pre-allocated -- no need for initialize
                      out_scalar_view_t t_vals  // pre-allocated -- no need for initialize
) {
  // allocate some memory for work for row pointers
  tempwork_row_view_t tmp_row_view(Kokkos::view_alloc(Kokkos::WithoutInitializing, "tmp_row_view"), num_cols + 1);

  // create the functor for tranpose.
  typedef TransposeMatrix<in_row_view_t, in_nnz_view_t, in_scalar_view_t, out_row_view_t, out_nnz_view_t,
                          out_scalar_view_t, tempwork_row_view_t, MyExecSpace>
      TransposeFunctor_t;

  typedef typename TransposeFunctor_t::team_count_policy_t count_tp_t;
  typedef typename TransposeFunctor_t::team_fill_policy_t fill_tp_t;

  typename in_row_view_t::non_const_value_type nnz = adj.extent(0);

  // determine vector lanes per thread
  int thread_size =
      kk_get_suggested_vector_size(num_rows, nnz, KokkosKernels::Impl::kk_get_exec_space_type<MyExecSpace>());

  // determine threads per team
  int team_size = kk_get_suggested_team_size(thread_size, KokkosKernels::Impl::kk_get_exec_space_type<MyExecSpace>());

  TransposeFunctor_t tm(num_rows, num_cols, xadj, adj, vals, t_xadj, t_adj, t_vals, tmp_row_view, true, team_size);

  Kokkos::parallel_for("KokkosSparse::Impl::transpose_matrix::S0",
                       count_tp_t((num_rows + team_size - 1) / team_size, team_size, thread_size), tm);

  KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<MyExecSpace>(num_cols + 1, t_xadj);

  Kokkos::deep_copy(tmp_row_view, t_xadj);

  Kokkos::parallel_for("KokkosSparse::Impl::transpose_matrix::S1",
                       fill_tp_t((num_rows + team_size - 1) / team_size, team_size, thread_size), tm);

  MyExecSpace().fence();
}

template <typename crsMat_t>
crsMat_t transpose_matrix(const crsMat_t &A) {
  // Allocate views and call the other version of transpose_matrix
  using c_rowmap_t  = typename crsMat_t::row_map_type;
  using c_entries_t = typename crsMat_t::index_type;
  using c_values_t  = typename crsMat_t::values_type;
  using rowmap_t    = typename crsMat_t::row_map_type::non_const_type;
  using entries_t   = typename crsMat_t::index_type::non_const_type;
  using values_t    = typename crsMat_t::values_type::non_const_type;
  rowmap_t AT_rowmap("Transpose rowmap", A.numCols() + 1);
  entries_t AT_entries(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Transpose entries"), A.nnz());
  values_t AT_values(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Transpose values"), A.nnz());
  transpose_matrix<c_rowmap_t, c_entries_t, c_values_t, rowmap_t, entries_t, values_t, rowmap_t,
                   typename crsMat_t::execution_space>(A.numRows(), A.numCols(), A.graph.row_map, A.graph.entries,
                                                       A.values, AT_rowmap, AT_entries, AT_values);
  // And construct the transpose crsMat_t
  return crsMat_t("Transpose", A.numCols(), A.numRows(), A.nnz(), AT_values, AT_rowmap, AT_entries);
}

template <typename in_row_view_t, typename in_nnz_view_t, typename out_row_view_t, typename out_nnz_view_t,
          typename tempwork_row_view_t, typename MyExecSpace>
void transpose_graph(typename in_nnz_view_t::non_const_value_type num_rows,
                     typename in_nnz_view_t::non_const_value_type num_cols, in_row_view_t xadj, in_nnz_view_t adj,
                     out_row_view_t t_xadj,  // pre-allocated -- initialized with 0
                     out_nnz_view_t t_adj    // pre-allocated -- no need for initialize
) {
  // allocate some memory for work for row pointers
  tempwork_row_view_t tmp_row_view(Kokkos::view_alloc(Kokkos::WithoutInitializing, "tmp_row_view"), num_cols + 1);

  in_nnz_view_t tmp1;
  out_nnz_view_t tmp2;

  // create the functor for tranpose.
  typedef TransposeMatrix<in_row_view_t, in_nnz_view_t, in_nnz_view_t, out_row_view_t, out_nnz_view_t, out_nnz_view_t,
                          tempwork_row_view_t, MyExecSpace>
      TransposeFunctor_t;

  typedef typename TransposeFunctor_t::team_count_policy_t count_tp_t;
  typedef typename TransposeFunctor_t::team_fill_policy_t fill_tp_t;

  typename in_row_view_t::non_const_value_type nnz = adj.extent(0);

  // determine vector lanes per thread
  int thread_size =
      kk_get_suggested_vector_size(num_rows, nnz, KokkosKernels::Impl::kk_get_exec_space_type<MyExecSpace>());

  // determine threads per team
  int team_size = kk_get_suggested_team_size(thread_size, KokkosKernels::Impl::kk_get_exec_space_type<MyExecSpace>());

  TransposeFunctor_t tm(num_rows, num_cols, xadj, adj, tmp1, t_xadj, t_adj, tmp2, tmp_row_view, false, team_size);

  Kokkos::parallel_for("KokkosKernels::Impl::transpose_graph::S0",
                       count_tp_t((num_rows + team_size - 1) / team_size, team_size, thread_size), tm);

  KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<MyExecSpace>(num_cols + 1, t_xadj);

  Kokkos::deep_copy(tmp_row_view, t_xadj);

  Kokkos::parallel_for("KokkosKernels::Impl::transpose_graph::S1",
                       fill_tp_t((num_rows + team_size - 1) / team_size, team_size, thread_size), tm);

  MyExecSpace().fence();
}

template <typename in_row_view_t, typename in_nnz_view_t, typename in_scalar_view_t, typename out_row_view_t,
          typename out_nnz_view_t, typename out_scalar_view_t>
struct TransposeBsrMatrix {
  using ordinal_type = typename in_nnz_view_t::non_const_value_type;
  using size_type    = typename in_row_view_t::non_const_value_type;

  int block_size;
  in_row_view_t Arow_map;
  in_nnz_view_t Aentries;
  in_scalar_view_t Avalues;
  out_row_view_t tArow_map;    // allocated
  out_nnz_view_t tAentries;    // allocated
  out_scalar_view_t tAvalues;  // allocated

  TransposeBsrMatrix(const int blockSize, in_row_view_t row_mapA, in_nnz_view_t entriesA, in_scalar_view_t valuesA,
                     out_row_view_t row_mapAt, out_nnz_view_t entriesAt, out_scalar_view_t valuesAt)
      : block_size(blockSize),
        Arow_map(row_mapA),
        Aentries(entriesA),
        Avalues(valuesA),
        tArow_map(row_mapAt),
        tAentries(entriesAt),
        tAvalues(valuesAt){};

  KOKKOS_INLINE_FUNCTION
  void operator()(const int tArowIdx) const {
    // Loop over entries in row
    for (size_type tAentryIdx = tArow_map(tArowIdx); tAentryIdx < tArow_map(tArowIdx + 1); ++tAentryIdx) {
      ordinal_type tAcolIdx = tAentries(tAentryIdx);

      // we have block tA(tArowIdx, tAcolIdx) starting at tAvalues(entryIdx)
      // we need to find AentryIdx corresponding to A(tAcolIdx, tArowIdx)
      size_type AentryIdx;
      for (AentryIdx = Arow_map(tAcolIdx); AentryIdx < Arow_map(tAcolIdx + 1); ++AentryIdx) {
        if (tArowIdx == Aentries(AentryIdx)) break;
      }

      // we loop over block_size*block_size Avalues starting at AentryIdx
      // and store them into tAvalues in transpose order starting at tAentryIdx
      for (int i = 0; i < block_size; ++i) {
        for (int j = 0; j < block_size; ++j) {
          tAvalues(tAentryIdx * block_size * block_size + i * block_size + j) =
              Avalues(AentryIdx * block_size * block_size + j * block_size + i);
        }
      }
    }
  }
};  // TransposeBsrMatrix

template <typename in_row_view_t, typename in_nnz_view_t, typename in_scalar_view_t, typename out_row_view_t,
          typename out_nnz_view_t, typename out_scalar_view_t, typename MyExecSpace>
void transpose_bsr_matrix(typename in_nnz_view_t::non_const_value_type num_rows,
                          typename in_nnz_view_t::non_const_value_type num_cols, const int block_size,
                          in_row_view_t xadj, in_nnz_view_t adj, in_scalar_view_t vals,
                          out_row_view_t t_xadj,    // pre-allocated -- initialized with 0
                          out_nnz_view_t t_adj,     // pre-allocated -- no need for initialize
                          out_scalar_view_t t_vals  // pre-allocated -- no need for initialize
) {
  using TransposeBsrFunctor_type = TransposeBsrMatrix<in_row_view_t, in_nnz_view_t, in_scalar_view_t, out_row_view_t,
                                                      out_nnz_view_t, out_scalar_view_t>;

  // Step 1: call transpose_graph of bsr matrix
  transpose_graph<in_row_view_t, in_nnz_view_t, out_row_view_t, out_nnz_view_t, out_row_view_t, MyExecSpace>(
      num_rows, num_cols, xadj, adj, t_xadj, t_adj);

  // Step 2: transpose the values of A
  Kokkos::RangePolicy<MyExecSpace> my_policy(0, num_cols);
  TransposeBsrFunctor_type my_functor(block_size, xadj, adj, vals, t_xadj, t_adj, t_vals);

  Kokkos::parallel_for(my_policy, my_functor);
  MyExecSpace().fence();
}

template <typename bsrMat_t>
bsrMat_t transpose_bsr_matrix(const bsrMat_t &A) {
  // Allocate views and call the other version of transpose_matrix
  using c_rowmap_t  = typename bsrMat_t::row_map_type;
  using c_entries_t = typename bsrMat_t::index_type;
  using c_values_t  = typename bsrMat_t::values_type;
  using rowmap_t    = typename bsrMat_t::row_map_type::non_const_type;
  using entries_t   = typename bsrMat_t::index_type::non_const_type;
  using values_t    = typename bsrMat_t::values_type::non_const_type;

  rowmap_t AT_rowmap("Transpose rowmap", A.numCols() + 1);
  entries_t AT_entries(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Transpose entries"), A.nnz());
  values_t AT_values(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Transpose values"),
                     A.nnz() * A.blockDim() * A.blockDim());
  transpose_bsr_matrix<c_rowmap_t, c_entries_t, c_values_t, rowmap_t, entries_t, values_t,
                       typename bsrMat_t::execution_space>(A.numRows(), A.numCols(), A.blockDim(), A.graph.row_map,
                                                           A.graph.entries, A.values, AT_rowmap, AT_entries, AT_values);
  // And construct the transpose crsMat_t
  return bsrMat_t("Transpose", A.numCols(), A.numRows(), A.nnz(), AT_values, AT_rowmap, AT_entries, A.blockDim());
}

template <typename forward_map_type, typename reverse_map_type>
struct Fill_Reverse_Scale_Functor {
  struct CountTag {};
  struct FillTag {};

  typedef struct CountTag CountTag;
  typedef struct FillTag FillTag;

  typedef typename forward_map_type::value_type forward_type;
  typedef typename reverse_map_type::value_type reverse_type;
  forward_map_type forward_map;
  reverse_map_type reverse_map_xadj;
  reverse_map_type reverse_map_adj;

  const reverse_type multiply_shift_for_scale;
  const reverse_type division_shift_for_bucket;

  Fill_Reverse_Scale_Functor(forward_map_type forward_map_, reverse_map_type reverse_map_xadj_,
                             reverse_map_type reverse_map_adj_, reverse_type multiply_shift_for_scale_,
                             reverse_type division_shift_for_bucket_)
      : forward_map(forward_map_),
        reverse_map_xadj(reverse_map_xadj_),
        reverse_map_adj(reverse_map_adj_),
        multiply_shift_for_scale(multiply_shift_for_scale_),
        division_shift_for_bucket(division_shift_for_bucket_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const CountTag &, const size_t &ii) const {
    forward_type fm = forward_map[ii];
    fm              = fm << multiply_shift_for_scale;
    fm += ii >> division_shift_for_bucket;
    typedef typename std::remove_reference<decltype(reverse_map_xadj(0))>::type atomic_incr_type;
    Kokkos::atomic_fetch_add(&(reverse_map_xadj(fm)), atomic_incr_type(1));
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const FillTag &, const size_t &ii) const {
    forward_type fm = forward_map[ii];

    fm = fm << multiply_shift_for_scale;
    fm += ii >> division_shift_for_bucket;
    typedef typename std::remove_reference<decltype(reverse_map_xadj(0))>::type atomic_incr_type;
    const reverse_type future_index = Kokkos::atomic_fetch_add(&(reverse_map_xadj(fm)), atomic_incr_type(1));
    reverse_map_adj(future_index)   = ii;
  }
};

template <typename from_view_t, typename to_view_t>
struct StridedCopy1 {
  const from_view_t from;
  to_view_t to;
  const size_t stride;
  StridedCopy1(const from_view_t from_, to_view_t to_, size_t stride_) : from(from_), to(to_), stride(stride_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &ii) const { to[ii] = from[(ii)*stride]; }
};

template <typename forward_map_type, typename reverse_map_type>
struct Reverse_Map_Functor {
  struct CountTag {};
  struct FillTag {};

  typedef struct CountTag CountTag;
  typedef struct FillTag FillTag;

  typedef typename forward_map_type::value_type forward_type;
  typedef typename reverse_map_type::value_type reverse_type;
  forward_map_type forward_map;
  reverse_map_type reverse_map_xadj;
  reverse_map_type reverse_map_adj;

  Reverse_Map_Functor(forward_map_type forward_map_, reverse_map_type reverse_map_xadj_,
                      reverse_map_type reverse_map_adj_)
      : forward_map(forward_map_), reverse_map_xadj(reverse_map_xadj_), reverse_map_adj(reverse_map_adj_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const CountTag &, const size_t &ii) const {
    forward_type fm = forward_map[ii];
    typedef typename std::remove_reference<decltype(reverse_map_xadj(0))>::type atomic_incr_type;
    Kokkos::atomic_fetch_add(&(reverse_map_xadj(fm)), atomic_incr_type(1));
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const FillTag &, const size_t &ii) const {
    forward_type c = forward_map[ii];
    typedef typename std::remove_reference<decltype(reverse_map_xadj(0))>::type atomic_incr_type;
    const reverse_type future_index = Kokkos::atomic_fetch_add(&(reverse_map_xadj(c)), atomic_incr_type(1));
    reverse_map_adj(future_index)   = ii;
  }
};

/// \brief Utility function to obtain a reverse map given a map.
/// Input is a map with the number of elements within the map.
/// forward_map[c] = i, where c is a forward element and forward_map has a size
/// of num_forward_elements. i is the value that c is mapped in the forward map,
/// and the range of that is num_reverse_elements. Output is the
/// reverse_map_xadj and reverse_map_adj such that, all c, forward_map[c] = i,
/// will appear in reverse_map_adj[ reverse_map_xadj[i]: reverse_map_xadj[i+1])

/// \param num_forward_elements the number of elements in the forward map,
///        the size of the forward map.
/// \param num_reverse_elements the number of elements that
///        forward map is mapped to. It is the value of max i.
/// \param forward_map input forward_map, where forward_map[c] = i.
/// \param reverse_map_xadj
///          reverse map xadj, that is it will hold the beginning and
///          end indices on reverse_map_adj such that all values mapped
///          to i will be [reverse_map_xadj[i]: reverse_map_xadj[i+1])
//           its size will be num_reverse_elements + 1.
///          NO NEED TO INITIALIZE.
/// \param reverse_map_adj reverse map adj, holds the values of reverse
///        maps. Its size is num_forward_elements.
template <typename forward_array_type, typename reverse_array_type,
          typename MyExecSpace>
void kk_create_reverse_map(const typename reverse_array_type::value_type &num_forward_elements,  // num_vertices
                           const typename forward_array_type::value_type &num_reverse_elements,  // num_colors

                           const forward_array_type &forward_map,        // vertex to colors
                           const reverse_array_type &reverse_map_xadj,   // colors to vertex xadj
                           const reverse_array_type &reverse_map_adj) {  // colros to vertex adj

  typedef typename reverse_array_type::value_type lno_t;
  typedef typename forward_array_type::value_type reverse_lno_t;

  const lno_t MINIMUM_TO_ATOMIC = 128;

  // typedef Kokkos::TeamPolicy<CountTag, MyExecSpace> team_count_policy_t ;
  // typedef Kokkos::TeamPolicy<FillTag, MyExecSpace> team_fill_policy_t ;

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;

  // IF There are very few reverse elements, atomics are likely to create
  // contention.
  if (num_reverse_elements < MINIMUM_TO_ATOMIC) {
    const lno_t scale_size               = 1024;
    const lno_t multiply_shift_for_scale = 10;

    // there will be 1024 buckets
    const lno_t division_shift_for_bucket = lno_t(ceil(log(double(num_forward_elements) / scale_size) / log(2)));

    // coloring indices are base-1. we end up using not using element 1.
    const reverse_lno_t tmp_reverse_size = (num_reverse_elements + 1) << multiply_shift_for_scale;

    typename reverse_array_type::non_const_type tmp_color_xadj("TMP_REVERSE_XADJ", tmp_reverse_size + 1);

    typedef Fill_Reverse_Scale_Functor<forward_array_type, reverse_array_type> frsf;
    typedef typename frsf::CountTag cnt_tag;
    typedef typename frsf::FillTag fill_tag;
    typedef Kokkos::RangePolicy<cnt_tag, MyExecSpace> my_cnt_exec_space;
    typedef Kokkos::RangePolicy<fill_tag, MyExecSpace> my_fill_exec_space;

    frsf frm(forward_map, tmp_color_xadj, reverse_map_adj, multiply_shift_for_scale, division_shift_for_bucket);

    Kokkos::parallel_for("KokkosKernels::Common::CreateReverseMap::NonAtomic::S0",
                         my_cnt_exec_space(0, num_forward_elements), frm);
    MyExecSpace().fence();

    // kk_inclusive_parallel_prefix_sum<reverse_array_type,
    // MyExecSpace>(tmp_reverse_size + 1, tmp_color_xadj);
    KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<MyExecSpace>(tmp_reverse_size + 1, tmp_color_xadj);
    MyExecSpace().fence();

    Kokkos::parallel_for(
        "KokkosKernels::Common::CreateReverseMap::NonAtomic::S1", my_exec_space(0, num_reverse_elements + 1),
        StridedCopy1<reverse_array_type, reverse_array_type>(tmp_color_xadj, reverse_map_xadj, scale_size));
    MyExecSpace().fence();
    Kokkos::parallel_for("KokkosKernels::Common::CreateReverseMap::NonAtomic::S2",
                         my_fill_exec_space(0, num_forward_elements), frm);
    MyExecSpace().fence();
  } else
  // atomic implementation.
  {
    reverse_array_type tmp_color_xadj("TMP_REVERSE_XADJ", num_reverse_elements + 1);

    typedef Reverse_Map_Functor<forward_array_type, reverse_array_type> rmp_functor_type;
    typedef typename rmp_functor_type::CountTag cnt_tag;
    typedef typename rmp_functor_type::FillTag fill_tag;
    typedef Kokkos::RangePolicy<cnt_tag, MyExecSpace> my_cnt_exec_space;
    typedef Kokkos::RangePolicy<fill_tag, MyExecSpace> my_fill_exec_space;

    rmp_functor_type frm(forward_map, tmp_color_xadj, reverse_map_adj);

    Kokkos::parallel_for("KokkosKernels::Common::CreateReverseMap::Atomic::S0",
                         my_cnt_exec_space(0, num_forward_elements), frm);
    MyExecSpace().fence();

    // kk_inclusive_parallel_prefix_sum<reverse_array_type,
    // MyExecSpace>(num_reverse_elements + 1, reverse_map_xadj);
    KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<MyExecSpace>(num_reverse_elements + 1, tmp_color_xadj);
    MyExecSpace().fence();

    Kokkos::deep_copy(reverse_map_xadj, tmp_color_xadj);
    MyExecSpace().fence();

    Kokkos::parallel_for("KokkosKernels::Common::CreateReverseMap::Atomic::S1",
                         my_fill_exec_space(0, num_forward_elements), frm);
    MyExecSpace().fence();
  }
}

template <typename in_row_view_t, typename in_nnz_view_t, typename in_color_view_t, typename team_member>
struct ColorChecker {
  typedef typename in_row_view_t::value_type size_type;
  typedef typename in_nnz_view_t::value_type lno_t;
  typedef typename in_color_view_t::value_type color_t;
  lno_t num_rows;
  in_row_view_t xadj;
  in_nnz_view_t adj;
  in_color_view_t color_view;
  lno_t team_row_chunk_size;

  ColorChecker(lno_t num_rows_, in_row_view_t xadj_, in_nnz_view_t adj_, in_color_view_t color_view_, lno_t chunk_size)
      : num_rows(num_rows_), xadj(xadj_), adj(adj_), color_view(color_view_), team_row_chunk_size(chunk_size) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member &teamMember, size_t &num_conflicts) const {
    // get the range of rows for team.
    const lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, num_rows);

    size_t nf = 0;
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
        [&](const lno_t &row_index, size_t &team_num_conf) {
          color_t my_color          = color_view(row_index);
          const size_type col_begin = xadj[row_index];
          const size_type col_end   = xadj[row_index + 1];
          const lno_t left_work     = col_end - col_begin;

          size_t conf1 = 0;
          Kokkos::parallel_reduce(
              Kokkos::ThreadVectorRange(teamMember, left_work),
              [&](lno_t i, size_t &valueToUpdate) {
                const size_type adjind = i + col_begin;
                const lno_t colIndex   = adj[adjind];
                if (colIndex != row_index) {
                  color_t second_color = color_view(colIndex);
                  if (second_color == my_color) valueToUpdate += 1;
                }
              },
              conf1);
          team_num_conf += conf1;
        },
        nf);
    num_conflicts += nf;
  }
};

/// \brief given a graph and a coloring function returns true or false if
///        distance-1 coloring is valid or not.
///
/// \param num_rows num rows in input graph
/// \param xadj     row pointers of the input graph
/// \param adj      column indices of the input graphw
/// \param v_colors The colors at each vertex in the graph.
template <typename in_row_view_t, typename in_nnz_view_t, typename in_color_view_t, typename MyExecSpace>
inline size_t kk_is_d1_coloring_valid(typename in_nnz_view_t::non_const_value_type num_rows,
                                      typename in_nnz_view_t::non_const_value_type /*num_cols*/, in_row_view_t xadj,
                                      in_nnz_view_t adj, in_color_view_t v_colors) {
  KokkosKernels::Impl::ExecSpaceType my_exec_space = KokkosKernels::Impl::kk_get_exec_space_type<MyExecSpace>();
  int vector_size         = kk_get_suggested_vector_size(num_rows, adj.extent(0), my_exec_space);
  int suggested_team_size = kk_get_suggested_team_size(vector_size, my_exec_space);
  ;
  typename in_nnz_view_t::non_const_value_type team_work_chunk_size = suggested_team_size;
  typedef Kokkos::TeamPolicy<MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic>> dynamic_team_policy;
  typedef typename dynamic_team_policy::member_type team_member_t;

  struct ColorChecker<in_row_view_t, in_nnz_view_t, in_color_view_t, team_member_t> cc(num_rows, xadj, adj, v_colors,
                                                                                       team_work_chunk_size);
  size_t num_conf = 0;
  Kokkos::parallel_reduce("KokkosKernels::Common::IsD1ColoringValid",
                          dynamic_team_policy(num_rows / team_work_chunk_size + 1, suggested_team_size, vector_size),
                          cc, num_conf);

  MyExecSpace().fence();
  return num_conf;
}

template <typename Reducer, typename ordinal_t, typename rowmap_t>
struct MinMaxDegreeFunctor {
  using ReducerVal = typename Reducer::value_type;
  MinMaxDegreeFunctor(const rowmap_t &rowmap_) : rowmap(rowmap_) {}
  KOKKOS_INLINE_FUNCTION void operator()(ordinal_t i, ReducerVal &lminmax) const {
    ordinal_t deg = rowmap(i + 1) - rowmap(i);
    if (deg < lminmax.min_val) lminmax.min_val = deg;
    if (deg > lminmax.max_val) lminmax.max_val = deg;
  }
  rowmap_t rowmap;
};

template <typename Reducer, typename ordinal_t, typename rowmap_t>
struct MaxDegreeFunctor {
  using ReducerVal = typename Reducer::value_type;
  MaxDegreeFunctor(const rowmap_t &rowmap_) : rowmap(rowmap_) {}
  KOKKOS_INLINE_FUNCTION void operator()(ordinal_t i, ReducerVal &lmax) const {
    ordinal_t deg = rowmap(i + 1) - rowmap(i);
    if (deg > lmax) lmax = deg;
  }
  rowmap_t rowmap;
};

template <typename device_t, typename ordinal_t, typename rowmap_t>
ordinal_t graph_max_degree(const rowmap_t &rowmap) {
  using Reducer   = Kokkos::Max<ordinal_t>;
  ordinal_t nrows = rowmap.extent(0);
  if (nrows) nrows--;
  if (nrows == 0) return 0;
  ordinal_t val;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<typename device_t::execution_space>(0, nrows),
                          MaxDegreeFunctor<Reducer, ordinal_t, rowmap_t>(rowmap), Reducer(val));
  return val;
}

template <typename execution_space, typename rowmap_t>
typename rowmap_t::non_const_value_type graph_max_degree(const execution_space &exec, const rowmap_t &rowmap) {
  using Offset  = typename rowmap_t::non_const_value_type;
  using Reducer = Kokkos::Max<Offset>;
  Offset nrows  = rowmap.extent(0);
  if (nrows) nrows--;
  if (nrows == 0) return 0;
  Offset val;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<execution_space>(exec, 0, nrows),
                          MaxDegreeFunctor<Reducer, Offset, rowmap_t>(rowmap), Reducer(val));
  return val;
}

template <typename device_t, typename ordinal_t, typename rowmap_t>
void graph_min_max_degree(const rowmap_t &rowmap, ordinal_t &min_degree, ordinal_t &max_degree) {
  using Reducer   = Kokkos::MinMax<ordinal_t>;
  ordinal_t nrows = rowmap.extent(0);
  if (nrows) nrows--;
  if (nrows == 0) {
    min_degree = 0;
    max_degree = 0;
    return;
  }
  typename Reducer::value_type result;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<typename device_t::execution_space>(0, nrows),
                          MinMaxDegreeFunctor<Reducer, ordinal_t, rowmap_t>(rowmap), Reducer(result));
  min_degree = result.min_val;
  max_degree = result.max_val;
}

template <typename size_type, typename lno_t, typename ExecutionSpace, typename scalar_t = double>
struct LowerTriangularMatrix {
  struct CountTag {};
  struct FillTag {};

  typedef struct CountTag CountTag;
  typedef struct FillTag FillTag;

  typedef Kokkos::TeamPolicy<CountTag, ExecutionSpace> team_count_policy_t;
  typedef Kokkos::TeamPolicy<FillTag, ExecutionSpace> team_fill_policy_t;

  typedef Kokkos::TeamPolicy<CountTag, ExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>> dynamic_team_count_policy_t;
  typedef Kokkos::TeamPolicy<FillTag, ExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>> dynamic_team_fill_policy_t;

  typedef typename team_count_policy_t::member_type team_count_member_t;
  typedef typename team_fill_policy_t::member_type team_fill_member_t;

  lno_t num_rows;
  const size_type *xadj;
  const lno_t *adj;
  const scalar_t *in_vals;
  const lno_t *permutation;

  size_type *t_xadj;  // allocated
  lno_t *t_adj;       // allocated
  scalar_t *t_vals;

  const lno_t team_work_size;
  const KokkosKernels::Impl::ExecSpaceType exec_space;
  const bool is_lower;
  const bool incl_diag;

  LowerTriangularMatrix(const lno_t num_rows_, const size_type *xadj_, const lno_t *adj_, const scalar_t *in_vals_,
                        const lno_t *permutation_, size_type *t_xadj_, lno_t *t_adj_, scalar_t *out_vals_,
                        const lno_t team_row_work_size_, bool is_lower_ = true, bool incl_diag_ = false)
      : num_rows(num_rows_),
        xadj(xadj_),
        adj(adj_),
        in_vals(in_vals_),
        permutation(permutation_),
        t_xadj(t_xadj_),
        t_adj(t_adj_),
        t_vals(out_vals_),
        team_work_size(team_row_work_size_),
        exec_space(KokkosKernels::Impl::kk_get_exec_space_type<ExecutionSpace>()),
        is_lower(is_lower_),
        incl_diag(incl_diag_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const CountTag &, const team_count_member_t &teamMember) const {
    const lno_t team_row_begin = teamMember.league_rank() * team_work_size;
    const lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, num_rows);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const lno_t &row_index) {
                           lno_t row_perm = row_index;
                           if (permutation != NULL) {
                             row_perm = permutation[row_perm];
                           }

                           const size_type col_begin = xadj[row_index];
                           const size_type col_end   = xadj[row_index + 1];
                           const lno_t left_work     = col_end - col_begin;
                           lno_t lower_row_size      = 0;
                           Kokkos::parallel_reduce(
                               Kokkos::ThreadVectorRange(teamMember, left_work),
                               [&](lno_t i, lno_t &rowsize_) {
                                 const size_type adjind = i + col_begin;
                                 lno_t colIndex         = adj[adjind];
                                 if (permutation != NULL) {
                                   colIndex = permutation[colIndex];
                                 }
                                 if ((is_lower && row_perm > colIndex) || (!is_lower && row_perm < colIndex) ||
                                     (incl_diag && row_perm == colIndex)) {
                                   rowsize_ += 1;
                                 }
                               },
                               lower_row_size);

                           t_xadj[row_index] = lower_row_size;
                         });
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const FillTag &, const team_fill_member_t &teamMember) const {
    const lno_t team_row_begin = teamMember.league_rank() * team_work_size;
    const lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, num_rows);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const lno_t &row_index) {
                           lno_t row_perm = row_index;
                           if (permutation != NULL) {
                             row_perm = permutation[row_perm];
                           }

                           const size_type col_begin  = xadj[row_index];
                           const size_type col_end    = xadj[row_index + 1];
                           const lno_t read_left_work = col_end - col_begin;

                           const size_type write_begin = t_xadj[row_index];
                           const size_type write_end   = t_xadj[row_index + 1];
                           const lno_t write_left_work = write_end - write_begin;

                           // TODO: Write GPU (vector-level) version here:
                           /*
                           if(kk_is_gpu_exec_space<ExecutionSpace>())
                           {
                             Kokkos::parallel_for(
                                 Kokkos::ThreadVectorRange(teamMember, read_left_work),
                                 [&] (lno_t i) {
                               const size_type adjind = i + col_begin;
                               const lno_t colIndex = adj[adjind];
                             });
                           }
                           else
                           ...
                           */

                           for (lno_t r = 0, w = 0; r < read_left_work && w < write_left_work; ++r) {
                             const size_type adjind = r + col_begin;
                             const lno_t colIndex   = adj[adjind];
                             lno_t colperm          = colIndex;
                             if (permutation != NULL) {
                               colperm = permutation[colIndex];
                             }
                             if ((is_lower && row_perm > colperm) || (!is_lower && row_perm < colperm) ||
                                 (incl_diag && row_perm == colperm)) {
                               if (in_vals != NULL) {
                                 t_vals[write_begin + w] = in_vals[adjind];
                               }
                               t_adj[write_begin + w++] = colIndex;
                             }
                           }
                         });
  }
};
template <typename size_type, typename lno_t, typename ExecutionSpace>
void kk_get_lower_triangle_count_parallel(const lno_t nv, const size_type ne, const size_type *in_xadj,
                                          const lno_t *in_adj, size_type *out_xadj, const lno_t *new_indices = NULL,
                                          bool use_dynamic_scheduling = false, int chunksize = 4, bool is_lower = true,
                                          bool incl_diag = false) {
  const int vector_size =
      kk_get_suggested_vector_size(nv, ne, KokkosKernels::Impl::kk_get_exec_space_type<ExecutionSpace>());
  const int suggested_team_size =
      kk_get_suggested_team_size(vector_size, KokkosKernels::Impl::kk_get_exec_space_type<ExecutionSpace>());
  const int team_work_chunk_size = suggested_team_size * chunksize;
  typedef LowerTriangularMatrix<size_type, lno_t, ExecutionSpace> ltm_t;

  ltm_t ltm(nv, in_xadj, in_adj, NULL, new_indices, out_xadj, NULL, NULL, team_work_chunk_size, is_lower, incl_diag);

  typedef typename ltm_t::team_count_policy_t count_tp_t;
  typedef typename ltm_t::dynamic_team_count_policy_t d_count_tp_t;

  if (use_dynamic_scheduling) {
    Kokkos::parallel_for("KokkosKernels::Common::GetLowerTriangleCount::DynamicSchedule",
                         d_count_tp_t(nv / team_work_chunk_size + 1, suggested_team_size, vector_size), ltm);
  } else {
    Kokkos::parallel_for("KokkosKernels::Common::GetLowerTriangleCount::StaticSchedule",
                         count_tp_t(nv / team_work_chunk_size + 1, suggested_team_size, vector_size), ltm);
  }
  ExecutionSpace().fence();
}

template <typename size_type, typename lno_t>
void kk_sort_by_row_size_sequential(const lno_t nv, const size_type *in_xadj, lno_t *new_indices,
                                    int sort_decreasing_order = 1) {
  std::vector<lno_t> begins(nv);
  std::vector<lno_t> nexts(nv);
  for (lno_t i = 0; i < nv; ++i) {
    nexts[i] = begins[i] = -1;
  }

  for (lno_t i = 0; i < nv; ++i) {
    lno_t row_size   = in_xadj[i + 1] - in_xadj[i];
    nexts[i]         = begins[row_size];
    begins[row_size] = i;
  }
  if (sort_decreasing_order == 1) {
    lno_t new_index     = nv;
    const lno_t row_end = -1;
    for (lno_t i = 0; i < nv; ++i) {
      lno_t row = begins[i];
      while (row != row_end) {
        new_indices[row] = --new_index;
        row              = nexts[row];
      }
    }
  } else if (sort_decreasing_order == 2) {
    lno_t new_index_top    = nv;
    lno_t new_index_bottom = 0;
    const lno_t row_end    = -1;
    bool is_even           = true;
    for (lno_t i = nv - 1;; --i) {
      lno_t row = begins[i];
      while (row != row_end) {
        if (is_even) {
          new_indices[row] = --new_index_top;
        } else {
          new_indices[row] = new_index_bottom++;
        }
        is_even = !is_even;
        row     = nexts[row];
      }
      if (i == 0) break;
    }
  } else {
    lno_t new_index     = 0;
    const lno_t row_end = -1;
    for (lno_t i = 0; i < nv; ++i) {
      lno_t row = begins[i];
      while (row != row_end) {
        new_indices[row] = new_index++;
        row              = nexts[row];
      }
    }
  }
}
#ifdef KOKKOSKERNELS_HAVE_PARALLEL_GNUSORT
template <typename size_type, typename lno_t, typename ExecutionSpace>
void kk_sort_by_row_size_parallel(const lno_t nv, const size_type *in_xadj, lno_t *new_indices,
                                  int sort_decreasing_order = 1, int num_threads = 1) {
  typedef Kokkos::RangePolicy<ExecutionSpace> my_exec_space;

  struct SortItem {
    lno_t id;
    lno_t size;
    bool operator<(const SortItem &a) const { return this->size > a.size; }
  };

  std::vector<SortItem> vnum_elements(nv);
  SortItem *num_elements = &(vnum_elements[0]);

  Kokkos::parallel_for(
      "KokkosKernels::Common::SortByRowSize::S0", my_exec_space(0, nv), KOKKOS_LAMBDA(const lno_t &row) {
        lno_t row_size         = in_xadj[row + 1] - in_xadj[row];
        num_elements[row].size = row_size;
        num_elements[row].id   = row;
      });
  __gnu_parallel::sort(&(num_elements[0]), &(num_elements[0]) + nv, std::less<struct SortItem>());

  if (sort_decreasing_order == 1) {
    Kokkos::parallel_for(
        "KokkosKernels::Common::SortByRowSize::S1", my_exec_space(0, nv),
        KOKKOS_LAMBDA(const lno_t &row) { new_indices[num_elements[row].id] = row; });
  } else if (sort_decreasing_order == 0) {
    Kokkos::parallel_for(
        "KokkosKernels::Common::SortByRowSize::S2", my_exec_space(0, nv),
        KOKKOS_LAMBDA(const lno_t &row) { new_indices[num_elements[row].id] = nv - row - 1; });
  } else {
    Kokkos::parallel_for(
        "KokkosKernels::Common::SortByRowSize::S3", my_exec_space(0, nv), KOKKOS_LAMBDA(const lno_t &row) {
          if (row & 1) {
            new_indices[num_elements[row].id] = nv - (row + 1) / 2;
          } else {
            new_indices[num_elements[row].id] = row / 2;
          }
        });
  }
}
#endif

#ifdef KOKKOSKERNELS_HAVE_PARALLEL_GNUSORT
template <typename size_type, typename lno_t, typename ExecutionSpace>
void kk_sort_by_row_size(const lno_t nv, const size_type *in_xadj, lno_t *new_indices, int sort_decreasing_order = 1,
                         int num_threads = 64) {
  std::cout << "Parallel Sort" << std::endl;
  kk_sort_by_row_size_parallel<size_type, lno_t, ExecutionSpace>(nv, in_xadj, new_indices, sort_decreasing_order,
                                                                 num_threads);
}
#else
template <typename size_type, typename lno_t, typename ExecutionSpace>
void kk_sort_by_row_size(const lno_t nv, const size_type *in_xadj, lno_t *new_indices, int sort_decreasing_order = 1,
                         int /*num_threads*/ = 64) {
  std::cout << "Sequential Sort" << std::endl;
  kk_sort_by_row_size_sequential(nv, in_xadj, new_indices, sort_decreasing_order);
}
#endif

template <typename size_type, typename lno_t, typename ExecutionSpace, typename scalar_t = double>
void kk_get_lower_triangle_fill_parallel(const lno_t nv, const size_type ne, const size_type *in_xadj,
                                         const lno_t *in_adj, const scalar_t *in_vals, size_type *out_xadj,
                                         lno_t *out_adj, scalar_t *out_vals, const lno_t *new_indices = NULL,
                                         bool use_dynamic_scheduling = false, bool chunksize = 4, bool is_lower = true,
                                         bool incl_diag = false) {
  const int vector_size =
      kk_get_suggested_vector_size(nv, ne, KokkosKernels::Impl::kk_get_exec_space_type<ExecutionSpace>());
  const int suggested_team_size =
      kk_get_suggested_team_size(vector_size, KokkosKernels::Impl::kk_get_exec_space_type<ExecutionSpace>());
  const int team_work_chunk_size = suggested_team_size * chunksize;

  typedef LowerTriangularMatrix<size_type, lno_t, ExecutionSpace, scalar_t> ltm_t;
  ltm_t ltm(nv, in_xadj, in_adj, in_vals, new_indices, out_xadj, out_adj, out_vals, team_work_chunk_size, is_lower,
            incl_diag);

  typedef typename ltm_t::team_fill_policy_t fill_p_t;
  typedef typename ltm_t::dynamic_team_fill_policy_t d_fill_p_t;

  if (use_dynamic_scheduling) {
    Kokkos::parallel_for("KokkosKernels::Common::GetLowerTriangleFill::DynamicSchedule",
                         d_fill_p_t(nv / team_work_chunk_size + 1, suggested_team_size, vector_size), ltm);
  } else {
    Kokkos::parallel_for("KokkosKernels::Common::GetLowerTriangleFill::StaticSchedule",
                         fill_p_t(nv / team_work_chunk_size + 1, suggested_team_size, vector_size), ltm);
  }
  ExecutionSpace().fence();
}

template <typename size_type, typename lno_t, typename ExecutionSpace>
void kk_get_lower_triangle_count(const lno_t nv, const size_type ne, const size_type *in_xadj, const lno_t *in_adj,
                                 size_type *out_xadj, const lno_t *new_indices = NULL,
                                 bool use_dynamic_scheduling = false, bool chunksize = 4, bool is_lower = true,
                                 bool incl_diag = false) {
  // Kokkos::Timer timer1;

  kk_get_lower_triangle_count_parallel<size_type, lno_t, ExecutionSpace>(
      nv, ne, in_xadj, in_adj, out_xadj, new_indices, use_dynamic_scheduling, chunksize, is_lower, incl_diag);
  // double count = timer1.seconds();
  // std::cout << "lower count time:" << count<< std::endl;
}
template <typename size_type, typename lno_t, typename scalar_t, typename ExecutionSpace>
void kk_get_lower_triangle_fill(lno_t nv, size_type ne, const size_type *in_xadj, const lno_t *in_adj,
                                const scalar_t *in_vals, size_type *out_xadj, lno_t *out_adj, scalar_t *out_vals,
                                const lno_t *new_indices = NULL, bool use_dynamic_scheduling = false,
                                bool chunksize = 4, bool is_lower = true, bool incl_diag = false) {
  // Kokkos::Timer timer1;

  kk_get_lower_triangle_fill_parallel<size_type, lno_t, ExecutionSpace, scalar_t>(
      nv, ne, in_xadj, in_adj, in_vals, out_xadj, out_adj, out_vals, new_indices, use_dynamic_scheduling, chunksize,
      is_lower, incl_diag);

  // double fill = timer1.seconds();
  // std::cout << "lower fill time:" << fill<< std::endl;
}

template <typename crstmat_t>
crstmat_t kk_get_lower_triangle(crstmat_t in_crs_matrix, typename crstmat_t::index_type::value_type *new_indices = NULL,
                                bool use_dynamic_scheduling = false, bool chunksize = 4, bool is_lower = true,
                                bool incl_diag = false) {
  typedef typename crstmat_t::execution_space exec_space;
  typedef typename crstmat_t::StaticCrsGraphType graph_t;
  typedef typename crstmat_t::row_map_type::non_const_type row_map_view_t;
  typedef typename crstmat_t::index_type::non_const_type cols_view_t;
  typedef typename crstmat_t::values_type::non_const_type values_view_t;
  // typedef typename crstmat_t::row_map_type::const_type const_row_map_view_t;
  // typedef typename crstmat_t::index_type::const_type   const_cols_view_t;
  // typedef typename crstmat_t::values_type::const_type const_values_view_t;

  typedef typename row_map_view_t::non_const_value_type size_type;
  typedef typename cols_view_t::non_const_value_type lno_t;
  typedef typename values_view_t::non_const_value_type scalar_t;

  lno_t nr = in_crs_matrix.numRows();

  const scalar_t *vals    = in_crs_matrix.values.data();
  const size_type *rowmap = in_crs_matrix.graph.row_map.data();
  const lno_t *entries    = in_crs_matrix.graph.entries.data();
  const size_type ne      = in_crs_matrix.graph.entries.extent(0);

  row_map_view_t new_row_map(Kokkos::view_alloc(Kokkos::WithoutInitializing, "LL"), nr + 1);
  kk_get_lower_triangle_count<size_type, lno_t, exec_space>(nr, ne, rowmap, entries, new_row_map.data(), new_indices,
                                                            use_dynamic_scheduling, chunksize, is_lower, incl_diag);

  KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<exec_space>(nr + 1, new_row_map);
  exec_space().fence();

  auto ll_size   = Kokkos::subview(new_row_map, nr);
  auto h_ll_size = Kokkos::create_mirror_view(ll_size);
  Kokkos::deep_copy(h_ll_size, ll_size);
  size_type ll_nnz_size = h_ll_size();

  // cols_view_t new_entries ("LL", ll_nnz_size);
  // values_view_t new_values ("LL", ll_nnz_size);
  cols_view_t new_entries(Kokkos::view_alloc(Kokkos::WithoutInitializing, "LL"), ll_nnz_size);
  values_view_t new_values(Kokkos::view_alloc(Kokkos::WithoutInitializing, "LL"), ll_nnz_size);

  kk_get_lower_triangle_fill<size_type, lno_t, scalar_t, exec_space>(
      nr, ne, rowmap, entries, vals, new_row_map.data(), new_entries.data(), new_values.data(), new_indices,
      use_dynamic_scheduling, chunksize, is_lower, incl_diag);

  graph_t g(new_entries, new_row_map);
  crstmat_t new_ll_mtx("lower triangle", in_crs_matrix.numCols(), new_values, g);
  return new_ll_mtx;
}

template <typename row_map_view_t, typename cols_view_t, typename values_view_t, typename out_row_map_view_t,
          typename out_cols_view_t, typename out_values_view_t, typename new_indices_t, typename exec_space>
void kk_get_lower_triangle(typename cols_view_t::non_const_value_type nr, row_map_view_t in_rowmap,
                           cols_view_t in_entries, values_view_t in_values, out_row_map_view_t &out_rowmap,
                           out_cols_view_t &out_entries, out_values_view_t &out_values, new_indices_t &new_indices,
                           bool use_dynamic_scheduling = false, bool chunksize = 4, bool is_lower = true,
                           bool incl_diag = false) {
  // typedef typename row_map_view_t::const_type const_row_map_view_t;
  // typedef typename cols_view_t::const_type   const_cols_view_t;
  // typedef typename values_view_t::const_type const_values_view_t;

  typedef typename row_map_view_t::non_const_value_type size_type;
  typedef typename cols_view_t::non_const_value_type lno_t;
  typedef typename values_view_t::non_const_value_type scalar_t;

  const scalar_t *vals    = in_values.data();
  const size_type *rowmap = in_rowmap.data();
  const lno_t *entries    = in_entries.data();
  const size_type ne      = in_entries.extent(0);

  out_rowmap = out_row_map_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "LL"), nr + 1);
  kk_get_lower_triangle_count<size_type, lno_t, exec_space>(nr, ne, rowmap, entries, out_rowmap.data(),
                                                            new_indices.data(), use_dynamic_scheduling, chunksize,
                                                            is_lower, incl_diag);

  KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<exec_space>(nr + 1, out_rowmap);
  exec_space().fence();

  auto ll_size   = Kokkos::subview(out_rowmap, nr);
  auto h_ll_size = Kokkos::create_mirror_view(ll_size);
  Kokkos::deep_copy(h_ll_size, ll_size);
  size_type ll_nnz_size = h_ll_size();

  // cols_view_t new_entries ("LL", ll_nnz_size);
  // values_view_t new_values ("LL", ll_nnz_size);
  out_entries = out_cols_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "LL"), ll_nnz_size);

  if (in_values.data() != NULL)
    out_values = out_values_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "LL"), ll_nnz_size);

  kk_get_lower_triangle_fill<size_type, lno_t, scalar_t, exec_space>(
      nr, ne, rowmap, entries, vals, out_rowmap.data(), out_entries.data(), out_values.data(), new_indices.data(),
      use_dynamic_scheduling, chunksize, is_lower, incl_diag);
}

template <typename row_map_view_t, typename cols_view_t, typename out_row_map_view_t, typename out_cols_view_t,
          typename exec_space>
void kk_create_incidence_tranpose_matrix_from_lower_triangle(typename cols_view_t::non_const_value_type nr,
                                                             row_map_view_t in_rowmap, cols_view_t in_entries,
                                                             out_row_map_view_t &out_rowmap,
                                                             out_cols_view_t &out_entries,
                                                             bool /*use_dynamic_scheduling */ = false,
                                                             bool /*chunksize*/               = 4) {
  // typedef typename row_map_view_t::const_type const_row_map_view_t;
  // typedef typename cols_view_t::const_type   const_cols_view_t;

  typedef typename row_map_view_t::non_const_value_type size_type;
  typedef typename cols_view_t::non_const_value_type lno_t;

  // const size_type *rowmap = in_rowmap.data();
  // const lno_t *entries= in_entries.data();
  const size_type ne = in_entries.extent(0);
  out_rowmap         = out_row_map_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "LL"), ne + 1);
  // const lno_t nr = in_rowmap.extent(0) - 1;
  typedef Kokkos::RangePolicy<exec_space> my_exec_space;

  Kokkos::parallel_for(
      "KokkosKernels::Common::CreateIncidenceTransposeMatrixFromLowerTriangle::"
      "S0",
      my_exec_space(0, ne + 1), KOKKOS_LAMBDA(const lno_t &i) { out_rowmap[i] = i * 2; });

  // typedef Kokkos::TeamPolicy<exec_space> team_policy_t;
  // int vector_size = 2;
  // team_policy_t(ne)
  // nv  / team_work_chunk_size + 1 , suggested_team_size, vector_size

  out_entries = out_cols_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "LL"), 2 * ne);

  // TODO MAKE IT WITH TEAMS.
  Kokkos::parallel_for(
      "KokkosKernels::Common::CreateIncidenceTransposeMatrixFromLowerTriangle::"
      "S1",
      my_exec_space(0, nr), KOKKOS_LAMBDA(const size_type &row) {
        size_type begin = in_rowmap(row);
        lno_t row_size  = in_rowmap(row + 1) - begin;
        for (int i = 0; i < row_size; ++i) {
          size_type edge_ind        = i + begin;
          lno_t col                 = in_entries(edge_ind);
          edge_ind                  = edge_ind * 2;
          out_entries[edge_ind]     = row;
          out_entries[edge_ind + 1] = col;
        }
      });
}

template <typename row_map_view_t, typename cols_view_t, typename out_row_map_view_t, typename out_cols_view_t,
          typename permutation_view_t, typename exec_space>
void kk_create_incidence_matrix_from_original_matrix(typename cols_view_t::non_const_value_type nr,
                                                     row_map_view_t in_rowmap, cols_view_t in_entries,
                                                     out_row_map_view_t &out_rowmap, out_cols_view_t &out_entries,
                                                     permutation_view_t permutation,
                                                     bool use_dynamic_scheduling = false, bool chunksize = 4) {
  // typedef typename row_map_view_t::const_type const_row_map_view_t;
  // typedef typename cols_view_t::const_type   const_cols_view_t;

  typedef typename row_map_view_t::non_const_value_type size_type;
  typedef typename cols_view_t::non_const_value_type lno_t;
  typedef Kokkos::RangePolicy<exec_space> my_exec_space;
  lno_t *perm        = permutation.data();
  const size_type ne = in_entries.extent(0);

  out_rowmap  = out_row_map_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "out_rowmap"), nr + 1);
  out_entries = out_cols_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "out_cols_view"), ne);

  // todo: need to try both true and false
  bool sort_decreasing_order = true;
  // find the size of rows at upper triangular.
  // this gives the size of each column in lower triangluar.
  kk_get_lower_triangle_count<size_type, lno_t, exec_space>(nr, ne, in_rowmap.data(), in_entries.data(),
                                                            out_rowmap.data(), permutation.data(),
                                                            use_dynamic_scheduling, chunksize, sort_decreasing_order);
  exec_space().fence();
  KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<exec_space>(nr + 1, out_rowmap);

  // kk_print_1Dview(out_rowmap, false, 20);

  out_row_map_view_t out_rowmap_copy(Kokkos::view_alloc(Kokkos::WithoutInitializing, "tmp"), nr + 1);
  // out_rowmap = out_row_map_view_t("LL", nr+1);
  Kokkos::parallel_for(
      "KokkosKernels::Common::"
      "CreateIncidenceTransposeMatrixFromOriginalTriangle::S0",
      my_exec_space(0, nr + 1), KOKKOS_LAMBDA(const lno_t &i) { out_rowmap_copy[i] = in_rowmap[i]; });

  if (sort_decreasing_order) {
    Kokkos::parallel_for(
        "KokkosKernels::Common::"
        "CreateIncidenceTransposeMatrixFromOriginalTriangle::S1",
        my_exec_space(0, nr), KOKKOS_LAMBDA(const size_type &row) {
          size_type begin = in_rowmap(row);
          lno_t row_size  = in_rowmap(row + 1) - begin;

          lno_t row_perm = row;
          if (perm) row_perm = perm[row];
          // std::cout << "row:" << row << " rowperm:" << row_perm << std::endl;
          size_type used_edge_index = out_rowmap[row];
          lno_t used_count          = 0;
          for (int i = 0; i < row_size; ++i) {
            size_type edge_ind = i + begin;
            lno_t col          = in_entries[edge_ind];

            lno_t col_perm = col;
            if (perm) col_perm = perm[col];
            if (row_perm > col_perm) {
              typedef typename std::remove_reference<decltype(out_rowmap_copy[0])>::type atomic_incr_type;
              size_type row_write_index    = Kokkos::atomic_fetch_add(&(out_rowmap_copy[row]), atomic_incr_type(1));
              size_type col_write_index    = Kokkos::atomic_fetch_add(&(out_rowmap_copy[col]), atomic_incr_type(1));
              out_entries[row_write_index] = used_edge_index + used_count;
              out_entries[col_write_index] = used_edge_index + used_count;
              ++used_count;
            }
          }
        });

  } else {
    Kokkos::parallel_for(
        "KokkosKernels::Common::"
        "CreateIncidenceTransposeMatrixFromOriginalTriangle::S2",
        my_exec_space(0, nr), KOKKOS_LAMBDA(const size_type &row) {
          size_type begin = in_rowmap(row);
          lno_t row_size  = in_rowmap(row + 1) - begin;

          lno_t row_perm = row;
          if (perm) row_perm = perm[row];
          // std::cout << "row:" << row << " rowperm:" << row_perm << std::endl;
          size_type used_edge_index = out_rowmap[row];
          lno_t used_count          = 0;
          for (int i = 0; i < row_size; ++i) {
            size_type edge_ind = i + begin;
            lno_t col          = in_entries[edge_ind];

            lno_t col_perm = col;
            if (perm) col_perm = perm[col];
            if (row_perm < col_perm) {
              typedef typename std::remove_reference<decltype(out_rowmap_copy[0])>::type atomic_incr_type;
              size_type row_write_index    = Kokkos::atomic_fetch_add(&(out_rowmap_copy[row]), atomic_incr_type(1));
              size_type col_write_index    = Kokkos::atomic_fetch_add(&(out_rowmap_copy[col]), atomic_incr_type(1));
              out_entries[row_write_index] = used_edge_index + used_count;
              out_entries[col_write_index] = used_edge_index + used_count;
              ++used_count;
            }
          }
        });
  }

  // out_rowmap = out_row_map_view_t("LL", nr+1);
  Kokkos::parallel_for(
      "KokkosKernels::Common::"
      "CreateIncidenceTransposeMatrixFromOriginalTriangle::S3",
      my_exec_space(0, nr + 1), KOKKOS_LAMBDA(const lno_t &i) { out_rowmap[i] = in_rowmap[i]; });
}

template <typename view_type>
struct ReduceLargerRowCount {
  view_type rowmap;
  typename view_type::const_value_type threshold;

  ReduceLargerRowCount(view_type view_to_reduce_, typename view_type::const_value_type threshold_)
      : rowmap(view_to_reduce_), threshold(threshold_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i, typename view_type::non_const_value_type &sum_reduction) const {
    if (rowmap(i + 1) - rowmap(i) > threshold) {
      sum_reduction += 1;
    }
  }
};

template <typename view_type, typename MyExecSpace>
void kk_reduce_numrows_larger_than_threshold(const MyExecSpace &my_exec_space, size_t num_elements,
                                             view_type view_to_reduce, typename view_type::const_value_type threshold,
                                             typename view_type::non_const_value_type &sum_reduction) {
  typedef Kokkos::RangePolicy<MyExecSpace> range_policy_t;
  Kokkos::parallel_reduce("KokkosKernels::Common::ReduceNumRowsLargerThanThreshold",
                          range_policy_t(my_exec_space, 0, num_elements),
                          ReduceLargerRowCount<view_type>(view_to_reduce, threshold), sum_reduction);
}

template <typename view_type, typename MyExecSpace>
void kk_reduce_numrows_larger_than_threshold(size_t num_elements, view_type view_to_reduce,
                                             typename view_type::const_value_type threshold,
                                             typename view_type::non_const_value_type &sum_reduction) {
  MyExecSpace my_exec_space;
  kk_reduce_numrows_larger_than_threshold(my_exec_space, num_elements, view_to_reduce, threshold, sum_reduction);
}

// Note: "block" in member name means it's block internal - otherwise it
// addresses sparse rows/columns (whole blocks) within whole matrix.
template <typename lno_t, typename size_type>
class RowIndexBase {
 public:
  KOKKOS_INLINE_FUNCTION
  RowIndexBase(const lno_t block_size_, const lno_t row_begin_, const lno_t row_end_)
      : block_size(block_size_), row_begin(row_begin_), row_end(row_end_), row_size(row_end_ - row_begin_) {
    row_off = row_begin_ * block_mtx_size();
  }

  KOKKOS_INLINE_FUNCTION
  lno_t begin() { return row_begin; }

  KOKKOS_INLINE_FUNCTION
  lno_t end() { return row_end; }

  KOKKOS_INLINE_FUNCTION
  lno_t size() { return row_size; }

  KOKKOS_INLINE_FUNCTION
  size_type block_mtx_size() { return static_cast<size_type>(block_size) * static_cast<size_type>(block_size); }

  KOKKOS_INLINE_FUNCTION
  size_type row_offset() { return row_off; }

 protected:
  lno_t block_size;  // = num_block_cols = num_block_rows
  lno_t row_begin;
  lno_t row_end;

 private:  // cache
  size_type row_off;
  lno_t row_size;
};

/* The only use of this is in Sparse Gauss Seidel, which is only implemented
   for BSR and CRS, which are identical when block size is 1

*/
template <SparseMatrixFormat /* format */, typename lno_t, typename size_type>
class MatrixRowIndex;

/* CWP August 11 2022
   This is pretty much the old BlockCRS one, but
   Should be able to create a specialized version of this for CRS because block
   size is 1
*/
template <typename lno_t, typename size_type>
class MatrixRowIndex<SparseMatrixFormat::CRS, lno_t, size_type> : public RowIndexBase<lno_t, size_type> {
 public:
  using Base = RowIndexBase<lno_t, size_type>;

  KOKKOS_INLINE_FUNCTION
  MatrixRowIndex(const lno_t block_size_, const lno_t row_begin_, const lno_t row_end_)
      : Base(block_size_, row_begin_, row_end_) {}

  KOKKOS_INLINE_FUNCTION
  size_type block(const lno_t col_idx) { return Base::row_offset() + col_idx * Base::block_size; }

  KOKKOS_INLINE_FUNCTION
  size_type block_stride() { return Base::size() * Base::block_size; }

  KOKKOS_INLINE_FUNCTION
  size_type value(const lno_t col_idx, const lno_t block_row, const lno_t block_col) {
    return block(col_idx) + block_row * block_stride() + block_col;
  }
};

template <typename lno_t, typename size_type>
class MatrixRowIndex<SparseMatrixFormat::BSR, lno_t, size_type> : public RowIndexBase<lno_t, size_type> {
 public:
  using Base = RowIndexBase<lno_t, size_type>;

  KOKKOS_INLINE_FUNCTION
  MatrixRowIndex(const lno_t block_size_, const lno_t row_begin_, const lno_t row_end_)
      : Base(block_size_, row_begin_, row_end_) {}

  KOKKOS_INLINE_FUNCTION
  size_type block(const lno_t col_idx) { return Base::row_offset() + col_idx * Base::block_mtx_size(); }

  KOKKOS_INLINE_FUNCTION
  size_type block_stride() { return Base::block_size; }

  KOKKOS_INLINE_FUNCTION
  size_type value(const lno_t col_idx, const lno_t block_row, const lno_t block_col) {
    return block(col_idx) + block_row * block_stride() + block_col;
  }
};

template <typename mtx_t>
struct MatrixTraits;

template <typename scalar_t, typename lno_t, typename device, typename mem_traits, typename size_type>
struct MatrixTraits<KokkosSparse::CrsMatrix<scalar_t, lno_t, device, mem_traits, size_type>> {
  static constexpr auto format = KokkosSparse::SparseMatrixFormat::CRS;
};

template <typename scalar_t, typename lno_t, typename device, typename mem_traits, typename size_type>
struct MatrixTraits<KokkosSparse::Experimental::BsrMatrix<scalar_t, lno_t, device, mem_traits, size_type>> {
  static constexpr auto format = KokkosSparse::SparseMatrixFormat::BSR;
};

template <SparseMatrixFormat /* outFormat */>
struct MatrixConverter;

template <>
struct MatrixConverter<SparseMatrixFormat::BSR> {
  template <typename scalar_t, typename lno_t, typename size_type, typename device,
            typename bsrMtx_t = KokkosSparse::Experimental::BsrMatrix<scalar_t, lno_t, device, void, size_type>>
  static bsrMtx_t from_bsr_formated_point_crsmatrix(
      const KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> &mtx, lno_t block_size) {
    return bsrMtx_t(mtx, block_size);
  }
};

template <typename Rowmap, typename Entries>
struct CountEntriesFallingEdges {
  using size_type = typename Rowmap::non_const_value_type;

  CountEntriesFallingEdges(const Entries &entries_) : entries(entries_) {}

  KOKKOS_INLINE_FUNCTION void operator()(size_type i, size_type &numFallingEdges) const {
    if (entries(i) > entries(i + 1)) numFallingEdges++;
  }

  Entries entries;
};

template <typename Rowmap, typename Entries>
struct CountRowBoundaryFallingEdges {
  using size_type    = typename Rowmap::non_const_value_type;
  using ordinal_type = typename Entries::non_const_value_type;

  CountRowBoundaryFallingEdges(const Rowmap &rowmap_, const Entries &entries_) : rowmap(rowmap_), entries(entries_) {}

  KOKKOS_INLINE_FUNCTION void operator()(ordinal_type i, size_type &numBoundaryFallingEdges) const {
    // Comparing the entries at end of row i, and beginning of row i+1
    size_type rowBegin = rowmap(i);
    size_type rowEnd   = rowmap(i + 1);
    // But skip this row if empty (meaning there is no last entry), because
    // there would be double-counting.
    if (rowBegin == rowEnd) return;
    // rowEnd is also the beginning of the next (nonempty) row.
    // But if it points the end of all entries, skip this row because it's the
    // last nonempty row.
    if (rowEnd == size_type(entries.extent(0))) return;
    if (entries(rowEnd - 1) > entries(rowEnd)) numBoundaryFallingEdges++;
  }

  Rowmap rowmap;
  Entries entries;
};

// Efficient test for whether a StaticCrsGraph has sorted rows
//(parallel and not affected by imbalanced rows).
// Unmerged/repeated entries in a row are still considered sorted.
template <typename Rowmap, typename Entries>
bool isCrsGraphSorted(const Rowmap &rowmap, const Entries &entries) {
  using size_type    = typename Rowmap::non_const_value_type;
  using ordinal_type = typename Entries::non_const_value_type;
  using exec_space   = typename Entries::execution_space;
  size_type nnz      = entries.extent(0);
  // Catch case of zero rows, and zero-length rowmap
  if (rowmap.extent(0) == size_t(0)) return true;
  // Eliminate cases with zero rows/cols/entries.
  // This also eliminates cases where row_map is extent 0.
  if (nnz == size_type(0)) return true;
  ordinal_type numRows = rowmap.extent(0) - 1;
  // A "falling edge" is where entry i is greater than entry i+1.
  // A graph is unsorted if and only if there is a falling edge where i, i+1 are
  // in the same row. So count the total falling edges, and then subtract the
  // falling edges which cross row boundaries.
  size_type totalFallingEdges = 0;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<exec_space>(0, nnz - 1),
                          CountEntriesFallingEdges<Rowmap, Entries>(entries), totalFallingEdges);
  size_type rowBoundaryFallingEdges = 0;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<exec_space>(0, numRows - 1),
                          CountRowBoundaryFallingEdges<Rowmap, Entries>(rowmap, entries), rowBoundaryFallingEdges);
  return totalFallingEdges == rowBoundaryFallingEdges;
}

template <typename Values, typename Mag, typename Offset>
struct CountDroppedEntriesFunctor {
  using Scalar = typename Values::non_const_value_type;
  CountDroppedEntriesFunctor(const Values &values_, Mag tol_) : values(values_), tol(tol_) {}

  KOKKOS_INLINE_FUNCTION void operator()(int64_t i, Offset &lcount) const {
    if (Kokkos::ArithTraits<Scalar>::abs(values(i)) <= tol) lcount++;
  }

  Values values;
  Mag tol;
};

template <typename Bitset, typename Rowmap>
struct MarkFinalRowEntries {
  MarkFinalRowEntries(const Bitset &rowEndMarkers_, const Rowmap &rowmap_)
      : rowEndMarkers(rowEndMarkers_), rowmap(rowmap_) {}

  KOKKOS_INLINE_FUNCTION void operator()(int64_t i) const {
    auto index = rowmap(i);
    if (index) rowEndMarkers.set(index - 1);
  }

  Bitset rowEndMarkers;
  Rowmap rowmap;
};

template <typename Offset>
struct DropEntriesScanner {
  KOKKOS_DEFAULTED_FUNCTION DropEntriesScanner() = default;
  KOKKOS_INLINE_FUNCTION DropEntriesScanner(Offset i_out_, Offset row_) : i_out(i_out_), row(row_) {}

  KOKKOS_INLINE_FUNCTION void operator+=(const DropEntriesScanner<Offset> &rhs) {
    i_out += rhs.i_out;
    row += rhs.row;
  }

  Offset i_out;  // The index to write in output entries/values
  Offset row;    // The row index (ignoring rows which were empty in input)
};

template <typename Bitset, typename RowmapIn, typename EntriesIn, typename ValuesIn, typename RowmapOut,
          typename EntriesOut, typename ValuesOut, typename Mag>
struct DropEntriesFunctor {
  using Offset = typename RowmapIn::non_const_value_type;
  using Scalar = typename ValuesIn::non_const_value_type;

  DropEntriesFunctor(const Bitset &rowEndMarkers_, const RowmapIn &rowmapIn_, const EntriesIn &entriesIn_,
                     const ValuesIn &valuesIn_, const RowmapOut &compactRowmapOut_, const EntriesOut &entriesOut_,
                     const ValuesOut &valuesOut_, Mag tol_)
      : rowEndMarkers(rowEndMarkers_),
        rowmapIn(rowmapIn_),
        entriesIn(entriesIn_),
        valuesIn(valuesIn_),
        compactRowmapOut(compactRowmapOut_),
        entriesOut(entriesOut_),
        valuesOut(valuesOut_),
        tol(tol_) {}

  KOKKOS_INLINE_FUNCTION void operator()(int64_t i_in, DropEntriesScanner<Offset> &scanval, bool finalPass) const {
    // i_in is the index of the input entry being processed
    // i_out (if finalPass == true) is the index of where that same entry goes
    // in the filtered matrix
    bool filter   = Kokkos::ArithTraits<Scalar>::abs(valuesIn(i_in)) <= tol;
    bool isRowEnd = rowEndMarkers.test(i_in);
    if (finalPass) {
      if (!filter) {
        // Keeping this entry, so copy it to the output.
        entriesOut(scanval.i_out) = entriesIn(i_in);
        valuesOut(scanval.i_out)  = valuesIn(i_in);
      }
      if (isRowEnd) {
        // Entry i_in was the last in its row of the input matrix.
        // We now know where that filtered row ends, so mark it in
        // compactRowmapOut.
        compactRowmapOut(scanval.row + 1) = scanval.i_out + (filter ? 0 : 1);
      }
      // Also, make one thread responsible for initializing first compact rowmap
      // entry
      if (i_in == 0) compactRowmapOut(0) = 0;
    }
    if (!filter) scanval.i_out++;
    if (isRowEnd) scanval.row++;
  }

  Bitset rowEndMarkers;
  RowmapIn rowmapIn;
  EntriesIn entriesIn;
  ValuesIn valuesIn;
  RowmapOut compactRowmapOut;
  EntriesOut entriesOut;
  ValuesOut valuesOut;
  Mag tol;
};

template <typename RowmapIn, typename RowmapOut, typename Ordinal>
struct ExpandRowmapFunctor {
  using Offset = typename RowmapIn::non_const_value_type;

  ExpandRowmapFunctor(const RowmapIn &rowmapIn_, const RowmapOut &compactRowmapOut_, const RowmapOut &rowmapOut_)
      : rowmapIn(rowmapIn_), compactRowmapOut(compactRowmapOut_), rowmapOut(rowmapOut_) {}

  KOKKOS_INLINE_FUNCTION void operator()(Ordinal row, Ordinal &compactRow, bool finalPass) const {
    if (finalPass) {
      rowmapOut(row) = compactRowmapOut(compactRow);
    }
    if (row + 1 < rowmapIn.extent_int(0) && rowmapIn(row + 1) != rowmapIn(row)) compactRow++;
  }

  RowmapIn rowmapIn;
  RowmapOut compactRowmapOut;
  RowmapOut rowmapOut;
};

// Given a CrsMatrix A, filter out all entries Aij where |Aij| <= tol.
// If there are no entries to remove, A is returned.
// Otherwise a new matrix is returned.
template <typename Matrix>
Matrix removeCrsMatrixZeros(const Matrix &A,
                            typename Kokkos::ArithTraits<typename Matrix::value_type>::mag_type tol = 0) {
  using Ordinal   = typename Matrix::non_const_ordinal_type;
  using Offset    = typename Matrix::non_const_size_type;
  using Device    = typename Matrix::device_type;
  using ExecSpace = typename Device::execution_space;
  using Mag       = decltype(tol);
  using RangePol  = Kokkos::RangePolicy<ExecSpace>;
  // First, count the number of entries to remove
  Offset entriesToRemove;
  Kokkos::parallel_reduce(RangePol(0, A.nnz()),
                          CountDroppedEntriesFunctor<typename Matrix::values_type, Mag, Offset>(A.values, tol),
                          entriesToRemove);
  if (entriesToRemove == Offset(0)) {
    // The matrix has no zeros to remove, so just return it as-is
    return A;
  }
  // Actually have to make the new matrix with (near-)zeros removed.
  // To help construct the new rowmap, for each original entry record whether
  // it's at the end of its row.
  Kokkos::Bitset<Device> rowEndMarkersNonconst(A.nnz());
  Kokkos::parallel_for(RangePol(0, A.graph.row_map.extent(0)),
                       MarkFinalRowEntries(rowEndMarkersNonconst, A.graph.row_map));
  Offset filteredNNZ = A.nnz() - entriesToRemove;
  typename Matrix::values_type::non_const_type filteredValues(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Afiltered values"), filteredNNZ);
  typename Matrix::index_type::non_const_type filteredEntries(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Afiltered entries"), filteredNNZ);
  typename Matrix::row_map_type::non_const_type compactFilteredRowmap(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Afiltered rowmap (compact)"), A.numRows() + 1);
  typename Matrix::row_map_type::non_const_type filteredRowmap(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Afiltered rowmap"), A.numRows() + 1);
  // Using a parallel scan, compact the non-filtered entries and partially fill
  // in the rowmap (only marking row begins for rows which were originally
  // non-empty) The rest can be filled in with a max-scan.
  Kokkos::ConstBitset<Device> rowEndMarkers(rowEndMarkersNonconst);
  Kokkos::parallel_scan(RangePol(0, A.nnz()),
                        DropEntriesFunctor(rowEndMarkers, A.graph.row_map, A.graph.entries, A.values,
                                           compactFilteredRowmap, filteredEntries, filteredValues, tol));
  Kokkos::parallel_scan(
      RangePol(0, A.numRows() + 1),
      ExpandRowmapFunctor<typename Matrix::row_map_type, typename Matrix::row_map_type::non_const_type, Ordinal>(
          A.graph.row_map, compactFilteredRowmap, filteredRowmap));
  ExecSpace().fence();
  return Matrix("A filtered", A.numRows(), A.numCols(), filteredNNZ, filteredValues, filteredRowmap, filteredEntries);
}

template <typename Rowmap, typename Entries, typename Values>
void validateCrsMatrix(int m, int n, const Rowmap &rowmapIn, const Entries &entriesIn, const Values &valuesIn) {
  auto rowmap  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rowmapIn);
  auto entries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), entriesIn);
  auto values  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), valuesIn);
  size_t nnz   = entries.extent(0);
  if (nnz != values.extent(0)) throw std::runtime_error("Matrix entries/values views have different lengths");
  if ((m == 0 && rowmap.extent(0) > size_t(1)) || (rowmap.extent(0) != size_t(m + 1)))
    throw std::runtime_error("Matrix rowmap has wrong length");
  if (m && nnz != rowmap(m)) throw std::runtime_error("Matrix rowmap final entry doesn't match nnz");
  for (int i = 0; i < m; i++) {
    if (rowmap(i) > rowmap(i + 1)) throw std::runtime_error("Matrix rowmap not ascending");
  }
  for (size_t i = 0; i < size_t(nnz); i++) {
    if (entries(i) >= n) throw std::runtime_error("Matrix entry out of bounds");
  }
}

/**
 * @brief Count the non-zeros of a sub-block in a CRS matrix and find the first
 * and last column indices at each row of the sub-block. This is a host function
 * used by the kk_extract_diagonal_blocks_crsmatrix_sequential()
 */
template <typename row_map_type, typename entries_type, typename ordinal_type, typename size_type,
          typename offset_view1d_type>
void kk_find_nnz_first_last_indices_subblock_crsmatrix_sequential(
    const row_map_type &A_row_map, const entries_type &A_entries, const ordinal_type &blk_row_start,
    const ordinal_type &blk_col_start, const ordinal_type &blk_nrows, const ordinal_type &blk_ncols, size_type &blk_nnz,
    offset_view1d_type &first_indices, offset_view1d_type &last_indices) {
  // Rowmap of i-th row-oriented sub-matrix
  auto A_row_map_sub = Kokkos::subview(A_row_map, Kokkos::make_pair(blk_row_start, blk_row_start + blk_nrows + 1));

  blk_nnz = 0;

  for (ordinal_type j = 0; j < blk_nrows; j++) {  // loop through each row
    size_type k1 = A_row_map_sub(j);
    size_type k2 = A_row_map_sub(j + 1);
    size_type k;
    // Assume column indices are sorted in ascending order
    // Find the position of the start column in the row
    for (k = k1; k < k2; k++) {
      ordinal_type col = A_entries(k);
      if (col >= blk_col_start) {
        break;
      }
    }
    first_indices(j) = k;
    // Find the position of the last column in the row
    for (k = k2 - 1; k >= k1; k--) {
      ordinal_type col = A_entries(k);
      if (col < blk_col_start + blk_ncols) {
        break;
      }
    }
    last_indices(j) = k;
    blk_nnz += (last_indices(j) - first_indices(j) + 1);
  }
}

/**
 * @brief Extract a CRS sub-block from a CRS matrix
 * This is a host function used by the
 * kk_extract_diagonal_blocks_crsmatrix_sequential()
 */
template <typename entries_type, typename values_type, typename ordinal_type, typename size_type,
          typename offset_view1d_type, typename out_row_map_type, typename out_entries_type, typename out_values_type>
void kk_extract_subblock_crsmatrix_sequential(const entries_type &A_entries, const values_type &A_values,
                                              const ordinal_type &blk_col_start, const ordinal_type &blk_nrows,
                                              const size_type &blk_nnz, const offset_view1d_type &first_indices,
                                              const offset_view1d_type &last_indices, out_row_map_type &blk_row_map,
                                              out_entries_type &blk_entries, out_values_type &blk_values) {
  // - create out_row_map
  // - copy A_entries to out_entries and update out_entries with local column
  // indices
  // - copy A_values to out_values
  size_type first_ = 0;
  for (ordinal_type j = 0; j < blk_nrows; j++) {  // loop through each row
    size_type nnz  = last_indices(j) - first_indices(j) + 1;
    blk_row_map(j) = first_;
    for (size_type k = 0; k < nnz; k++) {
      blk_entries(first_ + k) = A_entries(first_indices(j) + k) - blk_col_start;
      blk_values(first_ + k)  = A_values(first_indices(j) + k);
    }
    first_ += nnz;
  }
  blk_row_map(blk_nrows) = blk_nnz;  // last element
}

/**
 * @brief Extract the diagonal blocks out of a crs matrix.
 * This is a blocking function that runs on the host.
 *
 * @tparam crsMat_t The type of the CRS matrix.
 * @param A [in] The square CrsMatrix. It is expected that column indices are
 * in ascending order
 * @param UseRCMReordering [in] Boolean indicating whether applying (true) RCM
 * reordering to diagonal blocks or not (false) (default: false)
 * @param DiagBlk_v [out] The vector of the extracted the CRS diagonal blocks
 * (1 <= the number of diagonal blocks <= A_nrows)
 * @return a vector of lists of vertices in RCM order (a list per a diagonal
 * block) if UseRCMReordering is true, or an empty vector if UseRCMReordering is
 * false
 *
 * Usage Example:
 *    perm = kk_extract_diagonal_blocks_crsmatrix_sequential(A_in, diagBlk_out,
 * UseRCMReordering);
 */
template <typename crsMat_t>
std::vector<typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type>
kk_extract_diagonal_blocks_crsmatrix_sequential(const crsMat_t &A, std::vector<crsMat_t> &DiagBlk_v,
                                                bool UseRCMReordering = false) {
  using row_map_type                = typename crsMat_t::row_map_type;
  using entries_type                = typename crsMat_t::index_type;
  using values_type                 = typename crsMat_t::values_type;
  using graph_t                     = typename crsMat_t::StaticCrsGraphType;
  using out_row_map_type            = typename graph_t::row_map_type::non_const_type;
  using out_entries_type            = typename graph_t::entries_type::non_const_type;
  using out_values_type             = typename crsMat_t::values_type::non_const_type;
  using out_row_map_hostmirror_type = typename out_row_map_type::HostMirror;
  using out_entries_hostmirror_type = typename out_entries_type::HostMirror;
  using out_values_hostmirror_type  = typename out_values_type::HostMirror;

  using ordinal_type       = typename crsMat_t::non_const_ordinal_type;
  using size_type          = typename crsMat_t::non_const_size_type;
  using value_type         = typename crsMat_t::non_const_value_type;
  using offset_view1d_type = Kokkos::View<size_type *, Kokkos::LayoutLeft, Kokkos::HostSpace>;

  row_map_type A_row_map = A.graph.row_map;
  entries_type A_entries = A.graph.entries;
  values_type A_values   = A.values;

  auto A_row_map_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_row_map);
  auto A_entries_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_entries);
  auto A_values_h  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_values);

  ordinal_type A_nrows  = static_cast<ordinal_type>(A.numRows());
  ordinal_type A_ncols  = static_cast<ordinal_type>(A.numCols());
  ordinal_type n_blocks = static_cast<ordinal_type>(DiagBlk_v.size());

  if (A_nrows != A_ncols) {
    std::ostringstream os;
    os << "The diagonal block extraction only works with square matrices -- "
          "matrix A: "
       << A_nrows << " x " << A_ncols;
    throw std::runtime_error(os.str());
  }

  std::vector<out_entries_type> perm_v;
  std::vector<out_entries_hostmirror_type> perm_h_v;

  if (n_blocks == 1) {
    // One block case: simply shallow copy A to DiagBlk_v[0]
    // Note: always not applying RCM reordering, for now
    DiagBlk_v[0] = crsMat_t(A);
  } else {
    // n_blocks > 1
    if (A_nrows == 0) {
      // Degenerate case: A is an empty matrix
      for (ordinal_type i = 0; i < n_blocks; i++) {
        DiagBlk_v[i] = crsMat_t();
      }
    } else {
      // A_nrows >= 1
      if ((n_blocks < 1) || (A_nrows < n_blocks)) {
        std::ostringstream os;
        os << "The number of diagonal blocks (" << n_blocks
           << ") should be >=1 and <= the number of rows of the matrix A (" << A_nrows << ")";
        throw std::runtime_error(os.str());
      }

      ordinal_type rows_per_block = ((A_nrows % n_blocks) == 0) ? (A_nrows / n_blocks) : (A_nrows / n_blocks + 1);

      if (UseRCMReordering) {
        perm_v.resize(n_blocks);
        perm_h_v.resize(n_blocks);
      }

      ordinal_type blk_row_start = 0;     // first row index of i-th diagonal block
      ordinal_type blk_col_start = 0;     // first col index of i-th diagonal block
      ordinal_type blk_nrows, blk_ncols;  // Nrows, Ncols of i-th diagonal block

      for (ordinal_type i = 0; i < n_blocks; i++) {
        blk_nrows = rows_per_block;
        if ((blk_row_start + rows_per_block) > A_nrows) {
          blk_nrows = A_nrows - blk_row_start;
        }
        blk_col_start = blk_row_start;
        blk_ncols     = blk_nrows;

        // First round: count i-th non-zeros or size of entries_v[i] and find
        // the first and last column indices at each row
        size_type blk_nnz = 0;
        offset_view1d_type first(Kokkos::view_alloc(Kokkos::WithoutInitializing, "first"),
                                 blk_nrows);  // first position per row
        offset_view1d_type last(Kokkos::view_alloc(Kokkos::WithoutInitializing, "last"),
                                blk_nrows);  // last position per row

        kk_find_nnz_first_last_indices_subblock_crsmatrix_sequential(
            A_row_map_h, A_entries_h, blk_row_start, blk_col_start, blk_nrows, blk_ncols, blk_nnz, first, last);

        // Second round: extract
        out_row_map_type row_map(Kokkos::view_alloc(Kokkos::WithoutInitializing, "row_map"), blk_nrows + 1);
        out_entries_type entries(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entries"), blk_nnz);
        out_values_type values(Kokkos::view_alloc(Kokkos::WithoutInitializing, "values"), blk_nnz);
        out_row_map_hostmirror_type row_map_h(Kokkos::view_alloc(Kokkos::WithoutInitializing, "row_map_h"),
                                              blk_nrows + 1);
        out_entries_hostmirror_type entries_h(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entries_h"), blk_nnz);
        out_values_hostmirror_type values_h(Kokkos::view_alloc(Kokkos::WithoutInitializing, "values_h"), blk_nnz);

        kk_extract_subblock_crsmatrix_sequential(A_entries_h, A_values_h, blk_col_start, blk_nrows, blk_nnz, first,
                                                 last, row_map_h, entries_h, values_h);

        if (!UseRCMReordering) {
          Kokkos::deep_copy(row_map, row_map_h);
          Kokkos::deep_copy(entries, entries_h);
          Kokkos::deep_copy(values, values_h);
        } else {
          perm_h_v[i] = KokkosGraph::Experimental::graph_rcm<Kokkos::DefaultHostExecutionSpace>(row_map_h, entries_h);
          perm_v[i] =
              out_entries_type(Kokkos::view_alloc(Kokkos::WithoutInitializing, "perm_v"), perm_h_v[i].extent(0));

          out_row_map_hostmirror_type row_map_perm_h(Kokkos::view_alloc(Kokkos::WithoutInitializing, "row_map_perm_h"),
                                                     blk_nrows + 1);
          out_entries_hostmirror_type entries_perm_h(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entries_perm_h"),
                                                     blk_nnz);
          out_values_hostmirror_type values_perm_h(Kokkos::view_alloc(Kokkos::WithoutInitializing, "values_perm_h"),
                                                   blk_nnz);

          out_entries_hostmirror_type reverseperm_h(Kokkos::view_alloc(Kokkos::WithoutInitializing, "reverseperm_h"),
                                                    blk_nrows);
          for (ordinal_type ii = 0; ii < blk_nrows; ii++) reverseperm_h(perm_h_v[i](ii)) = ii;

          std::map<ordinal_type, value_type> colIdx_Value_rcm;

          // Loop through each row of the reordered matrix
          size_type cnt = 0;
          for (ordinal_type ii = 0; ii < blk_nrows; ii++) {
            colIdx_Value_rcm.clear();
            // ii: reordered index
            ordinal_type origRow = reverseperm_h(ii);  // get the original row idx of the reordered row idx, ii
            for (size_type j = row_map_h(origRow); j < row_map_h(origRow + 1); j++) {
              ordinal_type origEi = entries_h(j);
              value_type origV    = values_h(j);
              ordinal_type Ei     = perm_h_v[i](origEi);  // get the reordered col idx of the
                                                          // original col idx, origEi
              colIdx_Value_rcm[Ei] = origV;
            }
            row_map_perm_h(ii) = cnt;
            for (typename std::map<ordinal_type, value_type>::iterator it = colIdx_Value_rcm.begin();
                 it != colIdx_Value_rcm.end(); ++it) {
              entries_perm_h(cnt) = it->first;
              values_perm_h(cnt)  = it->second;
              cnt++;
            }
          }
          row_map_perm_h(blk_nrows) = cnt;

          Kokkos::deep_copy(row_map, row_map_perm_h);
          Kokkos::deep_copy(entries, entries_perm_h);
          Kokkos::deep_copy(values, values_perm_h);
          Kokkos::deep_copy(perm_v[i], perm_h_v[i]);
        }

        DiagBlk_v[i] = crsMat_t("CrsMatrix", blk_nrows, blk_ncols, blk_nnz, values, row_map, entries);

        blk_row_start += blk_nrows;
      }  // for (ordinal_type i = 0; i < n_blocks; i++)
    }    // A_nrows >= 1
  }      // n_blocks > 1
  return perm_v;
}

}  // namespace Impl

using Impl::isCrsGraphSorted;
using Impl::removeCrsMatrixZeros;

}  // namespace KokkosSparse

#endif
