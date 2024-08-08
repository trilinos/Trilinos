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

#include "KokkosGraph_Distance2Color.hpp"

namespace KokkosSparse {

namespace Impl {

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename a_row_view_t__, typename a_nnz_view_t__, typename a_scalar_view_t__, typename b_row_view_t__,
          typename b_nnz_view_t__, typename b_scalar_view_t__, typename c_row_view_t__, typename c_nnz_view_t__,
          typename c_scalar_view_t__>
struct KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::NumericCCOLOR {
  nnz_lno_t numrows;
  nnz_lno_t numcols;

  a_row_view_t__ row_mapA;
  a_nnz_view_t__ entriesA;
  a_scalar_view_t__ valuesA;

  b_row_view_t__ row_mapB;
  b_nnz_view_t__ entriesB;
  b_scalar_view_t__ valuesB;

  c_row_view_t__ rowmapC;
  c_nnz_view_t__ entriesC;
  c_scalar_view_t__ valuesC;
  nnz_lno_t *pEntriesC;
  scalar_t *pVals;

  scalar_temp_work_view_t denseAccumulator;  // initially all zeroes
  scalar_t *pdenseAccumulator;               // initially all zeroes

  // bool_temp_view_t denseAccumulatorFlags; //initially all false.
  // bool * pdenseAccumulatorFlags; //initially all false.

  nnz_lno_t team_work_size;

  nnz_lno_t color_begin;
  nnz_lno_t color_end;
  nnz_lno_persistent_work_view_t color_adj;
  nnz_lno_persistent_work_view_t vertex_colors;

  nnz_lno_t consecutive_chunk_size;
  nnz_lno_t consecutive_all_color_chunk_size;

  nnz_lno_t chunk_divison;
  nnz_lno_t chunk_and;
  NumericCCOLOR(

      nnz_lno_t m_, nnz_lno_t k_,

      a_row_view_t__ row_mapA_, a_nnz_view_t__ entriesA_, a_scalar_view_t__ valuesA_,

      b_row_view_t__ row_mapB_, b_nnz_view_t__ entriesB_, b_scalar_view_t__ valuesB_,

      c_row_view_t__ rowmapC_, c_nnz_view_t__ entriesC_, c_scalar_view_t__ valuesC_,
      scalar_temp_work_view_t denseAccumulator_,  // initially all zeroes
      // bool_temp_view_t denseAccumulatorFlags_, //initially all false.

      nnz_lno_t team_row_work_size_)
      : numrows(m_),
        numcols(k_),
        row_mapA(row_mapA_),
        entriesA(entriesA_),
        valuesA(valuesA_),

        row_mapB(row_mapB_),
        entriesB(entriesB_),
        valuesB(valuesB_),

        rowmapC(rowmapC_),
        entriesC(entriesC_),
        valuesC(valuesC_),
        pEntriesC(entriesC_.data()),
        pVals(valuesC.data()),
        denseAccumulator(denseAccumulator_),
        pdenseAccumulator(denseAccumulator_.data()),
        // denseAccumulatorFlags (denseAccumulatorFlags_),
        // pdenseAccumulatorFlags(denseAccumulatorFlags_.data()),
        team_work_size(team_row_work_size_),
        color_begin(0),
        color_end(0),
        color_adj(),
        vertex_colors(),
        consecutive_chunk_size(numcols),
        consecutive_all_color_chunk_size(numcols),
        chunk_divison(0),
        chunk_and(0) {}

  // one color at a time.
  KOKKOS_INLINE_FUNCTION
  void operator()(const Numeric1Tag &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size + color_begin;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, color_end);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &color_index_index) {
                           nnz_lno_t row_index = color_adj(color_index_index);

                           // nnz_lno_t color = vertex_colors(row_index);
                           scalar_t *mydenseAccumulator = pdenseAccumulator;  // + numcols * color;
                           const size_type col_begin    = row_mapA[row_index];
                           const nnz_lno_t left_work    = nnz_lno_t(row_mapA[row_index + 1] - col_begin);
                           const size_type c_row_begin  = rowmapC[row_index];
                           nnz_lno_t *my_entries        = pEntriesC + c_row_begin;
                           for (nnz_lno_t colind = 0; colind < left_work; ++colind) {
                             size_type a_col = colind + col_begin;
                             nnz_lno_t rowB  = entriesA[a_col];
                             scalar_t valA   = valuesA[a_col];

                             size_type rowBegin   = row_mapB(rowB);
                             nnz_lno_t left_work_ = row_mapB(rowB + 1) - rowBegin;

                             Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, left_work_), [&](nnz_lno_t i) {
                               const size_type adjind = i + rowBegin;
                               nnz_lno_t acc_index    = entriesB[adjind];
                               scalar_t b_val         = valuesB[adjind] * valA;
                               mydenseAccumulator[acc_index] += b_val;
                             });
                           }
                           scalar_t *my_vals = pVals + c_row_begin;

                           nnz_lno_t row_size = rowmapC[row_index + 1] - c_row_begin;
                           Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, row_size), [&](nnz_lno_t i) {
                             nnz_lno_t acc_index           = my_entries[i];
                             my_vals[i]                    = mydenseAccumulator[acc_index];
                             mydenseAccumulator[acc_index] = 0;
                           });
                         });
  }

  // multi-color minimized writes
  KOKKOS_INLINE_FUNCTION
  void operator()(const Numeric2Tag &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size + color_begin;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, color_end);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &color_index_index) {
                           nnz_lno_t row_index = color_adj(color_index_index);

                           nnz_lno_t color              = vertex_colors(row_index);
                           scalar_t *mydenseAccumulator = pdenseAccumulator + numcols * color;

                           const size_type col_begin   = row_mapA[row_index];
                           const nnz_lno_t left_work   = nnz_lno_t(row_mapA[row_index + 1] - col_begin);
                           const size_type c_row_begin = rowmapC[row_index];
                           nnz_lno_t *my_entries       = pEntriesC + c_row_begin;
                           for (nnz_lno_t colind = 0; colind < left_work; ++colind) {
                             size_type a_col = colind + col_begin;
                             nnz_lno_t rowB  = entriesA[a_col];
                             scalar_t valA   = valuesA[a_col];

                             size_type rowBegin   = row_mapB(rowB);
                             nnz_lno_t left_work_ = row_mapB(rowB + 1) - rowBegin;

                             Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, left_work_), [&](nnz_lno_t i) {
                               const size_type adjind = i + rowBegin;
                               nnz_lno_t acc_index    = entriesB[adjind];
                               scalar_t b_val         = valuesB[adjind] * valA;
                               mydenseAccumulator[acc_index] += b_val;
                             });
                           }
                           scalar_t *my_vals = pVals + c_row_begin;

                           nnz_lno_t row_size = rowmapC[row_index + 1] - c_row_begin;
                           Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, row_size), [&](nnz_lno_t i) {
                             nnz_lno_t acc_index           = my_entries[i];
                             my_vals[i]                    = mydenseAccumulator[acc_index];
                             mydenseAccumulator[acc_index] = 0;
                           });
                         });
  }

  // multi-color minimized reads
  KOKKOS_INLINE_FUNCTION
  void operator()(const Numeric3Tag &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size + color_begin;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, color_end);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &color_index_index) {
                           nnz_lno_t row_index = color_adj(color_index_index);

                           nnz_lno_t color              = vertex_colors(row_index);
                           scalar_t *mydenseAccumulator = pdenseAccumulator + numcols * color;
                           const size_type col_begin    = row_mapA[row_index];
                           const nnz_lno_t left_work    = nnz_lno_t(row_mapA[row_index + 1] - col_begin);
                           const size_type c_row_begin  = rowmapC[row_index];
                           nnz_lno_t *my_entries        = pEntriesC + c_row_begin;
                           for (nnz_lno_t colind = 0; colind < left_work; ++colind) {
                             size_type a_col = colind + col_begin;
                             nnz_lno_t rowB  = entriesA[a_col];
                             scalar_t valA   = valuesA[a_col];

                             size_type rowBegin   = row_mapB(rowB);
                             nnz_lno_t left_work_ = row_mapB(rowB + 1) - rowBegin;

                             Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, left_work_), [&](nnz_lno_t i) {
                               const size_type adjind = i + rowBegin;
                               nnz_lno_t col_ind      = entriesB[adjind];
                               nnz_lno_t acc_index    = col_ind;
                               // nnz_lno_t acc_index = (col_ind  >> chunk_divison) *
                               // (consecutive_all_color_chunk_size) + (color <<
                               // chunk_divison)+ (col_ind & chunk_and);
                               scalar_t b_val = valuesB[adjind] * valA;
                               mydenseAccumulator[acc_index] += b_val;
                             });
                           }
                           scalar_t *my_vals = pVals + c_row_begin;

                           nnz_lno_t row_size = rowmapC[row_index + 1] - c_row_begin;
                           Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, row_size), [&](nnz_lno_t i) {
                             nnz_lno_t col_ind   = my_entries[i];
                             nnz_lno_t acc_index = col_ind;

                             // nnz_lno_t acc_index = (col_ind  >>
                             // chunk_divison) *
                             // (consecutive_all_color_chunk_size) + (color
                             // << chunk_divison)+ (col_ind & chunk_and);
                             my_vals[i]                    = mydenseAccumulator[acc_index];
                             mydenseAccumulator[acc_index] = 0;
                           });
                         });
  }

  size_t team_shmem_size(int team_size) const { return team_size * sizeof(nnz_lno_t) * 8; }
};

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
void KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                  b_lno_nnz_view_t_,
                  b_scalar_nnz_view_t_>::KokkosSPGEMM_numeric_color(c_row_view_t rowmapC_, c_lno_nnz_view_t entriesC_,
                                                                    c_scalar_nnz_view_t valuesC_,
                                                                    SPGEMMAlgorithm spgemm_algorithm_) {
  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\tCOLOR MODE" << std::endl;
  }
  nnz_lno_temp_work_view_t entryIndicesC_ = this->handle->get_spgemm_handle()->get_c_column_indices();

  KokkosKernels::Impl::kk_copy_vector<nnz_lno_temp_work_view_t, c_lno_nnz_view_t, MyExecSpace>(
      entryIndicesC_.extent(0), entryIndicesC_, entriesC_);

  // KokkosKernels::Impl::ExecSpaceType my_exec_space =
  //    KokkosKernels::Impl::get_exec_space_type<MyExecSpace>();

  nnz_lno_t brows = row_mapB.extent(0) - 1;
  size_type bnnz  = valsB.extent(0);
  // get vector size, team size.
  int suggested_vector_size = this->handle->get_suggested_vector_size(brows, bnnz);
  int suggested_team_size   = this->handle->get_suggested_team_size(suggested_vector_size);

  // get color vertices
  nnz_lno_t num_colors, num_multi_colors, num_used_colors;
  nnz_lno_persistent_work_host_view_t color_xadj;
  nnz_lno_persistent_work_view_t color_adj, vertex_colors;
  this->handle->get_spgemm_handle()->get_color_xadj(num_colors, color_xadj, color_adj, vertex_colors, num_multi_colors,
                                                    num_used_colors);

  const nnz_lno_t block_size    = 64;
  const nnz_lno_t shift_divisor = 6;
  scalar_temp_work_view_t denseAccumulator("Scalar Accumulator", ((this->b_col_cnt + block_size) * num_multi_colors));

  // bool_temp_view_t denseAccumulatorFlags ("Accumulator flags",
  // ((this->k* 1.5)  * num_multi_colors));

  NumericCCOLOR<const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t, const_b_lno_row_view_t,
                const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t, c_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t>
      sc(a_row_cnt, b_col_cnt, row_mapA, entriesA, valsA,

         row_mapB, entriesB, valsB,

         rowmapC_, entriesC_, valuesC_,

         denseAccumulator,
         // denseAccumulatorFlags,
         -1);

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\tCOLORING-num_multi_colors:" << num_multi_colors << " num_used_colors:" << num_used_colors
              << std::endl;
  }
  sc.color_adj                        = color_adj;
  sc.vertex_colors                    = vertex_colors;
  sc.chunk_divison                    = shift_divisor;
  sc.chunk_and                        = block_size - 1;
  sc.consecutive_chunk_size           = block_size;
  sc.consecutive_all_color_chunk_size = sc.consecutive_chunk_size * num_multi_colors;

  Kokkos::Timer timer1;
  for (nnz_lno_t i = 0; i < num_used_colors;) {
    nnz_lno_t color_begin = color_xadj(i);
    nnz_lno_t lastcolor   = i + 1;
    if (spgemm_algorithm_ == SPGEMM_KK_MULTICOLOR2) {
      lastcolor = KOKKOSKERNELS_MACRO_MIN(i + num_multi_colors, num_used_colors);
      i += num_multi_colors;
    } else {
      ++i;
    }

    nnz_lno_t color_end = color_xadj(lastcolor);
    sc.color_begin      = color_begin;
    sc.color_end        = color_end;

    nnz_lno_t team_row_chunk_size =
        this->handle->get_team_work_size(suggested_team_size, concurrency, color_end - color_begin);
    sc.team_work_size = team_row_chunk_size;

    if (use_dynamic_schedule) {
      switch (spgemm_algorithm_) {
        default:
        case SPGEMM_KK_COLOR:
          Kokkos::parallel_for(dynamic_team_numeric1_policy_t((color_end - color_begin) / team_row_chunk_size + 1,
                                                              suggested_team_size, suggested_vector_size),
                               sc);
          break;
        case SPGEMM_KK_MULTICOLOR2:
          Kokkos::parallel_for(dynamic_team_numeric2_policy_t((color_end - color_begin) / team_row_chunk_size + 1,
                                                              suggested_team_size, suggested_vector_size),
                               sc);
          break;
        case SPGEMM_KK_MULTICOLOR:
          Kokkos::parallel_for(dynamic_team_numeric3_policy_t((color_end - color_begin) / team_row_chunk_size + 1,
                                                              suggested_team_size, suggested_vector_size),
                               sc);
          break;
      }
    } else {
      switch (spgemm_algorithm_) {
        default:
        case SPGEMM_KK_COLOR:
          Kokkos::parallel_for(team_numeric1_policy_t((color_end - color_begin) / team_row_chunk_size + 1,
                                                      suggested_team_size, suggested_vector_size),
                               sc);
          break;
        case SPGEMM_KK_MULTICOLOR2:
          Kokkos::parallel_for(team_numeric2_policy_t((color_end - color_begin) / team_row_chunk_size + 1,
                                                      suggested_team_size, suggested_vector_size),
                               sc);
          break;
        case SPGEMM_KK_MULTICOLOR:
          Kokkos::parallel_for(team_numeric3_policy_t((color_end - color_begin) / team_row_chunk_size + 1,
                                                      suggested_team_size, suggested_vector_size),
                               sc);
          break;
      }
    }
    MyExecSpace().fence();
  }

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\tNumeric TIME:" << timer1.seconds() << std::endl;
  }
}

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename c_row_view_t, typename c_nnz_view_t>
void KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                  b_lno_nnz_view_t_,
                  b_scalar_nnz_view_t_>::d2_color_c_matrix(c_row_view_t rowmapC, c_nnz_view_t entryIndicesC_,

                                                           nnz_lno_t &original_num_colors,
                                                           nnz_lno_persistent_work_host_view_t &h_color_xadj,
                                                           nnz_lno_persistent_work_view_t &color_adj,
                                                           nnz_lno_persistent_work_view_t &vertex_colors_to_store,

                                                           nnz_lno_t &num_colors_in_one_step,
                                                           nnz_lno_t &num_multi_color_steps,
                                                           SPGEMMAlgorithm spgemm_algorithm_) {
  nnz_lno_persistent_work_view_t color_xadj;

  size_type c_nnz_size = this->handle->get_spgemm_handle()->get_c_nnz();
  ;

  // first we need to transpose the C graph.
  // allocate memory for that.
  row_lno_temp_work_view_t transpose_col_xadj;  // ("transpose_col_xadj", b_col_cnt + 1);
  nnz_lno_temp_work_view_t transpose_col_adj;   // (Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                                // "tmp_row_view"), c_nnz_size);

  // KokkosKernels::Impl::ExecSpaceType my_exec_space =
  // this->handle->get_handle_exec_space();
  int suggested_vector_size     = this->handle->get_suggested_vector_size(rowmapC.extent(0) - 1, c_nnz_size);
  int suggested_team_size       = this->handle->get_suggested_team_size(suggested_vector_size);
  nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size, concurrency, a_row_cnt);

  Kokkos::Timer timer1;
  if (this->handle->get_spgemm_handle()->coloring_input_file == "") {
    transpose_col_xadj = row_lno_temp_work_view_t("transpose_col_xadj", b_col_cnt + 1);
    transpose_col_adj =
        nnz_lno_temp_work_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "tmp_row_view"), c_nnz_size);

    KokkosKernels::Impl::transpose_graph<c_row_view_t, c_nnz_view_t, row_lno_temp_work_view_t, nnz_lno_temp_work_view_t,
                                         row_lno_temp_work_view_t, MyExecSpace>(
        a_row_cnt, b_col_cnt, rowmapC, entryIndicesC_, transpose_col_xadj, transpose_col_adj, suggested_vector_size,
        suggested_team_size, team_row_chunk_size, use_dynamic_schedule);

    MyExecSpace().fence();
  }
  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\tTranspose Time:" << timer1.seconds() << std::endl;
  }

  {
    timer1.reset();

    this->handle->create_graph_coloring_handle();

    typename HandleType::GraphColoringHandleType::color_view_t vertex_color_view;

    if (this->handle->get_spgemm_handle()->coloring_input_file == "") {
      // for now only sequential one exists.
      // find distance-2 graph coloring

      auto gchD2 = handle->get_distance2_graph_coloring_handle();

      KokkosGraph::Experimental::graph_compute_distance2_color<HandleType, c_row_view_t, c_nnz_view_t,
                                                               row_lno_temp_work_view_t, nnz_lno_temp_work_view_t>(
          this->handle, a_row_cnt, b_col_cnt, rowmapC, entryIndicesC_, transpose_col_xadj, transpose_col_adj);

      original_num_colors = handle->get_graph_coloring_handle()->get_num_colors();

      if (KOKKOSKERNELS_VERBOSE) {
        std::cout << "\t\tNum colors:" << handle->get_graph_coloring_handle()->get_num_colors()
                  << " coloring time:" << timer1.seconds() << std::endl;
      }
      vertex_color_view = handle->get_graph_coloring_handle()->get_vertex_colors();

      if (this->handle->get_spgemm_handle()->coloring_output_file != "") {
        KokkosKernels::Impl::kk_write_1Dview_to_file(vertex_color_view,
                                                     this->handle->get_spgemm_handle()->coloring_output_file.c_str());
      }
    } else {
      vertex_color_view =
          typename HandleType::GraphColoringHandleType::color_view_t("vertex colors from file", a_row_cnt);
      KokkosKernels::Impl::kk_read_1Dview_from_file(vertex_color_view,
                                                    this->handle->get_spgemm_handle()->coloring_input_file.c_str());
      KokkosKernels::Impl::view_reduce_max<typename HandleType::GraphColoringHandleType::color_view_t, MyExecSpace>(
          a_row_cnt, vertex_color_view, original_num_colors);
      MyExecSpace().fence();

      // KokkosKernels::Impl::kk_print_1Dview(vertex_color_view);
    }
    num_multi_color_steps  = original_num_colors;
    num_colors_in_one_step = 1;

    vertex_colors_to_store = nnz_lno_persistent_work_view_t(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "persistent_color_view"), a_row_cnt);

    if (KOKKOSKERNELS_VERBOSE) {
      // create histogram and print it.
      nnz_lno_temp_work_view_t histogram("histogram", original_num_colors + 1);
      MyExecSpace().fence();
      timer1.reset();
      KokkosKernels::Impl::kk_get_histogram<typename HandleType::GraphColoringHandleType::color_view_t,
                                            nnz_lno_temp_work_view_t, MyExecSpace>(a_row_cnt, vertex_color_view,
                                                                                   histogram);
      std::cout << "\t\tHistogram"
                << " time:" << timer1.seconds() << std::endl
                << "\t\t";
      KokkosKernels::Impl::kk_print_1Dview(histogram);
    }

    {
      // TODO: this should have been as below.
      // nnz_lno_temp_work_view_t tmp_color_view = vertex_color_view;
      typename HandleType::GraphColoringHandleType::color_view_t tmp_color_view = vertex_color_view;
      // if the algorithm is spgemm, then we will have multiple colors per
      // iteration.
      if (spgemm_algorithm_ == SPGEMM_KK_MULTICOLOR) {
        // tmp_color_view = nnz_lno_temp_work_view_t(
        // Kokkos::view_alloc(Kokkos::WithoutInitializing, "tmp_color_view"),
        // a_row_cnt);
        tmp_color_view = typename HandleType::GraphColoringHandleType::color_view_t(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "tmp_color_view"), a_row_cnt);

        // upper bound is the output size for dense acumulators.
        num_colors_in_one_step = c_nnz_size / this->b_col_cnt;

        // scale if provided.
        double scale_ = this->handle->get_spgemm_handle()->get_multi_color_scale();

        num_colors_in_one_step = scale_ * num_colors_in_one_step;

        // at the end of this tmp_color_view holds the colors that correspond
        // the step colors.
        // that is if num_multi_colors is 32, first 32 is 1, next 32 is 2 and so
        // on.
        if (num_colors_in_one_step > 1) {
          float scale_factor = 1.0 / num_colors_in_one_step;

          // get the sets multicolors. color(i) / num_multi_colors + 1 is the
          // new color.
          KokkosKernels::Impl::kk_a_times_x_plus_b<
              typename HandleType::GraphColoringHandleType::color_view_t  // nnz_lno_temp_work_view_t
              ,
              typename HandleType::GraphColoringHandleType::color_view_t, float, float, MyExecSpace>(
              a_row_cnt, tmp_color_view, vertex_color_view, scale_factor, 0);
          num_multi_color_steps = original_num_colors / num_colors_in_one_step;
          if (original_num_colors % num_colors_in_one_step) ++num_multi_color_steps;
        } else {
          num_colors_in_one_step = 1;
        }
      } else {
        KokkosKernels::Impl::kk_a_times_x_plus_b<
            typename HandleType::GraphColoringHandleType::color_view_t  // nnz_lno_temp_work_view_t,
            ,
            typename HandleType::GraphColoringHandleType::color_view_t  // nnz_lno_temp_work_view_t
            ,
            int, int, MyExecSpace>(a_row_cnt, tmp_color_view, tmp_color_view, 1, -1);
      }

      if (spgemm_algorithm_ == SPGEMM_KK_MULTICOLOR2) {
        num_multi_color_steps  = original_num_colors;
        num_colors_in_one_step = c_nnz_size / this->b_col_cnt;
        double scale_          = this->handle->get_spgemm_handle()->get_multi_color_scale();
        num_colors_in_one_step = scale_ * num_colors_in_one_step;
      }

      // with the modular operation, we find the colors within a step.
      // that is if num_multi_colors is 32, then color 32 will be 0, 33 will be
      // 1, 34 will be 2.
      // it will hold their color within their multicolor step.
      KokkosKernels::Impl::kk_modular_view<nnz_lno_persistent_work_view_t,
                                           typename HandleType::GraphColoringHandleType::color_view_t, MyExecSpace>(
          a_row_cnt, vertex_colors_to_store, vertex_color_view, num_colors_in_one_step);
      timer1.reset();

      // allocate color xadj and adj arrays.
      color_xadj = nnz_lno_persistent_work_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Reverse xadj"),
                                                  num_multi_color_steps + 1);
      color_adj =
          nnz_lno_persistent_work_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Reverse xadj"), a_row_cnt);

      // create reverse map from colors.
      KokkosKernels::Impl::kk_create_reverse_map<
          typename HandleType::GraphColoringHandleType::color_view_t  // nnz_lno_temp_work_view_t
          ,
          nnz_lno_persistent_work_view_t, MyExecSpace>(a_row_cnt, num_multi_color_steps, tmp_color_view, color_xadj,
                                                       color_adj);
      MyExecSpace().fence();

      if (KOKKOSKERNELS_VERBOSE) {
        std::cout << "\t\tReverse Map Create Time:" << timer1.seconds() << std::endl;
      }
      h_color_xadj = Kokkos::create_mirror_view(color_xadj);
      Kokkos::deep_copy(h_color_xadj, color_xadj);
      MyExecSpace().fence();
    }
    this->handle->destroy_graph_coloring_handle();
  }
}

}  // namespace Impl
}  // namespace KokkosSparse
