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

namespace KokkosSparse {
namespace Impl {
#ifdef KOKKOSKERNELS_ANALYZE_MEMORYACCESS
template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
struct KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::Cache {
  const size_type cache_line_count;
  std::vector<size_type> cache_lines;
  std::vector<char> cache_line_owners;
  std::vector<int> cache_inserting_hyperthread;

  size_type *pcache_lines;
  char *pcache_line_owners;
  int *pcache_inserting_hyperthread;

  Cache(size_type cache_line_count_)
      : cache_line_count(cache_line_count_),
        cache_lines(cache_line_count_, size_type(0)),
        cache_line_owners(cache_line_count_, char(0)),
        cache_inserting_hyperthread(cache_line_count_, int(0)) {
    pcache_lines                 = &(cache_lines[0]);
    pcache_line_owners           = &(cache_line_owners[0]);
    pcache_inserting_hyperthread = &(cache_inserting_hyperthread[0]);
  }

  void clear_cache() {
    for (size_type i = 0; i < cache_line_count; ++i) {
      pcache_lines[i]                 = 0;
      pcache_line_owners[i]           = 0;
      pcache_inserting_hyperthread[i] = 0;
    }
  }
  bool is_in_cache(char owner, size_type line_index, size_type insertion_index) {
    if (pcache_line_owners[insertion_index] == owner && pcache_lines[insertion_index] == line_index) {
      return true;
    } else {
      return false;
    }
  }
  void insert_to_cache(char owner, size_type line_index, size_type insertion_index, int inserting_hyperthread = 0) {
    pcache_line_owners[insertion_index]           = owner;
    pcache_lines[insertion_index]                 = line_index;
    pcache_inserting_hyperthread[insertion_index] = inserting_hyperthread;
  }

  Cache &operator=(const Cache &a) {
    cache_line_count             = a.cache_line_count;
    pcache_lines                 = a.pcache_lines;
    pcache_line_owners           = a.pcache_line_owners;
    pcache_inserting_hyperthread = a.pcache_inserting_hyperthread;
    return *this;
  }
};

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename a_row_view_t, typename a_nnz_view_t, typename b_row_view_t, typename b_nnz_view_t,
          typename c_row_view_t>
struct KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::FlopsPerRow {
  nnz_lno_t m;
  a_row_view_t row_mapA;
  a_nnz_view_t entriesA;

  b_row_view_t row_mapB;
  b_nnz_view_t entriesB;

  c_row_view_t rough_row_mapC;

  const size_type min_val;
  nnz_lno_t team_row_chunk_size;

  row_lno_temp_work_view_t c_comp_a_net_index;
  row_lno_temp_work_view_t c_comp_b_net_index;
  nnz_lno_temp_work_view_t c_comp_row_index;
  nnz_lno_temp_work_view_t c_comp_col_index;

  FlopsPerRow(nnz_lno_t m_, a_row_view_t row_mapA_, a_nnz_view_t entriesA_, b_row_view_t row_mapB_,
              b_nnz_view_t entriesB_, c_row_view_t rough_row_mapC_, nnz_lno_t team_row_chunk_size_)
      : m(m_),
        row_mapA(row_mapA_),
        entriesA(entriesA_),
        row_mapB(row_mapB_),
        entriesB(entriesB_),
        rough_row_mapC(rough_row_mapC_),
        min_val(((std::numeric_limits<size_type>::lowest()))),
        team_row_chunk_size(team_row_chunk_size_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const CountTag &, const team_member_t &teamMember, size_t &overal_max) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, m);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_index) {
                           // nnz_lno_t row_index = teamMember.league_rank()  *
                           // teamMember.team_size()+ teamMember.team_rank(); check ii is out of
                           // range. if it is, just return. if (row_index >= m) return;
                           const size_type col_begin = row_mapA[row_index];
                           const size_type col_end   = row_mapA[row_index + 1];
                           const nnz_lno_t left_work = col_end - col_begin;

                           size_type max_num_results_in_row = 0;

                           Kokkos::parallel_reduce(
                               Kokkos::ThreadVectorRange(teamMember, left_work),
                               [&](nnz_lno_t i, size_type &valueToUpdate) {
                                 const size_type adjind   = i + col_begin;
                                 const nnz_lno_t colIndex = entriesA[adjind];
                                 valueToUpdate += row_mapB[colIndex + 1] - row_mapB[colIndex];
                                 // valueToUpdate += row_mapB [colIndex] - oldrow_mapB[colIndex];
                               },
                               max_num_results_in_row);
                           overal_max += max_num_results_in_row;
                           rough_row_mapC(row_index) = max_num_results_in_row;
                         });
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const FillTag &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, m);
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&](const nnz_lno_t &row_index) {
          // nnz_lno_t row_index = teamMember.league_rank()  *
          // teamMember.team_size()+ teamMember.team_rank(); check ii is out of
          // range. if it is, just return. if (row_index >= m) return;
          const size_type col_begin = row_mapA[row_index];
          const size_type col_end   = row_mapA[row_index + 1];
          const nnz_lno_t left_work = col_end - col_begin;

          size_type mult_index = rough_row_mapC(row_index);

          // size_type max_num_results_in_row = 0;

          for (nnz_lno_t i = 0; i < left_work; ++i) {
            const size_type adjind   = i + col_begin;
            const nnz_lno_t colIndex = entriesA[adjind];

            const size_type browbegin  = row_mapB[colIndex];
            const size_type browend    = row_mapB[colIndex + 1];
            const nnz_lno_t bleft_work = browend - browbegin;

            Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, bleft_work), [&](nnz_lno_t j) {
              size_type bcolind = browbegin + j;
              nnz_lno_t bcol    = entriesB(bcolind);

              if (mult_index + j >= rough_row_mapC(row_index + 1))
                std::cout << "mult_index:" << mult_index << " j:" << j
                          << " rough_row_mapC(row_index + 1):" << rough_row_mapC(row_index + 1) << std::endl;
              c_comp_a_net_index(mult_index + j) = adjind;
              c_comp_b_net_index(mult_index + j) = bcolind;
              c_comp_row_index(mult_index + j)   = row_index;
              c_comp_col_index(mult_index + j)   = bcol;
            });
            mult_index += bleft_work;
          }
        });
  }
};

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
void KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                  b_lno_nnz_view_t_,
                  b_scalar_nnz_view_t_>::create_read_write_hg(size_t &overall_flops,
                                                              row_lno_temp_work_view_t &c_flop_rowmap,
                                                              row_lno_temp_work_view_t &c_comp_a_net_index,
                                                              row_lno_temp_work_view_t &c_comp_b_net_index,
                                                              nnz_lno_temp_work_view_t &c_comp_row_index,
                                                              nnz_lno_temp_work_view_t &c_comp_col_index) {
  overall_flops = 0;
  c_flop_rowmap = row_lno_temp_work_view_t("flops per row", a_row_cnt + 1);

  // KokkosKernels::Impl::ExecSpaceType my_exec_space =
  // KokkosKernels::Impl::get_exec_space_type<MyExecSpace>();
  int suggested_vector_size     = this->handle->get_suggested_vector_size(a_row_cnt, entriesA.extent(0));
  int suggested_team_size       = this->handle->get_suggested_team_size(suggested_vector_size);
  nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size, concurrency, a_row_cnt);

  FlopsPerRow<const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_b_lno_row_view_t, const_b_lno_nnz_view_t,
              row_lno_temp_work_view_t>
      pcnnnz(a_row_cnt, row_mapA, entriesA, row_mapB, entriesB, c_flop_rowmap, team_row_chunk_size);

  // calculate how many flops per row is performed
  Kokkos::parallel_reduce(
      team_count_policy_t(a_row_cnt / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size), pcnnnz,
      overall_flops);
  MyExecSpace().fence();

  // do a parallel prefix sum
  KokkosKernels::Impl::exclusive_parallel_prefix_sum<row_lno_temp_work_view_t, MyExecSpace>(a_row_cnt + 1,
                                                                                            c_flop_rowmap);
  MyExecSpace().fence();

  std::cout << "overall_flops:" << overall_flops << std::endl;

  // hypergraph creation. each computation is a vertex. There are overall_flops
  // many computation vertices. Each of them is connected to nonzero of a,
  // nonzero of b and nonzero of c.
  c_comp_a_net_index = row_lno_temp_work_view_t("a_net", overall_flops);
  c_comp_b_net_index = row_lno_temp_work_view_t("b_net", overall_flops);
  c_comp_row_index   = nnz_lno_temp_work_view_t("row ", overall_flops);
  c_comp_col_index   = nnz_lno_temp_work_view_t("col ", overall_flops);

  pcnnnz.c_comp_a_net_index = c_comp_a_net_index;
  pcnnnz.c_comp_b_net_index = c_comp_b_net_index;
  pcnnnz.c_comp_row_index   = c_comp_row_index;
  pcnnnz.c_comp_col_index   = c_comp_col_index;

  // fill the hypergraph values.
  // indices of nnzs for a and b nets, for c nets, the row and column index.
  Kokkos::parallel_for(
      team_fill_policy_t(a_row_cnt / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size), pcnnnz);
  MyExecSpace().fence();
}

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename c_row_view_t>
void KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                  b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::print_read_write_cost(c_row_view_t rowmapC) {
  size_t overall_flops = 0;
  row_lno_temp_work_view_t c_flop_rowmap;  //("flops per row", a_row_cnt + 1);

  row_lno_temp_work_view_t c_comp_a_net_index;  //("a_net", overall_flops);
  row_lno_temp_work_view_t c_comp_b_net_index;  //("b_net", overall_flops);
  nnz_lno_temp_work_view_t c_comp_row_index;    //("row ", overall_flops);
  nnz_lno_temp_work_view_t c_comp_col_index;    //("col ", overall_flops);

  // create the hypergraph.
  this->create_read_write_hg(overall_flops, c_flop_rowmap, c_comp_a_net_index, c_comp_b_net_index, c_comp_row_index,
                             c_comp_col_index);

  int write_type                   = 0;  // 0 -- KKMEM, 1-KKSPEED, 2- KKCOLOR 3-KKMULTICOLOR 4-KKMULTICOLOR2
  SPGEMMAlgorithm spgemm_algorithm = this->handle->get_spgemm_handle()->get_algorithm_type();

  if (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_COLOR) {
    write_type = 2;
  } else if (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR) {
    write_type = 3;
  } else if (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2) {
    write_type = 4;
  } else if (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_SPEED) {
    write_type = 1;
  }

  typename row_lno_temp_work_view_t::HostMirror h_c_flop_rowmap      = Kokkos::create_mirror_view(c_flop_rowmap);
  typename row_lno_temp_work_view_t::HostMirror h_c_comp_a_net_index = Kokkos::create_mirror_view(c_comp_a_net_index);
  typename row_lno_temp_work_view_t::HostMirror h_c_comp_b_net_index = Kokkos::create_mirror_view(c_comp_b_net_index);
  typename nnz_lno_temp_work_view_t::HostMirror h_c_comp_row_index   = Kokkos::create_mirror_view(c_comp_row_index);
  typename nnz_lno_temp_work_view_t::HostMirror h_c_comp_col_index   = Kokkos::create_mirror_view(c_comp_col_index);

  Kokkos::deep_copy(h_c_flop_rowmap, c_flop_rowmap);
  Kokkos::deep_copy(h_c_comp_a_net_index, c_comp_a_net_index);
  Kokkos::deep_copy(h_c_comp_b_net_index, c_comp_b_net_index);
  Kokkos::deep_copy(h_c_comp_row_index, c_comp_row_index);
  Kokkos::deep_copy(h_c_comp_col_index, c_comp_col_index);

  // nnz_lno_t num_parallel_colors = 1;
  bool isGPU = false;
  // int num_threads=68;
  int vectorlane      = 1;
  int cache_line_size = 64;
  int data_size       = 8;
  int cache_l1_size   = 32 * 1024;

  nnz_lno_t num_colors = 0, num_multi_colors = 0, num_used_colors_steps = 0;
  nnz_lno_persistent_work_host_view_t color_xadj;
  nnz_lno_persistent_work_view_t color_adj, vertex_colors;

  this->handle->get_spgemm_handle()->get_color_xadj(num_colors, color_xadj, color_adj, vertex_colors, num_multi_colors,
                                                    num_used_colors_steps);

  nnz_lno_persistent_work_host_view_t h_color_adj;
  nnz_lno_persistent_work_host_view_t h_vertex_colors;

  if (write_type == 2 || write_type == 3) {
    num_used_colors_steps = num_colors / num_multi_colors;
    if (num_colors % num_multi_colors) num_used_colors_steps++;
    num_multi_colors = 1;
    KokkosKernels::Impl::print_1Dview(vertex_colors);
    KokkosKernels::Impl::print_1Dview(color_adj);
    KokkosKernels::Impl::print_1Dview(color_xadj);
    h_color_adj     = Kokkos::create_mirror_view(color_adj);
    h_vertex_colors = Kokkos::create_mirror_view(vertex_colors);
    Kokkos::deep_copy(h_color_adj, color_adj);
    Kokkos::deep_copy(h_vertex_colors, vertex_colors);
  } else if (write_type == 4) {
    // num_used_colors_steps  = num_colors / num_multi_colors;
    // if (num_colors % num_multi_colors) num_used_colors_steps++;
    num_used_colors_steps = num_colors;
    KokkosKernels::Impl::print_1Dview(vertex_colors);
    KokkosKernels::Impl::print_1Dview(color_adj);
    KokkosKernels::Impl::print_1Dview(color_xadj);
    h_color_adj     = Kokkos::create_mirror_view(color_adj);
    h_vertex_colors = Kokkos::create_mirror_view(vertex_colors);
    Kokkos::deep_copy(h_color_adj, color_adj);
    Kokkos::deep_copy(h_vertex_colors, vertex_colors);
  } else {
    num_used_colors_steps = 1;
    num_multi_colors      = 1;
    num_colors            = 1;
    color_xadj            = nnz_lno_persistent_work_host_view_t("tmp", 2);
    color_xadj(0)         = 0;
    color_xadj(1)         = a_row_cnt;

    h_vertex_colors = nnz_lno_persistent_work_host_view_t("tmp", a_row_cnt);
    h_color_adj     = nnz_lno_persistent_work_host_view_t("tt", a_row_cnt);

    for (int i = 0; i < a_row_cnt; ++i) h_color_adj(i) = i;
    for (int i = 0; i < a_row_cnt; ++i) h_vertex_colors(i) = 0;
  }
  std::cout << "num_colors:" << num_colors << " num_multi_colors:" << num_multi_colors
            << " num_used_colors:" << num_used_colors_steps << std::endl;

  typename c_row_view_t::HostMirror h_c_rowmap = Kokkos::create_mirror_view(rowmapC);
  Kokkos::deep_copy(h_c_rowmap, rowmapC);

  /*
  isGPU = true;
  vectorlane = get_suggested_vector__size(b_row_cnt, entriesB.extent(0),
  KokkosKernels::Impl::Exec_CUDA); cache_line_size = 32 * 8; cache_size =
  cache_line_size * 7;
  */

  int cores[]        = {64, 64, 64, 64, 64, 64};
  int hyperthreads[] = {1, 2, 2, 4, 4, 4};
  int team_size[]    = {1, 1, 2, 4, 2, 1};
  for (int k = 5; k >= 0; --k) {
    int num_cores                = cores[k];
    int num_hyperthreads_in_core = hyperthreads[k];
    int hyper_threads_in_team    = team_size[k];
    this->read_write_cost(

        num_colors, num_multi_colors, num_used_colors_steps, isGPU, num_cores, num_hyperthreads_in_core,
        hyper_threads_in_team, vectorlane, cache_line_size, data_size, cache_l1_size,

        color_xadj, h_color_adj, h_vertex_colors,

        overall_flops, h_c_flop_rowmap, h_c_comp_a_net_index, h_c_comp_b_net_index, h_c_comp_row_index,
        h_c_comp_col_index, h_c_rowmap, write_type);
  }
}

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename c_row_view_t>
void KokkosSPGEMM<
    HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_, b_lno_nnz_view_t_,
    b_scalar_nnz_view_t_>::read_write_cost(nnz_lno_t num_colors, nnz_lno_t num_multi_colors,
                                           nnz_lno_t num_parallel_colors, bool isGPU, int num_cores,

                                           nnz_lno_t num_hyperthreads_in_core, nnz_lno_t hyper_threads_in_team,

                                           int vectorlane, const int cache_line_size, const int data_size,
                                           const int cache_size,

                                           nnz_lno_persistent_work_host_view_t color_xadj,
                                           typename nnz_lno_persistent_work_view_t::HostMirror color_adj,
                                           typename nnz_lno_persistent_work_view_t::HostMirror vertex_colors,

                                           size_t overall_flops,
                                           typename row_lno_temp_work_view_t::HostMirror c_flop_rowmap,
                                           typename row_lno_temp_work_view_t::HostMirror c_comp_a_net_index,
                                           typename row_lno_temp_work_view_t::HostMirror c_comp_b_net_index,
                                           typename nnz_lno_temp_work_view_t::HostMirror c_comp_row_index,
                                           typename nnz_lno_temp_work_view_t::HostMirror c_comp_col_index,
                                           c_row_view_t rowmapC,
                                           int write_type  // 0 -- KKMEM, 1-KKSPEED, 2- KKCOLOR 3-KKMULTICOLOR
                                                           // 4-KKMULTICOLOR2
) {
  std::cout << "num_colors:" << num_colors << " num_multi_colors:" << num_multi_colors
            << " num_parallel_colors:" << num_parallel_colors << std::endl;

  const int cache_1_line_word_count = cache_line_size / data_size;
  const int cache_1_line_count      = cache_size / cache_line_size;

  size_t a_line_size = entriesA.extent(0) / cache_line_size + 1;
  size_t b_line_size = entriesB.extent(0) / cache_line_size + 1;

  typename nnz_lno_temp_work_view_t::HostMirror tester("t", overall_flops);
  typename c_row_view_t::HostMirror h_rowmapC = Kokkos::create_mirror_view(rowmapC);
  Kokkos::deep_copy(h_rowmapC, rowmapC);

  size_t overall_a_l1_missread = 0;
  size_t overall_b_l1_missread = 0;
  size_t overall_a_l1_reuse    = 0;
  size_t overall_b_l1_reuse    = 0;

  size_t overall_c_l1_misswrite = 0;
  size_t overall_c_l1_reuse     = 0;

  // std::cout << "I:" << i << " color_begin:" << color_begin <<  " color_end:"
  // << color_end << " percore:" << percore << std::endl;
  int num_threads = 1;
#if defined(KOKKOS_ENABLE_OPENMP)
#pragma omp parallel
  { num_threads = omp_get_num_threads(); }
#endif

  // std::cout << "num_threads:" << num_threads << std::endl;
  const int num_teams_in_core = num_hyperthreads_in_core / hyper_threads_in_team;

  std::vector<std::vector<size_type> > t_team_worksize(num_threads);
  std::vector<std::vector<size_type> > t_team_begins(num_threads);
  std::vector<std::vector<size_type> > t_team_ends(num_threads);
  std::vector<std::vector<nnz_lno_t> > t_team_row_col_indices(num_threads);
  std::vector<Cache *> t_team_caches(num_threads);

  for (int tid = 0; tid < num_threads; ++tid) {
    t_team_worksize[tid].resize(num_teams_in_core);
    t_team_begins[tid].resize(num_teams_in_core);
    t_team_ends[tid].resize(num_teams_in_core);

    t_team_caches[tid] = new Cache(cache_1_line_count);

    t_team_row_col_indices[tid].resize(b_col_cnt, -1);
  }

  // std::cout << "num_parallel_colors:" << num_parallel_colors << std::endl;
  for (nnz_lno_t i = 0; i < num_parallel_colors; i += num_multi_colors) {
    nnz_lno_t color_upperbound = KOKKOSKERNELS_MACRO_MIN(num_parallel_colors, i + num_multi_colors);
    std::cout << "i:" << i << " color_upperbound:" << color_upperbound << " num_parallel_colors:" << num_parallel_colors
              << " num_multi_colors:" << num_multi_colors << std::endl;
    ;

    nnz_lno_t color_begin = color_xadj(i);
    nnz_lno_t color_end   = color_xadj(color_upperbound);
    // std::cout << "i:" << i << " color_begin:" << color_begin  << "
    // color_end:" << color_end << std::endl;
    nnz_lno_t percore = 32 / vectorlane;
    if (!isGPU)
      percore = (color_end - color_begin) / num_cores + 1;
    else {
      num_cores = a_row_cnt / (percore) + 1;
    }

#pragma omp parallel for
    for (nnz_lno_t j = color_begin; j < color_end; j += percore) {
      size_t amissreadcount_l1 = 0;
      size_t areusecount_l1    = 0;
      size_t bmissreadcount_l1 = 0;
      size_t breusecount_l1    = 0;

      size_t cmisswritecount_l1 = 0;
      size_t creusecount_l1     = 0;

      const size_t upperbound = KOKKOSKERNELS_MACRO_MIN(color_end, j + percore);

      const size_t num_work_for_core = upperbound - j;

      // std::vector <size_type> num_work_per_hyperthread (num_hyperthreads,
      // num_work_for_core / num_hyperthreads); std::vector <size_type>
      // team_index(num_hyperthreads);

      /*
       *
      Cache L1_cache (cache_1_line_count);
      std::vector <size_type> team_worksize(num_teams_in_core, num_work_for_core
      / num_teams_in_core); std::vector <size_type>
      team_begins(num_teams_in_core, 0); std::vector <size_type>
      team_ends(num_teams_in_core, 0);
      */
      int mytid = 0;
#if defined(KOKKOS_ENABLE_OPENMP)
      mytid = omp_get_thread_num();
#endif
      Cache L1_cache                        = *(t_team_caches[mytid]);
      std::vector<size_type> &team_worksize = t_team_worksize[mytid];
      std::vector<size_type> &team_begins   = t_team_begins[mytid];
      std::vector<size_type> &team_ends     = t_team_ends[mytid];
      L1_cache.clear_cache();

      std::vector<nnz_lno_t> &rows_col_indices = t_team_row_col_indices[mytid];

      for (int ws = 0; ws < num_teams_in_core; ++ws) {
        team_worksize[ws] = num_work_for_core / num_teams_in_core;
        team_begins[ws] = team_ends[ws] = 0;
      }
      // std::vector <bool> teams_are_done(num_teams_in_core, false);
      int not_all_done = num_teams_in_core;

      for (int ht = 0; ht < num_teams_in_core; ht++) {
        if (ht < int(num_work_for_core % num_teams_in_core)) {
          team_worksize[ht]++;
        }
      }

      team_begins[0] = j;
      team_ends[0]   = team_worksize[0] + j;

      for (int ht = 1; ht < num_teams_in_core; ht++) {
        team_begins[ht] = team_ends[ht - 1];
        team_ends[ht]   = team_begins[ht] + team_worksize[ht];

        if (team_begins[ht] == team_ends[ht]) --not_all_done;

        // std::cout << "ht:" << ht << " team_begins[ht]:" << team_begins[ht] <<
        // " team_ends[ht]:" << team_ends[ht] << " upperbound: " << upperbound<<
        // std::endl;
      }

      while (not_all_done) {
        int hyper_thread_ind = 0;
        for (int k = 0; k < num_teams_in_core; k++) {
          // std::cout << "k:"<< k << " team_begins[k]:" << team_begins[k] <<  "
          // team_ends[k]:" << team_ends[k] << " not_all_done:" << not_all_done
          // <<  std::endl;

          for (nnz_lno_t z = 0; z < hyper_threads_in_team; z++) {
            ++hyper_thread_ind;
            if (team_begins[k] < team_ends[k]) {
              nnz_lno_t zz = team_begins[k]++;
              if (!(team_begins[k] < team_ends[k])) {
                not_all_done--;
              }
              nnz_lno_t row                     = color_adj(zz);
              const size_t row_flop_begin_index = c_flop_rowmap(row);
              const size_t row_flop_end_index   = c_flop_rowmap(row + 1);

              nnz_lno_t num_cols = 0;

              for (size_t z = row_flop_begin_index; z < row_flop_end_index; ++z) {
                // std::cout << "z:" << z << std::endl;
                size_t areadind = c_comp_a_net_index(z) / cache_1_line_word_count;
                size_t breadind = c_comp_b_net_index(z) / cache_1_line_word_count;

                ++tester(z);
                if (L1_cache.is_in_cache('a', areadind, areadind % cache_1_line_count)) {
                  areusecount_l1++;
                } else {
                  amissreadcount_l1++;
                  L1_cache.insert_to_cache('a', areadind, areadind % cache_1_line_count);
                }
                if (L1_cache.is_in_cache('b', breadind, (a_line_size + breadind) % cache_1_line_count)) {
                  breusecount_l1++;
                } else {
                  bmissreadcount_l1++;
                  L1_cache.insert_to_cache('b', breadind, (a_line_size + breadind) % cache_1_line_count);
                }

                nnz_lno_t result_col = c_comp_col_index(z);
                // if (row == 0) std::cout << c_comp_col_index(z) << std::endl;

                switch (write_type) {
                  case 0: {
                    nnz_lno_t colpos = 0;
                    bool found       = false;
                    for (colpos = 0; colpos < num_cols; ++colpos) {
                      if (rows_col_indices[colpos] == result_col) {
                        found = true;
                        break;
                      }
                    }
                    if (!found) {
                      num_cols++;
                      rows_col_indices[colpos] = result_col;
                    }
                    if (L1_cache.is_in_cache(
                            'c', (h_rowmapC(row) + colpos) / cache_1_line_word_count,
                            (a_line_size + b_line_size + (h_rowmapC(row) + colpos) / cache_1_line_word_count) %
                                cache_1_line_count)) {
                      creusecount_l1++;
                    } else {
                      cmisswritecount_l1++;
                      L1_cache.insert_to_cache(
                          'c', (h_rowmapC(row) + colpos) / cache_1_line_word_count,
                          (a_line_size + b_line_size + (colpos + h_rowmapC(row)) / cache_1_line_word_count) %
                              cache_1_line_count);
                    }
                  } break;
                  case 2: hyper_thread_ind = 1; [[fallthrough]];
                  case 1: {
                    nnz_lno_t result_col_ind = result_col + (hyper_thread_ind - 1) * b_col_cnt;
                    if (L1_cache.is_in_cache('c', (result_col_ind) / cache_1_line_word_count,
                                             (a_line_size + b_line_size + (result_col_ind) / cache_1_line_word_count) %
                                                 cache_1_line_count)) {
                      /*
                      if (L1_cache.cache_inserting_hyperthread[(
                          a_line_size + b_line_size + (result_col_ind)  /
                      cache_1_line_word_count) % cache_1_line_count] !=
                      hyper_thread_ind){ std::cout << "hyper_thread_ind:" <<
                      hyper_thread_ind
                                  << " result_col:" << result_col
                                  << " result_col_ind:" << result_col_ind
                                  << " is in the cache inserted by:" <<
                      L1_cache.cache_inserting_hyperthread[( a_line_size +
                      b_line_size + (result_col_ind)  / cache_1_line_word_count)
                      % cache_1_line_count] << std::endl;
                      }

                      */
                      creusecount_l1++;
                    } else {
                      cmisswritecount_l1++;
                      L1_cache.insert_to_cache(
                          'c', (result_col_ind) / cache_1_line_word_count,
                          (a_line_size + b_line_size + (result_col_ind) / cache_1_line_word_count) % cache_1_line_count,
                          hyper_thread_ind);
                    }
                  } break;
                  case 3:
                  case 4: {
                    nnz_lno_t row_color = vertex_colors(row);
                    // if (row_color < 0 || row_color >= num_multi_colors)
                    // std::cout << "row:" << row << " rowcol:" << row_color <<
                    // std::endl;
                    result_col = result_col + row_color * b_col_cnt;
                    if (L1_cache.is_in_cache('c', (result_col) / cache_1_line_word_count,
                                             (a_line_size + b_line_size + (result_col) / cache_1_line_word_count) %
                                                 cache_1_line_count)) {
                      creusecount_l1++;
                    } else {
                      cmisswritecount_l1++;
                      L1_cache.insert_to_cache(
                          'c', (result_col) / cache_1_line_word_count,
                          (a_line_size + b_line_size + (result_col) / cache_1_line_word_count) % cache_1_line_count);
                    }
                  } break;
                }
              }

              for (size_type z = 0; z < num_cols; ++z) {
                rows_col_indices[z] = -1;
              }

              if (write_type != 0) {
                cmisswritecount_l1 += ceil((h_rowmapC(row + 1) - h_rowmapC(row)) / double(cache_1_line_word_count));
                creusecount_l1 += row_flop_end_index - row_flop_begin_index -
                                  ceil((h_rowmapC(row + 1) - h_rowmapC(row)) / double(cache_1_line_word_count));
              }
            }
          }
        }
      }

      Kokkos::atomic_fetch_add(&overall_c_l1_misswrite, cmisswritecount_l1);
      // overall_c_l1_misswrite += cmisswritecount_l1;
      Kokkos::atomic_fetch_add(&overall_c_l1_reuse, creusecount_l1);
      // overall_c_l1_reuse += creusecount_l1;
      Kokkos::atomic_fetch_add(&overall_a_l1_missread, amissreadcount_l1);
      // overall_a_l1_missread += amissreadcount_l1;
      Kokkos::atomic_fetch_add(&overall_b_l1_missread, bmissreadcount_l1);
      // overall_b_l1_missread += bmissreadcount_l1;
      Kokkos::atomic_fetch_add(&overall_a_l1_reuse, areusecount_l1);
      // overall_a_l1_reuse += areusecount_l1;
      Kokkos::atomic_fetch_add(&overall_b_l1_reuse, breusecount_l1);
      // overall_b_l1_reuse += breusecount_l1;
    }
  }

  for (int tid = 0; tid < num_threads; ++tid) {
    delete t_team_caches[tid];
  }
  // std::cout << "write_type:" << write_type << std::endl;
  std::string algo = "KKMEM";
  switch (write_type) {
    case 0: algo = "KKMEM"; break;
    case 1: algo = "KKSPEED"; break;
    case 2: algo = "KKCOLOR"; break;
    case 3: algo = "KKMULTICOLOR"; break;
    case 4: algo = "KKMULTICOLOR2"; break;
  }
  // 0 -- KKMEM, 1-KKSPEED, 2- KKCOLOR 3-KKMULTICOLOR 4-KKMULTICOLOR2
  std::cout << algo << " numthreads:" << num_cores * num_hyperthreads_in_core << " teamsize:" << hyper_threads_in_team
            << " overall_flops: " << overall_flops << " flops per row: " << overall_flops / double(a_row_cnt)
            << " a_read_perrow:" << overall_a_l1_missread / double(a_row_cnt)
            << " b_read_perrow:" << overall_b_l1_missread / double(a_row_cnt)
            << " reads per row: " << (overall_a_l1_missread + overall_b_l1_missread) / double(a_row_cnt)
            << " writes per row: " << (overall_c_l1_misswrite) / double(a_row_cnt) << std::endl;

  std::cout << algo << " numthreads:" << num_cores * num_hyperthreads_in_core << " teamsize:" << hyper_threads_in_team
            << " overall_flops: " << overall_flops
            << " a_read_perflop:" << overall_a_l1_missread / double(overall_flops)
            << " b_read_perflop:" << overall_b_l1_missread / double(overall_flops)
            << " total_read_per_flop:" << (overall_a_l1_missread + overall_b_l1_missread) / double(overall_flops)
            << " write per flop:" << overall_c_l1_misswrite / double(overall_flops) << std::endl;

  std::cout << "\n\ta_read_pernnz:" << overall_a_l1_missread / double(rowmapC(a_row_cnt))
            << "\n\tb_read_pernnz:" << overall_b_l1_missread / double(rowmapC(a_row_cnt))
            << "\n\ttotal_read_per_nnz:" << (overall_a_l1_missread + overall_b_l1_missread) / double(rowmapC(a_row_cnt))
            << "\n\twrite per nnz:" << overall_c_l1_misswrite / double(rowmapC(a_row_cnt)) << std::endl;

  std::cout << "\topt_a_read:" << entriesA.extent(0) / cache_1_line_word_count
            << " opt_b_read:" << entriesB.extent(0) / cache_1_line_word_count << " opt_total_read:"
            << entriesA.extent(0) / cache_1_line_word_count + entriesB.extent(0) / cache_1_line_word_count << std::endl;

  std::cout << "\taverage_per_thread a_read:" << overall_a_l1_missread / num_cores
            << " average_per_thread b_read:" << overall_b_l1_missread / num_cores
            << " average_per_thread total_read:" << (overall_a_l1_missread + overall_b_l1_missread) / num_cores
            << std::endl;

  std::cout << "\toverall_a_read:" << overall_a_l1_missread << " overall_b_read:" << overall_b_l1_missread
            << " total_read:" << overall_a_l1_missread + overall_b_l1_missread << std::endl;

  std::cout << "\toverall_a_reuse:" << overall_a_l1_reuse << " overall_b_reuse:" << overall_b_l1_reuse
            << " total_reuse:" << overall_a_l1_reuse + overall_b_l1_reuse << std::endl;

  std::cout << "\toverall_c_write:" << overall_c_l1_misswrite << " overall_c_reuse:" << overall_c_l1_reuse
            << " opt_c_write:" << double(rowmapC(a_row_cnt)) / cache_1_line_word_count << std::endl;
}
#endif

}  // namespace Impl
}  // namespace KokkosSparse
