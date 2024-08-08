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

#pragma once
// exclude from Cuda builds without lambdas enabled
#if !defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_CUDA_LAMBDA)
#include <list>
#include <limits>
#include <Kokkos_Core.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spgemm.hpp"
#include "KokkosSparse_SortCrs.hpp"
#include "KokkosKernels_HashmapAccumulator.hpp"
#include "KokkosKernels_Uniform_Initialized_MemoryPool.hpp"
#include "KokkosGraph_CoarsenHeuristics.hpp"

namespace KokkosSparse {

namespace Impl {

template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t>
struct SortLowDegreeCrsMatrixFunctor {
  using size_type  = typename rowmap_t::non_const_value_type;
  using lno_t      = typename entries_t::non_const_value_type;
  using scalar_t   = typename values_t::non_const_value_type;
  using team_mem   = typename Kokkos::TeamPolicy<execution_space>::member_type;
  using value_type = lno_t;

  SortLowDegreeCrsMatrixFunctor(bool usingRangePol, const rowmap_t& _rowmap, const entries_t& _entries,
                                const values_t& _values, const lno_t _degreeLimit)
      : rowmap(_rowmap), entries(_entries), values(_values), degreeLimit(_degreeLimit) {
    if (usingRangePol) {
      entriesAux = entries_t(Kokkos::ViewAllocateWithoutInitializing("Entries aux"), entries.extent(0));
      valuesAux  = values_t(Kokkos::ViewAllocateWithoutInitializing("Values aux"), values.extent(0));
    }
    // otherwise, aux arrays won't be allocated (sorting in place)
  }

  KOKKOS_INLINE_FUNCTION void operator()(const lno_t i, value_type& reducer) const {
    size_type rowStart = rowmap(i);
    size_type rowEnd   = rowmap(i + 1);
    lno_t rowNum       = rowEnd - rowStart;
    if (rowNum > degreeLimit) {
      reducer++;
      return;
    }
    // Radix sort requires unsigned keys for comparison
    using unsigned_lno_t = typename std::make_unsigned<lno_t>::type;
    KokkosKernels::SerialRadixSort2<lno_t, unsigned_lno_t, scalar_t>(
        (unsigned_lno_t*)entries.data() + rowStart, (unsigned_lno_t*)entriesAux.data() + rowStart,
        values.data() + rowStart, valuesAux.data() + rowStart, rowNum);
  }

  KOKKOS_INLINE_FUNCTION void operator()(const team_mem t, value_type& reducer) const {
    size_type i        = t.league_rank();
    size_type rowStart = rowmap(i);
    size_type rowEnd   = rowmap(i + 1);
    lno_t rowNum       = rowEnd - rowStart;
    if (rowNum > degreeLimit) {
      Kokkos::single(Kokkos::PerTeam(t), [&]() { reducer++; });
      return;
    }
    KokkosKernels::TeamBitonicSort2<lno_t, lno_t, scalar_t, team_mem>(entries.data() + rowStart,
                                                                      values.data() + rowStart, rowNum, t);
  }

  rowmap_t rowmap;
  entries_t entries;
  entries_t entriesAux;
  values_t values;
  values_t valuesAux;
  lno_t degreeLimit;
};

}  // namespace Impl

// Sort a CRS matrix: within each row, sort entries ascending by column.
// At the same time, permute the values.
// Only modifies rows below the degreeLimit
template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t>
typename entries_t::non_const_value_type sort_low_degree_rows_crs_matrix(
    const rowmap_t& rowmap, const entries_t& entries, const values_t& values,
    const typename entries_t::non_const_value_type degreeLimit) {
  using lno_t    = typename entries_t::non_const_value_type;
  using team_pol = Kokkos::TeamPolicy<execution_space>;
  bool useRadix  = !KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>();
  Impl::SortLowDegreeCrsMatrixFunctor<execution_space, rowmap_t, entries_t, values_t> funct(useRadix, rowmap, entries,
                                                                                            values, degreeLimit);
  lno_t numRows   = rowmap.extent(0) ? rowmap.extent(0) - 1 : 0;
  lno_t notSorted = 0;
  if (useRadix) {
    Kokkos::parallel_reduce("sort_crs_matrix", Kokkos::RangePolicy<execution_space>(0, numRows), funct, notSorted);
  } else {
    // Try to get teamsize to be largest power of 2 not greater than avg entries
    // per row
    // TODO (probably important for performnce): add thread-level sort also, and
    // use that for small avg degree. But this works for now. probably important
    // for this particular use case of only low-degree rows
    int teamSize = 1;
    lno_t avgDeg = 0;
    if (numRows) avgDeg = (entries.extent(0) + numRows - 1) / numRows;
    if (avgDeg > degreeLimit) {
      avgDeg = degreeLimit;
    }
    while (teamSize * 2 * 2 <= avgDeg) {
      teamSize *= 2;
    }
    team_pol temp(numRows, teamSize);
    teamSize = std::min(teamSize, temp.team_size_max(funct, Kokkos::ParallelReduceTag()));
    Kokkos::parallel_reduce("sort_crs_matrix", team_pol(numRows, teamSize), funct, notSorted);
  }
  return notSorted;
}

}  // namespace KokkosSparse

namespace KokkosGraph {

namespace Experimental {

// this class is not meant to be instantiated
// think of it like a templated namespace
template <class crsMat>
class coarse_builder {
 public:
  // define internal types
  using matrix_t              = crsMat;
  using exec_space            = typename matrix_t::execution_space;
  using mem_space             = typename matrix_t::memory_space;
  using Device                = typename matrix_t::device_type;
  using ordinal_t             = typename matrix_t::ordinal_type;
  using edge_offset_t         = typename matrix_t::size_type;
  using scalar_t              = typename matrix_t::value_type;
  using vtx_view_t            = Kokkos::View<ordinal_t*, Device>;
  using wgt_view_t            = Kokkos::View<scalar_t*, Device>;
  using edge_view_t           = Kokkos::View<edge_offset_t*, Device>;
  using edge_subview_t        = Kokkos::View<edge_offset_t, Device>;
  using graph_type            = typename matrix_t::staticcrsgraph_type;
  using policy_t              = Kokkos::RangePolicy<exec_space>;
  using dyn_policy_t          = Kokkos::RangePolicy<Kokkos::Schedule<Kokkos::Dynamic>, exec_space>;
  using team_policy_t         = Kokkos::TeamPolicy<exec_space>;
  using dyn_team_policy_t     = Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Dynamic>, exec_space>;
  using member                = typename team_policy_t::member_type;
  using spgemm_kernel_handle  = KokkosKernels::Experimental::KokkosKernelsHandle<edge_offset_t, ordinal_t, scalar_t,
                                                                                exec_space, mem_space, mem_space>;
  using uniform_memory_pool_t = KokkosKernels::Impl::UniformMemoryPool<exec_space, ordinal_t>;
  using mapper_t              = coarsen_heuristics<crsMat>;
  static constexpr ordinal_t get_null_val() {
    // this value must line up with the null value used by the hashmap
    // accumulator
    if (std::is_signed<ordinal_t>::value) {
      return -1;
    } else {
      return std::numeric_limits<ordinal_t>::max();
    }
  }
  static constexpr ordinal_t ORD_MAX = get_null_val();
  static constexpr bool is_host_space =
      std::is_same<typename exec_space::memory_space, typename Kokkos::DefaultHostExecutionSpace::memory_space>::value;
  static constexpr bool scal_eq_ord = std::is_same<ordinal_t, scalar_t>::value;
  // contains matrix and vertex weights corresponding to current level
  // interp matrix maps previous level to this level
  struct coarse_level_triple {
    matrix_t mtx;
    vtx_view_t vtx_wgts;
    matrix_t interp_mtx;
    int level;
    bool uniform_weights;
  };

  // define behavior-controlling enums
  enum Heuristic { HECv1, Match, MtMetis, MIS2, GOSHv1, GOSHv2 };
  enum Builder { Sort, Hashmap, Hybrid, Spgemm, Spgemm_transpose_first };

  struct coarsen_handle {
    // internal parameters and data
    // default heuristic is HEC
    Heuristic h = HECv1;
    // default builder is Hybrid
    Builder b = Hybrid;
    std::list<coarse_level_triple> results;
    ordinal_t coarse_vtx_cutoff = 50;
    ordinal_t min_allowed_vtx   = 10;
    unsigned int max_levels     = 200;
    size_t max_mem_allowed      = 536870912;
  };

  // determine if dynamic scheduling should be used
  static bool should_use_dyn(const ordinal_t n, const Kokkos::View<const edge_offset_t*, Device> work, int t_count) {
    bool use_dyn      = false;
    edge_offset_t max = 0;
    edge_offset_t min = std::numeric_limits<edge_offset_t>::max();
    if (is_host_space) {
      ordinal_t static_size = (n + t_count) / t_count;
      for (ordinal_t i = 0; i < t_count; i++) {
        ordinal_t start = i * static_size;
        ordinal_t end   = start + static_size;
        if (start > n) start = n;
        if (end > n) end = n;
        edge_offset_t size = work(end) - work(start);
        if (size > max) {
          max = size;
        }
        if (size < min) {
          min = size;
        }
      }
      if (n > 500000 && max > 5 * min) {
        use_dyn = true;
      }
    }
    return use_dyn;
  }

  // build the course graph according to ((B^T A) B) or (B^T (A B)), where B is
  // aggregator matrix
  static coarse_level_triple build_coarse_graph_spgemm(coarsen_handle& handle, const coarse_level_triple level,
                                                       const matrix_t interp_mtx) {
    vtx_view_t f_vtx_w = level.vtx_wgts;
    matrix_t g         = level.mtx;
    if (!KokkosSparse::Impl::isCrsGraphSorted(g.graph.row_map, g.graph.entries)) KokkosSparse::sort_crs_matrix(g);

    ordinal_t n  = g.numRows();
    ordinal_t nc = interp_mtx.numCols();

    matrix_t interp_transpose = KokkosSparse::Impl::transpose_matrix(interp_mtx);
    KokkosSparse::sort_crs_matrix(interp_transpose);

    spgemm_kernel_handle kh;
    kh.set_team_work_size(64);
    kh.set_dynamic_scheduling(true);

    vtx_view_t adj_coarse;
    wgt_view_t wgt_coarse;
    edge_view_t row_map_coarse;

    if (handle.b == Spgemm_transpose_first) {
      kh.create_spgemm_handle();
      edge_view_t row_map_p1("rows_partial", nc + 1);
      KokkosSparse::Experimental::spgemm_symbolic(&kh, nc, n, n, interp_transpose.graph.row_map,
                                                  interp_transpose.graph.entries, false, g.graph.row_map,
                                                  g.graph.entries, false, row_map_p1);

      // partial-result matrix
      vtx_view_t entries_p1("adjacencies_partial", kh.get_spgemm_handle()->get_c_nnz());
      wgt_view_t values_p1("weights_partial", kh.get_spgemm_handle()->get_c_nnz());

      KokkosSparse::Experimental::spgemm_numeric(
          &kh, nc, n, n, interp_transpose.graph.row_map, interp_transpose.graph.entries, interp_transpose.values, false,
          g.graph.row_map, g.graph.entries, g.values, false, row_map_p1, entries_p1, values_p1);
      kh.destroy_spgemm_handle();

      row_map_coarse = edge_view_t("rows_coarse", nc + 1);
      kh.create_spgemm_handle();
      KokkosSparse::Experimental::spgemm_symbolic(&kh, nc, n, nc, row_map_p1, entries_p1, false,
                                                  interp_mtx.graph.row_map, interp_mtx.graph.entries, false,
                                                  row_map_coarse);
      // coarse-graph adjacency matrix
      adj_coarse = vtx_view_t("adjacencies_coarse", kh.get_spgemm_handle()->get_c_nnz());
      wgt_coarse = wgt_view_t("weights_coarse", kh.get_spgemm_handle()->get_c_nnz());

      KokkosSparse::Experimental::spgemm_numeric(&kh, nc, n, nc, row_map_p1, entries_p1, values_p1, false,
                                                 interp_mtx.graph.row_map, interp_mtx.graph.entries, interp_mtx.values,
                                                 false, row_map_coarse, adj_coarse, wgt_coarse);
      kh.destroy_spgemm_handle();
    } else {
      edge_view_t row_map_p1("rows_partial", n + 1);
      kh.create_spgemm_handle();
      KokkosSparse::Experimental::spgemm_symbolic(&kh, n, n, nc, g.graph.row_map, g.graph.entries, false,
                                                  interp_mtx.graph.row_map, interp_mtx.graph.entries, false,
                                                  row_map_p1);

      // partial-result matrix
      vtx_view_t entries_p1("adjacencies_partial", kh.get_spgemm_handle()->get_c_nnz());
      wgt_view_t values_p1("weights_partial", kh.get_spgemm_handle()->get_c_nnz());

      KokkosSparse::Experimental::spgemm_numeric(&kh, n, n, nc, g.graph.row_map, g.graph.entries, g.values, false,
                                                 interp_mtx.graph.row_map, interp_mtx.graph.entries, interp_mtx.values,
                                                 false, row_map_p1, entries_p1, values_p1);
      kh.destroy_spgemm_handle();

      row_map_coarse = edge_view_t("rows_coarse", nc + 1);
      kh.create_spgemm_handle();
      KokkosSparse::Experimental::spgemm_symbolic(&kh, nc, n, nc, interp_transpose.graph.row_map,
                                                  interp_transpose.graph.entries, false, row_map_p1, entries_p1, false,
                                                  row_map_coarse);
      // coarse-graph adjacency matrix
      adj_coarse = vtx_view_t("adjacencies_coarse", kh.get_spgemm_handle()->get_c_nnz());
      wgt_coarse = wgt_view_t("weights_coarse", kh.get_spgemm_handle()->get_c_nnz());

      KokkosSparse::Experimental::spgemm_numeric(
          &kh, nc, n, nc, interp_transpose.graph.row_map, interp_transpose.graph.entries, interp_transpose.values,
          false, row_map_p1, entries_p1, values_p1, false, row_map_coarse, adj_coarse, wgt_coarse);
      kh.destroy_spgemm_handle();
    }

    // now we must remove self-loop edges
    edge_view_t nonLoops("nonLoop", nc);

    // gonna reuse this to count non-self loop edges
    Kokkos::parallel_for(
        policy_t(0, nc), KOKKOS_LAMBDA(ordinal_t i) { nonLoops(i) = 0; });

    Kokkos::parallel_for(
        policy_t(0, nc), KOKKOS_LAMBDA(ordinal_t u) {
          for (edge_offset_t j = row_map_coarse(u); j < row_map_coarse(u + 1); j++) {
            if (adj_coarse(j) != u) {
              nonLoops(u)++;
            }
          }
        });

    edge_view_t row_map_nonloop("nonloop row map", nc + 1);

    Kokkos::parallel_scan(
        policy_t(0, nc), KOKKOS_LAMBDA(const ordinal_t i, edge_offset_t& update, const bool final) {
          const edge_offset_t val_i = nonLoops(i);
          update += val_i;
          if (final) {
            row_map_nonloop(i + 1) = update;
          }
        });

    edge_subview_t rmn_subview = Kokkos::subview(row_map_nonloop, nc);
    edge_offset_t rmn          = 0;
    Kokkos::deep_copy(rmn, rmn_subview);

    vtx_view_t entries_nonloop("nonloop entries", rmn);
    wgt_view_t values_nonloop("nonloop values", rmn);

    Kokkos::parallel_for(
        policy_t(0, nc), KOKKOS_LAMBDA(const ordinal_t i) { nonLoops(i) = 0; });

    Kokkos::parallel_for(
        policy_t(0, nc), KOKKOS_LAMBDA(const ordinal_t u) {
          for (edge_offset_t j = row_map_coarse(u); j < row_map_coarse(u + 1); j++) {
            if (adj_coarse(j) != u) {
              edge_offset_t offset    = row_map_nonloop(u) + nonLoops(u)++;
              entries_nonloop(offset) = adj_coarse(j);
              values_nonloop(offset)  = wgt_coarse(j);
            }
          }
        });
    // done removing self-loop edges

    kh.destroy_spgemm_handle();

    graph_type gc_graph(entries_nonloop, row_map_nonloop);
    matrix_t gc("gc", nc, values_nonloop, gc_graph);

    vtx_view_t c_vtx_w("coarse vtx weights", interp_mtx.numCols());
    Kokkos::parallel_for(
        "compute coarse vtx wgts", policy_t(0, n), KOKKOS_LAMBDA(const ordinal_t i) {
          ordinal_t u = interp_mtx.graph.entries(i);
          Kokkos::atomic_add(&c_vtx_w(u), f_vtx_w(i));
        });

    coarse_level_triple next_level;
    next_level.mtx             = gc;
    next_level.vtx_wgts        = c_vtx_w;
    next_level.level           = level.level + 1;
    next_level.interp_mtx      = interp_mtx;
    next_level.uniform_weights = false;
    return next_level;
  }

  struct prefix_sum {
    vtx_view_t input;
    edge_view_t output;

    prefix_sum(vtx_view_t _input, edge_view_t _output) : input(_input), output(_output) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const ordinal_t i, edge_offset_t& update, const bool final) const {
      const edge_offset_t val_i = input(i);
      update += val_i;
      if (final) {
        output(i + 1) = update;
      }
    }
  };

  struct functorDedupeLowDegreeAfterSort {
    // compiler may get confused what the reduction type is without this
    typedef edge_offset_t value_type;

    edge_view_t row_map;
    vtx_view_t entries, entriesOut;
    wgt_view_t wgts, wgtsOut;
    vtx_view_t dedupe_edge_count;
    ordinal_t degreeLimit;

    functorDedupeLowDegreeAfterSort(edge_view_t _row_map, vtx_view_t _entries, vtx_view_t _entriesOut, wgt_view_t _wgts,
                                    wgt_view_t _wgtsOut, vtx_view_t _dedupe_edge_count, ordinal_t _degreeLimit_)
        : row_map(_row_map),
          entries(_entries),
          entriesOut(_entriesOut),
          wgts(_wgts),
          wgtsOut(_wgtsOut),
          dedupe_edge_count(_dedupe_edge_count),
          degreeLimit(_degreeLimit_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member& thread, edge_offset_t& thread_sum) const {
      ordinal_t u         = thread.league_rank();
      edge_offset_t start = row_map(u);
      edge_offset_t end   = row_map(u + 1);
      ordinal_t degree    = end - start;
      if (degree > degreeLimit) {
        return;
      }
      Kokkos::parallel_scan(Kokkos::TeamThreadRange(thread, start, end),
                            [&](const edge_offset_t& i, edge_offset_t& update, const bool final) {
                              if (i == start) {
                                update += 1;
                              } else if (entries(i) != entries(i - 1)) {
                                update += 1;
                              }
                              if (final) {
                                entriesOut(start + update - 1) = entries(i);
                                // requires that wgtsOut be initialized to 0
                                Kokkos::atomic_add(&wgtsOut(start + update - 1), wgts(i));
                                if (i + 1 == end) {
                                  dedupe_edge_count(u) = update;
                                }
                              }
                            });
      Kokkos::parallel_for(Kokkos::TeamThreadRange(thread, start, start + dedupe_edge_count(u)),
                           [&](const edge_offset_t& i) {
                             entries(i) = entriesOut(i);
                             wgts(i)    = wgtsOut(i);
                           });
      Kokkos::single(Kokkos::PerTeam(thread), [&]() { thread_sum += dedupe_edge_count(u); });
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const ordinal_t& u, edge_offset_t& thread_sum) const {
      ordinal_t offset = row_map(u);
      ordinal_t last   = ORD_MAX;
      ordinal_t degree = row_map(u + 1) - row_map(u);
      if (degree > degreeLimit) {
        return;
      }
      for (edge_offset_t i = row_map(u); i < row_map(u + 1); i++) {
        if (last != entries(i)) {
          entriesOut(offset) = entries(i);
          wgtsOut(offset)    = wgts(i);
          last               = entries(offset);
          offset++;
        } else {
          wgtsOut(offset - 1) += wgts(i);
        }
      }
      dedupe_edge_count(u) = offset - row_map(u);
      thread_sum += offset - row_map(u);
    }
  };

  struct functorDedupeAfterSort {
    // compiler may get confused what the reduction type is without this
    typedef edge_offset_t value_type;

    edge_view_t row_map;
    vtx_view_t entries, entriesOut;
    wgt_view_t wgts, wgtsOut;
    vtx_view_t dedupe_edge_count;

    functorDedupeAfterSort(edge_view_t _row_map, vtx_view_t _entries, vtx_view_t _entriesOut, wgt_view_t _wgts,
                           wgt_view_t _wgtsOut, vtx_view_t _dedupe_edge_count)
        : row_map(_row_map),
          entries(_entries),
          entriesOut(_entriesOut),
          wgts(_wgts),
          wgtsOut(_wgtsOut),
          dedupe_edge_count(_dedupe_edge_count) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member& thread, edge_offset_t& thread_sum) const {
      ordinal_t u         = thread.league_rank();
      edge_offset_t start = row_map(u);
      edge_offset_t end   = row_map(u + 1);
      Kokkos::parallel_scan(Kokkos::TeamThreadRange(thread, start, end),
                            [&](const edge_offset_t& i, edge_offset_t& update, const bool final) {
                              if (i == start) {
                                update += 1;
                              } else if (entries(i) != entries(i - 1)) {
                                update += 1;
                              }
                              if (final) {
                                entriesOut(start + update - 1) = entries(i);
                                // requires that wgtsOut be initialized to 0
                                Kokkos::atomic_add(&wgtsOut(start + update - 1), wgts(i));
                                if (i + 1 == end) {
                                  dedupe_edge_count(u) = update;
                                }
                              }
                            });
      Kokkos::single(Kokkos::PerTeam(thread), [&]() { thread_sum += dedupe_edge_count(u); });
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const ordinal_t& u, edge_offset_t& thread_sum) const {
      ordinal_t offset = row_map(u);
      ordinal_t last   = ORD_MAX;
      for (edge_offset_t i = row_map(u); i < row_map(u + 1); i++) {
        if (last != entries(i)) {
          entries(offset) = entries(i);
          wgtsOut(offset) = wgts(i);
          last            = entries(offset);
          offset++;
        } else {
          wgtsOut(offset - 1) += wgts(i);
        }
      }
      dedupe_edge_count(u) = offset - row_map(u);
      thread_sum += offset - row_map(u);
    }
  };

  struct functorCollapseDirectedToUndirected {
    const edge_view_t source_row_map;
    const edge_view_t target_row_map;
    const vtx_view_t source_edge_counts;
    vtx_view_t target_edge_counts;
    const vtx_view_t source_destinations;
    vtx_view_t target_destinations;
    const wgt_view_t source_wgts;
    wgt_view_t target_wgts;

    functorCollapseDirectedToUndirected(const edge_view_t _source_row_map, const edge_view_t _target_row_map,
                                        const vtx_view_t _source_edge_counts, vtx_view_t _target_edge_counts,
                                        const vtx_view_t _source_destinations, vtx_view_t _target_destinations,
                                        const wgt_view_t _source_wgts, wgt_view_t _target_wgts)
        : source_row_map(_source_row_map),
          target_row_map(_target_row_map),
          source_edge_counts(_source_edge_counts),
          target_edge_counts(_target_edge_counts),
          source_destinations(_source_destinations),
          target_destinations(_target_destinations),
          source_wgts(_source_wgts),
          target_wgts(_target_wgts) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member& thread) const {
      ordinal_t u                 = thread.league_rank();
      edge_offset_t u_origin      = source_row_map(u);
      edge_offset_t u_dest_offset = target_row_map(u);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(thread, source_edge_counts(u)), [&](const edge_offset_t u_idx) {
        ordinal_t v                 = source_destinations(u_origin + u_idx);
        scalar_t wgt                = source_wgts(u_origin + u_idx);
        edge_offset_t v_dest_offset = target_row_map(v);
        edge_offset_t v_dest        = v_dest_offset + Kokkos::atomic_fetch_add(&target_edge_counts(v), 1);
        edge_offset_t u_dest        = u_dest_offset + Kokkos::atomic_fetch_add(&target_edge_counts(u), 1);

        target_destinations(u_dest) = v;
        target_wgts(u_dest)         = wgt;
        target_destinations(v_dest) = u;
        target_wgts(v_dest)         = wgt;
      });
    }
  };

  struct functorHashmapAccumulator {
    // compiler may get confused what the reduction type is without this
    typedef ordinal_t value_type;
    edge_view_t row_map;
    vtx_view_t entries_in, entries_out;
    wgt_view_t wgts_in, wgts_out;
    vtx_view_t dedupe_edge_count;
    uniform_memory_pool_t memory_pool;
    const ordinal_t hash_size;
    const ordinal_t max_hash_entries;
    vtx_view_t remaining;
    bool use_out;

    functorHashmapAccumulator(edge_view_t _row_map, vtx_view_t _entries_in, vtx_view_t _entries_out,
                              wgt_view_t _wgts_in, wgt_view_t _wgts_out, vtx_view_t _dedupe_edge_count,
                              uniform_memory_pool_t _memory_pool, const ordinal_t _hash_size,
                              const ordinal_t _max_hash_entries, vtx_view_t _remaining, bool _use_out)
        : row_map(_row_map),
          entries_in(_entries_in),
          entries_out(_entries_out),
          wgts_in(_wgts_in),
          wgts_out(_wgts_out),
          dedupe_edge_count(_dedupe_edge_count),
          memory_pool(_memory_pool),
          hash_size(_hash_size),
          max_hash_entries(_max_hash_entries),
          remaining(_remaining),
          use_out(_use_out) {}

    KOKKOS_INLINE_FUNCTION
    ordinal_t get_thread_id(const ordinal_t row_index) const {
#if defined(KOKKOS_ENABLE_SERIAL)
      if (std::is_same<exec_space, Kokkos::Serial>::value) return 0;
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
      if (std::is_same<exec_space, Kokkos::OpenMP>::value) return Kokkos::OpenMP::impl_hardware_thread_id();
#endif
#if defined(KOKKOS_ENABLE_THREADS)
      if (std::is_same<exec_space, Kokkos::Threads>::value) return Kokkos::Threads::impl_hardware_thread_id();
#endif
      return row_index;
    }

    // reduces to find total number of rows that were too large
    KOKKOS_INLINE_FUNCTION
    void operator()(const ordinal_t& rem_idx, ordinal_t& thread_sum) const {
      ordinal_t idx = remaining(rem_idx);
      typedef ordinal_t hash_size_type;
      typedef ordinal_t hash_key_type;
      typedef scalar_t hash_value_type;

      // can't do this row at current hashmap size
      ordinal_t hash_entries = row_map(idx + 1) - row_map(idx);
      if (hash_entries >= max_hash_entries) {
        thread_sum++;
        return;
      }
      volatile ordinal_t* ptr_temp = nullptr;
      ordinal_t t_id               = rem_idx;
      // need to use the hardware thread id if the pool type is
      // OneThread2OneChunk
      t_id = get_thread_id(t_id);
      while (nullptr == ptr_temp) {
        ptr_temp = (volatile ordinal_t*)(memory_pool.allocate_chunk(t_id));
      }
      if (ptr_temp == nullptr) {
        return;
      }
      ordinal_t* ptr_memory_pool_chunk = (ordinal_t*)(ptr_temp);

      // These are updated by Hashmap_Accumulator insert functions.
      ordinal_t* used_hash_size = (ordinal_t*)(ptr_temp);
      ptr_temp++;
      ordinal_t* used_hash_count = (ordinal_t*)(ptr_temp);
      ptr_temp++;
      *used_hash_size  = 0;
      *used_hash_count = 0;

      // hash function is hash_size-1 (note: hash_size must be a power of 2)
      ordinal_t hash_func_pow2 = hash_size - 1;

      // Set pointer to hash indices
      ordinal_t* used_hash_indices = (ordinal_t*)(ptr_temp);
      ptr_temp += hash_size;

      // Set pointer to hash begins
      ordinal_t* hash_begins = (ordinal_t*)(ptr_temp);
      ptr_temp += hash_size;

      // Set pointer to hash nexts
      ordinal_t* hash_nexts = (ordinal_t*)(ptr_temp);

      // Set pointer to hash keys
      ordinal_t* keys = (ordinal_t*)entries_out.data() + row_map(idx);

      // Set pointer to hash values
      scalar_t* values = (scalar_t*)wgts_out.data() + row_map(idx);

      KokkosKernels::Experimental::HashmapAccumulator<hash_size_type, hash_key_type, hash_value_type,
                                                      KokkosKernels::Experimental::HashOpType::bitwiseAnd>
          hash_map(hash_size, hash_func_pow2, hash_begins, hash_nexts, keys, values);

      for (edge_offset_t i = row_map(idx); i < row_map(idx + 1); i++) {
        ordinal_t key  = entries_in(i);
        scalar_t value = wgts_in(i);
        hash_map.sequential_insert_into_hash_mergeAdd_TrackHashes(key, value, used_hash_size, used_hash_count,
                                                                  used_hash_indices);
      };

      // Reset the Begins values to -1 before releasing the memory pool chunk.
      // If you don't do this the next thread that grabs this memory chunk will
      // not work properly.
      for (ordinal_t i = 0; i < *used_hash_count; i++) {
        ordinal_t dirty_hash = used_hash_indices[i];

        hash_map.hash_begins[dirty_hash] = ORD_MAX;
      };

      // used_hash_size gives the number of entries, used_hash_count gives the
      // number of dirty hash values (I think)
      dedupe_edge_count(idx) = *used_hash_size;
      // Release the memory pool chunk back to the pool
      memory_pool.release_chunk(ptr_memory_pool_chunk);

    }  // operator()

    // reduces to find total number of rows that were too large
    KOKKOS_INLINE_FUNCTION
    void operator()(const member& thread, ordinal_t& thread_sum) const {
      ordinal_t idx = remaining(thread.league_rank());
      typedef ordinal_t hash_size_type;
      typedef ordinal_t hash_key_type;
      typedef scalar_t hash_value_type;

      // can't do this row at current hashmap size
      ordinal_t hash_entries = row_map(idx + 1) - row_map(idx);
      if (hash_entries >= max_hash_entries) {
        Kokkos::single(Kokkos::PerTeam(thread), [&]() { thread_sum++; });
        thread.team_barrier();
        return;
      }
      volatile ordinal_t* ptr_temp = nullptr;
      Kokkos::single(
          Kokkos::PerTeam(thread),
          [&](volatile ordinal_t*& ptr_write) {
            // Acquire a chunk from the memory pool using a spin-loop.
            ptr_write = nullptr;
            while (nullptr == ptr_write) {
              ptr_write = (volatile ordinal_t*)(memory_pool.allocate_chunk(thread.league_rank()));
            }
          },
          ptr_temp);
      thread.team_barrier();
      if (ptr_temp == nullptr) {
        return;
      }
      ordinal_t* ptr_memory_pool_chunk = (ordinal_t*)(ptr_temp);

      // These are updated by Hashmap_Accumulator insert functions.
      ordinal_t* used_hash_size = (ordinal_t*)(ptr_temp);
      ptr_temp++;
      ordinal_t* used_hash_count = (ordinal_t*)(ptr_temp);
      ptr_temp++;
      ordinal_t* write_idx = (ordinal_t*)(ptr_temp);
      ptr_temp++;
      Kokkos::single(Kokkos::PerTeam(thread), [&]() {
        *used_hash_size  = 0;
        *used_hash_count = 0;
        *write_idx       = 0;
      });

      // hash function is hash_size-1 (note: hash_size must be a power of 2)
      ordinal_t hash_func_pow2 = hash_size - 1;

      // Set pointer to hash indices
      ordinal_t* used_hash_indices = (ordinal_t*)(ptr_temp);
      ptr_temp += hash_size;

      // Set pointer to hash begins
      ordinal_t* hash_begins = (ordinal_t*)(ptr_temp);
      ptr_temp += hash_size;

      // Set pointer to hash nexts
      ordinal_t* hash_nexts = (ordinal_t*)(ptr_temp);
      ptr_temp += max_hash_entries;

      // Set pointer to hash keys
      ordinal_t* keys = (ordinal_t*)(ptr_temp);
      ptr_temp += max_hash_entries;

      // Set pointer to hash values
      scalar_t* values;
      if (use_out) {
        values = (scalar_t*)wgts_out.data() + row_map(idx);
      } else {
        values = (scalar_t*)(ptr_temp);
      }

      KokkosKernels::Experimental::HashmapAccumulator<hash_size_type, hash_key_type, hash_value_type,
                                                      KokkosKernels::Experimental::HashOpType::bitwiseAnd>
          hash_map(hash_size, hash_func_pow2, hash_begins, hash_nexts, keys, values);

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(thread, row_map(idx), row_map(idx + 1)),
                           [&](const edge_offset_t& i) {
                             ordinal_t key  = entries_in(i);
                             scalar_t value = wgts_in(i);
                             // duplicate keys may be inserted simultaneously, this causes
                             // problems we must handle later
                             int r = hash_map.vector_atomic_insert_into_hash_mergeAtomicAdd_TrackHashes(
                                 key, value, used_hash_size, used_hash_count, used_hash_indices);

                             // Check return code
                             if (r) {
                             }
                           });
      thread.team_barrier();

      // Reset the Begins values to -1 before releasing the memory pool chunk.
      // If you don't do this the next thread that grabs this memory chunk will
      // not work properly. Also merges values inside each linked list, because
      // there can be duplicate key insertions (these are hopefully rare or else
      // performance will suffer) This did not work as a TeamThreadRange, don't
      // know why (possibly issues with atomic addition on write_idx)
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(thread, (ordinal_t)0, *used_hash_count), [&](const ordinal_t& i) {
        ordinal_t dirty_hash = used_hash_indices[i];

        ordinal_t bucket = hash_begins[dirty_hash];

        // ascending-key bubble-sort the linked list
        // it really do be like that sometimes
        ordinal_t end_inner = ORD_MAX;
        while (end_inner != bucket) {
          ordinal_t last_idx = bucket;
          ordinal_t last_key = keys[last_idx];
          scalar_t last_val  = values[last_idx];
          bool is_sorted     = true;
          // bubble-up
          for (ordinal_t k = hash_nexts[bucket]; k != end_inner; k = hash_nexts[k]) {
            // swap
            if (keys[k] < last_key) {
              keys[last_idx]   = keys[k];
              values[last_idx] = values[k];
              keys[k]          = last_key;
              values[k]        = last_val;
              is_sorted        = false;
            }
            // increment last
            last_key = keys[k];
            last_val = values[k];
            last_idx = k;
          }
          end_inner = last_idx;
          if (is_sorted) {
            // end the outer loop
            end_inner = bucket;
          }
        }
        ordinal_t key  = keys[bucket];
        scalar_t val   = values[bucket];
        ordinal_t last = bucket;
        // merge linked list and write out
        for (ordinal_t j = hash_nexts[bucket]; j != ORD_MAX; j = hash_nexts[j]) {
          if (keys[j] == key) {
            val += values[j];
          } else {
            ordinal_t write_at    = row_map(idx) + Kokkos::atomic_fetch_add(write_idx, 1);
            entries_out(write_at) = key;
            if (use_out) {
              // reuse wgts_in as scratch space because we are overwriting
              // working memory if we use wgts_out
              wgts_in(write_at) = val;
            } else {
              wgts_out(write_at) = val;
            }
            key = keys[j];
            val = values[j];
          }
          hash_nexts[last] = ORD_MAX;
          last             = j;
        }
        hash_nexts[last] = ORD_MAX;
        // write out the final entry in linked list
        ordinal_t write_at    = row_map(idx) + Kokkos::atomic_fetch_add(write_idx, 1);
        entries_out(write_at) = key;
        if (use_out) {
          // reuse wgts_in as scratch space because we are overwriting
          // working memory if we use wgts_out
          wgts_in(write_at) = val;
        } else {
          wgts_out(write_at) = val;
        }
        hash_begins[dirty_hash] = ORD_MAX;
      });
      thread.team_barrier();
      // need to copy from wgts_in to wgts_out if we used wgts_in as scratch
      // space
      if (use_out) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(thread, (ordinal_t)0, *write_idx),
                             [&](const ordinal_t& i) { wgts_out(row_map(idx) + i) = wgts_in(row_map(idx) + i); });
      }

      Kokkos::single(Kokkos::PerTeam(thread), [&]() {
        // used_hash_size gives the number of entries, used_hash_count gives the
        // number of dirty hash values (I think)
        dedupe_edge_count(idx) = *write_idx;
        // Release the memory pool chunk back to the pool
        memory_pool.release_chunk(ptr_memory_pool_chunk);
      });

    }  // operator()

  };  // functorHashmapAccumulator

  static void getHashmapSizeAndCount(coarsen_handle& handle, const ordinal_t n, const ordinal_t remaining_count,
                                     vtx_view_t remaining, vtx_view_t edges_per_source, ordinal_t& hash_size,
                                     ordinal_t& max_entries, ordinal_t& mem_chunk_size, ordinal_t& mem_chunk_count) {
    ordinal_t avg_entries = 0;
    if (!is_host_space && static_cast<double>(remaining_count) / static_cast<double>(n) > 0.01) {
      Kokkos::parallel_reduce(
          "calc average among remaining", policy_t(0, remaining_count),
          KOKKOS_LAMBDA(const ordinal_t i, ordinal_t& thread_sum) {
            ordinal_t u      = remaining(i);
            ordinal_t degree = edges_per_source(u);
            thread_sum += degree;
          },
          avg_entries);
      // degrees are often skewed so we want to err on the side of bigger
      // hashmaps
      avg_entries = avg_entries * 2 / remaining_count;
      avg_entries++;
      if (avg_entries < 50) avg_entries = 50;
    } else {
      Kokkos::parallel_reduce(
          "calc max", policy_t(0, remaining_count),
          KOKKOS_LAMBDA(const ordinal_t i, ordinal_t& thread_max) {
            ordinal_t u      = remaining(i);
            ordinal_t degree = edges_per_source(u);
            if (degree > thread_max) {
              thread_max = degree;
            }
          },
          Kokkos::Max<ordinal_t, Kokkos::HostSpace>(avg_entries));
      // need precisely one larger than max, don't remember why atm
      avg_entries++;
    }

    // Set the hash_size as the next power of 2 bigger than hash_size_hint.
    // - hash_size must be a power of two since we use & rather than % (which is
    // slower) for computing the hash value for HashmapAccumulator.
    max_entries = avg_entries;
    hash_size   = 1;
    while (hash_size < max_entries) {
      hash_size *= 2;
    }

    // Determine memory chunk size for UniformMemoryPool
    mem_chunk_size = hash_size;         // for hash indices
    mem_chunk_size += hash_size;        // for hash begins
    mem_chunk_size += 3 * max_entries;  // for hash nexts, keys, and values (unless scalar_t
                                        // != ordinal_t, in which case memory is unused)
    mem_chunk_size += 10;               // for metadata
    mem_chunk_count = exec_space().concurrency();
    if (mem_chunk_count > remaining_count) {
      mem_chunk_count = remaining_count + 1;
    }

    if (!is_host_space) {
      // decrease number of mem_chunks to reduce memory usage if necessary
      size_t mem_needed =
          static_cast<size_t>(mem_chunk_count) * static_cast<size_t>(mem_chunk_size) * sizeof(ordinal_t);
      //~500MB
      size_t max_mem_allowed = handle.max_mem_allowed;
      if (mem_needed > max_mem_allowed) {
        size_t chunk_dif = mem_needed - max_mem_allowed;
        chunk_dif        = chunk_dif / (static_cast<size_t>(mem_chunk_size) * sizeof(ordinal_t));
        chunk_dif++;
        mem_chunk_count -= chunk_dif;
      }
    }
  }

  static void deduplicate_graph(coarsen_handle& handle, const ordinal_t n, const bool use_team,
                                vtx_view_t edges_per_source, vtx_view_t dest_by_source, wgt_view_t wgt_by_source,
                                const edge_view_t source_bucket_offset, edge_offset_t& gc_nedges) {
    if (handle.b == Hashmap || is_host_space) {
      ordinal_t remaining_count = n;
      vtx_view_t remaining("remaining vtx", n);
      Kokkos::parallel_for(
          policy_t(0, n), KOKKOS_LAMBDA(const ordinal_t i) { remaining(i) = i; });
      // deduplicate rows in phases starting with the small degree rows so we
      // can use small hashmaps increase the hashmap size each phase to the
      // necessary size for twice the average of remaining rows
      wgt_view_t wgt_out = wgt_by_source;
      if (!scal_eq_ord && !is_host_space) {
        // only necessary if teams might be used
        wgt_out = wgt_view_t("wgts out", wgt_by_source.extent(0));
      }
      do {
        // determine size for hashmap
        ordinal_t hash_size, max_entries, mem_chunk_size, mem_chunk_count;
        getHashmapSizeAndCount(handle, n, remaining_count, remaining, edges_per_source, hash_size, max_entries,
                               mem_chunk_size, mem_chunk_count);
        // Create Uniform Initialized Memory Pool
        KokkosKernels::Impl::PoolType pool_type = KokkosKernels::Impl::ManyThread2OneChunk;

        if (is_host_space) {
          pool_type = KokkosKernels::Impl::OneThread2OneChunk;
        }

        bool use_dyn = should_use_dyn(n, source_bucket_offset, mem_chunk_count);

        uniform_memory_pool_t memory_pool(mem_chunk_count, mem_chunk_size, ORD_MAX, pool_type);

        functorHashmapAccumulator hashmapAccumulator(source_bucket_offset, dest_by_source, dest_by_source,
                                                     wgt_by_source, wgt_out, edges_per_source, memory_pool, hash_size,
                                                     max_entries, remaining, !scal_eq_ord);

        ordinal_t old_remaining_count = remaining_count;
        if (!is_host_space && max_entries >= 128) {
          Kokkos::parallel_reduce("hashmap time", team_policy_t(old_remaining_count, 1, 64), hashmapAccumulator,
                                  remaining_count);
        } else {
          if (use_dyn) {
            Kokkos::parallel_reduce("hashmap time", dyn_policy_t(0, old_remaining_count, Kokkos::ChunkSize(128)),
                                    hashmapAccumulator, remaining_count);
          } else {
            Kokkos::parallel_reduce("hashmap time", policy_t(0, old_remaining_count), hashmapAccumulator,
                                    remaining_count);
          }
        }

        if (remaining_count > 0) {
          vtx_view_t new_remaining("new remaining vtx", remaining_count);

          Kokkos::parallel_scan(
              "move remaining vertices", policy_t(0, old_remaining_count),
              KOKKOS_LAMBDA(const ordinal_t i, ordinal_t& update, const bool final) {
                ordinal_t u = remaining(i);
                if (edges_per_source(u) >= max_entries) {
                  if (final) {
                    new_remaining(update) = u;
                  }
                  update++;
                }
              });

          remaining = new_remaining;
        }
      } while (remaining_count > 0);
      Kokkos::parallel_reduce(
          policy_t(0, n), KOKKOS_LAMBDA(const ordinal_t i, edge_offset_t& sum) { sum += edges_per_source(i); },
          gc_nedges);
      if (!scal_eq_ord && !is_host_space) {
        Kokkos::deep_copy(wgt_by_source, wgt_out);
      }
    } else if (handle.b == Sort) {
      // sort the (implicit) crs matrix
      KokkosSparse::sort_crs_matrix<exec_space, edge_view_t, vtx_view_t, wgt_view_t>(source_bucket_offset,
                                                                                     dest_by_source, wgt_by_source);

      // combine adjacent entries that are equal
      if (use_team) {
        // thread team version
        wgt_view_t wgts_out("wgts after dedupe", wgt_by_source.extent(0));
        vtx_view_t dest_out("dest after dedupe", dest_by_source.extent(0));
        functorDedupeAfterSort deduper(source_bucket_offset, dest_by_source, dest_out, wgt_by_source, wgts_out,
                                       edges_per_source);
        Kokkos::parallel_reduce("deduplicated sorted", team_policy_t(n, 64), deduper, gc_nedges);
        Kokkos::deep_copy(wgt_by_source, wgts_out);
        Kokkos::deep_copy(dest_by_source, dest_out);
      } else {
        // no thread team version
        functorDedupeAfterSort deduper(source_bucket_offset, dest_by_source, dest_by_source, wgt_by_source,
                                       wgt_by_source, edges_per_source);
        Kokkos::parallel_reduce("deduplicated sorted", policy_t(0, n), deduper, gc_nedges);
      }

    } else if (handle.b == Hybrid) {
      wgt_view_t wgt_out = wgt_by_source;
      if (!scal_eq_ord && !is_host_space) {
        // only necessary if teams might be used
        wgt_out = wgt_view_t("wgts out", wgt_by_source.extent(0));
      }
      ordinal_t limit = 128;
      // sort the (implicit) crs matrix, but only the low degree rows
      ordinal_t remaining_count =
          KokkosSparse::sort_low_degree_rows_crs_matrix<exec_space, edge_view_t, vtx_view_t, wgt_view_t>(
              source_bucket_offset, dest_by_source, wgt_by_source, limit);
      // combine adjacent entries that are equal
      {
        // no thread team version
        functorDedupeLowDegreeAfterSort deduper(source_bucket_offset, dest_by_source, dest_by_source, wgt_by_source,
                                                wgt_out, edges_per_source, limit);
        Kokkos::parallel_reduce("deduplicated sorted", policy_t(0, n), deduper, gc_nedges);
      }
      vtx_view_t remaining("remaining vtx", remaining_count);
      Kokkos::parallel_scan(
          "move remaining vertices", policy_t(0, n),
          KOKKOS_LAMBDA(const ordinal_t i, ordinal_t& update, const bool final) {
            if (edges_per_source(i) > limit) {
              if (final) {
                remaining(update) = i;
              }
              update++;
            }
          });
      // deduplicate rows in phases starting with the small degree rows so we
      // can use small hashmaps increase the hashmap size each phase to the
      // necessary size for twice the average of remaining rows
      while (remaining_count > 0) {
        // determine size for hashmap
        ordinal_t hash_size, max_entries, mem_chunk_size, mem_chunk_count;
        getHashmapSizeAndCount(handle, n, remaining_count, remaining, edges_per_source, hash_size, max_entries,
                               mem_chunk_size, mem_chunk_count);
        // Create Uniform Initialized Memory Pool
        KokkosKernels::Impl::PoolType pool_type = KokkosKernels::Impl::ManyThread2OneChunk;

        if (is_host_space) {
          pool_type = KokkosKernels::Impl::OneThread2OneChunk;
        }

        uniform_memory_pool_t memory_pool(mem_chunk_count, mem_chunk_size, ORD_MAX, pool_type);

        functorHashmapAccumulator hashmapAccumulator(source_bucket_offset, dest_by_source, dest_by_source,
                                                     wgt_by_source, wgt_out, edges_per_source, memory_pool, hash_size,
                                                     max_entries, remaining, !scal_eq_ord);

        ordinal_t old_remaining_count = remaining_count;
        if (!is_host_space && max_entries >= 128) {
          Kokkos::parallel_reduce("hashmap time", dyn_team_policy_t(old_remaining_count, 1, 64), hashmapAccumulator,
                                  remaining_count);
        } else {
          Kokkos::parallel_reduce("hashmap time", dyn_policy_t(0, old_remaining_count), hashmapAccumulator,
                                  remaining_count);
        }

        if (remaining_count > 0) {
          vtx_view_t new_remaining("new remaining vtx", remaining_count);

          Kokkos::parallel_scan(
              "move remaining vertices", policy_t(0, old_remaining_count),
              KOKKOS_LAMBDA(const ordinal_t i, ordinal_t& update, const bool final) {
                ordinal_t u = remaining(i);
                if (edges_per_source(u) >= max_entries) {
                  if (final) {
                    new_remaining(update) = u;
                  }
                  update++;
                }
              });

          remaining = new_remaining;
        }
      }
      gc_nedges = 0;
      Kokkos::parallel_reduce(
          policy_t(0, n), KOKKOS_LAMBDA(const ordinal_t i, edge_offset_t& sum) { sum += edges_per_source(i); },
          gc_nedges);
      if (!scal_eq_ord && !is_host_space) {
        Kokkos::deep_copy(wgt_by_source, wgt_out);
      }
    }
  }

  struct translationFunctor {
    matrix_t vcmap, g;
    vtx_view_t mapped_edges, edges_per_source;
    edge_view_t source_bucket_offset;
    vtx_view_t edges_out;
    wgt_view_t wgts_out;
    ordinal_t workLength;

    translationFunctor(matrix_t _vcmap, matrix_t _g, vtx_view_t _mapped_edges, vtx_view_t _edges_per_source,
                       edge_view_t _source_bucket_offset, vtx_view_t _edges_out, wgt_view_t _wgts_out)
        : vcmap(_vcmap),
          g(_g),
          mapped_edges(_mapped_edges),
          edges_per_source(_edges_per_source),
          source_bucket_offset(_source_bucket_offset),
          edges_out(_edges_out),
          wgts_out(_wgts_out),
          workLength(_g.numRows()) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member& t) const {
      ordinal_t i = t.league_rank() * t.team_size() + t.team_rank();
      if (i >= workLength) return;
      ordinal_t u         = vcmap.graph.entries(i);
      edge_offset_t start = g.graph.row_map(i);
      edge_offset_t end   = g.graph.row_map(i + 1);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(t, start, end), [&](const edge_offset_t idx) {
        ordinal_t v = mapped_edges(idx);
        if (u != v) {
          // fix this, inefficient
          edge_offset_t offset = Kokkos::atomic_fetch_add(&edges_per_source(u), 1);

          offset += source_bucket_offset(u);

          edges_out(offset) = v;
          wgts_out(offset)  = g.values(idx);
        }
      });
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const ordinal_t& i) const {
      ordinal_t u         = vcmap.graph.entries(i);
      edge_offset_t start = g.graph.row_map(i);
      edge_offset_t end   = g.graph.row_map(i + 1);
      for (edge_offset_t idx = start; idx < end; idx++) {
        ordinal_t v = mapped_edges(idx);
        if (u != v) {
          // fix this
          edge_offset_t offset = Kokkos::atomic_fetch_add(&edges_per_source(u), 1);

          offset += source_bucket_offset(u);

          edges_out(offset) = v;
          wgts_out(offset)  = g.values(idx);
        }
      }
    }
  };

  // optimized for regular distribution low degree rows
  static coarse_level_triple build_nonskew(coarsen_handle& handle, const matrix_t g, const matrix_t vcmap,
                                           vtx_view_t mapped_edges, vtx_view_t edges_per_source) {
    ordinal_t n  = g.numRows();
    ordinal_t nc = vcmap.numCols();
    edge_view_t source_bucket_offset("source_bucket_offsets", nc + 1);
    edge_offset_t gc_nedges = 0;

    Kokkos::parallel_scan("calc source offsets", policy_t(0, nc), prefix_sum(edges_per_source, source_bucket_offset));

    Kokkos::deep_copy(edges_per_source, static_cast<ordinal_t>(0));

    edge_subview_t sbo_subview   = Kokkos::subview(source_bucket_offset, nc);
    edge_offset_t nnz_pre_dedupe = 0;
    Kokkos::deep_copy(nnz_pre_dedupe, sbo_subview);

    vtx_view_t dest_by_source("dest_by_source", nnz_pre_dedupe);
    wgt_view_t wgt_by_source("wgt_by_source", nnz_pre_dedupe);

    // translates fine entries into coarse entries and writes into coarse rows
    translationFunctor translateF(vcmap, g, mapped_edges, edges_per_source, source_bucket_offset, dest_by_source,
                                  wgt_by_source);
    if (is_host_space) {
      bool use_dyn = should_use_dyn(n, g.graph.row_map, exec_space().concurrency());
      if (use_dyn) {
        Kokkos::parallel_for("move edges to coarse matrix", dyn_policy_t(0, n), translateF);
      } else {
        Kokkos::parallel_for("move edges to coarse matrix", policy_t(0, n), translateF);
      }
    } else {
      auto execSpaceEnum = KokkosKernels::Impl::kk_get_exec_space_type<exec_space>();
      int vectorLength   = KokkosKernels::Impl::kk_get_suggested_vector_size(n, g.nnz(), execSpaceEnum);
      team_policy_t dummy(1, 1, vectorLength);
      int teamSize = dummy.team_size_max(translateF, Kokkos::ParallelForTag());
      Kokkos::parallel_for("move edges to coarse matrix",
                           team_policy_t((n + teamSize - 1) / teamSize, teamSize, vectorLength), translateF);
    }

    deduplicate_graph(handle, nc, false, edges_per_source, dest_by_source, wgt_by_source, source_bucket_offset,
                      gc_nedges);

    edge_view_t source_offsets("source_offsets", nc + 1);

    Kokkos::parallel_scan("calc source offsets again", policy_t(0, nc), prefix_sum(edges_per_source, source_offsets));

    edge_subview_t edge_total_subview = Kokkos::subview(source_offsets, nc);
    Kokkos::deep_copy(gc_nedges, edge_total_subview);

    vtx_view_t dest_idx("dest_idx", gc_nedges);
    wgt_view_t wgts("wgts", gc_nedges);

    if (is_host_space) {
      bool use_dyn = should_use_dyn(nc, source_offsets, exec_space().concurrency());
      if (use_dyn) {
        Kokkos::parallel_for(
            "move deduped edges to new coarse matrix", dyn_policy_t(0, nc), KOKKOS_LAMBDA(const ordinal_t& u) {
              edge_offset_t start_origin = source_bucket_offset(u);
              edge_offset_t start_dest   = source_offsets(u);
              for (ordinal_t idx = 0; idx < edges_per_source(u); idx++) {
                dest_idx(start_dest + idx) = dest_by_source(start_origin + idx);
                wgts(start_dest + idx)     = wgt_by_source(start_origin + idx);
              }
            });
      } else {
        Kokkos::parallel_for(
            "move deduped edges to new coarse matrix", policy_t(0, nc), KOKKOS_LAMBDA(const ordinal_t& u) {
              edge_offset_t start_origin = source_bucket_offset(u);
              edge_offset_t start_dest   = source_offsets(u);
              for (ordinal_t idx = 0; idx < edges_per_source(u); idx++) {
                dest_idx(start_dest + idx) = dest_by_source(start_origin + idx);
                wgts(start_dest + idx)     = wgt_by_source(start_origin + idx);
              }
            });
      }
    } else {
      Kokkos::parallel_for(
          "move deduped edges to new coarse matrix", team_policy_t(nc, Kokkos::AUTO),
          KOKKOS_LAMBDA(const member& thread) {
            ordinal_t u                = thread.league_rank();
            edge_offset_t start_origin = source_bucket_offset(u);
            edge_offset_t start_dest   = source_offsets(u);
            Kokkos::parallel_for(Kokkos::TeamThreadRange(thread, edges_per_source(u)), [=](const ordinal_t idx) {
              dest_idx(start_dest + idx) = dest_by_source(start_origin + idx);
              wgts(start_dest + idx)     = wgt_by_source(start_origin + idx);
            });
          });
    }

    graph_type gc_graph(dest_idx, source_offsets);
    matrix_t gc("gc", nc, wgts, gc_graph);

    coarse_level_triple next_level;
    next_level.mtx = gc;
    return next_level;
  }

  // forms the explicit matrix created by symmetrizing the implicit matrix
  static matrix_t collapse_directed_to_undirected(const ordinal_t nc, const vtx_view_t source_edge_counts,
                                                  const edge_view_t source_row_map,
                                                  const vtx_view_t source_destinations, const wgt_view_t source_wgts) {
    vtx_view_t coarse_degree("coarse degree", nc);
    Kokkos::deep_copy(coarse_degree, source_edge_counts);

    Kokkos::parallel_for(
        "count directed edges owned by opposite endpoint", team_policy_t(nc, Kokkos::AUTO),
        KOKKOS_LAMBDA(const member& thread) {
          ordinal_t u         = thread.league_rank();
          edge_offset_t start = source_row_map(u);
          edge_offset_t end   = start + source_edge_counts(u);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(thread, start, end), [=](const edge_offset_t idx) {
            ordinal_t v = source_destinations(idx);
            // increment other vertex
            Kokkos::atomic_fetch_add(&coarse_degree(v), 1);
          });
        });

    edge_view_t target_row_map("target row map", nc + 1);

    Kokkos::parallel_scan("calc target row map", policy_t(0, nc), prefix_sum(coarse_degree, target_row_map));

    Kokkos::deep_copy(coarse_degree, static_cast<ordinal_t>(0));

    edge_offset_t coarse_edges_total         = 0;
    edge_subview_t coarse_edge_total_subview = Kokkos::subview(target_row_map, nc);
    Kokkos::deep_copy(coarse_edges_total, coarse_edge_total_subview);

    vtx_view_t dest_idx("dest_idx", coarse_edges_total);
    wgt_view_t wgts("wgts", coarse_edges_total);

    Kokkos::parallel_for(
        "move edges into correct size matrix", team_policy_t(nc, Kokkos::AUTO),
        functorCollapseDirectedToUndirected(source_row_map, target_row_map, source_edge_counts, coarse_degree,
                                            source_destinations, dest_idx, source_wgts, wgts));

    graph_type gc_graph(dest_idx, target_row_map);
    matrix_t gc("gc", nc, wgts, gc_graph);
    return gc;
  }

  // optimized for skewed degree distributions
  static coarse_level_triple build_skew(coarsen_handle& handle, const matrix_t g, const matrix_t vcmap,
                                        vtx_view_t mapped_edges, vtx_view_t degree_initial) {
    ordinal_t n             = g.numRows();
    ordinal_t nc            = vcmap.numCols();
    edge_offset_t gc_nedges = 0;

    vtx_view_t edges_per_source("edges_per_source", nc);

    // recount with edges only belonging to coarse vertex of smaller degree
    // matrix becomes directed
    Kokkos::parallel_for(
        "recount edges", team_policy_t(n, Kokkos::AUTO), KOKKOS_LAMBDA(const member& thread) {
          ordinal_t outer_idx         = thread.league_rank();
          ordinal_t u                 = vcmap.graph.entries(outer_idx);
          edge_offset_t start         = g.graph.row_map(outer_idx);
          edge_offset_t end           = g.graph.row_map(outer_idx + 1);
          ordinal_t nonLoopEdgesTotal = 0;
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange(thread, start, end),
              [=](const edge_offset_t idx, ordinal_t& local_sum) {
                ordinal_t v       = mapped_edges(idx);
                bool degree_less  = degree_initial(u) < degree_initial(v);
                bool degree_equal = degree_initial(u) == degree_initial(v);
                if (u != v && (degree_less || (degree_equal && u < v))) {
                  local_sum++;
                }
              },
              nonLoopEdgesTotal);
          Kokkos::single(Kokkos::PerTeam(thread),
                         [=]() { Kokkos::atomic_add(&edges_per_source(u), nonLoopEdgesTotal); });
        });

    edge_view_t source_bucket_offset("source_bucket_offsets", nc + 1);

    Kokkos::parallel_scan("calc source offsets", policy_t(0, nc), prefix_sum(edges_per_source, source_bucket_offset));
    edge_subview_t sbo_subview   = Kokkos::subview(source_bucket_offset, nc);
    edge_offset_t nnz_pre_dedupe = 0;
    Kokkos::deep_copy(nnz_pre_dedupe, sbo_subview);

    Kokkos::deep_copy(edges_per_source, static_cast<ordinal_t>(0));
    vtx_view_t dest_by_source("dest by source", nnz_pre_dedupe);
    wgt_view_t wgt_by_source("wgt by source", nnz_pre_dedupe);
    Kokkos::parallel_for(
        "combine fine rows", team_policy_t(n, Kokkos::AUTO), KOKKOS_LAMBDA(const member& thread) {
          ordinal_t outer_idx = thread.league_rank();
          ordinal_t u         = vcmap.graph.entries(outer_idx);
          edge_offset_t start = g.graph.row_map(outer_idx);
          edge_offset_t end   = g.graph.row_map(outer_idx + 1);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(thread, start, end), [=](const edge_offset_t idx) {
            ordinal_t v       = mapped_edges(idx);
            bool degree_less  = degree_initial(u) < degree_initial(v);
            bool degree_equal = degree_initial(u) == degree_initial(v);
            if (degree_less || (degree_equal && u < v)) {
              edge_offset_t offset = Kokkos::atomic_fetch_add(&edges_per_source(u), 1);

              offset += source_bucket_offset(u);

              dest_by_source(offset) = v;
              wgt_by_source(offset)  = g.values(idx);
            }
          });
        });
    gc_nedges = 0;

    deduplicate_graph(handle, nc, true, edges_per_source, dest_by_source, wgt_by_source, source_bucket_offset,
                      gc_nedges);

    // form the final coarse graph, which requires symmetrizing the matrix
    matrix_t gc =
        collapse_directed_to_undirected(nc, edges_per_source, source_bucket_offset, dest_by_source, wgt_by_source);

    coarse_level_triple next_level;
    next_level.mtx = gc;
    return next_level;
  }

  // optimized for very large row sizes caused by lots of duplicate entries
  // first translates each fine entry to its coarse vertex label
  // deduplicates within each fine row
  // combines fine rows into coarse rows
  // deduplicates within each coarse row
  static coarse_level_triple build_high_duplicity(coarsen_handle& handle, const matrix_t g, const matrix_t vcmap,
                                                  vtx_view_t mapped_edges, vtx_view_t degree_initial) {
    ordinal_t n             = g.numRows();
    ordinal_t nc            = vcmap.numCols();
    edge_offset_t gc_nedges = 0;

    vtx_view_t dedupe_count("dedupe count", n);
    edge_view_t row_map_copy("row map copy", n + 1);

    // recount fine row sizes with edges only belonging to fine vertex of coarse
    // vertex of smaller degree matrix becomes directed
    Kokkos::parallel_for(
        "recount edges", team_policy_t(n, Kokkos::AUTO), KOKKOS_LAMBDA(const member& thread) {
          ordinal_t outer_idx         = thread.league_rank();
          ordinal_t u                 = vcmap.graph.entries(outer_idx);
          edge_offset_t start         = g.graph.row_map(outer_idx);
          edge_offset_t end           = g.graph.row_map(outer_idx + 1);
          ordinal_t nonLoopEdgesTotal = 0;
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange(thread, start, end),
              [=](const edge_offset_t idx, ordinal_t& local_sum) {
                ordinal_t v       = mapped_edges(idx);
                bool degree_less  = degree_initial(u) < degree_initial(v);
                bool degree_equal = degree_initial(u) == degree_initial(v);
                if (u != v && (degree_less || (degree_equal && u < v))) {
                  local_sum++;
                }
              },
              nonLoopEdgesTotal);
          Kokkos::single(Kokkos::PerTeam(thread), [=]() { dedupe_count(outer_idx) = nonLoopEdgesTotal; });
        });

    Kokkos::parallel_scan("calc source offsets", policy_t(0, n), prefix_sum(dedupe_count, row_map_copy));
    // reset counters to 0
    Kokkos::deep_copy(dedupe_count, static_cast<ordinal_t>(0));

    edge_subview_t fine_recount_subview = Kokkos::subview(row_map_copy, n);
    edge_offset_t fine_recount          = 0;
    Kokkos::deep_copy(fine_recount, fine_recount_subview);

    vtx_view_t dest_fine("fine to coarse dests", fine_recount);
    wgt_view_t wgt_fine("fine to coarse wgts", fine_recount);

    // create a new directed version of the fine matrix
    Kokkos::parallel_for(
        "move edges to new matrix", team_policy_t(n, Kokkos::AUTO), KOKKOS_LAMBDA(const member& thread) {
          ordinal_t outer_idx = thread.league_rank();
          ordinal_t u         = vcmap.graph.entries(outer_idx);
          edge_offset_t start = g.graph.row_map(outer_idx);
          edge_offset_t end   = g.graph.row_map(outer_idx + 1);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(thread, start, end), [=](const edge_offset_t idx) {
            ordinal_t v       = mapped_edges(idx);
            bool degree_less  = degree_initial(u) < degree_initial(v);
            bool degree_equal = degree_initial(u) == degree_initial(v);
            if (u != v && (degree_less || (degree_equal && u < v))) {
              edge_offset_t offset = Kokkos::atomic_fetch_add(&dedupe_count(outer_idx), 1);

              offset += row_map_copy(outer_idx);

              dest_fine(offset) = v;
              wgt_fine(offset)  = g.values(idx);
            }
          });
        });
    //"delete" these views
    Kokkos::resize(mapped_edges, 0);

    // deduplicate coarse adjacencies within each fine row
    deduplicate_graph(handle, n, true, dedupe_count, dest_fine, wgt_fine, row_map_copy, gc_nedges);

    edge_view_t source_bucket_offset("source_bucket_offsets", nc + 1);
    vtx_view_t edges_per_source("edges_per_source", nc);

    Kokkos::parallel_for(
        "sum fine row sizes", policy_t(0, n), KOKKOS_LAMBDA(const ordinal_t i) {
          ordinal_t u = vcmap.graph.entries(i);
          Kokkos::atomic_fetch_add(&edges_per_source(u), dedupe_count(i));
        });
    Kokkos::parallel_scan("calc source offsets", policy_t(0, nc), prefix_sum(edges_per_source, source_bucket_offset));
    Kokkos::deep_copy(edges_per_source, static_cast<ordinal_t>(0));
    vtx_view_t dest_by_source("dest by source", gc_nedges);
    wgt_view_t wgt_by_source("wgt by source", gc_nedges);
    Kokkos::parallel_for(
        "combine deduped fine rows", team_policy_t(n, Kokkos::AUTO), KOKKOS_LAMBDA(const member& thread) {
          ordinal_t outer_idx = thread.league_rank();
          ordinal_t u         = vcmap.graph.entries(outer_idx);
          edge_offset_t start = row_map_copy(outer_idx);
          edge_offset_t end   = start + dedupe_count(outer_idx);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(thread, start, end), [=](const edge_offset_t idx) {
            ordinal_t v       = dest_fine(idx);
            bool degree_less  = degree_initial(u) < degree_initial(v);
            bool degree_equal = degree_initial(u) == degree_initial(v);
            if (degree_less || (degree_equal && u < v)) {
              edge_offset_t offset = Kokkos::atomic_fetch_add(&edges_per_source(u), 1);

              offset += source_bucket_offset(u);

              dest_by_source(offset) = v;
              wgt_by_source(offset)  = wgt_fine(idx);
            }
          });
        });
    gc_nedges = 0;
    Kokkos::resize(dest_fine, 0);
    Kokkos::resize(wgt_fine, 0);

    deduplicate_graph(handle, nc, true, edges_per_source, dest_by_source, wgt_by_source, source_bucket_offset,
                      gc_nedges);

    // form the final coarse graph, which requires symmetrizing the matrix
    matrix_t gc =
        collapse_directed_to_undirected(nc, edges_per_source, source_bucket_offset, dest_by_source, wgt_by_source);

    coarse_level_triple next_level;
    next_level.mtx = gc;
    return next_level;
  }

  struct countingFunctor {
    // counts adjancies for each coarse vertex
    // also calculates coarse vertex wgts
    matrix_t vcmap, g;
    vtx_view_t mapped_edges, degree_initial;
    vtx_view_t c_vtx_w, f_vtx_w;
    ordinal_t workLength;

    countingFunctor(matrix_t _vcmap, matrix_t _g, vtx_view_t _mapped_edges, vtx_view_t _degree_initial,
                    vtx_view_t _c_vtx_w, vtx_view_t _f_vtx_w)
        : vcmap(_vcmap),
          g(_g),
          mapped_edges(_mapped_edges),
          degree_initial(_degree_initial),
          c_vtx_w(_c_vtx_w),
          f_vtx_w(_f_vtx_w),
          workLength(_g.numRows()) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member& t) const {
      ordinal_t i = t.league_rank() * t.team_size() + t.team_rank();
      if (i >= workLength) return;
      ordinal_t u                 = vcmap.graph.entries(i);
      edge_offset_t start         = g.graph.row_map(i);
      edge_offset_t end           = g.graph.row_map(i + 1);
      ordinal_t nonLoopEdgesTotal = 0;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(t, start, end),
          [&](const edge_offset_t idx, ordinal_t& local_sum) {
            ordinal_t v       = vcmap.graph.entries(g.graph.entries(idx));
            mapped_edges(idx) = v;
            if (u != v) {
              local_sum++;
            }
          },
          nonLoopEdgesTotal);
      Kokkos::single(Kokkos::PerThread(t), [&]() {
        Kokkos::atomic_add(&degree_initial(u), nonLoopEdgesTotal);
        Kokkos::atomic_add(&c_vtx_w(u), f_vtx_w(i));
      });
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const ordinal_t& i) const {
      ordinal_t u                 = vcmap.graph.entries(i);
      edge_offset_t start         = g.graph.row_map(i);
      edge_offset_t end           = g.graph.row_map(i + 1);
      ordinal_t nonLoopEdgesTotal = 0;
      for (edge_offset_t idx = start; idx < end; idx++) {
        ordinal_t v       = vcmap.graph.entries(g.graph.entries(idx));
        mapped_edges(idx) = v;
        if (u != v) {
          nonLoopEdgesTotal++;
        }
      }
      Kokkos::atomic_add(&degree_initial(u), nonLoopEdgesTotal);
      Kokkos::atomic_add(&c_vtx_w(u), f_vtx_w(i));
    }
  };

  static coarse_level_triple build_coarse_graph(coarsen_handle& handle, const coarse_level_triple level,
                                                const matrix_t vcmap) {
    if (handle.b == Spgemm || handle.b == Spgemm_transpose_first) {
      return build_coarse_graph_spgemm(handle, level, vcmap);
    }

    matrix_t g   = level.mtx;
    ordinal_t n  = g.numRows();
    ordinal_t nc = vcmap.numCols();

    vtx_view_t mapped_edges("mapped edges", g.nnz());

    vtx_view_t degree_initial("edges_per_source", nc);
    vtx_view_t f_vtx_w = level.vtx_wgts;
    vtx_view_t c_vtx_w = vtx_view_t("coarse vertex weights", nc);

    // count non-self loop edges per coarse vertex
    // also computes coarse vertex weights
    countingFunctor countF(vcmap, g, mapped_edges, degree_initial, c_vtx_w, f_vtx_w);
    if (is_host_space) {
      Kokkos::parallel_for("count edges per coarse vertex (also compute coarse vertex weights)", policy_t(0, n),
                           countF);
    } else {
      auto execSpaceEnum = KokkosKernels::Impl::kk_get_exec_space_type<exec_space>();
      int vectorLength   = KokkosKernels::Impl::kk_get_suggested_vector_size(n, g.nnz(), execSpaceEnum);
      team_policy_t dummy(1, 1, vectorLength);
      int teamSize = dummy.team_size_max(countF, Kokkos::ParallelForTag());
      // count edges per vertex
      Kokkos::parallel_for("count edges per coarse vertex (also compute coarse vertex weights)",
                           team_policy_t((n + teamSize - 1) / teamSize, teamSize, vectorLength), countF);
    }

    // compute max row size and avg row size
    // use this to determine most efficient method for building coarse graph
    // (for load balance primarily)
    edge_offset_t total_unduped = 0;
    ordinal_t max_unduped       = 0;
    Kokkos::parallel_reduce(
        "find max", policy_t(0, nc),
        KOKKOS_LAMBDA(const ordinal_t i, ordinal_t& l_max) {
          if (l_max <= degree_initial(i)) {
            l_max = degree_initial(i);
          }
        },
        Kokkos::Max<ordinal_t, Kokkos::HostSpace>(max_unduped));
    Kokkos::parallel_reduce(
        "find total", policy_t(0, nc),
        KOKKOS_LAMBDA(const ordinal_t i, edge_offset_t& sum) { sum += degree_initial(i); }, total_unduped);
    ordinal_t avg_unduped = total_unduped / nc;

    coarse_level_triple next_level;
    // optimized subroutines for sufficiently irregular graphs or high average
    // adjacency rows don't do optimizations if running on CPU (the default host
    // space)
    if (avg_unduped > (nc / 4) && !is_host_space) {
      next_level = build_high_duplicity(handle, g, vcmap, mapped_edges, degree_initial);
    } else if (avg_unduped > 50 && (max_unduped / 10) > avg_unduped && !is_host_space) {
      next_level = build_skew(handle, g, vcmap, mapped_edges, degree_initial);
    } else {
      next_level = build_nonskew(handle, g, vcmap, mapped_edges, degree_initial);
    }

    next_level.vtx_wgts        = c_vtx_w;
    next_level.level           = level.level + 1;
    next_level.interp_mtx      = vcmap;
    next_level.uniform_weights = false;
    return next_level;
  }

  static matrix_t generate_coarse_mapping(coarsen_handle& handle, const matrix_t g, bool uniform_weights) {
    matrix_t interpolation_graph;
    int choice = 0;

    switch (handle.h) {
      case Match: choice = 0; break;
      case MtMetis: choice = 1; break;
      default: choice = 0;
    }

    switch (handle.h) {
      case HECv1: interpolation_graph = mapper_t::coarsen_HEC(g, uniform_weights); break;
      case Match:
      case MtMetis: interpolation_graph = mapper_t::coarsen_match(g, uniform_weights, choice); break;
      case MIS2: interpolation_graph = mapper_t::coarsen_mis_2(g); break;
      case GOSHv2: interpolation_graph = mapper_t::coarsen_GOSH_v2(g); break;
      case GOSHv1: interpolation_graph = mapper_t::coarsen_GOSH(g); break;
    }
    return interpolation_graph;
  }

  // we can support weighted vertices pretty easily, but we don't rn
  // this function can't return the generated list directly because of an NVCC
  // compiler bug caller must use the get_levels() method after calling this
  // function
  static void generate_coarse_graphs(coarsen_handle& handle, const matrix_t fine_g, bool uniform_weights = false) {
    ordinal_t fine_n                       = fine_g.numRows();
    std::list<coarse_level_triple>& levels = handle.results;
    levels.clear();
    coarse_level_triple finest;
    finest.mtx = fine_g;
    // 1-indexed, not zero indexed
    finest.level           = 1;
    finest.uniform_weights = uniform_weights;
    vtx_view_t vtx_weights("vertex weights", fine_n);
    Kokkos::deep_copy(vtx_weights, static_cast<ordinal_t>(1));
    finest.vtx_wgts = vtx_weights;
    levels.push_back(finest);
    while (levels.rbegin()->mtx.numRows() > handle.coarse_vtx_cutoff) {
      coarse_level_triple current_level = *levels.rbegin();

      matrix_t interp_graph = generate_coarse_mapping(handle, current_level.mtx, current_level.uniform_weights);

      if (interp_graph.numCols() < handle.min_allowed_vtx) {
        break;
      }

      coarse_level_triple next_level = build_coarse_graph(handle, current_level, interp_graph);

      levels.push_back(next_level);

      if (levels.size() > handle.max_levels) break;
    }
  }
};

}  // end namespace Experimental
}  // end namespace KokkosGraph
// exclude from Cuda builds without lambdas enabled
#endif
