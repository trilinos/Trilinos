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
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_UnorderedMap.hpp"
#include <iostream>
#include <limits>

#include "KokkosKernels_ExecSpaceUtils.hpp"
#include "KokkosKernels_SimpleUtils.hpp"
#include "KokkosKernels_PrintUtils.hpp"
#include "KokkosKernels_VectorUtils.hpp"

#ifndef _KOKKOSKERNELSUTILS_HPP
#define _KOKKOSKERNELSUTILS_HPP

namespace KokkosKernels {

namespace Impl {

template <typename ExecutionSpace>
ExecSpaceType get_exec_space_type() {
  return kk_get_exec_space_type<ExecutionSpace>();
}

inline int get_suggested_vector__size(size_t nr, size_t nnz, ExecSpaceType exec_space) {
  return kk_get_suggested_vector_size(nr, nnz, exec_space);
}

template <typename in_lno_view_t, typename out_lno_view_t, typename MyExecSpace>
void get_histogram(typename in_lno_view_t::size_type in_elements, in_lno_view_t in_view,
                   out_lno_view_t histogram /*must be initialized with 0s*/) {
  kk_get_histogram<in_lno_view_t, out_lno_view_t, MyExecSpace>(in_elements, in_view, histogram);
}

template <typename idx, typename ExecutionSpace>
void get_suggested_vector_size(int &suggested_vector_size_, idx nr, idx nnz) {
  suggested_vector_size_ = kk_get_suggested_vector_size(nr, nnz, get_exec_space_type<ExecutionSpace>());
}

// Get the best team size for the given functor.
// If it uses shared memory, the amount used must be available through
// f.team_shmem_size(n), not through the TeamPolicy. If this is how dynamic
// shared is set, just use AUTO for the team size.
template <typename team_policy_t, typename Functor, typename ParallelTag = Kokkos::ParallelForTag>
int get_suggested_team_size(Functor &f, int vector_size) {
  using execution_space = typename team_policy_t::traits::execution_space;
  if (kk_is_gpu_exec_space<execution_space>()) {
    team_policy_t temp(1, 1, vector_size);
    return temp.team_size_recommended(f, ParallelTag());
  } else
    return 1;
}

template <typename team_policy_t, typename Functor, typename ParallelTag = Kokkos::ParallelForTag>
int get_suggested_team_size(Functor &f, int vector_size, size_t sharedPerTeam, size_t sharedPerThread) {
  using execution_space = typename team_policy_t::traits::execution_space;
  if (kk_is_gpu_exec_space<execution_space>()) {
    team_policy_t temp = team_policy_t(1, 1, vector_size)
                             .set_scratch_size(0, Kokkos::PerTeam(sharedPerTeam), Kokkos::PerThread(sharedPerThread));
    return temp.team_size_recommended(f, ParallelTag());
  } else
    return 1;
}

template <typename idx_array_type, typename idx_edge_array_type, typename idx_out_edge_array_type, typename team_member>
struct FillSymmetricEdges {
  typedef typename idx_array_type::value_type idx;
  idx num_rows;
  idx nnz;
  idx_array_type xadj;
  idx_edge_array_type adj;

  idx_out_edge_array_type srcs;
  idx_out_edge_array_type dsts;

  FillSymmetricEdges(typename idx_array_type::value_type num_rows_, idx_array_type xadj_, idx_edge_array_type adj_,

                     idx_out_edge_array_type srcs_, idx_out_edge_array_type dsts_)
      : num_rows(num_rows_), nnz(adj_.extent(0)), xadj(xadj_), adj(adj_), srcs(srcs_), dsts(dsts_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member &teamMember) const {
    idx ii = teamMember.league_rank() * teamMember.team_size() + teamMember.team_rank();
    if (ii >= num_rows) return;
    idx row_begin = xadj[ii];
    idx row_end   = xadj[ii + 1];

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, row_end - row_begin), [&](idx i) {
      idx adjind   = i + row_begin;
      idx colIndex = adj[adjind];
      if (colIndex < num_rows) {
        srcs[adjind] = ii + 1;
        dsts[adjind] = colIndex + 1;
        if (colIndex != ii) {
          srcs[adjind + nnz] = colIndex + 1;
          dsts[adjind + nnz] = ii + 1;
        }
      }
    });
  }
};

template <typename in_lno_row_view_t, typename in_lno_nnz_view_t, typename hashmap_t, typename out_lno_row_view_t,
          typename team_member>
struct FillSymmetricEdgesHashMap {
  typedef typename in_lno_row_view_t::value_type idx;
  idx num_rows;
  idx nnz;
  in_lno_row_view_t xadj;
  in_lno_nnz_view_t adj;
  hashmap_t umap;
  out_lno_row_view_t pre_pps;
  bool lower_only;

  FillSymmetricEdgesHashMap(idx num_rows_, in_lno_row_view_t xadj_, in_lno_nnz_view_t adj_, hashmap_t hashmap_,
                            out_lno_row_view_t pre_pps_)
      : num_rows(num_rows_), nnz(adj_.extent(0)), xadj(xadj_), adj(adj_), umap(hashmap_), pre_pps(pre_pps_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member &teamMember /*, idx &nnz*/) const {
    typedef typename std::remove_reference<decltype(pre_pps(0))>::type atomic_incr_type;
    idx ii = teamMember.league_rank() * teamMember.team_size() + teamMember.team_rank();
    if (ii >= num_rows) {
      return;
    }
    idx row_begin = xadj[ii];
    idx row_end   = xadj[ii + 1];
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, row_end - row_begin), [&](idx i) {
      idx adjind   = i + row_begin;
      idx colIndex = adj[adjind];
      if (colIndex < num_rows) {
        if (colIndex < ii) {
          Kokkos::UnorderedMapInsertResult r = umap.insert(Kokkos::pair<idx, idx>(colIndex, ii));
          if (r.success()) {
            Kokkos::atomic_fetch_add(&(pre_pps(ii)), atomic_incr_type(1));

            Kokkos::atomic_fetch_add(&(pre_pps(colIndex)), atomic_incr_type(1));
          }
        } else if (colIndex > ii) {
          Kokkos::UnorderedMapInsertResult r = umap.insert(Kokkos::pair<idx, idx>(ii, colIndex));
          if (r.success()) {
            Kokkos::atomic_fetch_add(&(pre_pps(colIndex)), atomic_incr_type(1));

            Kokkos::atomic_fetch_add(&(pre_pps(ii)), atomic_incr_type(1));
          }
        } else {
          Kokkos::atomic_fetch_add(&(pre_pps(ii)), atomic_incr_type(1));
        }
      }
    });
  }
};

template <typename in_lno_row_view_t, typename in_lno_nnz_view_t, typename hashmap_t, typename out_lno_row_view_t,
          typename team_member>
struct FillSymmetricLowerEdgesHashMap {
  typedef typename in_lno_row_view_t::value_type idx;
  idx num_rows;
  idx nnz;
  in_lno_row_view_t xadj;
  in_lno_nnz_view_t adj;
  hashmap_t umap;
  out_lno_row_view_t pre_pps;

  FillSymmetricLowerEdgesHashMap(idx num_rows_, in_lno_row_view_t xadj_, in_lno_nnz_view_t adj_, hashmap_t hashmap_,
                                 out_lno_row_view_t pre_pps_, bool /* lower_only_ */ = false)
      : num_rows(num_rows_), nnz(adj_.extent(0)), xadj(xadj_), adj(adj_), umap(hashmap_), pre_pps(pre_pps_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member &teamMember /*, idx &nnz*/) const {
    typedef typename std::remove_reference<decltype(pre_pps(0))>::type atomic_incr_type;
    idx ii = teamMember.league_rank() * teamMember.team_size() + teamMember.team_rank();
    if (ii >= num_rows) {
      return;
    }
    idx row_begin = xadj[ii];
    idx row_end   = xadj[ii + 1];

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, row_end - row_begin), [&](idx i) {
      idx adjind   = i + row_begin;
      idx colIndex = adj[adjind];
      if (colIndex < num_rows) {
        if (colIndex < ii) {
          Kokkos::UnorderedMapInsertResult r = umap.insert(Kokkos::pair<idx, idx>(colIndex, ii));
          if (r.success()) {
            Kokkos::atomic_fetch_add(&(pre_pps(colIndex)), atomic_incr_type(1));
          }
        } else if (colIndex > ii) {
          Kokkos::UnorderedMapInsertResult r = umap.insert(Kokkos::pair<idx, idx>(ii, colIndex));
          if (r.success()) {
            Kokkos::atomic_fetch_add(&(pre_pps(ii)), atomic_incr_type(1));
          }
        }
      }
    });
  }
};

template <typename in_lno_row_view_t, typename in_lno_nnz_view_t, typename hashmap_t, typename out_lno_row_view_t,
          typename out_lno_nnz_view_t, typename team_member_t>
struct FillSymmetricCRS_HashMap {
  typedef typename in_lno_row_view_t::value_type idx;
  idx num_rows;
  idx nnz;
  in_lno_row_view_t xadj;
  in_lno_nnz_view_t adj;
  hashmap_t umap;
  out_lno_row_view_t pre_pps;
  out_lno_nnz_view_t sym_adj;

  FillSymmetricCRS_HashMap(idx num_rows_, in_lno_row_view_t xadj_, in_lno_nnz_view_t adj_, hashmap_t hashmap_,
                           out_lno_row_view_t pre_pps_, out_lno_nnz_view_t sym_adj_)
      : num_rows(num_rows_),
        nnz(adj_.extent(0)),
        xadj(xadj_),
        adj(adj_),
        umap(hashmap_),
        pre_pps(pre_pps_),
        sym_adj(sym_adj_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member_t &teamMember) const {
    typedef typename std::remove_reference<decltype(pre_pps(0))>::type atomic_incr_type;
    idx ii = teamMember.league_rank() * teamMember.team_size() + teamMember.team_rank();
    if (ii >= num_rows) {
      return;
    }
    idx row_begin = xadj[ii];
    idx row_end   = xadj[ii + 1];

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, row_end - row_begin), [&](idx i) {
      idx adjind   = i + row_begin;
      idx colIndex = adj[adjind];
      if (colIndex < num_rows) {
        if (colIndex < ii) {
          if (umap.insert(Kokkos::pair<idx, idx>(colIndex, ii)).success()) {
            idx cAdjInd      = Kokkos::atomic_fetch_add(&(pre_pps(colIndex)), atomic_incr_type(1));
            idx iAdjInd      = Kokkos::atomic_fetch_add(&(pre_pps(ii)), atomic_incr_type(1));
            sym_adj[cAdjInd] = ii;
            sym_adj[iAdjInd] = colIndex;
          }
        } else if (colIndex > ii) {
          if (umap.insert(Kokkos::pair<idx, idx>(ii, colIndex)).success()) {
            idx cAdjInd      = Kokkos::atomic_fetch_add(&(pre_pps(colIndex)), atomic_incr_type(1));
            idx iAdjInd      = Kokkos::atomic_fetch_add(&(pre_pps(ii)), atomic_incr_type(1));
            sym_adj[cAdjInd] = ii;
            sym_adj[iAdjInd] = colIndex;
          }
        } else {
          idx cAdjInd      = Kokkos::atomic_fetch_add(&(pre_pps(colIndex)), atomic_incr_type(1));
          sym_adj[cAdjInd] = ii;
        }
      }
    });
  }
};

template <typename in_lno_row_view_t, typename in_lno_nnz_view_t, typename hashmap_t, typename out_lno_nnz_view_t,
          typename out_lno_row_view_t, typename team_member_t>
struct FillSymmetricEdgeList_HashMap {
  typedef typename in_lno_row_view_t::value_type idx;
  idx num_rows;
  idx nnz;
  in_lno_row_view_t xadj;
  in_lno_nnz_view_t adj;
  hashmap_t umap;
  out_lno_nnz_view_t sym_src;
  out_lno_nnz_view_t sym_dst;
  out_lno_row_view_t pps;

  FillSymmetricEdgeList_HashMap(idx num_rows_, in_lno_row_view_t xadj_, in_lno_nnz_view_t adj_, hashmap_t hashmap_,
                                out_lno_nnz_view_t sym_src_, out_lno_nnz_view_t sym_dst_, out_lno_row_view_t pps_)
      : num_rows(num_rows_),
        nnz(adj_.extent(0)),
        xadj(xadj_),
        adj(adj_),
        umap(hashmap_),
        sym_src(sym_src_),
        sym_dst(sym_dst_),
        pps(pps_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member_t &teamMember) const {
    typedef typename std::remove_reference<decltype(pps(0))>::type atomic_incr_type;
    idx ii = teamMember.league_rank() * teamMember.team_size() + teamMember.team_rank();
    if (ii >= num_rows) {
      return;
    }
    idx row_begin = xadj[ii];
    idx row_end   = xadj[ii + 1];

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, row_end - row_begin), [&](idx i) {
      idx adjind   = i + row_begin;
      idx colIndex = adj[adjind];
      if (colIndex < num_rows) {
        if (colIndex < ii) {
          if (umap.insert(Kokkos::pair<idx, idx>(colIndex, ii)).success()) {
            idx cAdjInd      = Kokkos::atomic_fetch_add(&(pps(colIndex)), atomic_incr_type(1));
            sym_src[cAdjInd] = colIndex;
            sym_dst[cAdjInd] = ii;
          }
        } else if (colIndex > ii) {
          if (umap.insert(Kokkos::pair<idx, idx>(ii, colIndex)).success()) {
            idx cAdjInd      = Kokkos::atomic_fetch_add(&(pps(ii)), atomic_incr_type(1));
            sym_src[cAdjInd] = ii;
            sym_dst[cAdjInd] = colIndex;
          }
        }
      }
    });
  }
};

template <typename idx_array_type>
void print_1Dview(std::ostream &os, idx_array_type view, bool print_all = false, const char *sep = " ") {
  kk_print_1Dview(os, view, print_all, sep);
}

template <typename idx_array_type>
void print_1Dview(idx_array_type view, bool print_all = false) {
  kk_print_1Dview(view, print_all);
}

template <typename lno_t, typename memory_space>
void print_1Dpointer(const lno_t *pview, size_t size, bool print_all = false) {
  typedef Kokkos::View<const lno_t *, memory_space, Kokkos::MemoryUnmanaged> um_array_type;
  um_array_type view(pview, size);
  kk_print_1Dview(view, print_all);
}

template <typename forward_map_type, typename reverse_map_type>
struct Reverse_Map_Init {
  typedef typename forward_map_type::value_type forward_type;
  typedef typename reverse_map_type::value_type reverse_type;
  forward_map_type forward_map;
  reverse_map_type reverse_map_xadj;
  Reverse_Map_Init(forward_map_type forward_map_, reverse_map_type reverse_xadj_)
      : forward_map(forward_map_), reverse_map_xadj(reverse_xadj_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &ii) const {
    typedef typename std::remove_reference<decltype(reverse_map_xadj(0))>::type atomic_incr_type;
    forward_type fm = forward_map[ii];
    Kokkos::atomic_fetch_add(&(reverse_map_xadj(fm)), atomic_incr_type(1));
  }
};

template <typename forward_map_type, typename reverse_map_type>
struct Fill_Reverse_Map {
  typedef typename forward_map_type::value_type forward_type;
  typedef typename reverse_map_type::value_type reverse_type;
  forward_map_type forward_map;
  reverse_map_type reverse_map_xadj;
  reverse_map_type reverse_map_adj;

  Fill_Reverse_Map(forward_map_type forward_map_, reverse_map_type reverse_map_xadj_, reverse_map_type reverse_map_adj_)
      : forward_map(forward_map_), reverse_map_xadj(reverse_map_xadj_), reverse_map_adj(reverse_map_adj_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &ii) const {
    typedef typename std::remove_reference<decltype(reverse_map_xadj(0))>::type atomic_incr_type;
    forward_type c                  = forward_map[ii];
    const reverse_type future_index = Kokkos::atomic_fetch_add(&(reverse_map_xadj(c - 1)), atomic_incr_type(1));
    reverse_map_adj(future_index)   = ii;
  }
};

template <typename forward_array_type, typename MyExecSpace>
void inclusive_parallel_prefix_sum(MyExecSpace my_exec_space, typename forward_array_type::value_type num_elements,
                                   forward_array_type arr) {
  return kk_inclusive_parallel_prefix_sum(my_exec_space, num_elements, arr);
}

template <typename forward_array_type, typename MyExecSpace>
void inclusive_parallel_prefix_sum(typename forward_array_type::value_type num_elements, forward_array_type arr) {
  MyExecSpace my_exec_space;
  return inclusive_parallel_prefix_sum(my_exec_space, num_elements, arr);
}

template <typename forward_array_type, typename MyExecSpace>
void exclusive_parallel_prefix_sum(typename forward_array_type::value_type num_elements, forward_array_type arr) {
  kk_exclusive_parallel_prefix_sum<MyExecSpace>(num_elements, arr);
}

template <typename array_type>
struct PropogataMaxValstoZeros {
  typedef typename array_type::value_type idx;
  array_type array_sum;
  PropogataMaxValstoZeros(array_type arr_) : array_sum(arr_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, idx &update, const bool final) const {
    idx value = array_sum(ii);
    if (value != 0) {
      update = value;
    } else if (final) {
      array_sum(ii) = idx(update);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join(idx &update, const idx &input) const {
    if (input > update) update = input;
  }
};

template <typename out_array_t, typename in_array_t, typename scalar_1, typename scalar_2, typename MyExecSpace>
void a_times_x_plus_b(typename in_array_t::value_type num_elements, in_array_t out_arr, in_array_t in_arr, scalar_1 a,
                      scalar_2 b) {
  kk_a_times_x_plus_b<out_array_t, in_array_t, scalar_1, scalar_2, MyExecSpace>(num_elements, out_arr, in_arr, a, b);
}

template <typename out_array_type, typename in_array_type, typename MyExecSpace>
void modular_view(typename in_array_type::value_type num_elements, out_array_type out_arr, in_array_type in_arr,
                  int mod_factor_) {
  kk_modular_view<out_array_type, in_array_type, MyExecSpace>(num_elements, out_arr, in_arr, mod_factor_);
}

template <typename array_type>
struct LinearInitialization {
  typedef typename array_type::value_type idx;
  array_type array_sum;
  LinearInitialization(array_type arr_) : array_sum(arr_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii) const { array_sum(ii) = ii; }
};
template <typename array_type, typename MyExecSpace>
void linear_init(typename array_type::value_type num_elements, array_type arr) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_for("KokkosKernels::Common::LinearInit", my_exec_space(0, num_elements),
                       LinearInitialization<array_type>(arr));
}

template <typename forward_array_type, typename MyExecSpace>
void remove_zeros_in_xadj_vector(typename forward_array_type::value_type num_elements, forward_array_type arr) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_scan("KokkosKernels::Common::RemoveZerosInXadjVector", my_exec_space(0, num_elements),
                        PropogataMaxValstoZeros<forward_array_type>(arr));
}

template <typename forward_array_type, typename reverse_array_type>
struct FillReverseBegins {
  const forward_array_type &forward_map;  // vertex to colors
  reverse_array_type &reverse_map_xadj;   // colors to vertex xadj

  FillReverseBegins(const forward_array_type &forward_map_,  // vertex to colors
                    reverse_array_type &reverse_map_xadj_    // colors to vertex xadj
                    )
      : forward_map(forward_map_), reverse_map_xadj(reverse_map_xadj_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii) const {
    typename forward_array_type::value_type prev_col = forward_map(ii - 1);
    typename forward_array_type::value_type cur_col  = forward_map(ii);
    while (prev_col < cur_col) {
      prev_col += 1;
      forward_map(prev_col) = ii + 1;
    }
  }
};

template <typename forward_map_type, typename reverse_map_type>
struct Reverse_Map_Scale_Init {
  typedef typename forward_map_type::value_type forward_type;
  typedef typename reverse_map_type::value_type reverse_type;
  forward_map_type forward_map;
  reverse_map_type reverse_map_xadj;

  const reverse_type multiply_shift_for_scale;
  const reverse_type division_shift_for_bucket;

  Reverse_Map_Scale_Init(forward_map_type forward_map_, reverse_map_type reverse_xadj_,
                         reverse_type multiply_shift_for_scale_, reverse_type division_shift_for_bucket_)
      : forward_map(forward_map_),
        reverse_map_xadj(reverse_xadj_),
        multiply_shift_for_scale(multiply_shift_for_scale_),
        division_shift_for_bucket(division_shift_for_bucket_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &ii) const {
    typedef typename std::remove_reference<decltype(reverse_map_xadj(0))>::type atomic_incr_type;
    forward_type fm = forward_map[ii];
    fm              = fm << multiply_shift_for_scale;
    fm += ii >> division_shift_for_bucket;
    Kokkos::atomic_fetch_add(&(reverse_map_xadj(fm)), atomic_incr_type(1));
  }
};

template <typename forward_map_type, typename reverse_map_type>
struct Fill_Reverse_Scale_Map {
  typedef typename forward_map_type::value_type forward_type;
  typedef typename reverse_map_type::value_type reverse_type;
  forward_map_type forward_map;
  reverse_map_type reverse_map_xadj;
  reverse_map_type reverse_map_adj;

  const reverse_type multiply_shift_for_scale;
  const reverse_type division_shift_for_bucket;

  Fill_Reverse_Scale_Map(forward_map_type forward_map_, reverse_map_type reverse_map_xadj_,
                         reverse_map_type reverse_map_adj_, reverse_type multiply_shift_for_scale_,
                         reverse_type division_shift_for_bucket_)
      : forward_map(forward_map_),
        reverse_map_xadj(reverse_map_xadj_),
        reverse_map_adj(reverse_map_adj_),
        multiply_shift_for_scale(multiply_shift_for_scale_),
        division_shift_for_bucket(division_shift_for_bucket_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &ii) const {
    typedef typename std::remove_reference<decltype(reverse_map_xadj(0))>::type atomic_incr_type;
    forward_type fm = forward_map[ii];

    fm = fm << multiply_shift_for_scale;
    fm += ii >> division_shift_for_bucket;
    const reverse_type future_index = Kokkos::atomic_fetch_add(&(reverse_map_xadj(fm - 1)), atomic_incr_type(1));
    reverse_map_adj(future_index)   = ii;
  }
};

template <typename from_view_t, typename to_view_t>
struct StridedCopy {
  const from_view_t from;
  to_view_t to;
  const size_t stride;
  StridedCopy(const from_view_t from_, to_view_t to_, size_t stride_) : from(from_), to(to_), stride(stride_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &ii) const {
    // std::cout << "ii:" << ii << " ii * stride:" << ii * stride << std::endl;
    to[ii] = from[(ii + 1) * stride - 1];
  }
};

/**
 * \brief Utility function to obtain a reverse map given a map.
 * Input is a map with the number of elements within the map.
 * forward_map[c] = i, where c is a forward elements and forward_map has a size
 * of num_forward_elements. i is the value that c is mapped in the forward map,
 * and the range of that is num_reverse_elements. Output is the reverse_map_xadj
 * and reverse_map_adj such that, all c, forward_map[c] = i, will appear in
 * reverse_map_adj[ reverse_map_xadj[i]: reverse_map_xadj[i+1]) \param:
 * num_forward_elements: the number of elements in the forward map, the size of
 * the forward map. \param: num_reverse_elements: the number of elements that
 * forward map is mapped to. It is the value of max i. \param: forward_map:
 * input forward_map, where forward_map[c] = i. \param: reverse_map_xadj:
 * reverse map xadj, that is it will hold the beginning and end indices on
 * reverse_map_adj such that all values mapped to i will be [
 * reverse_map_xadj[i]: reverse_map_xadj[i+1]) its size will be
 * num_reverse_elements + 1. \param: reverse_map_adj: reverse map adj, holds the
 * values of reverse maps. Its size will be num_forward_elements.
 *
 */
template <typename forward_array_type, typename reverse_array_type, typename MyExecSpace>
void create_reverse_map(MyExecSpace my_exec_space,
                        const typename reverse_array_type::value_type &num_forward_elements,  // num_vertices
                        const typename forward_array_type::value_type &num_reverse_elements,  // num_colors

                        const forward_array_type &forward_map,  // vertex to colors
                        reverse_array_type &reverse_map_xadj,   // colors to vertex xadj
                        reverse_array_type &reverse_map_adj) {  // colros to vertex adj

  typedef typename reverse_array_type::value_type lno_t;
  typedef typename forward_array_type::value_type reverse_lno_t;

  const lno_t MINIMUM_TO_ATOMIC = 64;

  typedef Kokkos::RangePolicy<MyExecSpace> range_policy_t;
  reverse_map_xadj =
      reverse_array_type(Kokkos::view_alloc(my_exec_space, "Reverse Map Xadj"), num_reverse_elements + 1);
  reverse_map_adj = reverse_array_type(Kokkos::view_alloc(my_exec_space, Kokkos::WithoutInitializing, "REVERSE_ADJ"),
                                       num_forward_elements);

  if (num_reverse_elements < MINIMUM_TO_ATOMIC) {
    const lno_t scale_size                = 1024;
    const lno_t multiply_shift_for_scale  = 10;
    const lno_t division_shift_for_bucket = lno_t(ceil(log(double(num_forward_elements) / scale_size) / log(2)));
    // const lno_t bucket_range_size = pow(2, division_shift_for_bucket);

    // coloring indices are base-1. we end up using not using element 1.
    const reverse_lno_t tmp_reverse_size = (num_reverse_elements + 1) << multiply_shift_for_scale;

    reverse_array_type tmp_color_xadj(Kokkos::view_alloc(my_exec_space, "TMP_REVERSE_XADJ"), tmp_reverse_size + 1);

    Reverse_Map_Scale_Init<forward_array_type, reverse_array_type> rmi(
        forward_map, tmp_color_xadj, multiply_shift_for_scale, division_shift_for_bucket);
    Kokkos::parallel_for("KokkosKernels::Common::ReverseMapScaleInit",
                         range_policy_t(my_exec_space, 0, num_forward_elements), rmi);
    my_exec_space.fence();

    inclusive_parallel_prefix_sum<reverse_array_type, MyExecSpace>(my_exec_space, tmp_reverse_size + 1, tmp_color_xadj);
    my_exec_space.fence();

    Kokkos::parallel_for(
        "KokkosKernels::Common::StridedCopy", range_policy_t(my_exec_space, 0, num_reverse_elements + 1),
        StridedCopy<reverse_array_type, reverse_array_type>(tmp_color_xadj, reverse_map_xadj, scale_size));
    my_exec_space.fence();
    Fill_Reverse_Scale_Map<forward_array_type, reverse_array_type> frm(
        forward_map, tmp_color_xadj, reverse_map_adj, multiply_shift_for_scale, division_shift_for_bucket);
    Kokkos::parallel_for("KokkosKernels::Common::FillReverseMap",
                         range_policy_t(my_exec_space, 0, num_forward_elements), frm);
    my_exec_space.fence();
  } else
  // atomic implementation.
  {
    reverse_array_type tmp_color_xadj(
        Kokkos::view_alloc(my_exec_space, Kokkos::WithoutInitializing, "TMP_REVERSE_XADJ"), num_reverse_elements + 1);

    Reverse_Map_Init<forward_array_type, reverse_array_type> rmi(forward_map, reverse_map_xadj);

    Kokkos::parallel_for("KokkosKernels::Common::ReverseMapInit",
                         range_policy_t(my_exec_space, 0, num_forward_elements), rmi);
    my_exec_space.fence();
    // print_1Dview(reverse_map_xadj);

    inclusive_parallel_prefix_sum<reverse_array_type, MyExecSpace>(my_exec_space, num_reverse_elements + 1,
                                                                   reverse_map_xadj);
    Kokkos::deep_copy(my_exec_space, tmp_color_xadj, reverse_map_xadj);
    my_exec_space.fence();
    Fill_Reverse_Map<forward_array_type, reverse_array_type> frm(forward_map, tmp_color_xadj, reverse_map_adj);
    Kokkos::parallel_for("KokkosKernels::Common::FillReverseMap",
                         range_policy_t(my_exec_space, 0, num_forward_elements), frm);
    my_exec_space.fence();
  }
}

template <typename forward_array_type, typename reverse_array_type,
          typename MyExecSpace>
void create_reverse_map(const typename reverse_array_type::value_type &num_forward_elements,  // num_vertices
                        const typename forward_array_type::value_type &num_reverse_elements,  // num_colors

                        const forward_array_type &forward_map,  // vertex to colors
                        reverse_array_type &reverse_map_xadj,   // colors to vertex xadj
                        reverse_array_type &reverse_map_adj) {
  MyExecSpace my_exec_space;
  return create_reverse_map(my_exec_space, num_forward_elements, num_reverse_elements, forward_map, reverse_map_xadj,
                            reverse_map_adj);
}

template <typename value_array_type, typename out_value_array_type, typename idx_array_type>
struct PermuteVector {
  typedef typename idx_array_type::value_type idx;
  value_array_type old_vector;
  out_value_array_type new_vector;
  idx_array_type old_to_new_mapping;
  idx mapping_size;
  PermuteVector(value_array_type old_vector_, out_value_array_type new_vector_, idx_array_type old_to_new_mapping_)
      : old_vector(old_vector_),
        new_vector(new_vector_),
        old_to_new_mapping(old_to_new_mapping_),
        mapping_size(old_to_new_mapping_.extent(0)) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const idx &ii) const {
    idx mapping = ii;
    if (ii < mapping_size) mapping = old_to_new_mapping[ii];
    for (idx j = 0; j < static_cast<idx>(new_vector.extent(1)); j++) {
      new_vector.access(mapping, j) = old_vector.access(ii, j);
    }
  }
};

template <typename value_array_type, typename out_value_array_type, typename idx_array_type, typename MyExecSpace>
void permute_vector(MyExecSpace my_exec_space, typename idx_array_type::value_type num_elements,
                    idx_array_type &old_to_new_index_map, value_array_type &old_vector,
                    out_value_array_type &new_vector) {
  using range_policy_t = Kokkos::RangePolicy<MyExecSpace>;

  Kokkos::parallel_for("KokkosKernels::Common::PermuteVector", range_policy_t(my_exec_space, 0, num_elements),
                       PermuteVector<value_array_type, out_value_array_type, idx_array_type>(old_vector, new_vector,
                                                                                             old_to_new_index_map));
}

template <typename value_array_type, typename out_value_array_type, typename idx_array_type, typename MyExecSpace>
void permute_vector(typename idx_array_type::value_type num_elements, idx_array_type &old_to_new_index_map,
                    value_array_type &old_vector, out_value_array_type &new_vector) {
  permute_vector(MyExecSpace(), num_elements, old_to_new_index_map, old_vector, new_vector);
}

template <typename value_array_type, typename out_value_array_type, typename idx_array_type>
struct PermuteBlockVector {
  typedef typename idx_array_type::value_type idx;
  int block_size;
  value_array_type old_vector;
  out_value_array_type new_vector;
  idx_array_type old_to_new_mapping;
  idx mapping_size;
  PermuteBlockVector(int block_size_, value_array_type old_vector_, out_value_array_type new_vector_,
                     idx_array_type old_to_new_mapping_)
      : block_size(block_size_),
        old_vector(old_vector_),
        new_vector(new_vector_),
        old_to_new_mapping(old_to_new_mapping_),
        mapping_size(old_to_new_mapping_.extent(0)) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const idx &ii) const {
    idx mapping = ii;
    if (ii < mapping_size) mapping = old_to_new_mapping[ii];
    for (idx j = 0; j < static_cast<idx>(new_vector.extent(1)); j++) {
      for (int i = 0; i < block_size; ++i) {
        new_vector.access(mapping * block_size + i, j) = old_vector.access(ii * block_size + i, j);
      }
    }
  }
};

template <typename value_array_type, typename out_value_array_type, typename idx_array_type, typename MyExecSpace>
void permute_block_vector(MyExecSpace my_exec_space, typename idx_array_type::value_type num_elements, int block_size,
                          idx_array_type &old_to_new_index_map, value_array_type &old_vector,
                          out_value_array_type &new_vector) {
  using range_policy_t = Kokkos::RangePolicy<MyExecSpace>;
  Kokkos::parallel_for("KokkosKernels::Common::PermuteVector", range_policy_t(my_exec_space, 0, num_elements),
                       PermuteBlockVector<value_array_type, out_value_array_type, idx_array_type>(
                           block_size, old_vector, new_vector, old_to_new_index_map));
}

template <typename value_array_type, typename out_value_array_type, typename idx_array_type, typename MyExecSpace>
void permute_block_vector(typename idx_array_type::value_type num_elements, int block_size,
                          idx_array_type &old_to_new_index_map, value_array_type &old_vector,
                          out_value_array_type &new_vector) {
  permute_block_vector(MyExecSpace(), num_elements, block_size, old_to_new_index_map, old_vector, new_vector);
}

// TODO BMK: clean this up by removing 1st argument. It is unused but
// its name gives the impression that only num_elements of the vector are
// zeroed, when really it's always the whole thing.
template <class ExecSpaceIn, typename value_array_type>
void zero_vector(ExecSpaceIn &exec_space_in, typename value_array_type::value_type /* num_elements */,
                 value_array_type &vector) {
  typedef typename value_array_type::non_const_value_type val_type;
  Kokkos::deep_copy(exec_space_in, vector, Kokkos::ArithTraits<val_type>::zero());
  exec_space_in.fence();
}

template <typename value_array_type, typename MyExecSpace>
void zero_vector(typename value_array_type::value_type /* num_elements */, value_array_type &vector) {
  using ne_tmp_t  = typename value_array_type::value_type;
  ne_tmp_t ne_tmp = ne_tmp_t(0);
  MyExecSpace my_exec_space;
  zero_vector(my_exec_space, ne_tmp, vector);
}

template <typename v1, typename v2, typename v3>
struct MarkDuplicateSortedKeyValuePairs {
  v1 keys;
  v2 vals;
  v3 prefix_sum;
  typename v1::size_type overall_size;
  MarkDuplicateSortedKeyValuePairs(v1 keys_, v2 vals_, v3 prefix_sum_, typename v1::size_type overall_size_)
      : keys(keys_), vals(vals_), prefix_sum(prefix_sum_), overall_size(overall_size_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i, typename v3::value_type &num_result) const {
    typename v1::value_type my_key = keys(i);
    typename v2::value_type my_val = vals(i);

    if ((my_key != 0 && my_val != 0) && ((i + 1 >= overall_size) || (my_key != keys(i + 1) || my_val != vals(i + 1)))) {
      prefix_sum(i) = 1;
      num_result += 1;
    }
  }
};

template <typename v1, typename v2, typename v3, typename v4, typename v5>
struct FillSymmetricCSR {
  v1 keys;
  v2 vals;
  v3 prefix_sum;
  typename v3::size_type array_size;
  v4 out_xadj;
  v5 out_adj;
  FillSymmetricCSR(v1 keys_, v2 vals_, v3 prefix_sum_, typename v3::size_type array_size_, v4 out_xadj_, v5 out_adj_)
      : keys(keys_),
        vals(vals_),
        prefix_sum(prefix_sum_),
        array_size(array_size_),
        out_xadj(out_xadj_),
        out_adj(out_adj_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i) const {
    typename v3::value_type my_pos = prefix_sum(i);

    if (i + 1 >= array_size) {
      typename v2::value_type my_val = vals(i);
      typename v1::value_type my_key = keys(i);
      out_adj(my_pos)                = my_val - 1;
      out_xadj(my_key)               = my_pos + 1;
    } else {
      typename v3::value_type next_pos = prefix_sum(i + 1);
      if (my_pos != next_pos) {
        typename v2::value_type my_val   = vals(i);
        typename v1::value_type my_key   = keys(i);
        typename v1::value_type next_key = keys(i + 1);
        out_adj(my_pos)                  = my_val - 1;
        if (my_key != next_key) {
          out_xadj(my_key) = my_pos + 1;
        }
      }
    }
  }
};

template <typename in_lno_row_view_t, typename in_lno_nnz_view_t, typename out_lno_nnz_view_t, typename MyExecSpace>
void symmetrize_and_get_lower_diagonal_edge_list(typename in_lno_nnz_view_t::value_type num_rows_to_symmetrize,
                                                 in_lno_row_view_t xadj, in_lno_nnz_view_t adj,
                                                 out_lno_nnz_view_t &sym_srcs, out_lno_nnz_view_t &sym_dsts_) {
  typedef typename in_lno_row_view_t::non_const_value_type idx;

  idx nnz = adj.extent(0);

  // idx_out_edge_array_type tmp_srcs("tmpsrc", nnz * 2);
  // idx_out_edge_array_type tmp_dsts("tmpdst",nnz * 2);

  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy;
  typedef typename team_policy::member_type team_member_t;

  // typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;

  // TODO: Should change this to temporary memory space?
  typedef Kokkos::UnorderedMap<Kokkos::pair<idx, idx>, void, MyExecSpace> hashmap_t;

  out_lno_nnz_view_t pre_pps_("pre_pps", num_rows_to_symmetrize + 1);

  idx num_symmetric_edges = 0;
  {
    hashmap_t umap(nnz);
    umap.clear();
    umap.end_erase();
    FillSymmetricLowerEdgesHashMap<in_lno_row_view_t, in_lno_nnz_view_t, hashmap_t, out_lno_nnz_view_t, team_member_t>
        fse(num_rows_to_symmetrize, xadj, adj, umap, pre_pps_);

    int teamSizeMax = 0;
    int vector_size = 0;

    get_suggested_vector_size<idx, MyExecSpace>(vector_size, xadj.extent(0) - 1, nnz);

    teamSizeMax = get_suggested_team_size<team_policy>(fse, vector_size);
    // std::cout << "max_allowed_team_size:" << max_allowed_team_size << " vs:"
    // << vector_size << " tsm:" << teamSizeMax<< std::endl;

    team_policy pol((num_rows_to_symmetrize + teamSizeMax - 1) / teamSizeMax, teamSizeMax, vector_size);
    Kokkos::parallel_for("KokkosKernels::Common::SymmetrizeAndGetLowerDiagonalEdgeList::S0", pol,
                         fse /*, num_symmetric_edges*/);
    MyExecSpace().fence();
  }

  if (num_rows_to_symmetrize > 0)
    exclusive_parallel_prefix_sum<out_lno_nnz_view_t, MyExecSpace>(num_rows_to_symmetrize + 1, pre_pps_);
  MyExecSpace().fence();

  auto d_sym_edge_size = Kokkos::subview(pre_pps_, num_rows_to_symmetrize);
  auto h_sym_edge_size = Kokkos::create_mirror_view(d_sym_edge_size);
  Kokkos::deep_copy(h_sym_edge_size, d_sym_edge_size);
  num_symmetric_edges = h_sym_edge_size();
  /*
  typename out_lno_nnz_view_t::HostMirror h_sym_edge_size =
  Kokkos::create_mirror_view (pre_pps_);

  Kokkos::deep_copy (h_sym_edge_size , pre_pps_);
  num_symmetric_edges = h_sym_edge_size(h_sym_edge_size.extent(0) - 1);
  */

  sym_srcs  = out_lno_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "sym_srcs"), num_symmetric_edges);
  sym_dsts_ = out_lno_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "sym_dsts_"), num_symmetric_edges);
  MyExecSpace().fence();
  {
    hashmap_t umap(nnz);
    FillSymmetricEdgeList_HashMap<in_lno_row_view_t, in_lno_nnz_view_t, hashmap_t, out_lno_nnz_view_t,
                                  out_lno_nnz_view_t, team_member_t>
        FSCH(num_rows_to_symmetrize, xadj, adj, umap, sym_srcs, sym_dsts_, pre_pps_);

    int teamSizeMax = 0;
    int vector_size = 0;

    get_suggested_vector_size<idx, MyExecSpace>(vector_size, xadj.extent(0) - 1, nnz);

    teamSizeMax = get_suggested_team_size<team_policy>(FSCH, vector_size);

    team_policy pol((num_rows_to_symmetrize + teamSizeMax - 1) / teamSizeMax, teamSizeMax, vector_size);
    Kokkos::parallel_for("KokkosKernels::Common::SymmetrizeAndGetLowerDiagonalEdgeList::S1", pol, FSCH);
    MyExecSpace().fence();
  }
}

template <typename in_lno_row_view_t, typename in_lno_nnz_view_t, typename out_lno_row_view_t,
          typename out_lno_nnz_view_t, typename MyExecSpace>
void symmetrize_graph_symbolic_hashmap(typename in_lno_row_view_t::value_type num_rows_to_symmetrize,
                                       in_lno_row_view_t xadj, in_lno_nnz_view_t adj, out_lno_row_view_t &sym_xadj,
                                       out_lno_nnz_view_t &sym_adj) {
  typedef typename in_lno_row_view_t::non_const_value_type idx;

  idx nnz = adj.extent(0);

  // idx_out_edge_array_type tmp_srcs("tmpsrc", nnz * 2);
  // idx_out_edge_array_type tmp_dsts("tmpdst",nnz * 2);

  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy;
  typedef typename team_policy::member_type team_member_t;

  // typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;

  // TODO: Should change this to temporary memory space?
  typedef Kokkos::UnorderedMap<Kokkos::pair<idx, idx>, void, MyExecSpace> hashmap_t;

  out_lno_row_view_t pre_pps_("pre_pps", num_rows_to_symmetrize + 1);

  idx num_symmetric_edges = 0;
  {
    hashmap_t umap(nnz);
    umap.clear();
    umap.end_erase();
    FillSymmetricEdgesHashMap<in_lno_row_view_t, in_lno_nnz_view_t, hashmap_t, out_lno_row_view_t, team_member_t> fse(
        num_rows_to_symmetrize, xadj, adj, umap, pre_pps_);

    int teamSizeMax = 0;
    int vector_size = 0;

    get_suggested_vector_size<idx, MyExecSpace>(vector_size, xadj.extent(0) - 1, nnz);

    teamSizeMax = get_suggested_team_size<team_policy>(fse, vector_size);

    team_policy pol((num_rows_to_symmetrize + teamSizeMax - 1) / teamSizeMax, teamSizeMax, vector_size);
    Kokkos::parallel_for("KokkosKernels::Common::SymmetrizeGraphSymbolicHashMap::S0", pol,
                         fse /*, num_symmetric_edges*/);
    MyExecSpace().fence();
  }

  if (num_rows_to_symmetrize > 0)
    exclusive_parallel_prefix_sum<out_lno_row_view_t, MyExecSpace>(num_rows_to_symmetrize + 1, pre_pps_);
  MyExecSpace().fence();

  // out_lno_row_view_t d_sym_edge_size = Kokkos::subview(pre_pps_,
  // num_rows_to_symmetrize, num_rows_to_symmetrize );
  typename out_lno_row_view_t::HostMirror h_sym_edge_size = Kokkos::create_mirror_view(pre_pps_);

  Kokkos::deep_copy(h_sym_edge_size, pre_pps_);
  num_symmetric_edges = h_sym_edge_size(h_sym_edge_size.extent(0) - 1);

  sym_adj = out_lno_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "sym_adj"), num_symmetric_edges);
  MyExecSpace().fence();
  sym_xadj =
      out_lno_row_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "sym_xadj"), num_rows_to_symmetrize + 1);
  Kokkos::deep_copy(sym_xadj, pre_pps_);
  {
    hashmap_t umap(nnz);
    FillSymmetricCRS_HashMap<in_lno_row_view_t, in_lno_nnz_view_t, hashmap_t, out_lno_row_view_t, out_lno_nnz_view_t,
                             team_member_t>
        FSCH(num_rows_to_symmetrize, xadj, adj, umap, pre_pps_, sym_adj);

    int teamSizeMax = 0;
    int vector_size = 0;

    get_suggested_vector_size<idx, MyExecSpace>(vector_size, xadj.extent(0) - 1, nnz);

    teamSizeMax = get_suggested_team_size<team_policy>(FSCH, vector_size);

    team_policy pol((num_rows_to_symmetrize + teamSizeMax - 1) / teamSizeMax, teamSizeMax, vector_size);
    Kokkos::parallel_for("KokkosKernels::Common::SymmetrizeGraphSymbolicHashMap::S1", pol, FSCH);
    MyExecSpace().fence();
  }

  MyExecSpace().fence();
}

template <typename from_vector, typename to_vector, typename MyExecSpace>
void copy_vector(size_t num_elements, from_vector from, to_vector to) {
  kk_copy_vector<from_vector, to_vector, MyExecSpace>(num_elements, from, to);
}

template <typename from_vector, typename to_vector>
struct CopyView {
  from_vector from;
  to_vector to;

  CopyView(from_vector &from_, to_vector to_) : from(from_), to(to_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i) const { to(i) = from(i); }
};
template <typename from_vector, typename to_vector, typename MyExecSpace>
void copy_view(size_t num_elements, from_vector from, to_vector to) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_for("KokkosKernels::Common::CopyView", my_exec_space(0, num_elements),
                       CopyView<from_vector, to_vector>(from, to));
}

template <typename from_view>
void safe_device_to_host_deep_copy(size_t num_elements, from_view from, typename from_view::HostMirror to) {
  typedef typename from_view::value_type scalar_t;
  typedef typename from_view::device_type device_t;

  typedef Kokkos::View<scalar_t *, device_t> unstrided_from_view_t;
  unstrided_from_view_t unstrided_from("unstrided", num_elements);

  copy_view<from_view, unstrided_from_view_t, typename device_t::execution_space>(num_elements, from, unstrided_from);

  Kokkos::fence();

  typedef typename unstrided_from_view_t::HostMirror host_unstrided_from_view_t;
  host_unstrided_from_view_t h_unstrided_from = Kokkos::create_mirror_view(unstrided_from);

  Kokkos::deep_copy(h_unstrided_from, unstrided_from);
  Kokkos::fence();

  copy_view<host_unstrided_from_view_t, typename from_view::HostMirror,
            typename host_unstrided_from_view_t::device_type::execution_space>(num_elements, h_unstrided_from, to);

  Kokkos::fence();
}

template <typename to_view>
void safe_host_to_device_deep_copy(size_t num_elements, typename to_view::HostMirror from, to_view to) {
  typedef typename to_view::value_type scalar_t;
  typedef typename to_view::device_type device_t;

  typedef typename to_view::HostMirror::device_type h_device_t;

  typedef Kokkos::View<scalar_t *, h_device_t> host_unstrided_view_t;
  typedef Kokkos::View<scalar_t *, device_t> device_unstrided_view_t;

  host_unstrided_view_t host_unstrided_from("unstrided", num_elements);
  device_unstrided_view_t device_unstrided_to("unstrided", num_elements);

  copy_view<typename to_view::HostMirror, host_unstrided_view_t, typename h_device_t::execution_space>(
      num_elements, from, host_unstrided_from);

  Kokkos::fence();
  Kokkos::deep_copy(device_unstrided_to, host_unstrided_from);
  Kokkos::fence();

  copy_view<device_unstrided_view_t, to_view, typename device_t::execution_space>(num_elements, device_unstrided_to,
                                                                                  to);

  Kokkos::fence();
}

template <typename view_type>
struct ReduceSumFunctor {
  view_type view_to_reduce;

  ReduceSumFunctor(view_type view_to_reduce_) : view_to_reduce(view_to_reduce_) {}

  void operator()(const size_t &i, typename view_type::non_const_value_type &sum_reduction) const {
    sum_reduction += view_to_reduce(i);
  }
};

template <typename view_type, typename MyExecSpace>
void view_reduce_sum(size_t num_elements, view_type view_to_reduce,
                     typename view_type::non_const_value_type &sum_reduction) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_reduce("KokkosKernels::Common::ViewReduceSum", my_exec_space(0, num_elements),
                          ReduceSumFunctor<view_type>(view_to_reduce), sum_reduction);
}

template <typename view_type, typename MyExecSpace>
void view_reduce_max(size_t num_elements, view_type view_to_reduce,
                     typename view_type::non_const_value_type &max_reduction) {
  kk_view_reduce_max<view_type, MyExecSpace>(num_elements, view_to_reduce, max_reduction);
}

template <typename view_type, typename MyExecSpace>
void view_reduce_max(const MyExecSpace &exec, size_t num_elements, view_type view_to_reduce,
                     typename view_type::non_const_value_type &max_reduction) {
  kk_view_reduce_max<view_type, MyExecSpace>(exec, num_elements, view_to_reduce, max_reduction);
}

template <typename size_type>
struct ReduceRowSizeFunctor {
  const size_type *rowmap_view_begins;
  const size_type *rowmap_view_ends;
  const size_type min_val;
  ReduceRowSizeFunctor(const size_type *rb, const size_type *re)
      : rowmap_view_begins(rb), rowmap_view_ends(re), min_val(0) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i, size_type &max_reduction) const {
    size_type val = rowmap_view_ends[i] - rowmap_view_begins[i];
    if (max_reduction < val) {
      max_reduction = val;
    }
  }
  KOKKOS_INLINE_FUNCTION
  void join(size_type &dst, const size_type &src) const {
    if (dst < src) {
      dst = src;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(size_type &dst) const {
    // The identity under max is -Inf.
    // Kokkos does not come with a portable way to access
    // floating -point Inf and NaN. Trilinos does , however;
    // see Kokkos :: ArithTraits in the Tpetra package.
    dst = min_val;
  }
};

// view has num_rows+1 elements.
template <typename size_type, typename MyExecSpace>
void kk_view_reduce_max_row_size(MyExecSpace my_exec_space, const size_t num_rows, const size_type *rowmap_view_begins,
                                 const size_type *rowmap_view_ends, size_type &max_row_size) {
  typedef Kokkos::RangePolicy<MyExecSpace> range_policy_t;
  Kokkos::parallel_reduce("KokkosKernels::Common::ViewReduceMaxRowSize", range_policy_t(my_exec_space, 0, num_rows),
                          ReduceRowSizeFunctor<size_type>(rowmap_view_begins, rowmap_view_ends), max_row_size);
}

// view has num_rows+1 elements.
template <typename size_type, typename MyExecSpace>
void kk_view_reduce_max_row_size(const size_t num_rows, const size_type *rowmap_view_begins,
                                 const size_type *rowmap_view_ends, size_type &max_row_size) {
  return kk_view_reduce_max_row_size(MyExecSpace(), num_rows, rowmap_view_begins, rowmap_view_ends, max_row_size);
}

template <typename view_type>
struct ReduceMaxRowFunctor {
  view_type rowmap_view;
  typedef typename view_type::non_const_value_type value_type;
  const value_type min_val;
  ReduceMaxRowFunctor(view_type rowmap_view_) : rowmap_view(rowmap_view_), min_val(0) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i, value_type &max_reduction) const {
    value_type val = rowmap_view(i + 1) - rowmap_view(i);
    if (max_reduction < val) {
      max_reduction = val;
    }
  }
  KOKKOS_INLINE_FUNCTION
  void join(value_type &dst, const value_type &src) const {
    if (dst < src) {
      dst = src;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type &dst) const {
    // The identity under max is -Inf.
    // Kokkos does not come with a portable way to access
    // floating -point Inf and NaN. Trilinos does , however;
    // see Kokkos :: ArithTraits in the Tpetra package.
    dst = min_val;
  }
};

// view has num_rows+1 elements.
template <typename view_type, typename MyExecSpace>
void view_reduce_maxsizerow(size_t num_rows, view_type rowmap_view,
                            typename view_type::non_const_value_type &max_reduction) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_reduce("KokkosKernels::Common::ViewReduceMaxSizeRow", my_exec_space(0, num_rows),
                          ReduceMaxRowFunctor<view_type>(rowmap_view), max_reduction);
}

template <typename view_type1, typename view_type2>
struct IsEqualFunctor {
  view_type1 view1;
  view_type2 view2;

  IsEqualFunctor(view_type1 view1_, view_type2 view2_) : view1(view1_), view2(view2_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i, int &is_equal) const {
    if (view1(i) != view2(i)) {
      // std::cout << "i:" << i << "view1:" << view1(i) << " view2:" << view2(i)
      // << std::endl; printf("i:%d v1:")
      is_equal = 0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join(int &dst, const int &src) const { dst = dst & src; }
  KOKKOS_INLINE_FUNCTION
  void init(int &dst) const { dst = 1; }
};
template <typename view_type1, typename view_type2, typename MyExecSpace>
bool isSame(size_t num_elements, view_type1 view1, view_type2 view2) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  int issame = 1;
  Kokkos::parallel_reduce("KokkosKernels::Common::isSame", my_exec_space(0, num_elements),
                          IsEqualFunctor<view_type1, view_type2>(view1, view2), issame);
  MyExecSpace().fence();
  return issame;
}

template <typename a_view_t, typename b_view_t, typename size_type>
struct MaxHeap {
  a_view_t heap_keys;
  b_view_t heap_values;
  size_type max_size;
  size_type current_size;

  MaxHeap(a_view_t heap_keys_, b_view_t heap_values_, size_type max_size_)
      : heap_keys(heap_keys_), heap_values(heap_values_), max_size(max_size_), current_size(0) {}

  KOKKOS_INLINE_FUNCTION
  void insert(typename a_view_t::value_type &key, typename b_view_t::value_type &val) {
    for (size_type i = 0; i < current_size; ++i) {
      if (key == heap_keys(i)) {
        heap_values(i) = heap_values(i) & val;
        return;
      }
    }
    heap_keys(current_size)     = key;
    heap_values(current_size++) = val;
  }
};

template <typename in_view_t, typename MyExecSpace>
struct InitScalar {
  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy_t;
  typedef typename team_policy_t::member_type team_member_t;

  typedef typename in_view_t::non_const_value_type nnz_lno_t;
  typedef typename in_view_t::size_type size_type;

  in_view_t view_to_init;
  size_type num_elements;
  size_type team_row_chunk_size;
  nnz_lno_t init_val;

  InitScalar(size_type num_elements_, in_view_t view_to_init_, size_type chunk_size_, nnz_lno_t init_val_)
      : num_elements(num_elements_),
        view_to_init(view_to_init_),
        team_row_chunk_size(chunk_size_),
        init_val(init_val_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member_t &teamMember) const {
    // const nnz_lno_t row_index = teamMember.league_rank() *
    // team_row_chunk_size;

    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, num_elements);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_ind) { view_to_init[row_ind] = init_val; });
  }
};
template <typename in_row_view_t, typename MyExecSpace>
void init_view_withscalar(typename in_row_view_t::size_type num_elements, in_row_view_t arr,
                          typename in_row_view_t::size_type team_size,
                          typename in_row_view_t::non_const_value_type init_val) {
  typename in_row_view_t::size_type chunk_size = num_elements / team_size;
  typedef InitScalar<in_row_view_t, MyExecSpace> InitScalar_t;
  InitScalar_t tm(num_elements, arr, chunk_size, init_val);
  typedef typename InitScalar_t::team_policy_t tcp_t;
  int vector_size = 1;

  Kokkos::Timer timer1;
  Kokkos::parallel_for("KokkosKernels::Common::InitViewWithScalar",
                       tcp_t(num_elements / chunk_size + 1, team_size, vector_size), tm);
  MyExecSpace().fence();
}

template <typename scalar_t, int N>
struct array_sum_reduce {
  static_assert(N <= 8, "array_sum_reduce has only been tested up to N=8");
  using ValueType = array_sum_reduce<scalar_t, N>;
  // Workaround for https://github.com/kokkos/kokkos/issues/5860
  static constexpr int N_internal =
      ((N == 3 || N == 5 || N == 7) && std::is_same<scalar_t, Kokkos::Experimental::half_t>::value &&
       sizeof(Kokkos::Experimental::half_t) == 2)
          ? (N + 1)
          : N;

  scalar_t data[N_internal];
  KOKKOS_INLINE_FUNCTION
  array_sum_reduce() {
    // Initialize all the elements, even those at index >= N (prevent valgrind
    // warnings, etc.)
    for (int i = 0; i < N_internal; i++) data[i] = scalar_t();
  }

  KOKKOS_INLINE_FUNCTION  // add operator
      array_sum_reduce &
      operator+=(const ValueType &src) {
    // Don't bother summing elements >= N though as they will never be used
    for (int i = 0; i < N; i++) data[i] += src.data[i];
    return *this;
  }
};

template <typename T, typename InPtr>
KOKKOS_INLINE_FUNCTION T *alignPtrTo(InPtr *p) {
  // ugly but computationally free and the "right" way to do this in C++
  const std::uintptr_t ptrVal = reinterpret_cast<std::uintptr_t>(p);
  // ptrVal + (align - 1) lands inside the next valid aligned scalar_t,
  // and the mask produces the start of that scalar_t.
  const std::uintptr_t ptrValNew = (ptrVal + alignof(T) - 1) & (~(alignof(T) - 1));
  return reinterpret_cast<T *>(reinterpret_cast<char *>(const_cast<std::remove_cv_t<InPtr> *>(p)) +
                               (ptrValNew - ptrVal));
}

}  // namespace Impl
}  // namespace KokkosKernels

// Define the identity for array_sum_reduce
namespace Kokkos {
template <typename scalar_t, int N>
struct reduction_identity<KokkosKernels::Impl::array_sum_reduce<scalar_t, N>> {
  typedef KokkosKernels::Impl::array_sum_reduce<scalar_t, N> T;
  KOKKOS_FORCEINLINE_FUNCTION static T sum() {
    // default constructor default-initializes each element (this should always
    // be 0)
    return T();
  }
};
}  // namespace Kokkos

#endif
