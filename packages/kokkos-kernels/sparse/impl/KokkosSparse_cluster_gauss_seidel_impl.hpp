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

#ifndef _KOKKOSCGSIMP_HPP
#define _KOKKOSCGSIMP_HPP

#include "KokkosKernels_Utils.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Bitset.hpp>
#include <Kokkos_Sort.hpp>
#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosKernels_BitUtils.hpp"
#include "KokkosKernels_SimpleUtils.hpp"
#include "KokkosSparse_partitioning_impl.hpp"
#include "KokkosGraph_MIS2.hpp"
#include "KokkosGraph_ExplicitCoarsening.hpp"

namespace KokkosSparse {
namespace Impl {

template <typename HandleType, typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_>
class ClusterGaussSeidel {
 public:
  typedef lno_row_view_t_ in_lno_row_view_t;
  typedef lno_nnz_view_t_ in_lno_nnz_view_t;
  typedef scalar_nnz_view_t_ in_scalar_nnz_view_t;

  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;

  typedef typename in_lno_row_view_t::non_const_value_type row_lno_t;

  typedef typename HandleType::size_type size_type;
  typedef typename HandleType::nnz_lno_t nnz_lno_t;
  typedef typename HandleType::nnz_scalar_t nnz_scalar_t;

  static_assert(std::is_same<size_type, typename in_lno_row_view_t::non_const_value_type>::value,
                "ClusterGaussSeidel: Handle's size_type does not match input rowmap's "
                "element type.");
  static_assert(std::is_same<nnz_lno_t, typename in_lno_nnz_view_t::non_const_value_type>::value,
                "ClusterGaussSeidel: Handle's nnz_lno_t does not match input entries's "
                "element type.");

  typedef typename in_lno_row_view_t::const_type const_lno_row_view_t;
  typedef typename in_lno_row_view_t::non_const_type non_const_lno_row_view_t;

  typedef typename lno_nnz_view_t_::const_type const_lno_nnz_view_t;
  typedef typename lno_nnz_view_t_::non_const_type non_const_lno_nnz_view_t;

  typedef typename scalar_nnz_view_t_::const_type const_scalar_nnz_view_t;
  typedef typename scalar_nnz_view_t_::non_const_type non_const_scalar_nnz_view_t;

  typedef typename HandleType::row_lno_temp_work_view_t row_lno_temp_work_view_t;
  typedef typename HandleType::row_lno_persistent_work_view_t row_lno_persistent_work_view_t;
  typedef
      typename HandleType::row_lno_persistent_work_host_view_t row_lno_persistent_work_host_view_t;  // Host view type

  typedef typename HandleType::nnz_lno_temp_work_view_t nnz_lno_temp_work_view_t;
  typedef typename HandleType::nnz_lno_persistent_work_view_t nnz_lno_persistent_work_view_t;
  typedef
      typename HandleType::nnz_lno_persistent_work_host_view_t nnz_lno_persistent_work_host_view_t;  // Host view type

  typedef typename HandleType::scalar_temp_work_view_t scalar_temp_work_view_t;
  typedef typename HandleType::scalar_persistent_work_view_t scalar_persistent_work_view_t;

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  typedef nnz_lno_t color_t;
  typedef Kokkos::View<color_t*, MyTempMemorySpace> color_view_t;
  typedef Kokkos::Bitset<MyExecSpace> bitset_t;
  typedef Kokkos::ConstBitset<MyExecSpace> const_bitset_t;

  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy_t;
  typedef typename team_policy_t::member_type team_member_t;

 private:
  HandleType* handle;

  // Get the specialized ClusterGaussSeidel handle from the main handle
  typename HandleType::ClusterGaussSeidelHandleType* get_gs_handle() {
    auto* gsHandle = dynamic_cast<typename HandleType::ClusterGaussSeidelHandleType*>(this->handle->get_gs_handle());
    if (!gsHandle) {
      throw std::runtime_error(
          "ClusterGaussSeidel: GS handle has not been created, or is set up "
          "for Point GS.");
    }
    return gsHandle;
  }

  nnz_lno_t num_rows, num_cols;

  const_lno_row_view_t row_map;
  const_lno_nnz_view_t entries;
  const_scalar_nnz_view_t values;

  const_scalar_nnz_view_t given_inverse_diagonal;

  bool have_diagonal_given;
  bool is_symmetric;

  static constexpr nnz_lno_t apply_batch_size = 16;

 public:
  struct PSGS_ForwardTag {};
  struct PSGS_BackwardTag {};

  template <typename x_value_array_type, typename y_value_array_type>
  struct PSGS {
    // CSR storage of the matrix.
    const_lno_row_view_t _xadj;
    const_lno_nnz_view_t _adj;
    const_scalar_nnz_view_t _adj_vals;

    // Input/output vectors, as in Ax = y
    x_value_array_type _Xvector;
    y_value_array_type _Yvector;
    nnz_lno_persistent_work_view_t _color_adj;
    nnz_lno_persistent_work_view_t _cluster_offsets;
    nnz_lno_persistent_work_view_t _cluster_verts;
    scalar_persistent_work_view_t _inverse_diagonal;
    nnz_scalar_t _omega;

    nnz_lno_t _color_set_begin;
    nnz_lno_t _color_set_end;
    bool _forward_direction;

    PSGS(const_lno_row_view_t xadj_, const_lno_nnz_view_t adj_, const_scalar_nnz_view_t adj_vals_,
         x_value_array_type Xvector_, y_value_array_type Yvector_, nnz_lno_persistent_work_view_t color_adj_,
         nnz_lno_persistent_work_view_t cluster_offsets_, nnz_lno_persistent_work_view_t cluster_verts_,
         nnz_scalar_t omega_, scalar_persistent_work_view_t inverse_diagonal_)
        : _xadj(xadj_),
          _adj(adj_),
          _adj_vals(adj_vals_),
          _Xvector(Xvector_),
          _Yvector(Yvector_),
          _color_adj(color_adj_),
          _cluster_offsets(cluster_offsets_),
          _cluster_verts(cluster_verts_),
          _inverse_diagonal(inverse_diagonal_),
          _omega(omega_),
          _color_set_begin(0),
          _color_set_end(0),
          _forward_direction(true) {}

    KOKKOS_FORCEINLINE_FUNCTION
    void rowApply(nnz_scalar_t* sum, const nnz_lno_t row) const {
      size_type row_begin = _xadj(row);
      size_type row_end   = _xadj(row + 1);
      nnz_lno_t num_vecs  = _Xvector.extent(1);
      for (nnz_lno_t batch_start = 0; batch_start < num_vecs; batch_start += apply_batch_size) {
        nnz_lno_t this_batch_size = apply_batch_size;
        if (batch_start + this_batch_size >= num_vecs) this_batch_size = num_vecs - batch_start;
        // the current batch of columns given by: batch_start, this_batch_size
        for (nnz_lno_t i = 0; i < this_batch_size; i++) sum[i] = _Yvector(row, batch_start + i);
        for (size_type adjind = row_begin; adjind < row_end; ++adjind) {
          nnz_lno_t col    = _adj(adjind);
          nnz_scalar_t val = _adj_vals(adjind);
          for (nnz_lno_t i = 0; i < this_batch_size; i++) sum[i] -= val * _Xvector(col, batch_start + i);
        }
        nnz_scalar_t invDiagonalVal = _inverse_diagonal(row);
        for (nnz_lno_t i = 0; i < this_batch_size; i++)
          _Xvector(row, batch_start + i) += _omega * sum[i] * invDiagonalVal;
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const PSGS_ForwardTag, const nnz_lno_t ii) const {
      // color_adj(ii) is a cluster in the current color set.
      nnz_lno_t clusterColorsetIndex = _color_set_begin + ii;
      nnz_lno_t cluster              = _color_adj(clusterColorsetIndex);
      nnz_scalar_t sum[apply_batch_size];
      for (nnz_lno_t j = _cluster_offsets(cluster); j < _cluster_offsets(cluster + 1); j++) {
        rowApply(sum, _cluster_verts(j));
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const PSGS_BackwardTag, const nnz_lno_t ii) const {
      // color_adj(ii) is a cluster in the current color set.
      nnz_lno_t clusterColorsetIndex = _color_set_end - 1 - ii;
      nnz_lno_t cluster              = _color_adj(clusterColorsetIndex);
      nnz_scalar_t sum[apply_batch_size];
      for (nnz_lno_t j = _cluster_offsets(cluster + 1); j > _cluster_offsets(cluster); j--) {
        rowApply(sum, _cluster_verts(j - 1));
      }
    }
  };

  template <typename x_value_array_type, typename y_value_array_type>
  struct Team_PSGS {
    // CSR storage of the matrix
    const_lno_row_view_t _xadj;
    const_lno_nnz_view_t _adj;
    const_scalar_nnz_view_t _adj_vals;

    // X,Y vectors, as in Ax = y
    x_value_array_type _Xvector;
    y_value_array_type _Yvector;
    nnz_lno_t _color_set_begin;
    nnz_lno_t _color_set_end;
    nnz_lno_persistent_work_view_t _color_adj;
    nnz_lno_persistent_work_view_t _cluster_offsets;
    nnz_lno_persistent_work_view_t _cluster_verts;

    //_clusters_per_team tries to reach the same total work per
    // team by dividing the handle's heuristic get_
    nnz_lno_t _clusters_per_team;

    scalar_persistent_work_view_t _inverse_diagonal;

    bool _is_backward;

    nnz_scalar_t _omega;

    Team_PSGS(const_lno_row_view_t xadj_, const_lno_nnz_view_t adj_, const_scalar_nnz_view_t adj_vals_,
              x_value_array_type Xvector_, y_value_array_type Yvector_, nnz_lno_t color_set_begin_,
              nnz_lno_t color_set_end_, nnz_lno_persistent_work_view_t color_adj_,
              nnz_lno_persistent_work_view_t cluster_offsets_, nnz_lno_persistent_work_view_t cluster_verts_,
              scalar_persistent_work_view_t inverse_diagonal_, nnz_lno_t clusters_per_team_,
              nnz_scalar_t omega_ = Kokkos::ArithTraits<nnz_scalar_t>::one())
        : _xadj(xadj_),
          _adj(adj_),
          _adj_vals(adj_vals_),
          _Xvector(Xvector_),
          _Yvector(Yvector_),
          _color_set_begin(color_set_begin_),
          _color_set_end(color_set_end_),
          _color_adj(color_adj_),
          _cluster_offsets(cluster_offsets_),
          _cluster_verts(cluster_verts_),
          _clusters_per_team(clusters_per_team_),
          _inverse_diagonal(inverse_diagonal_),
          _is_backward(false),
          _omega(omega_) {}

    template <int N>
    KOKKOS_INLINE_FUNCTION void runColBatch(const team_member_t& teamMember, nnz_lno_t row, nnz_lno_t colStart) const {
      typedef KokkosKernels::Impl::array_sum_reduce<nnz_scalar_t, N> reducer;
      size_type row_begin = _xadj(row);
      size_type row_end   = _xadj(row + 1);
      reducer sum;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(teamMember, row_end - row_begin),
          [&](size_type i, reducer& lsum) {
            size_type adjind   = row_begin + i;
            nnz_lno_t colIndex = _adj(adjind);
            nnz_scalar_t val   = _adj_vals(adjind);
            for (int j = 0; j < N; j++) lsum.data[j] += val * _Xvector(colIndex, colStart + j);
          },
          sum);
      Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
        nnz_scalar_t invDiagonalVal = _inverse_diagonal(row);
        for (int i = 0; i < N; i++) {
          _Xvector(row, colStart + i) += _omega * (_Yvector(row, colStart + i) - sum.data[i]) * invDiagonalVal;
        }
      });
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t& teamMember) const {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, _clusters_per_team), [&](const nnz_lno_t work) {
        nnz_lno_t ii = _color_set_begin + (teamMember.league_rank() * _clusters_per_team) + work;
        if (ii >= _color_set_end) return;
        nnz_lno_t cluster      = _color_adj(ii);
        nnz_lno_t clusterBegin = _cluster_offsets(cluster);
        nnz_lno_t clusterEnd   = _cluster_offsets(cluster + 1);
        for (nnz_lno_t jcount = 0; jcount < clusterEnd - clusterBegin; jcount++) {
          nnz_lno_t j        = _is_backward ? (clusterEnd - 1 - jcount) : clusterBegin + jcount;
          nnz_lno_t row      = _cluster_verts(j);
          nnz_lno_t num_vecs = _Xvector.extent(1);
          for (nnz_lno_t batch_start = 0; batch_start < num_vecs;) {
            switch (num_vecs - batch_start) {
#define COL_BATCH_CASE(n)                         \
  case n:                                         \
    runColBatch<n>(teamMember, row, batch_start); \
    batch_start += n;                             \
    break;
              COL_BATCH_CASE(1)
              COL_BATCH_CASE(2)
              COL_BATCH_CASE(3)
#undef COL_BATCH_CASE
              default: runColBatch<4>(teamMember, row, batch_start); batch_start += 4;
            }
          }
        }
      });
    }
  };

  /**
   * \brief constructor
   */

  ClusterGaussSeidel(HandleType* handle_, nnz_lno_t num_rows_, nnz_lno_t num_cols_, const_lno_row_view_t row_map_,
                     const_lno_nnz_view_t entries_, const_scalar_nnz_view_t values_)
      : handle(handle_),
        num_rows(num_rows_),
        num_cols(num_cols_),
        row_map(row_map_),
        entries(entries_),
        values(values_),
        have_diagonal_given(false),
        is_symmetric(true) {}

  ClusterGaussSeidel(HandleType* handle_, nnz_lno_t num_rows_, nnz_lno_t num_cols_, const_lno_row_view_t row_map_,
                     const_lno_nnz_view_t entries_, bool is_symmetric_ = true)
      : handle(handle_),
        num_rows(num_rows_),
        num_cols(num_cols_),
        row_map(row_map_),
        entries(entries_),
        values(),
        have_diagonal_given(false),
        is_symmetric(is_symmetric_) {}

  /**
   * \brief constructor
   */
  ClusterGaussSeidel(HandleType* handle_, nnz_lno_t num_rows_, nnz_lno_t num_cols_, const_lno_row_view_t row_map_,
                     const_lno_nnz_view_t entries_, const_scalar_nnz_view_t values_, bool is_symmetric_)
      : handle(handle_),
        num_rows(num_rows_),
        num_cols(num_cols_),
        row_map(row_map_),
        entries(entries_),
        values(values_),
        have_diagonal_given(false),
        is_symmetric(is_symmetric_) {}

  ClusterGaussSeidel(HandleType* handle_, nnz_lno_t num_rows_, nnz_lno_t num_cols_, const_lno_row_view_t row_map_,
                     const_lno_nnz_view_t entries_, const_scalar_nnz_view_t values_,
                     const_scalar_nnz_view_t given_inverse_diagonal_, bool is_symmetric_)
      : handle(handle_),
        num_rows(num_rows_),
        num_cols(num_cols_),
        row_map(row_map_),
        entries(entries_),
        values(values_),
        given_inverse_diagonal(given_inverse_diagonal_),
        have_diagonal_given(true),
        is_symmetric(is_symmetric_) {}

  // Functor to swap the numbers of two colors,
  // so that the last cluster has the last color.
  // Except, doesn't touch the color of the last cluster,
  // since that value is needed the entire time this is running.
  struct ClusterColorRelabelFunctor {
    typedef typename HandleType::GraphColoringHandleType GCHandle;
    typedef typename GCHandle::color_view_t ColorView;
    typedef Kokkos::View<row_lno_t*, MyTempMemorySpace> RowmapView;
    typedef Kokkos::View<nnz_lno_t*, MyTempMemorySpace> EntriesView;
    ClusterColorRelabelFunctor(ColorView& colors_, color_t numClusterColors_, nnz_lno_t numClusters_)
        : colors(colors_), numClusterColors(numClusterColors_), numClusters(numClusters_) {}

    KOKKOS_INLINE_FUNCTION void operator()(const size_type i) const {
      if (colors(i) == numClusterColors)
        colors(i) = colors(numClusters - 1);
      else if (colors(i) == colors(numClusters - 1))
        colors(i) = numClusterColors;
    }

    ColorView colors;
    color_t numClusterColors;
    nnz_lno_t numClusters;
  };

  // Relabel the last cluster, after running ClusterColorRelabelFunctor.
  // Call with a one-element range policy.
  struct RelabelLastColorFunctor {
    typedef typename HandleType::GraphColoringHandleType GCHandle;
    typedef typename GCHandle::color_view_t ColorView;

    RelabelLastColorFunctor(ColorView& colors_, color_t numClusterColors_, nnz_lno_t numClusters_)
        : colors(colors_), numClusterColors(numClusterColors_), numClusters(numClusters_) {}

    KOKKOS_INLINE_FUNCTION void operator()(const size_type) const { colors(numClusters - 1) = numClusterColors; }

    ColorView colors;
    color_t numClusterColors;
    nnz_lno_t numClusters;
  };

  struct ClusterToVertexColoring {
    typedef typename HandleType::GraphColoringHandleType GCHandle;
    typedef typename GCHandle::color_view_t ColorView;

    ClusterToVertexColoring(ColorView& clusterColors_, ColorView& vertexColors_, nnz_lno_t numRows_,
                            nnz_lno_t numClusters_, nnz_lno_t clusterSize_)
        : clusterColors(clusterColors_),
          vertexColors(vertexColors_),
          numRows(numRows_),
          numClusters(numClusters_),
          clusterSize(clusterSize_) {}

    KOKKOS_INLINE_FUNCTION void operator()(const size_type i) const {
      size_type cluster       = i / clusterSize;
      size_type clusterOffset = i - cluster * clusterSize;
      vertexColors(i)         = ((clusterColors(cluster) - 1) * clusterSize) + clusterOffset + 1;
    }

    ColorView clusterColors;
    ColorView vertexColors;
    nnz_lno_t numRows;
    nnz_lno_t numClusters;
    nnz_lno_t clusterSize;
  };

  // Assign cluster labels to vertices, given that the vertices are naturally
  // ordered so that contiguous groups of vertices form decent clusters.
  template <typename View>
  struct NopVertClusteringFunctor {
    NopVertClusteringFunctor(View& vertClusters_, nnz_lno_t clusterSize_)
        : vertClusters(vertClusters_), numRows(vertClusters.extent(0)), clusterSize(clusterSize_) {}
    KOKKOS_INLINE_FUNCTION void operator()(const nnz_lno_t i) const { vertClusters(i) = i / clusterSize; }
    View vertClusters;
    nnz_lno_t numRows;
    nnz_lno_t clusterSize;
  };

  template <typename View>
  struct ReorderedClusteringFunctor {
    ReorderedClusteringFunctor(View& vertClusters_, View& ordering_, nnz_lno_t clusterSize_)
        : vertClusters(vertClusters_),
          ordering(ordering_),
          numRows(vertClusters.extent(0)),
          clusterSize(clusterSize_) {}
    KOKKOS_INLINE_FUNCTION void operator()(const nnz_lno_t i) const { vertClusters(i) = ordering(i) / clusterSize; }
    View vertClusters;
    View ordering;
    nnz_lno_t numRows;
    nnz_lno_t clusterSize;
  };

  void initialize_symbolic() {
    using nnz_view_t    = nnz_lno_persistent_work_view_t;
    using in_rowmap_t   = const_lno_row_view_t;
    using in_colinds_t  = const_lno_nnz_view_t;
    using rowmap_t      = Kokkos::View<size_type*, MyTempMemorySpace>;
    using colinds_t     = Kokkos::View<nnz_lno_t*, MyTempMemorySpace>;
    using raw_rowmap_t  = Kokkos::View<const size_type*, MyTempMemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using raw_colinds_t = Kokkos::View<const nnz_lno_t*, MyTempMemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    auto gsHandle       = get_gs_handle();
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    Kokkos::Timer timer;
#endif
    // sym_xadj/sym_adj is only used here for clustering.
    // Create them as non-const, unmanaged views to avoid
    // duplicating a bunch of code between the
    // symmetric and non-symmetric input cases.
    rowmap_t sym_xadj;
    colinds_t sym_adj;
    if (!this->is_symmetric) {
      KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap<in_rowmap_t, in_colinds_t, rowmap_t, colinds_t,
                                                             MyExecSpace>(num_rows, this->row_map, this->entries,
                                                                          sym_xadj, sym_adj);
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
      MyExecSpace().fence();
      std::cout << "SYMMETRIZING TIME: " << timer.seconds() << std::endl;
      timer.reset();
#endif
    }
    // Now that a symmetric graph is available, build the cluster graph (also
    // symmetric)
    nnz_lno_t clusterSize = gsHandle->get_cluster_size();
    nnz_lno_t numClusters = (num_rows + clusterSize - 1) / clusterSize;
    raw_rowmap_t raw_sym_xadj;
    raw_colinds_t raw_sym_adj;
    if (this->is_symmetric) {
      raw_sym_xadj = raw_rowmap_t(this->row_map.data(), this->row_map.extent(0));
      raw_sym_adj  = raw_colinds_t(this->entries.data(), this->entries.extent(0));
    } else {
      raw_sym_xadj = raw_rowmap_t(sym_xadj.data(), sym_xadj.extent(0));
      raw_sym_adj  = raw_colinds_t(sym_adj.data(), sym_adj.extent(0));
    }
    nnz_view_t vertClusters;
    auto clusterAlgo = gsHandle->get_clustering_algo();
    if (clusterAlgo == CLUSTER_DEFAULT) clusterAlgo = CLUSTER_MIS2;
    switch (clusterAlgo) {
      case CLUSTER_MIS2: {
        vertClusters = KokkosGraph::graph_mis2_aggregate<MyExecSpace, raw_rowmap_t, raw_colinds_t, nnz_view_t>(
            raw_sym_xadj, raw_sym_adj, numClusters);
        break;
      }
      case CLUSTER_BALLOON: {
        BalloonClustering<HandleType, raw_rowmap_t, raw_colinds_t> balloon(num_rows, raw_sym_xadj, raw_sym_adj);
        vertClusters = balloon.run(clusterSize);
        break;
      }
      case CLUSTER_DEFAULT: {
        throw std::logic_error("Logic to choose default clustering algorithm is incorrect");
      }
      default: throw std::runtime_error("Clustering algo " + std::to_string((int)clusterAlgo) + " is not implemented");
    }
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    MyExecSpace().fence();
    std::cout << "Graph clustering: " << timer.seconds() << '\n';
    timer.reset();
#endif
    rowmap_t clusterRowmap;
    colinds_t clusterEntries;
    nnz_view_t clusterOffsets;
    nnz_view_t clusterVerts;
    KokkosGraph::Experimental::graph_explicit_coarsen_with_inverse_map<Kokkos::Device<MyExecSpace, MyTempMemorySpace>,
                                                                       raw_rowmap_t, raw_colinds_t, nnz_view_t,
                                                                       rowmap_t, colinds_t, nnz_view_t>(
        raw_sym_xadj, raw_sym_adj, vertClusters, numClusters, clusterRowmap, clusterEntries, clusterOffsets,
        clusterVerts, false);
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    MyExecSpace().fence();
    std::cout << "Building explicit cluster graph: " << timer.seconds() << '\n';
    timer.reset();
#endif
#if KOKKOSSPARSE_IMPL_PRINTDEBUG
    {
      auto clusterRowmapHost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), clusterRowmap);
      auto clusterEntriesHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), clusterEntries);
      puts("Cluster graph (cluster #, and neighbors):");
      for (nnz_lno_t i = 0; i < numClusters; i++) {
        printf("%d: ", (int)i);
        for (nnz_lno_t j = clusterRowmapHost(i); j < clusterRowmapHost(i + 1); j++) {
          printf("%d ", (int)clusterEntriesHost(j));
        }
        putchar('\n');
      }
      printf("\n\n\n");
    }
#endif
    // Get the coloring of the cluster graph.
    typename HandleType::GraphColoringHandleType::color_view_t colors;
    color_t numColors;
#if KOKKOSSPARSE_IMPL_RUNSEQUENTIAL
    numColors = numClusters;
    std::cout << "SEQUENTIAL CGS: numColors = numClusters = " << numClusters << '\n';
    typename HandleType::GraphColoringHandleType::color_view_t::HostMirror h_colors =
        Kokkos::create_mirror_view(colors);
    for (int i = 0; i < numClusters; ++i) {
      h_colors(i) = i + 1;
    }
    Kokkos::deep_copy(colors, h_colors);
#else
    // Create a handle that uses nnz_lno_t as the size_type, since the cluster
    // graph should never be larger than 2^31 entries.
    HandleType kh;
    kh.create_graph_coloring_handle(gsHandle->get_coloring_algorithm());
    KokkosGraph::Experimental::graph_color_symbolic(&kh, numClusters, numClusters, clusterRowmap, clusterEntries);
    // retrieve colors
    auto coloringHandle = kh.get_graph_coloring_handle();
    colors              = coloringHandle->get_vertex_colors();
    numColors           = coloringHandle->get_num_colors();
    kh.destroy_graph_coloring_handle();
#endif
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    MyExecSpace().fence();
    std::cout << "Coloring: " << timer.seconds() << '\n';
    timer.reset();
#endif
    nnz_lno_persistent_work_view_t color_xadj;
    nnz_lno_persistent_work_view_t color_adj;
    KokkosKernels::Impl::create_reverse_map<typename HandleType::GraphColoringHandleType::color_view_t,
                                            nnz_lno_persistent_work_view_t, MyExecSpace>(numClusters, numColors, colors,
                                                                                         color_xadj, color_adj);
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    MyExecSpace().fence();
    std::cout << "CREATE_REVERSE_MAP:" << timer.seconds() << std::endl;
    timer.reset();
#endif
    nnz_lno_persistent_work_host_view_t color_xadj_host(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Color xadj"),
                                                        color_xadj.extent(0));
    Kokkos::deep_copy(color_xadj_host, color_xadj);
    gsHandle->set_color_xadj(color_xadj_host);
    gsHandle->set_color_adj(color_adj);
    gsHandle->set_num_colors(numColors);
    gsHandle->set_cluster_xadj(clusterOffsets);
    gsHandle->set_cluster_adj(clusterVerts);
    gsHandle->set_call_symbolic(true);
  }

  struct Get_Matrix_Diagonals {
    const_lno_row_view_t _xadj;
    const_lno_nnz_view_t _adj;          // CSR storage of the graph.
    const_scalar_nnz_view_t _adj_vals;  // CSR storage of the graph.
    scalar_persistent_work_view_t _diagonals;

    nnz_lno_t num_total_rows;
    nnz_lno_t rows_per_team;

    nnz_scalar_t one;

    Get_Matrix_Diagonals(const_lno_row_view_t xadj_, const_lno_nnz_view_t adj_, const_scalar_nnz_view_t adj_vals_,
                         scalar_persistent_work_view_t diagonals_, nnz_lno_t num_total_rows_, nnz_lno_t rows_per_team_)
        : _xadj(xadj_),
          _adj(adj_),
          _adj_vals(adj_vals_),
          _diagonals(diagonals_),
          num_total_rows(num_total_rows_),
          rows_per_team(rows_per_team_),
          one(Kokkos::ArithTraits<nnz_scalar_t>::one()) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t row_id) const {
      size_type row_begin = _xadj(row_id);
      size_type row_end   = _xadj(row_id + 1);
      for (size_type j = row_begin; j < row_end; j++) {
        nnz_lno_t column_id = _adj(j);
        if (column_id == row_id) {
          nnz_scalar_t val   = _adj_vals(j);
          _diagonals(row_id) = one / val;
          break;
        }
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t& team) const {
      const nnz_lno_t i_begin = team.league_rank() * rows_per_team;
      const nnz_lno_t i_end   = i_begin + rows_per_team <= num_total_rows ? i_begin + rows_per_team : num_total_rows;
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, i_begin, i_end), [&](const nnz_lno_t row_id) {
        size_type row_begin = _xadj(row_id);
        size_type row_end   = _xadj(row_id + 1);
        nnz_lno_t row_size  = row_end - row_begin;

        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, row_size), [&](const nnz_lno_t col_ind) {
          size_type val_index = col_ind + row_begin;
          nnz_lno_t column_id = _adj(val_index);
          if (column_id == row_id) {
            _diagonals(row_id) = one / _adj_vals(val_index);
            return;
          }
        });
      });
    }
  };

  void initialize_numeric() {
    auto gsHandle = get_gs_handle();
    if (!gsHandle->is_symbolic_called()) {
      this->initialize_symbolic();
    }
    // Timer for whole numeric
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    Kokkos::Timer timer;
#endif
    size_type nnz = this->entries.extent(0);

    int suggested_vector_size = this->handle->get_suggested_vector_size(num_rows, nnz);
    int suggested_team_size   = this->handle->get_suggested_team_size(suggested_vector_size);

    scalar_persistent_work_view_t inverse_diagonal(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Aii^-1"), num_rows);
    nnz_lno_t rows_per_team =
        this->handle->get_team_work_size(suggested_team_size, MyExecSpace().concurrency(), num_rows);

    if (have_diagonal_given) {
      Kokkos::deep_copy(inverse_diagonal, this->given_inverse_diagonal);
    } else {
      // extract inverse diagonal from matrix
      Get_Matrix_Diagonals gmd(this->row_map, this->entries, this->values, inverse_diagonal, num_rows, rows_per_team);
      if (gsHandle->use_teams()) {
        Kokkos::parallel_for(
            "KokkosSparse::GaussSeidel::team_get_matrix_diagonals",
            team_policy_t((num_rows + rows_per_team - 1) / rows_per_team, suggested_team_size, suggested_vector_size),
            gmd);
      } else {
        Kokkos::parallel_for("KokkosSparse::GaussSeidel::get_matrix_diagonals", my_exec_space(0, num_rows), gmd);
      }
    }
    gsHandle->set_inverse_diagonal(inverse_diagonal);
    gsHandle->set_call_numeric(true);
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    MyExecSpace().fence();
    std::cout << "NUMERIC:" << timer.seconds() << std::endl;
#endif
  }

  template <typename x_value_array_type, typename y_value_array_type>
  void apply(x_value_array_type x_lhs_output_vec, y_value_array_type y_rhs_input_vec, bool init_zero_x_vector = false,
             int numIter = 1, nnz_scalar_t omega = Kokkos::ArithTraits<nnz_scalar_t>::one(), bool apply_forward = true,
             bool apply_backward = true, bool /*update_y_vector*/ = true) {
    auto gsHandle = get_gs_handle();

    size_type nnz                                    = entries.extent(0);
    nnz_lno_persistent_work_view_t color_adj         = gsHandle->get_color_adj();
    nnz_lno_persistent_work_host_view_t h_color_xadj = gsHandle->get_color_xadj();

    color_t numColors = gsHandle->get_num_colors();

    if (init_zero_x_vector) {
      KokkosKernels::Impl::zero_vector<x_value_array_type, MyExecSpace>(num_cols, x_lhs_output_vec);
    }

    scalar_persistent_work_view_t inverse_diagonal = gsHandle->get_inverse_diagonal();

    if (gsHandle->use_teams()) {
      int suggested_vector_size = this->handle->get_suggested_vector_size(num_rows, nnz);
      int suggested_team_size   = this->handle->get_suggested_team_size(suggested_vector_size);

      nnz_lno_t rows_per_team =
          this->handle->get_team_work_size(suggested_team_size, MyExecSpace().concurrency(), num_rows);
      // Get clusters per team. Round down to favor finer granularity, since
      // this is sensitive to load imbalance
      nnz_lno_t clusters_per_team = rows_per_team / gsHandle->get_cluster_size();
      if (clusters_per_team == 0) clusters_per_team = 1;

      Team_PSGS<x_value_array_type, y_value_array_type> gs(
          this->row_map, this->entries, this->values, x_lhs_output_vec, y_rhs_input_vec, 0,
          0,  // color set range is set right before launch for each color set
          color_adj, gsHandle->get_cluster_xadj(), gsHandle->get_cluster_adj(), inverse_diagonal, clusters_per_team,
          omega);

      this->IterativeTeamPSGS(gs, numColors, h_color_xadj, suggested_team_size, suggested_vector_size, numIter,
                              apply_forward, apply_backward);
    } else {
      PSGS<x_value_array_type, y_value_array_type> gs(this->row_map, this->entries, this->values, x_lhs_output_vec,
                                                      y_rhs_input_vec, color_adj, gsHandle->get_cluster_xadj(),
                                                      gsHandle->get_cluster_adj(), omega, inverse_diagonal);

      this->IterativePSGS(gs, numColors, h_color_xadj, numIter, apply_forward, apply_backward);
    }
  }

  template <typename TPSGS>
  void IterativeTeamPSGS(TPSGS& gs, color_t numColors, nnz_lno_persistent_work_host_view_t h_color_xadj,
                         nnz_lno_t team_size, nnz_lno_t vec_size, int num_iteration, bool apply_forward,
                         bool apply_backward) {
    for (int i = 0; i < num_iteration; ++i)
      this->DoTeamPSGS(gs, numColors, h_color_xadj, team_size, vec_size, apply_forward, apply_backward);
  }

  template <typename TPSGS>
  void DoTeamPSGS(TPSGS& gs, color_t numColors, nnz_lno_persistent_work_host_view_t h_color_xadj, nnz_lno_t team_size,
                  nnz_lno_t vec_size, bool apply_forward, bool apply_backward) {
    if (apply_forward) {
      gs._is_backward = false;
      for (color_t i = 0; i < numColors; ++i) {
        nnz_lno_t color_index_begin = h_color_xadj(i);
        nnz_lno_t color_index_end   = h_color_xadj(i + 1);
        int overall_work            = color_index_end - color_index_begin;  // /256 + 1;
        gs._color_set_begin         = color_index_begin;
        gs._color_set_end           = color_index_end;
        Kokkos::parallel_for(
            "KokkosSparse::GaussSeidel::Team_PSGS::forward",
            team_policy_t((overall_work + gs._clusters_per_team - 1) / gs._clusters_per_team, team_size, vec_size), gs);
      }
    }
    if (apply_backward) {
      gs._is_backward = true;
      if (numColors > 0)
        for (color_t i = numColors - 1;; --i) {
          nnz_lno_t color_index_begin = h_color_xadj(i);
          nnz_lno_t color_index_end   = h_color_xadj(i + 1);
          nnz_lno_t overall_work      = color_index_end - color_index_begin;  // /256 + 1;
          gs._color_set_begin         = color_index_begin;
          gs._color_set_end           = color_index_end;
          Kokkos::parallel_for(
              "KokkosSparse::GaussSeidel::Team_PSGS::forward",
              team_policy_t((overall_work + gs._clusters_per_team - 1) / gs._clusters_per_team, team_size, vec_size),
              gs);
          if (i == 0) {
            break;
          }
        }
    }
  }

  template <typename PSGS>
  void IterativePSGS(PSGS& gs, color_t numColors, nnz_lno_persistent_work_host_view_t h_color_xadj, int num_iteration,
                     bool apply_forward, bool apply_backward) {
    for (int i = 0; i < num_iteration; ++i) {
      this->DoPSGS(gs, numColors, h_color_xadj, apply_forward, apply_backward);
    }
  }

  template <typename PSGS>
  void DoPSGS(PSGS& gs, color_t numColors, nnz_lno_persistent_work_host_view_t h_color_xadj, bool apply_forward,
              bool apply_backward) {
    if (apply_forward) {
      for (color_t i = 0; i < numColors; ++i) {
        nnz_lno_t color_index_begin = h_color_xadj(i);
        nnz_lno_t color_index_end   = h_color_xadj(i + 1);
        gs._color_set_begin         = color_index_begin;
        gs._color_set_end           = color_index_end;
        Kokkos::parallel_for("KokkosSparse::GaussSeidel::PSGS::forward",
                             Kokkos::RangePolicy<MyExecSpace, PSGS_ForwardTag>(0, color_index_end - color_index_begin),
                             gs);
      }
    }
    if (apply_backward && numColors) {
      for (size_type i = numColors - 1;; --i) {
        nnz_lno_t color_index_begin = h_color_xadj(i);
        nnz_lno_t color_index_end   = h_color_xadj(i + 1);
        gs._color_set_begin         = color_index_begin;
        gs._color_set_end           = color_index_end;
        Kokkos::parallel_for("KokkosSparse::GaussSeidel::PSGS::backward",
                             Kokkos::RangePolicy<MyExecSpace, PSGS_BackwardTag>(0, color_index_end - color_index_begin),
                             gs);
        if (i == 0) {
          break;
        }
      }
    }
  }
};  // class ClusterGaussSeidel
}  // namespace Impl
}  // namespace KokkosSparse

#endif
