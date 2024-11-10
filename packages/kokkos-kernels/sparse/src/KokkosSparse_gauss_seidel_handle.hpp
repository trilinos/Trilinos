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

#include <Kokkos_Core.hpp>
#include <KokkosKernels_Utils.hpp>
// needed for two-stage/classical GS
#include <KokkosSparse_CrsMatrix.hpp>
// needed for the set of available coloring algorithms
#include <KokkosGraph_Distance1ColorHandle.hpp>

#ifndef _GAUSSSEIDELHANDLE_HPP
#define _GAUSSSEIDELHANDLE_HPP
// #define VERBOSE

namespace KokkosSparse {

enum GSAlgorithm { GS_DEFAULT, GS_PERMUTED, GS_TEAM, GS_CLUSTER, GS_TWOSTAGE };
enum GSDirection { GS_FORWARD, GS_BACKWARD, GS_SYMMETRIC };
enum struct ClusteringAlgorithm { CLUSTER_DEFAULT, CLUSTER_MIS2, CLUSTER_BALLOON, NUM_CLUSTERING_ALGORITHMS };

static constexpr ClusteringAlgorithm CLUSTER_DEFAULT           = ClusteringAlgorithm::CLUSTER_DEFAULT;
static constexpr ClusteringAlgorithm CLUSTER_MIS2              = ClusteringAlgorithm::CLUSTER_MIS2;
static constexpr ClusteringAlgorithm CLUSTER_BALLOON           = ClusteringAlgorithm::CLUSTER_BALLOON;
static constexpr ClusteringAlgorithm NUM_CLUSTERING_ALGORITHMS = ClusteringAlgorithm::NUM_CLUSTERING_ALGORITHMS;

inline const char *getClusterAlgoName(ClusteringAlgorithm ca) {
  switch (ca) {
    case CLUSTER_BALLOON: return "Balloon";
    case CLUSTER_MIS2: return "MIS(2)";
    default:;
  }
  return "INVALID CLUSTERING ALGORITHM";
}

template <class size_type_, class lno_t_, class scalar_t_, class ExecutionSpace, class TemporaryMemorySpace,
          class PersistentMemorySpace>
class GaussSeidelHandle {
 public:
  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;

  typedef typename std::remove_const<size_type_>::type size_type;
  typedef const size_type const_size_type;

  typedef typename std::remove_const<lno_t_>::type nnz_lno_t;
  typedef const nnz_lno_t const_nnz_lno_t;

  typedef typename std::remove_const<scalar_t_>::type nnz_scalar_t;
  typedef const nnz_scalar_t const_nnz_scalar_t;

  typedef typename Kokkos::View<size_type *, HandleTempMemorySpace> row_lno_temp_work_view_t;
  typedef typename Kokkos::View<size_type *, HandlePersistentMemorySpace> row_lno_persistent_work_view_t;
  typedef typename row_lno_persistent_work_view_t::HostMirror row_lno_persistent_work_host_view_t;  // Host view type

  typedef typename Kokkos::View<nnz_scalar_t *, HandleTempMemorySpace> scalar_temp_work_view_t;
  typedef typename Kokkos::View<nnz_scalar_t *, HandlePersistentMemorySpace> scalar_persistent_work_view_t;
  typedef typename scalar_persistent_work_view_t::HostMirror scalar_persistent_work_host_view_t;  // Host view type

  typedef typename Kokkos::View<nnz_lno_t *, HandleTempMemorySpace> nnz_lno_temp_work_view_t;
  typedef typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace> nnz_lno_persistent_work_view_t;
  typedef typename nnz_lno_persistent_work_view_t::HostMirror nnz_lno_persistent_work_host_view_t;  // Host view type

 protected:
  HandleExecSpace execution_space;
  int num_streams;

  GSAlgorithm algorithm_type;

  nnz_lno_persistent_work_host_view_t color_xadj;
  nnz_lno_persistent_work_view_t color_adj;
  nnz_lno_t numColors;

  bool called_symbolic;
  bool called_numeric;

  int suggested_vector_size;
  int suggested_team_size;

 public:
  /**
   * \brief Default constructor.
   */
  GaussSeidelHandle(GSAlgorithm gs)
      : execution_space(HandleExecSpace()),
        num_streams(1),
        algorithm_type(gs),
        color_xadj(),
        color_adj(),
        numColors(0),
        called_symbolic(false),
        called_numeric(false),
        suggested_vector_size(0),
        suggested_team_size(0) {}

  GaussSeidelHandle(HandleExecSpace handle_exec_space, int n_streams, GSAlgorithm gs)
      : execution_space(handle_exec_space),
        num_streams(n_streams),
        algorithm_type(gs),
        color_xadj(),
        color_adj(),
        numColors(0),
        called_symbolic(false),
        called_numeric(false),
        suggested_vector_size(0),
        suggested_team_size(0) {}

  virtual ~GaussSeidelHandle() = default;

  // getters
  int get_num_streams() const { return num_streams; }

  HandleExecSpace get_execution_space() const { return this->execution_space; }

  GSAlgorithm get_algorithm_type() const { return this->algorithm_type; }

  nnz_lno_persistent_work_host_view_t get_color_xadj() const { return this->color_xadj; }
  nnz_lno_persistent_work_view_t get_color_adj() const { return this->color_adj; }
  nnz_lno_t get_num_colors() const { return this->numColors; }

  bool is_symbolic_called() const { return this->called_symbolic; }
  bool is_numeric_called() const { return this->called_numeric; }

  template <class ExecSpaceIn>
  void set_execution_space(const ExecSpaceIn exec_space_in) {
    static bool is_set = false;
    if (!is_set) {
      static_assert(std::is_same<ExecSpaceIn, HandleExecSpace>::value,
                    "The type of exec_space_in should be the same as "
                    "GaussSeidelHandle::HandleExecSpace");
      this->execution_space = exec_space_in;
    } else {
      if (exec_space_in != this->execution_space)
        throw std::runtime_error(
            "Gauss Seidel cannot be called on different execution spaces "
            "without multiple handles. Please create a new handle via "
            "create_gs_handle.\n");
    }
    is_set = true;
  }

  void set_algorithm_type(const GSAlgorithm sgs_algo) {
    this->algorithm_type  = sgs_algo;
    this->called_symbolic = false;
  }

  void set_call_symbolic(bool call = true) { this->called_symbolic = call; }
  void set_call_numeric(bool call = true) { this->called_numeric = call; }

  void set_color_xadj(const nnz_lno_persistent_work_host_view_t &color_xadj_) { this->color_xadj = color_xadj_; }
  void set_color_adj(const nnz_lno_persistent_work_view_t &color_adj_) { this->color_adj = color_adj_; }
  void set_num_colors(const nnz_lno_t &numColors_) { this->numColors = numColors_; }

  void vector_team_size(int max_allowed_team_size,
                        int &suggested_vector_size_,  // output
                        int &suggested_team_size_,    // output
                        size_type nr, size_type nnz) {
    if (this->suggested_team_size && this->suggested_vector_size) {
      suggested_vector_size_ = this->suggested_vector_size;
      suggested_team_size_   = this->suggested_team_size;
      return;
    } else {
      KokkosKernels::Impl::get_suggested_vector_size<size_type, ExecutionSpace>(suggested_vector_size_, nr, nnz);
      KokkosKernels::Impl::get_suggested_team_size<ExecutionSpace>(max_allowed_team_size, suggested_vector_size_,
                                                                   suggested_team_size_);
      this->suggested_team_size   = suggested_vector_size_;
      this->suggested_vector_size = suggested_vector_size_;
    }
  }
};

template <class size_type_, class lno_t_, class scalar_t_, class ExecutionSpace, class TemporaryMemorySpace,
          class PersistentMemorySpace>
class PointGaussSeidelHandle : public GaussSeidelHandle<size_type_, lno_t_, scalar_t_, ExecutionSpace,
                                                        TemporaryMemorySpace, PersistentMemorySpace> {
 public:
  typedef GaussSeidelHandle<size_type_, lno_t_, scalar_t_, ExecutionSpace, TemporaryMemorySpace, PersistentMemorySpace>
      GSHandle;
  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;

  typedef typename std::remove_const<size_type_>::type size_type;
  typedef const size_type const_size_type;

  typedef typename std::remove_const<lno_t_>::type nnz_lno_t;
  typedef const nnz_lno_t const_nnz_lno_t;

  typedef typename std::remove_const<scalar_t_>::type nnz_scalar_t;
  typedef const nnz_scalar_t const_nnz_scalar_t;

  typedef typename Kokkos::View<size_type *, HandleTempMemorySpace> row_lno_temp_work_view_t;
  typedef typename Kokkos::View<size_type *, HandlePersistentMemorySpace> row_lno_persistent_work_view_t;
  typedef typename row_lno_persistent_work_view_t::HostMirror row_lno_persistent_work_host_view_t;  // Host view type

  typedef typename Kokkos::View<nnz_scalar_t *, HandleTempMemorySpace> scalar_temp_work_view_t;
  typedef typename Kokkos::View<nnz_scalar_t *, HandlePersistentMemorySpace> scalar_persistent_work_view_t;
  typedef typename Kokkos::View<nnz_scalar_t **, default_layout, HandlePersistentMemorySpace>
      scalar_persistent_work_view2d_t;
  typedef typename scalar_persistent_work_view_t::HostMirror scalar_persistent_work_host_view_t;  // Host view type

  typedef typename Kokkos::View<nnz_lno_t *, HandleTempMemorySpace> nnz_lno_temp_work_view_t;
  typedef typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace> nnz_lno_persistent_work_view_t;
  typedef typename nnz_lno_persistent_work_view_t::HostMirror nnz_lno_persistent_work_host_view_t;  // Host view type

 private:
  row_lno_persistent_work_view_t permuted_xadj;
  nnz_lno_persistent_work_view_t permuted_adj;
  scalar_persistent_work_view_t permuted_adj_vals;
  nnz_lno_persistent_work_view_t old_to_new_map;

  scalar_persistent_work_view2d_t permuted_y_vector;
  scalar_persistent_work_view2d_t permuted_x_vector;

  scalar_persistent_work_view_t permuted_inverse_diagonal;
  nnz_lno_t block_size;  // this is for block sgs

  nnz_lno_t num_values_in_l1, num_values_in_l2, num_big_rows;
  size_t level_1_mem, level_2_mem;

  // Option set by user: rows with at least this many nonzeros are handled by a
  // separate kernel
  nnz_lno_t long_row_threshold;
  // Number of long rows per color set. They are all grouped at the end of each
  // color set.
  nnz_lno_persistent_work_host_view_t long_rows_per_color;
  // Maximum row length in each color set.
  nnz_lno_persistent_work_host_view_t max_row_length_per_color;
  // Temporary space for matvec over long rows - size is only max num long rows
  // in a color.
  scalar_persistent_work_view_t long_row_x;

  // Coloring algorithm to use
  KokkosGraph::ColoringAlgorithm coloring_algo;

 public:
  /**
   * \brief Default constructor.
   */
  PointGaussSeidelHandle(GSHandle gs_handle,
                         KokkosGraph::ColoringAlgorithm coloring_algo_ = KokkosGraph::COLORING_DEFAULT)
      : GSHandle(gs_handle),
        permuted_xadj(),
        permuted_adj(),
        permuted_adj_vals(),
        old_to_new_map(),
        permuted_y_vector(),
        permuted_x_vector(),
        permuted_inverse_diagonal(),
        block_size(1),
        num_values_in_l1(-1),
        num_values_in_l2(-1),
        num_big_rows(0),
        level_1_mem(0),
        level_2_mem(0),
        long_row_threshold(0),
        coloring_algo(coloring_algo_) {
    if (gs_handle.get_algorithm_type() == GS_DEFAULT) this->choose_default_algorithm();
  }

  PointGaussSeidelHandle(GSAlgorithm gs                                = GS_DEFAULT,
                         KokkosGraph::ColoringAlgorithm coloring_algo_ = KokkosGraph::COLORING_DEFAULT)
      : PointGaussSeidelHandle(GSHandle(gs), coloring_algo_) {}

  PointGaussSeidelHandle(HandleExecSpace handle_exec_space, int n_streams, GSAlgorithm gs = GS_DEFAULT,
                         KokkosGraph::ColoringAlgorithm coloring_algo_ = KokkosGraph::COLORING_DEFAULT)
      : PointGaussSeidelHandle(GSHandle(handle_exec_space, n_streams, gs), coloring_algo_) {}

  void set_block_size(nnz_lno_t bs) { this->block_size = bs; }
  nnz_lno_t get_block_size() const { return this->block_size; }

  void choose_default_algorithm() {
    if (KokkosKernels::Impl::kk_is_gpu_exec_space<ExecutionSpace>())
      this->algorithm_type = GS_TEAM;
    else
      this->algorithm_type = GS_PERMUTED;
  }

  KokkosGraph::ColoringAlgorithm get_coloring_algorithm() const { return this->coloring_algo; }
  void set_coloring_algorithm(KokkosGraph::ColoringAlgorithm algo) { this->coloring_algo = algo; }

  ~PointGaussSeidelHandle() = default;

  // getters
  row_lno_persistent_work_view_t get_new_xadj() const { return this->permuted_xadj; }
  nnz_lno_persistent_work_view_t get_new_adj() const { return this->permuted_adj; }
  scalar_persistent_work_view_t get_new_adj_val() const { return this->permuted_adj_vals; }
  nnz_lno_persistent_work_view_t get_old_to_new_map() const { return this->old_to_new_map; }

  // setters
  void set_algorithm_type(const GSAlgorithm &sgs_algo) { this->algorithm_type = sgs_algo; }

  void set_call_symbolic(bool call = true) { this->called_symbolic = call; }
  void set_call_numeric(bool call = true) { this->called_numeric = call; }

  void set_num_colors(const nnz_lno_t &numColors_) { this->numColors = numColors_; }

  void set_new_xadj(const row_lno_persistent_work_view_t &xadj_) { this->permuted_xadj = xadj_; }
  void set_new_adj(const nnz_lno_persistent_work_view_t &adj_) { this->permuted_adj = adj_; }
  void set_new_adj_val(const scalar_persistent_work_view_t &adj_vals_) { this->permuted_adj_vals = adj_vals_; }
  void set_old_to_new_map(const nnz_lno_persistent_work_view_t &old_to_new_map_) {
    this->old_to_new_map = old_to_new_map_;
  }
  void set_permuted_inverse_diagonal(const scalar_persistent_work_view_t permuted_inverse_diagonal_) {
    this->permuted_inverse_diagonal = permuted_inverse_diagonal_;
  }

  scalar_persistent_work_view_t get_permuted_inverse_diagonal() const { return this->permuted_inverse_diagonal; }

  void set_level_1_mem(size_t _level_1_mem) { this->level_1_mem = _level_1_mem; }
  void set_level_2_mem(size_t _level_2_mem) { this->level_2_mem = _level_2_mem; }

  void set_num_values_in_l1(nnz_lno_t _num_values_in_l1) { this->num_values_in_l1 = _num_values_in_l1; }
  void set_num_values_in_l2(nnz_lno_t _num_values_in_l2) { this->num_values_in_l2 = _num_values_in_l2; }

  void set_num_big_rows(nnz_lno_t _big_rows) { this->num_big_rows = _big_rows; }

  size_t get_level_1_mem() const { return this->level_1_mem; }
  size_t get_level_2_mem() const { return this->level_2_mem; }

  nnz_lno_t get_num_values_in_l1() const { return this->num_values_in_l1; }
  nnz_lno_t get_num_values_in_l2() const { return this->num_values_in_l2; }
  nnz_lno_t get_num_big_rows() const { return this->num_big_rows; }

  nnz_lno_t get_long_row_threshold() const { return long_row_threshold; }

  void set_long_row_threshold(nnz_lno_t lrt) { long_row_threshold = lrt; }

  nnz_lno_persistent_work_host_view_t get_long_rows_per_color() const { return long_rows_per_color; }

  void set_long_rows_per_color(const nnz_lno_persistent_work_host_view_t &long_rows_per_color_) {
    long_rows_per_color = long_rows_per_color_;
  }

  nnz_lno_persistent_work_host_view_t get_max_row_length_per_color() const { return max_row_length_per_color; }

  void set_max_row_length_per_color(const nnz_lno_persistent_work_host_view_t &max_row_length_per_color_) {
    max_row_length_per_color = max_row_length_per_color_;
  }

  scalar_persistent_work_view_t get_long_row_x() const { return long_row_x; }

  void set_long_row_x(const scalar_persistent_work_view_t &long_row_x_) { long_row_x = long_row_x_; }

  void allocate_x_y_vectors(nnz_lno_t num_rows, nnz_lno_t num_cols, nnz_lno_t num_vecs) {
    if (permuted_y_vector.extent(0) != size_t(num_rows) || permuted_y_vector.extent(1) != size_t(num_vecs)) {
      permuted_y_vector = scalar_persistent_work_view2d_t("PERMUTED Y VECTOR", num_rows, num_vecs);
    }
    if (permuted_x_vector.extent(0) != size_t(num_cols) || permuted_x_vector.extent(1) != size_t(num_vecs)) {
      permuted_x_vector = scalar_persistent_work_view2d_t("PERMUTED X VECTOR", num_cols, num_vecs);
    }
  }

  scalar_persistent_work_view2d_t get_permuted_y_vector() const { return this->permuted_y_vector; }
  scalar_persistent_work_view2d_t get_permuted_x_vector() const { return this->permuted_x_vector; }

  void vector_team_size(int max_allowed_team_size, int &suggested_vector_size_, int &suggested_team_size_, size_type nr,
                        size_type nnz) {
    // suggested_team_size_ =  this->suggested_team_size = 1;
    // suggested_vector_size_=this->suggested_vector_size = 1;
    // return;
    if (this->suggested_team_size && this->suggested_vector_size) {
      suggested_vector_size_ = this->suggested_vector_size;
      suggested_team_size_   = this->suggested_team_size;
      return;
    } else {
      KokkosKernels::Impl::get_suggested_vector_size<size_type, ExecutionSpace>(suggested_vector_size_, nr, nnz);
      KokkosKernels::Impl::get_suggested_team_size<ExecutionSpace>(max_allowed_team_size, suggested_vector_size_,
                                                                   suggested_team_size_);
      this->suggested_team_size   = suggested_vector_size_;
      this->suggested_vector_size = suggested_vector_size_;
    }
  }
};

template <class size_type_, class lno_t_, class scalar_t_, class ExecutionSpace, class TemporaryMemorySpace,
          class PersistentMemorySpace>
class ClusterGaussSeidelHandle : public GaussSeidelHandle<size_type_, lno_t_, scalar_t_, ExecutionSpace,
                                                          TemporaryMemorySpace, PersistentMemorySpace> {
 public:
  typedef GaussSeidelHandle<size_type_, lno_t_, scalar_t_, ExecutionSpace, TemporaryMemorySpace, PersistentMemorySpace>
      GSHandle;
  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;

  typedef typename std::remove_const<size_type_>::type size_type;
  typedef const size_type const_size_type;

  typedef typename std::remove_const<lno_t_>::type nnz_lno_t;
  typedef const nnz_lno_t const_nnz_lno_t;

  typedef typename std::remove_const<scalar_t_>::type nnz_scalar_t;
  typedef const nnz_scalar_t const_nnz_scalar_t;

  typedef typename Kokkos::View<size_type *, HandleTempMemorySpace> row_lno_temp_work_view_t;
  typedef typename Kokkos::View<size_type *, HandlePersistentMemorySpace> row_lno_persistent_work_view_t;
  typedef typename row_lno_persistent_work_view_t::HostMirror row_lno_persistent_work_host_view_t;  // Host view type

  typedef typename Kokkos::View<nnz_scalar_t *, HandleTempMemorySpace> scalar_temp_work_view_t;
  typedef typename Kokkos::View<nnz_scalar_t *, HandlePersistentMemorySpace> scalar_persistent_work_view_t;
  typedef typename scalar_persistent_work_view_t::HostMirror scalar_persistent_work_host_view_t;  // Host view type

  typedef typename Kokkos::View<nnz_lno_t *, HandleTempMemorySpace> nnz_lno_temp_work_view_t;
  typedef typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace> nnz_lno_persistent_work_view_t;
  typedef typename nnz_lno_persistent_work_view_t::HostMirror nnz_lno_persistent_work_host_view_t;  // Host view type

 private:
  ClusteringAlgorithm cluster_algo;

  // This is the user-configurable target cluster size.
  // Some clusters may be slightly larger or smaller,
  // but cluster_xadj always has the exact size of each.
  nnz_lno_t cluster_size;

  // Coloring algorithm to use on the coarsened graph
  KokkosGraph::ColoringAlgorithm coloring_algo;

  int suggested_vector_size;
  int suggested_team_size;

  scalar_persistent_work_view_t inverse_diagonal;

  // cluster_xadj and cluster_adj encode the vertices in each cluster
  nnz_lno_persistent_work_view_t cluster_xadj;
  nnz_lno_persistent_work_view_t cluster_adj;
  // vert_clusters(i) is the cluster that vertex i belongs to
  nnz_lno_persistent_work_view_t vert_clusters;

 public:
  /**
   * \brief Default constructor.
   */

  // Constructor for cluster-coloring based GS and SGS
  ClusterGaussSeidelHandle(ClusteringAlgorithm cluster_algo_, nnz_lno_t cluster_size_,
                           KokkosGraph::ColoringAlgorithm coloring_algo_)
      : GSHandle(GS_CLUSTER),
        cluster_algo(cluster_algo_),
        cluster_size(cluster_size_),
        coloring_algo(coloring_algo_),
        inverse_diagonal(),
        cluster_xadj(),
        cluster_adj(),
        vert_clusters() {}

  void set_cluster_size(nnz_lno_t cs) { this->cluster_size = cs; }
  nnz_lno_t get_cluster_size() const { return this->cluster_size; }

  KokkosGraph::ColoringAlgorithm get_coloring_algorithm() const { return this->coloring_algo; }
  void set_coloring_algorithm(KokkosGraph::ColoringAlgorithm algo) { this->coloring_algo = algo; }

  void set_vert_clusters(nnz_lno_persistent_work_view_t &vert_clusters_) { this->vert_clusters = vert_clusters_; }
  void set_cluster_xadj(nnz_lno_persistent_work_view_t &cluster_xadj_) { this->cluster_xadj = cluster_xadj_; }
  void set_cluster_adj(nnz_lno_persistent_work_view_t &cluster_adj_) { this->cluster_adj = cluster_adj_; }

  nnz_lno_persistent_work_view_t get_vert_clusters() const {
    if (!this->is_symbolic_called())
      throw std::runtime_error("vert_clusters does not exist until after symbolic setup.");
    return vert_clusters;
  }

  nnz_lno_persistent_work_view_t get_cluster_xadj() const {
    if (!this->is_symbolic_called())
      throw std::runtime_error("cluster_xadj does not exist until after symbolic setup.");
    return cluster_xadj;
  }

  nnz_lno_persistent_work_view_t get_cluster_adj() const {
    if (!this->is_symbolic_called()) throw std::runtime_error("cluster_adj does not exist until after symbolic setup.");
    return cluster_adj;
  }

  void set_inverse_diagonal(scalar_persistent_work_view_t &inv_diag) { this->inverse_diagonal = inv_diag; }

  scalar_persistent_work_view_t get_inverse_diagonal() const {
    if (!this->is_symbolic_called())
      throw std::runtime_error("inverse diagonal does not exist until after numeric setup.");
    return inverse_diagonal;
  }

  bool use_teams() const { return KokkosKernels::Impl::kk_is_gpu_exec_space<ExecutionSpace>(); }

  ~ClusterGaussSeidelHandle() = default;

  ClusteringAlgorithm get_clustering_algo() const { return this->cluster_algo; }
};

// -------------------------------------
// Handle for Two-stage/Classical GS
template <typename input_size_t, typename input_ordinal_t, typename input_scalar_t, class ExecutionSpace,
          class TemporaryMemorySpace, class PersistentMemorySpace>
class TwoStageGaussSeidelHandle
    : public GaussSeidelHandle<input_size_t, input_ordinal_t, input_scalar_t, ExecutionSpace, TemporaryMemorySpace,
                               PersistentMemorySpace> {
 public:
  using memory_space = typename ExecutionSpace::memory_space;
  using scalar_t     = typename std::remove_const<input_scalar_t>::type;
  using ordinal_t    = typename std::remove_const<input_ordinal_t>::type;
  using size_type    = typename std::remove_const<input_size_t>::type;

  using device_t = Kokkos::Device<ExecutionSpace, TemporaryMemorySpace>;
  using crsmat_t = KokkosSparse::CrsMatrix<scalar_t, ordinal_t, device_t, void, size_type>;
  using graph_t  = typename crsmat_t::StaticCrsGraphType;

  using input_row_map_view_t = typename graph_t::row_map_type;
  using input_entries_view_t = typename graph_t::entries_type;
  using input_values_view_t  = typename crsmat_t::values_type;

  using const_row_map_view_t = typename input_row_map_view_t::const_type;
  using row_map_view_t       = typename input_row_map_view_t::non_const_type;

  using const_entries_view_t = typename input_entries_view_t::const_type;
  using entries_view_t       = typename input_entries_view_t::non_const_type;

  using const_values_view_t = typename input_values_view_t::const_type;
  using values_view_t       = typename input_values_view_t::non_const_type;

  using const_ordinal_t = typename const_entries_view_t::value_type;
  using const_scalar_t  = typename const_values_view_t::value_type;

  using vector_view_t = Kokkos::View<scalar_t **, default_layout, device_t>;

  using GSHandle = GaussSeidelHandle<input_size_t, input_ordinal_t, input_scalar_t, ExecutionSpace,
                                     TemporaryMemorySpace, PersistentMemorySpace>;

  using HandleExecSpace = typename GSHandle::HandleExecSpace;

  /**
   * @brief Construct a new Two Stage Gauss Seidel Handle object
   *
   * @param gs_handle The GaussSeidel handle.
   */
  TwoStageGaussSeidelHandle(GSHandle gs_handle)
      : GSHandle(gs_handle),
        nrows(0),
        nrhs(1),
        direction(GS_SYMMETRIC),
        two_stage(true),
        compact_form(false),
        num_inner_sweeps(1),
        num_outer_sweeps(1) {
    const scalar_t one(1.0);
    inner_omega = one;
  }

  /**
   * @brief Construct a new Two Stage Gauss Seidel Handle object
   *
   */
  TwoStageGaussSeidelHandle() : TwoStageGaussSeidelHandle(GSHandle(GS_TWOSTAGE)) {}

  /**
   * @brief Construct a new Two Stage Gauss Seidel Handle object
   *
   * @param handle_exec_space The execution space instance
   * @param n_streams the number of streams
   */
  TwoStageGaussSeidelHandle(HandleExecSpace handle_exec_space, int n_streams)
      : TwoStageGaussSeidelHandle(GSHandle(handle_exec_space, n_streams, GS_TWOSTAGE)) {}

  // Sweep direction
  void setSweepDirection(GSDirection direction_) { this->direction = direction_; }
  GSDirection getSweepDirection() { return this->direction; }

  // specify whether to perform inner sweeps
  void setTwoStage(bool two_stage_) { this->two_stage = two_stage_; }
  bool isTwoStage() { return this->two_stage; }

  // specify whether to use compact form of recurrence
  void setCompactForm(bool compact_form_) { this->compact_form = compact_form_; }
  bool isCompactForm() { return this->compact_form; }

  // Number of outer sweeps
  void setNumOuterSweeps(int num_outer_sweeps_) { this->num_outer_sweeps = num_outer_sweeps_; }
  int getNumOuterSweeps() { return this->num_outer_sweeps; }

  // Number of inner sweeps
  void setNumInnerSweeps(int num_inner_sweeps_) { this->num_inner_sweeps = num_inner_sweeps_; }
  int getNumInnerSweeps() { return this->num_inner_sweeps; }

  // Inner damping factor
  void setInnerDampFactor(scalar_t inner_omega_) { this->inner_omega = inner_omega_; }
  scalar_t getInnerDampFactor() { return this->inner_omega; }

  // Workspaces
  // > diagonal (inverse)
  void setD(values_view_t D_) { this->D = D_; }
  values_view_t getD() { return this->D; }
  // > Lower part of diagonal block
  void setL(crsmat_t L) { this->crsmatL = L; }
  crsmat_t getL() { return this->crsmatL; }
  // > Upper part of diagonal block
  void setU(crsmat_t U) { this->crsmatU = U; }
  crsmat_t getU() { return this->crsmatU; }
  // > Complement of U
  void setLa(crsmat_t La) { this->crsmatLa = La; }
  crsmat_t getLa() { return this->crsmatLa; }
  // > Complement of L
  void setUa(crsmat_t Ua) { this->crsmatUa = Ua; }
  crsmat_t getUa() { return this->crsmatUa; }
  // > diagonal (not-inverse)
  void setDa(values_view_t Da_) { this->Da = Da_; }
  values_view_t getDa() { return this->Da; }

  void initVectors(int nrows_, int nrhs_) {
    if (this->nrows != nrows_ || this->nrhs != nrhs_) {
      this->localR = vector_view_t("temp", nrows_, nrhs_);
      this->localT = vector_view_t("temp", nrows_, nrhs_);
      this->localZ = vector_view_t("temp", nrows_, nrhs_);
      this->nrows  = nrows_;
      this->nrhs   = nrhs_;
    }
  }
  vector_view_t getVectorR() { return this->localR; }
  vector_view_t getVectorT() { return this->localT; }
  vector_view_t getVectorZ() { return this->localZ; }

 private:
  int nrows;
  int nrhs;

  // workspaces
  // > A = L + D + U
  values_view_t D;
  crsmat_t crsmatL;
  crsmat_t crsmatU;
  // > complements for compact form of recurrence
  //   where La = A - U and Ua = A - L
  values_view_t Da;
  crsmat_t crsmatLa;
  crsmat_t crsmatUa;

  // > residual vector for outer GS, Rk = B-A*Xk
  vector_view_t localR;
  // > workspace used for inner JR (for SpMV)
  vector_view_t localT;
  // > solultion correction from inner JR
  vector_view_t localZ;

  // solver parameters
  GSDirection direction;
  bool two_stage;
  bool compact_form;
  int num_inner_sweeps;
  int num_outer_sweeps;
  scalar_t inner_omega;
};
// -------------------------------------
}  // namespace KokkosSparse
#endif
