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

#include <fstream>
#include <ostream>

#include <KokkosKernels_Utils.hpp>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_Error.hpp>

#ifndef _GRAPHCOLORDISTANCE2HANDLE_HPP
#define _GRAPHCOLORDISTANCE2HANDLE_HPP

namespace KokkosGraph {

enum GraphColoringAlgorithmDistance2 {
  COLORING_D2_DEFAULT,    // Distance-2 Graph Coloring default algorithm
  COLORING_D2_SERIAL,     // Distance-2 Graph Coloring (SERIAL)
  COLORING_D2_VB,         // Distance-2 Graph Coloring Vertex Based
  COLORING_D2_VB_BIT,     // Distance-2 Graph Coloring Vertex Based BIT
  COLORING_D2_VB_BIT_EF,  // Distance-2 Graph Coloring Vertex Based BIT + Edge
                          // Filtering
  COLORING_D2_NB_BIT      // Distance-2 Graph Coloring Net Based BIT
};

template <class size_type_, class color_t_, class lno_t_, class ExecutionSpace, class TemporaryMemorySpace,
          class PersistentMemorySpace>
class GraphColorDistance2Handle {
 public:
  using HandleExecSpace                          = ExecutionSpace;
  using HandleTempMemorySpace                    = TemporaryMemorySpace;
  using HandlePersistentMemorySpace              = PersistentMemorySpace;
  using size_type                                = typename std::remove_const<size_type_>::type;
  using const_size_type                          = const size_type;
  using nnz_lno_type                             = typename std::remove_const<lno_t_>::type;
  using const_nnz_lno_type                       = const nnz_lno_type;
  using color_type                               = typename std::remove_const<color_t_>::type;
  using const_color_type                         = const color_type;
  using color_view_type                          = typename Kokkos::View<color_type*, HandlePersistentMemorySpace>;
  using color_view_array_layout                  = typename color_view_type::array_layout;
  using color_view_device_type                   = typename color_view_type::device_type;
  using color_view_memory_traits                 = typename color_view_type::memory_traits;
  using color_host_view_type                     = typename color_view_type::HostMirror;
  using size_type_temp_work_view_type            = typename Kokkos::View<size_type*, HandleTempMemorySpace>;
  using size_type_persistent_work_view_type      = typename Kokkos::View<size_type*, HandlePersistentMemorySpace>;
  using size_type_persistent_work_host_view_type = typename size_type_persistent_work_view_type::HostMirror;
  using nnz_lno_temp_work_view_type              = typename Kokkos::View<nnz_lno_type*, HandleTempMemorySpace>;
  using nnz_lno_persistent_work_view_type        = typename Kokkos::View<nnz_lno_type*, HandlePersistentMemorySpace>;
  using nnz_lno_persistent_work_host_view_type   = typename nnz_lno_persistent_work_view_type::HostMirror;
  using team_policy_type                         = Kokkos::TeamPolicy<HandleExecSpace>;
  using team_member_type                         = typename team_policy_type::member_type;
  using non_const_1d_size_type_view_type         = typename Kokkos::View<size_t*>;

 private:
  // Parameters
  GraphColoringAlgorithmDistance2 coloring_algorithm_type;  // Which algorithm type to use.

  bool verbose;  // verbosity flag
  bool tictoc;   // print time at every step

  bool vb_edge_filtering;  // whether to do edge filtering or not in vertex
                           // based algorithms.

  int vb_chunk_size;             // the (minimum) size of the consecutive works that a
                                 // thread will be assigned to.
  int max_number_of_iterations;  // maximum allowed number of phases that

  // STATISTICS
  double overall_coloring_time;         // The overall time taken to color the graph.
                                        // In the case of the iterative calls.
  double overall_coloring_time_phase1;  //
  double overall_coloring_time_phase2;  //
  double overall_coloring_time_phase3;  // Some timer accumulators for internal
                                        // phases.
  double overall_coloring_time_phase4;  //
  double overall_coloring_time_phase5;  //
  double coloring_time;                 // the time that it took to color the graph

  bool use_vtx_list;
  nnz_lno_temp_work_view_type vertex_list;
  size_type vertex_list_size;

  int num_phases;  // Number of phases used by the coloring algorithm

  color_view_type vertex_colors;
  bool is_coloring_called_before;
  nnz_lno_type num_colors;

 public:
  /**
   * Default constructor
   */
  GraphColorDistance2Handle()
      : coloring_algorithm_type(COLORING_D2_DEFAULT),
        verbose(false),
        tictoc(false),
        vb_edge_filtering(false),
        vb_chunk_size(8),
        max_number_of_iterations(200),
        overall_coloring_time(0),
        overall_coloring_time_phase1(0),
        overall_coloring_time_phase2(0),
        overall_coloring_time_phase3(0),
        overall_coloring_time_phase4(0),
        overall_coloring_time_phase5(0),
        coloring_time(0),
        use_vtx_list(false),
        num_phases(0),
        vertex_colors(),
        is_coloring_called_before(false),
        num_colors(0) {
    this->choose_default_algorithm();
    this->set_defaults(this->coloring_algorithm_type);

    // Throw an error if PersistentMemSpace != TempMemSpace since we don't
    // support them being different (for now).
    if (!std::is_same<PersistentMemorySpace, TemporaryMemorySpace>::value) {
      std::string message =
          "Distance-2 Graph Coloring Handle does not currently support "
          "different mem spaces";
      KokkosKernels::Impl::throw_runtime_exception(message);
    }
  }

  /**
   * Changes the graph coloring algorithm.
   *
   * @param[in] col_algo Coloring algorithm, one of:
   *                     - COLORING_D2_DEFAULT
   *                     - COLORING_D2_SERIAL
   *                     - COLORING_D2_VB
   *                     - COLORING_D2_VB_BIT
   *                     - COLORING_D2_VB_BIT_EF
   *                     - COLORING_D2_NB_BIT
   *
   *  @param[in] set_default_parameters Whether or not to reset the default
   * parameters for the given algorithm. Default = true.
   *
   *  @return None
   */
  void set_algorithm(const GraphColoringAlgorithmDistance2& col_algo, bool set_default_parameters = true) {
    if (col_algo == COLORING_D2_DEFAULT) {
      this->choose_default_algorithm();
    } else {
      this->coloring_algorithm_type = col_algo;
    }
    if (set_default_parameters) {
      this->set_defaults(this->coloring_algorithm_type);
    }
  }

  /**
   * Chooses best algorithm based on the execution space.
   *
   * This chooses the best algorithm based on the execution space:
   * - COLORING_D2_SERIAL if the execution space is SERIAL (more work efficient
   * than NB_BIT)
   * - COLORING_D2_NB_BIT otherwise (fastest parallel algorithm)
   *
   */

  void choose_default_algorithm() {
    if (KokkosKernels::Impl::kk_get_exec_space_type<ExecutionSpace>() == KokkosKernels::Impl::Exec_SERIAL) {
      this->coloring_algorithm_type = COLORING_D2_SERIAL;
#ifdef VERBOSE
      std::cout << "Serial Execution Space, Default Algorithm: COLORING_D2_SERIAL\n";
#endif
    } else {
      this->coloring_algorithm_type = COLORING_D2_NB_BIT;
#ifdef VERBOSE
      std::cout << ExecutionSpace::name() << " Execution Space, Default Algorithm: COLORING_D2_NB_BIT\n";
#endif
    }
  }

  nnz_lno_type get_num_colors() {
    if (num_colors == 0)
      KokkosKernels::Impl::view_reduce_max<color_view_type, ExecutionSpace>(vertex_colors.extent(0), vertex_colors,
                                                                            num_colors);
    return num_colors;
  }

  /** \brief Sets Default Parameter settings for the given algorithm.
   */
  void set_defaults(const GraphColoringAlgorithmDistance2& col_algo) {
    switch (col_algo) {
      case COLORING_D2_SERIAL:
      case COLORING_D2_VB:
      case COLORING_D2_VB_BIT:
      case COLORING_D2_VB_BIT_EF:
      case COLORING_D2_NB_BIT:
        this->tictoc                   = false;
        this->vb_edge_filtering        = false;
        this->vb_chunk_size            = 8;
        this->max_number_of_iterations = 200;
        break;
      default: throw std::runtime_error("Unknown Distance-2 Graph Coloring Algorithm\n");
    }
  }

  /**
   * \brief Destructor
   */
  virtual ~GraphColorDistance2Handle(){};

  // getters and setters
  GraphColoringAlgorithmDistance2 get_coloring_algo_type() const { return this->coloring_algorithm_type; }

  bool get_verbose() const { return this->verbose; }
  double get_coloring_time() const { return this->coloring_time; }
  int get_max_number_of_iterations() const { return this->max_number_of_iterations; }
  int get_num_phases() const { return this->num_phases; }

  double get_overall_coloring_time() const { return this->overall_coloring_time; }
  double get_overall_coloring_time_phase1() const { return this->overall_coloring_time_phase1; }
  double get_overall_coloring_time_phase2() const { return this->overall_coloring_time_phase2; }
  double get_overall_coloring_time_phase3() const { return this->overall_coloring_time_phase3; }
  double get_overall_coloring_time_phase4() const { return this->overall_coloring_time_phase4; }
  double get_overall_coloring_time_phase5() const { return this->overall_coloring_time_phase5; }

  bool get_tictoc() const { return this->tictoc; }

  int get_vb_chunk_size() const { return this->vb_chunk_size; }

  bool get_vb_edge_filtering() const { return this->vb_edge_filtering; }

  color_view_type get_vertex_colors() const { return this->vertex_colors; }

  bool is_coloring_called() const { return this->is_coloring_called_before; }

  bool get_use_vtx_list() const { return this->use_vtx_list; }
  nnz_lno_temp_work_view_type get_vertex_list() const { return this->vertex_list; }
  size_type get_vertex_list_size() const { return this->vertex_list_size; }

  // setters
  void set_vertex_list(nnz_lno_temp_work_view_type vertex_list_, size_type vertex_list_size_) {
    this->vertex_list      = vertex_list_;
    this->vertex_list_size = vertex_list_size_;
    this->use_vtx_list     = true;
  }
  void set_coloring_called() { this->is_coloring_called_before = true; }

  void set_coloring_algo_type(const GraphColoringAlgorithmDistance2& col_algo) {
    this->coloring_algorithm_type = col_algo;
  }

  void set_verbose(const bool verbose_) { this->verbose = verbose_; }
  void set_coloring_time(const double& coloring_time_) { this->coloring_time = coloring_time_; }
  void set_max_number_of_iterations(const int& max_phases) { this->max_number_of_iterations = max_phases; }
  void set_num_phases(const int& num_phases_) { this->num_phases = num_phases_; }

  void add_to_overall_coloring_time(const double& coloring_time_) { this->overall_coloring_time += coloring_time_; }
  void add_to_overall_coloring_time_phase1(const double& coloring_time_) {
    this->overall_coloring_time_phase1 += coloring_time_;
  }
  void add_to_overall_coloring_time_phase2(const double& coloring_time_) {
    this->overall_coloring_time_phase2 += coloring_time_;
  }
  void add_to_overall_coloring_time_phase3(const double& coloring_time_) {
    this->overall_coloring_time_phase3 += coloring_time_;
  }
  void add_to_overall_coloring_time_phase4(const double& coloring_time_) {
    this->overall_coloring_time_phase4 += coloring_time_;
  }
  void add_to_overall_coloring_time_phase5(const double& coloring_time_) {
    this->overall_coloring_time_phase5 += coloring_time_;
  }

  void set_tictoc(const bool use_tictoc) { this->tictoc = use_tictoc; }

  void set_vb_chunk_size(const int& chunksize) { this->vb_chunk_size = chunksize; }

  void set_vb_edge_filtering(const bool& use_vb_edge_filtering) { this->vb_edge_filtering = use_vb_edge_filtering; }

  void set_vertex_colors(const color_view_type vertex_colors_) {
    this->vertex_colors             = vertex_colors_;
    this->is_coloring_called_before = true;
    this->num_colors                = 0;
  }

  /**
   * Print / write out the graph in a GraphVIZ format.
   *
   * - Nodes colored with 1-5 will be filled in with some colors.
   * - Nodes colored > 5 will be unfilled (i.e., white background).
   * - Nodes colored with 0 (i.e., uncolored) will be filled in with red
   *   and will have a dashed border line. Uncolored nodes indicate a failed
   *   coloring.
   *
   * @param[in] os use std::cout for output to STDOUT stream, or a ofstream
   * object (i.e., `std::ofstream os("G.dot", std::ofstream::out);`) to write to
   * a file.
   */
  template <typename kokkos_view_type, typename rowmap_type, typename entries_type>
  void dump_graphviz(std::ostream& os, const size_t num_verts, rowmap_type& rowmap, entries_type& entries,
                     kokkos_view_type& colors) const {
    using h_colors_type  = typename kokkos_view_type::HostMirror;
    using h_rowmap_type  = typename rowmap_type::HostMirror;
    using h_entries_type = typename entries_type::HostMirror;

    h_colors_type h_colors = Kokkos::create_mirror_view(colors);
    Kokkos::deep_copy(h_colors, colors);

    h_rowmap_type h_rowmap = Kokkos::create_mirror_view(rowmap);
    Kokkos::deep_copy(h_rowmap, rowmap);

    h_entries_type h_entries = Kokkos::create_mirror_view(entries);
    Kokkos::deep_copy(h_entries, entries);

    os << "Graph G" << std::endl
       << "{" << std::endl
       << "    rankdir=LR;" << std::endl
       << "    overlap=false;" << std::endl
       << "    splines=true;" << std::endl
       << "    maxiter=2000;" << std::endl
       << "    node [shape=Mrecord];" << std::endl
       << "    edge [len=2.0];" << std::endl
       << std::endl;

    for (size_t vid = 0; vid < num_verts; vid++) {
      // Set color scheme for colors 1-5
      std::string fillcolor = "";
      std::string penwidth  = "";
      std::string color     = "";
      std::string style     = "";
      std::string fontcolor = "";
      if (1 == h_colors(vid)) {
        fillcolor = ", fillcolor=\"#CC4C51\"";
        style     = ", style=\"filled\"";
      } else if (2 == h_colors(vid)) {
        fillcolor = ", fillcolor=\"#CEBA5A\"";
        style     = ", style=\"filled\"";
      } else if (3 == h_colors(vid)) {
        fillcolor = ", fillcolor=\"#838FA4\"";
        style     = ", style=\"filled\"";
      } else if (4 == h_colors(vid)) {
        fillcolor = ", fillcolor=\"#C3AD86\"";
        style     = ", style=\"filled\"";
      } else if (5 == h_colors(vid)) {
        fillcolor = ", fillcolor=\"#935C44\"";
        style     = ", style=\"filled\"";
      } else if (0 == h_colors(vid)) {
        style     = ", style=\"filled,dashed\"";
        fillcolor = ", fillcolor=\"red\"";
        fontcolor = ", fontcolor=\"black\"";
        color     = ", color=\"black\"";
        penwidth  = ", penwidth=\"2.0\"";
      }

      os << "    " << vid << " [ label=\"" << vid << "|" << h_colors(vid) << "\"" << style << fontcolor << color
         << fillcolor << penwidth << "];" << std::endl;

      // Add the node's edges
      for (size_t iadj = h_rowmap(vid); iadj < (size_t)h_rowmap(vid + 1); iadj++) {
        size_t vadj = h_entries(iadj);
        if (vadj >= vid) {
          os << "    " << vid << " -- " << vadj << ";" << std::endl;
        }
      }
      os << std::endl;
    }
    os << "}" << std::endl;
  }  // dump_graphviz (end)

  const char* getD2AlgorithmName() const {
    switch (coloring_algorithm_type) {
      case COLORING_D2_DEFAULT: return "COLORING_D2_DEFAULT";
      case COLORING_D2_SERIAL: return "COLORING_D2_SERIAL";
      case COLORING_D2_VB: return "COLORING_D2_VB";
      case COLORING_D2_VB_BIT: return "COLORING_D2_VB_BIT";
      case COLORING_D2_VB_BIT_EF: return "COLORING_D2_VB_BIT_EF";
      case COLORING_D2_NB_BIT: return "COLORING_D2_NB_BIT";
    }
    return "ERROR: unregistered algorithm";
  }
};  // class GraphColorDistance2Handle (end)

}  // namespace KokkosGraph

#endif  // _GRAPHCOLORHANDLE_HPP
