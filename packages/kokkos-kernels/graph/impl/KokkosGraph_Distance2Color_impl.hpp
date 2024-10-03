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

#ifndef _KOKKOSGRAPH_DISTANCE2COLOR_IMPL_HPP
#define _KOKKOSGRAPH_DISTANCE2COLOR_IMPL_HPP

#include <iomanip>
#include <stdexcept>
#include <vector>
#include <type_traits>

#include <Kokkos_Core.hpp>

#include <KokkosKernels_Uniform_Initialized_MemoryPool.hpp>
#include <KokkosKernels_HashmapAccumulator.hpp>
#include <KokkosKernels_BitUtils.hpp>

#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosGraph_Distance1ColorHandle.hpp"  // todo: remove this  (SCAFFOLDING - WCMCLEN)
#include "KokkosGraph_Distance2ColorHandle.hpp"
#include "KokkosKernels_Handle.hpp"

namespace KokkosGraph {

namespace Impl {

#define VB_D2_COLORING_FORBIDDEN_SIZE 64
#define VBBIT_D2_COLORING_FORBIDDEN_SIZE 64

/*!
 * \brief Distance-2 Graph Coloring class
 *
 * This class supports direct methods for distance-2 graph coloring.
 *
 * doing_bipartite == false: (xadj,adj) is expected to be a
 * symmetric graph and the algorithm will always check for distance-1 conflicts.
 * (t_xadj,t_adj) should alias (xadj,adj).
 *
 * doing_bipartite == true: (t_xadj,t_adj) must be the transpose of (xadj,adj).
 * Distance-1 conflicts will not be checked.
 *
 */
template <typename HandleType, typename rowmap_t, typename entries_t, bool doing_bipartite>
class GraphColorDistance2 {
  // Need mutable entries type for edge filtering
  using nc_entries_t = typename entries_t::non_const_type;

 public:
  using execution_space = typename HandleType::HandleExecSpace;
  using memory_space    = typename HandleType::HandleTempMemorySpace;
  using device_type     = Kokkos::Device<execution_space, memory_space>;
  using color_view_type = typename HandleType::color_view_type;
  using color_type      = typename HandleType::color_type;
  using size_type       = typename rowmap_t::non_const_value_type;
  using lno_t           = typename entries_t::non_const_value_type;
  // Temporary view of lno_t
  using lno_view_t = typename HandleType::nnz_lno_temp_work_view_type;
  // Single-element (0-dim) view of lno_t:
  using single_lno_view_t     = typename Kokkos::View<lno_t, memory_space>;
  using range_policy_type     = Kokkos::RangePolicy<execution_space>;
  using team_policy_type      = Kokkos::TeamPolicy<execution_space>;
  using team_member_type      = typename team_policy_type::member_type;
  using bool_view_t           = Kokkos::View<bool*, memory_space>;
  using bit_64_forbidden_type = uint64_t;
  using bitset_t              = Kokkos::Bitset<device_type>;
  using forbidden_view        = Kokkos::View<uint32_t*, device_type>;

 protected:
  lno_t nr;               // num_rows  (also #verts, the objects being colored)
  lno_t nc;               // num cols
  size_type ne;           // # edges
  rowmap_t xadj;          // rowmap, transpose of rowmap
  entries_t adj;          // entries, transpose of entries   (size = # edges)
  rowmap_t t_xadj;        // transpose rowmap (aliases xadj if !doing_bipartite)
  entries_t t_adj;        // transpose entries (aliases adj if !doing_bipartite)
  HandleType* gc_handle;  // pointer to the graph coloring handle

 private:
  int _chunkSize;  // the size of the minimum work unit assigned to threads.
                   // Changes the convergence on GPUs
  int _max_num_iterations;
  bool _ticToc;   // if true print info in each step
  bool _verbose;  // if true, print verbose information
  bool using_edge_filtering;
  static constexpr color_type UNPROCESSED = 0;
  static constexpr color_type UNCOLORABLE = ~UNPROCESSED;
  static constexpr color_type CONFLICTED  = UNCOLORABLE - 1;

 public:
  /**
   * \brief GraphColorDistance2 constructor.
   * \param nr_: number of vertices in the graph
   * \param row_map: the xadj array of the graph. Its size is nr_ +1
   * \param entries: adjacency array of the graph.
   * \param handle: GraphColoringHandle object that holds the specification
   * about the graph coloring, including parameters.
   */
  GraphColorDistance2(lno_t nr_, lno_t nc_, rowmap_t row_map, entries_t entries, rowmap_t t_row_map,
                      entries_t t_entries, HandleType* handle)
      : nr(nr_),
        nc(nc_),
        ne(entries.extent(0)),
        xadj(row_map),
        adj(entries),
        t_xadj(t_row_map),
        t_adj(t_entries),
        gc_handle(handle),
        _chunkSize(handle->get_vb_chunk_size()),
        _max_num_iterations(handle->get_max_number_of_iterations()),
        _ticToc(handle->get_verbose()),
        _verbose(handle->get_verbose()) {
    // Internal logic check: distance-2 coloring (non-bipartite) requires a
    // square graph
    if (!doing_bipartite && nr != nc) {
      throw std::runtime_error(
          "D2 INTERNAL ERROR: requested undirected d2 coloring but input graph "
          "is not square (nr_ != nc_)");
    }
  }

  /**
   *  \brief GraphColor destructor.
   */
  virtual ~GraphColorDistance2() {}

  // -----------------------------------------------------------------
  //
  // GraphColorDistance2::execute()
  //
  // -----------------------------------------------------------------

  /**
   * \brief Computes the distance-2 graph coloring.
   */
  void compute_distance2_color() {
    // Delegate to different coloring functions, depending on algorithm
    using_edge_filtering = false;
    color_view_type colors_out;
    if (gc_handle->get_vertex_colors().use_count() > 0) {
      colors_out = gc_handle->get_vertex_colors();
    } else {
      colors_out = color_view_type("Graph Colors", this->nr);
    }
    switch (this->gc_handle->get_coloring_algo_type()) {
      case COLORING_D2_VB_BIT_EF: using_edge_filtering = true; [[fallthrough]];
      case COLORING_D2_VB_BIT:
      case COLORING_D2_VB: compute_d2_coloring_vb(colors_out); break;
      case COLORING_D2_NB_BIT: compute_d2_coloring_nb(colors_out); break;
      case COLORING_D2_SERIAL: compute_d2_coloring_serial(colors_out); break;
      default:
        throw std::runtime_error(std::string("D2 coloring handle has invalid algorithm: ") +
                                 std::to_string((int)this->gc_handle->get_coloring_algo_type()));
    }
  }

  void compute_d2_coloring_vb(const color_view_type& colors_out) {
    // Data:
    // gc_handle = graph coloring handle
    // nr        = num_rows  (scalar)
    // nc        = num_cols  (scalar)
    // xadj      = row_map   (view 1 dimension - [num_verts+1] - entries index
    // into adj ) adj       = entries   (view 1 dimension - [num_edges]   -
    // adjacency list )
    if (this->_ticToc) {
      std::cout << "\tcolor_graph_d2 params:" << std::endl
                << "\t  algorithm                : " << this->gc_handle->getD2AlgorithmName() << std::endl
                << "\t  ticToc                   : " << this->_ticToc << std::endl
                << "\t  max_num_iterations       : " << this->_max_num_iterations << std::endl
                << "\t  chunkSize                : " << this->_chunkSize << std::endl
                << "\t  Edge Filtering Pass?     : " << (int)using_edge_filtering << std::endl
                << "\tgraph information:" << std::endl
                << "\t  nr                       : " << this->nr << std::endl
                << "\t  ne                       : " << this->ne << std::endl;
      /*
      // For extra verbosity if debugging...
      prettyPrint1DView(this->xadj,   ">>> xadj   ", 500);
      prettyPrint1DView(this->adj,    ">>>  adj   ", 500);
      prettyPrint1DView(this->t_xadj, ">>> t_xadj ", 500);
      prettyPrint1DView(this->t_adj,  ">>> t_adj  ", 500);
      */
    }

    // conflictlist - store conflicts that can happen when we're coloring in
    // parallel.
    lno_view_t current_vertexList(Kokkos::view_alloc(Kokkos::WithoutInitializing, "vertexList"), this->nr);

    lno_t current_vertexListLength = this->nr;

    if (this->gc_handle->get_use_vtx_list()) {
      // init conflict list from coloring handle
      current_vertexList       = this->gc_handle->get_vertex_list();
      current_vertexListLength = this->gc_handle->get_vertex_list_size();
    } else {
      // init conflictlist sequentially.
      Kokkos::parallel_for("InitList", range_policy_type(0, this->nr), functorInitList<lno_view_t>(current_vertexList));
    }
    // Next iteratons's conflictList
    lno_view_t next_iteration_recolorList(Kokkos::view_alloc(Kokkos::WithoutInitializing, "recolorList"), this->nr);

    // Size the next iteration conflictList
    single_lno_view_t next_iteration_recolorListLength("recolorListLength");

    lno_t numUncolored             = this->nr;
    lno_t numUncoloredPreviousIter = this->nr + 1;

    double time;
    double total_time = 0.0;
    Kokkos::Timer timer;

    int iter = 0;
    for (; (iter < _max_num_iterations) && (numUncolored > 0); iter++) {
      timer.reset();

      // Save the # of uncolored from the previous iteration
      numUncoloredPreviousIter = numUncolored;

      // ------------------------------------------
      // Do greedy color
      // ------------------------------------------
      if (using_edge_filtering) {
        // Temporary mutable copies of adj array
        // * This is required for edge-filtering passes to avoid
        //   side effects since edge filtering modifies the adj array.
        // * Allocate using lno_view_t (managed) but then access as an
        // entries_t,
        //   so that it has the same type as adj
        // * on the other hand, t_adj is not actually modified by EF functor
        lno_view_t adj_copy(Kokkos::view_alloc(Kokkos::WithoutInitializing, "adj copy"), this->ne);
        Kokkos::deep_copy(adj_copy, this->adj);
        this->colorGreedyEF(this->xadj, adj_copy, this->t_xadj, this->t_adj, colors_out);
      } else {
        this->colorGreedy(this->xadj, this->adj, this->t_xadj, this->t_adj, colors_out, current_vertexList,
                          current_vertexListLength);
      }

      execution_space().fence();

      if (this->_ticToc) {
        time = timer.seconds();
        total_time += time;
        std::cout << "\tIteration: " << iter << std::endl
                  << "\t  - Time speculative greedy phase : " << time << std::endl
                  << "\t  - Num Uncolored (greedy-color)  : " << numUncolored << std::endl;

        gc_handle->add_to_overall_coloring_time_phase1(time);

        timer.reset();
      }

      // ------------------------------------------
      // Find conflicts
      // ------------------------------------------
      bool swap_work_arrays = true;  // NOTE: swap_work_arrays can go away in
                                     // this example -- was only ever
                                     //       set false in the PPS code in the
                                     //       original D1 coloring...

      // NOTE: not using colorset algorithm in this so we don't include colorset
      // data
      numUncolored = this->findConflicts(swap_work_arrays, this->xadj, this->adj, this->t_xadj, this->t_adj, colors_out,
                                         current_vertexList, current_vertexListLength, next_iteration_recolorList,
                                         next_iteration_recolorListLength);

      execution_space().fence();

      if (_ticToc) {
        time = timer.seconds();
        total_time += time;
        std::cout << "\t  - Time conflict detection       : " << time << std::endl;
        std::cout << "\t  - Num Uncolored (conflicts)     : " << numUncolored << std::endl;
        gc_handle->add_to_overall_coloring_time_phase2(time);
        timer.reset();
      }

      // Swap the work arrays (for conflictlist)
      if (swap_work_arrays) {
        // Swap Work Arrays
        if (iter + 1 < this->_max_num_iterations) {
          lno_view_t temp            = current_vertexList;
          current_vertexList         = next_iteration_recolorList;
          next_iteration_recolorList = temp;

          current_vertexListLength         = numUncolored;
          next_iteration_recolorListLength = single_lno_view_t("recolorListLength");
        }
      }

      // Bail out if we didn't make any progress since we're probably stuck and
      // it's better to just clean up in serial.
      if (numUncolored == numUncoloredPreviousIter) break;

    }  // end for iterations && numUncolored > 0...

    // ------------------------------------------
    // clean up in serial (resolveConflictsSerial)
    // ------------------------------------------
    if (numUncolored > 0) {
      this->resolveConflictsSerial(this->xadj, this->adj, this->t_xadj, this->t_adj, colors_out, current_vertexList,
                                   current_vertexListLength);
    }

    execution_space().fence();

    if (_ticToc) {
      time = timer.seconds();
      total_time += time;
      std::cout << "\tTime serial conflict resolution   : " << time << std::endl;
      std::cout << "\tTotal time for coloring           : " << total_time << std::endl;
      gc_handle->add_to_overall_coloring_time_phase3(time);
    }

    // Save the number of phases and vertex colors to the graph coloring handle
    this->gc_handle->set_vertex_colors(colors_out);
    this->gc_handle->set_num_phases(iter);

  }  // color_graph_d2 (end)

  template <int batch>
  struct NB_Coloring {
    NB_Coloring(const lno_view_t& worklist_, const single_lno_view_t& worklen_, color_type colorBase_,
                const forbidden_view& forbidden_, color_view_type colors_, const rowmap_t& Vrowmap_,
                const entries_t& Vcolinds_, lno_t vertsPerThread_, lno_t numCols_)
        : worklist(worklist_),
          worklen(worklen_),
          colorBase(colorBase_),
          forbidden(forbidden_),
          colors(colors_),
          Vrowmap(Vrowmap_),
          Vcolinds(Vcolinds_),
          vertsPerThread(vertsPerThread_),
          numCols(numCols_) {}
    KOKKOS_INLINE_FUNCTION void operator()(const lno_t ti) const {
      for (lno_t i = ti * vertsPerThread; i < (ti + 1) * vertsPerThread; i++) {
        if (i >= worklen()) return;
        lno_t v = worklist(i);
        // compute forbidden for v
        unsigned forbid[batch] = {0};
        // union the forbidden of all incident c's
        size_type rowBegin = Vrowmap(v);
        size_type rowEnd   = Vrowmap(v + 1);
        if (!doing_bipartite) {
          // gather distance-1 forbidden colors
          for (int b = 0; b < batch; b++) forbid[b] = forbidden(v * batch + b);
        }
        // gather distance-2 forbidden colors
        for (size_type j = rowBegin; j < rowEnd; j++) {
          lno_t nei = Vcolinds(j);
          if (nei < numCols) {
            for (int b = 0; b < batch; b++) forbid[b] |= forbidden(nei * batch + b);
          }
        }
        // Find the first 0 bit in forbid
        color_type color = 0;
        int colorWord    = 0;
        int colorBit     = 0;
        for (int b = 0; b < batch; b++) {
          if (~forbid[b]) {
            // least_set_bit returns 1 for the least significant bit, so
            // subtracting 1
            colorWord = b;
            colorBit  = KokkosKernels::Impl::least_set_bit(~forbid[b]) - 1;
            color     = colorBase + 32 * b + colorBit;
            break;
          }
        }
        if (color && (colors(v) == 0 || colors(v) == CONFLICTED || colors(v) == UNCOLORABLE)) {
          // Color v
          colors(v) = color;
          if (!doing_bipartite) {
            // Update forbidden for v (preventing dist-1 conflicts)
            if (v < numCols) Kokkos::atomic_fetch_or(&forbidden(v * batch + colorWord), (uint32_t)1 << colorBit);
          }
          // Update forbidden for all of v's neighbors
          for (size_type j = rowBegin; j < rowEnd; j++) {
            lno_t nei = Vcolinds(j);
            if (nei < numCols) {
              // Update column forbidden
              Kokkos::atomic_fetch_or(&forbidden(nei * batch + colorWord), (uint32_t)1 << colorBit);
            }
          }
        } else if (colors(v) == 0 || colors(v) == CONFLICTED || colors(v) == UNCOLORABLE) {
          colors(v) = UNCOLORABLE;
        }
      }
    }

    lno_view_t worklist;
    single_lno_view_t worklen;
    color_type colorBase;
    forbidden_view forbidden;  // forbidden color bitset for columns
    color_view_type colors;
    rowmap_t Vrowmap;  // V <-> C graph (row v is the columns incident to v)
    entries_t Vcolinds;
    lno_t vertsPerThread;
    const lno_t numCols;
  };

  template <int batch>
  struct NB_Conflict {
    NB_Conflict(color_type colorBase_, const forbidden_view& forbidden_, const color_view_type& colors_,
                const rowmap_t& Crowmap_, const entries_t& Ccolinds_, lno_t numVerts_)
        : colorBase(colorBase_),
          forbidden(forbidden_),
          colors(colors_),
          Crowmap(Crowmap_),
          Ccolinds(Ccolinds_),
          numVerts(numVerts_) {}

    KOKKOS_INLINE_FUNCTION void operator()(const lno_t c) const {
      // Here, only processing 32 colors at a time.
      // forbidNei is the (minimum) neighbor ID where each forbidden color was
      // observed. This is why the whole batch can't be processed at once
      lno_t forbidNei[32];
      // Go over all the v neighbors, updating forbidden
      size_type rowBegin = Crowmap(c);
      size_type rowEnd   = Crowmap(c + 1);
      for (int b = 0; b < batch; b++) {
        unsigned forbid       = 0U;
        color_type batchBegin = colorBase + 32 * b;
        for (size_type j = rowBegin; j <= rowEnd; j++) {
          lno_t nei;
          if (j == rowEnd) {
            if (!doing_bipartite)  // note: compile-time branch
              nei = c;
            else
              break;
          } else
            nei = Ccolinds(j);
          if (nei >= numVerts) continue;
          color_type neiColor = colors(nei);
          int colorOffset     = neiColor - batchBegin;
          if (colorOffset >= 0 && colorOffset < 32) {
            // if this is the first time the color has been seen, register nei
            // in forbidNei
            unsigned mask = 1U << colorOffset;
            if (0 == (forbid & mask)) {
              // First time seeing this color
              forbidNei[colorOffset] = nei;
            } else {
              // Have seen this color before:
              // must uncolor either nei or forbidNei[colorOffset] (whichever
              // has higher ID)
              if (nei > forbidNei[colorOffset]) {
                // nei has a higher id than another neighbor of the same color,
                // so uncolor nei
                colors(nei) = CONFLICTED;
              } else if (nei < forbidNei[colorOffset]) {
                colors(forbidNei[colorOffset]) = CONFLICTED;
                forbidNei[colorOffset]         = nei;
              }
            }
            forbid |= mask;
          }
        }
      }
    }

    color_type colorBase;
    forbidden_view forbidden;  // forbidden color bitset for columns
    color_view_type colors;
    rowmap_t Crowmap;  // C <-> V graph (row c is the vertices incident to c)
    entries_t Ccolinds;
    const lno_t numVerts;
  };

  template <int batch>
  struct NB_RefreshForbidden {
    NB_RefreshForbidden(color_type colorBase_, const forbidden_view& forbidden_, const color_view_type& colors_,
                        const rowmap_t& Crowmap_, const entries_t& Ccolinds_, lno_t numVerts_)
        : colorBase(colorBase_),
          colorEnd(colorBase + 32 * batch),
          forbidden(forbidden_),
          colors(colors_),
          Crowmap(Crowmap_),
          Ccolinds(Ccolinds_),
          numVerts(numVerts_) {}

    KOKKOS_INLINE_FUNCTION void operator()(const lno_t c) const {
      // compute this in registers before storing to forbidden
      unsigned newForbid[batch] = {0};
      // Go over all the v neighbors, updating forbidden
      size_type rowBegin = Crowmap(c);
      size_type rowEnd   = Crowmap(c + 1);
      if (!doing_bipartite) {
        // first, add d-1 conflict
        color_type selfColor = colors(c);
        if (colorBase <= selfColor && selfColor < colorEnd) {
          int colorWord = (selfColor - colorBase) / 32;
          int colorBit  = (selfColor - colorBase) % 32;
          newForbid[colorWord] |= ((uint32_t)1 << colorBit);
        }
      }
      for (size_type i = rowBegin; i < rowEnd; i++) {
        lno_t nei = Ccolinds(i);
        if (nei >= numVerts) continue;
        color_type neiColor = colors(nei);
        if (colorBase <= neiColor && neiColor < colorEnd) {
          int colorWord = (neiColor - colorBase) / 32;
          int colorBit  = (neiColor - colorBase) % 32;
          newForbid[colorWord] |= ((uint32_t)1 << colorBit);
        }
      }
      for (int i = 0; i < batch; i++) forbidden(c * batch + i) = newForbid[i];
    }

    color_type colorBase;
    color_type colorEnd;
    forbidden_view forbidden;  // forbidden color bitset for columns
    color_view_type colors;
    rowmap_t Crowmap;  // C <-> V graph (row c is the vertices incident to c)
    entries_t Ccolinds;
    const lno_t numVerts;
  };

  struct NB_Worklist {
    NB_Worklist(const color_view_type colors_, const lno_view_t& worklist_, const single_lno_view_t& worklen_,
                lno_t nr_)
        : colors(colors_), worklist(worklist_), worklen(worklen_), nr(nr_) {}

    KOKKOS_INLINE_FUNCTION void operator()(const lno_t v, lno_t& lnum, bool finalPass) const {
      if (colors(v) == CONFLICTED) {
        if (finalPass) worklist(lnum) = v;
        lnum++;
      }
      if (finalPass && v == nr - 1) {
        // The very last thread in the kernel knows how many items are in the
        // next worklist
        worklen() = lnum;
      }
    }

    color_view_type colors;
    lno_view_t worklist;
    single_lno_view_t worklen;
    lno_t nr;
  };

  struct NB_UpdateBatch {
    NB_UpdateBatch(const color_view_type& colors_, const lno_view_t& worklist_, const single_lno_view_t& worklen_,
                   lno_t nr_)
        : colors(colors_), worklist(worklist_), worklen(worklen_), nr(nr_) {}

    KOKKOS_INLINE_FUNCTION void operator()(const lno_t v, lno_t& lnum, bool finalPass) const {
      if (colors(v) == UNCOLORABLE) {
        if (finalPass) worklist(lnum) = v;
        lnum++;
      }
      if (finalPass && v == nr - 1) {
        // The very last thread in the kernel knows the length of the new
        // worklist.
        worklen() = lnum;
      }
    }

    color_view_type colors;
    lno_view_t worklist;
    single_lno_view_t worklen;
    lno_t nr;
  };

  void compute_d2_coloring_nb(const color_view_type& colors_out) {
    // Member data used:
    // gc_handle    = graph coloring handle
    // nr           = #vertices
    // nc           = #columns
    // xadj/adj     = graph where rows are vertices, and adjacent columns are
    // listed t_xadj/t_adj = graph where rows are columns, and adjacent vertices
    // are listed.
    //                Allowed to alias xadj/adj if same.
    if (this->_ticToc) {
      std::cout << "\tcolor_symmetric_graph_d2 params:\n"
                << "\t\t#vertices : " << this->nr << '\n'
                << "\t\t#edges: " << this->ne << '\n';
    }
    // Initialize worklist with every vertex
    lno_view_t worklist("Worklist", this->nr);
    single_lno_view_t worklen("Worklist length");
    Kokkos::deep_copy(worklen, this->nr);

    // init conflictlist sequentially.
    Kokkos::parallel_for("InitList", range_policy_type(0, this->nr), functorInitList<lno_view_t>(worklist));

    // Estimate the number of colors that will be needed
    // The algorithm can't use more colors than the max distance-2 degree,
    // but it can use fewer.
    // Here, subtracting 1 to represent the self-edge
    lno_t avgDeg = adj.extent(0) / (xadj.extent(0) - 1) - 1;
    // This d-2 chromatic number estimate is based on the following assumptions:
    //  -the number of self-loops to v is just deg(v), and these can't cause
    //  conflicts -each node has about the same degree, so
    //   the total number of length-2 walks from v is about deg(v)^2
    //  -the constant was determined experimentally
    int estNumColors = 0.1 * (avgDeg * (avgDeg - 1));
    // 8 words (256 bits/colors) is the maximum allowed batch size
    int batch = 8;
    // but don't use more than the estimate
    for (int tryBatch = 1; tryBatch < 8; tryBatch *= 2) {
      if (estNumColors <= 32 * tryBatch) {
        batch = tryBatch;
        break;
      }
    }
    const lno_t numVerts = this->nr;
    const lno_t numCols  = this->nc;
    // note: relying on forbidden and colors_out being initialized to 0
    forbidden_view forbidden("Forbidden", batch * numCols);
    int iter = 0;
    Kokkos::Timer timer;
    lno_t currentWork    = this->nr;
    batch                = 1;
    double colorTime     = 0;
    double conflictTime  = 0;
    double forbiddenTime = 0;
    double worklistTime  = 0;
    for (color_type colorBase = 1;; colorBase += 32 * batch) {
      // Until the worklist is completely empty, run the functor specialization
      // for batch size
      while (currentWork) {
        lno_t vertsPerThread = 1;
        lno_t workBatches    = (currentWork + vertsPerThread - 1) / vertsPerThread;
        timer.reset();
        // if still using this color set, refresh forbidden.
        // This avoids using too many colors, by relying on forbidden from
        // before previous conflict resolution (which is now stale). Refreshing
        // forbidden before conflict resolution ensures that previously-colored
        // vertices do not get recolored.
        switch (batch) {
          case 1:
            Kokkos::parallel_for(
                "NB D2 Forbidden", range_policy_type(0, numCols),
                NB_RefreshForbidden<1>(colorBase, forbidden, colors_out, this->t_xadj, this->t_adj, numVerts));
            break;
          case 2:
            Kokkos::parallel_for(
                "NB D2 Forbidden", range_policy_type(0, numCols),
                NB_RefreshForbidden<2>(colorBase, forbidden, colors_out, this->t_xadj, this->t_adj, numVerts));
            break;
          case 4:
            Kokkos::parallel_for(
                "NB D2 Forbidden", range_policy_type(0, numCols),
                NB_RefreshForbidden<4>(colorBase, forbidden, colors_out, this->t_xadj, this->t_adj, numVerts));
            break;
          case 8:
            Kokkos::parallel_for(
                "NB D2 Forbidden", range_policy_type(0, numCols),
                NB_RefreshForbidden<8>(colorBase, forbidden, colors_out, this->t_xadj, this->t_adj, numVerts));
            break;
          default:;
        }
        forbiddenTime += timer.seconds();
        timer.reset();
        switch (batch) {
          case 1:
            timer.reset();
            Kokkos::parallel_for("NB D2 Coloring", range_policy_type(0, workBatches),
                                 NB_Coloring<1>(worklist, worklen, colorBase, forbidden, colors_out, this->xadj,
                                                this->adj, vertsPerThread, numCols));
            colorTime += timer.seconds();
            timer.reset();
            Kokkos::parallel_for("NB D2 Conflict Resolution", range_policy_type(0, numCols),
                                 NB_Conflict<1>(colorBase, forbidden, colors_out, this->t_xadj, this->t_adj, numVerts));
            conflictTime += timer.seconds();
            break;
          case 2:
            timer.reset();
            Kokkos::parallel_for("NB D2 Coloring", range_policy_type(0, workBatches),
                                 NB_Coloring<2>(worklist, worklen, colorBase, forbidden, colors_out, this->xadj,
                                                this->adj, vertsPerThread, numCols));
            colorTime += timer.seconds();
            timer.reset();
            Kokkos::parallel_for("NB D2 Conflict Resolution", range_policy_type(0, numCols),
                                 NB_Conflict<2>(colorBase, forbidden, colors_out, this->t_xadj, this->t_adj, numVerts));
            conflictTime += timer.seconds();
            break;
          case 4:
            timer.reset();
            Kokkos::parallel_for("NB D2 Coloring", range_policy_type(0, workBatches),
                                 NB_Coloring<4>(worklist, worklen, colorBase, forbidden, colors_out, this->xadj,
                                                this->adj, vertsPerThread, numCols));
            colorTime += timer.seconds();
            timer.reset();
            Kokkos::parallel_for("NB D2 Conflict Resolution", range_policy_type(0, numCols),
                                 NB_Conflict<4>(colorBase, forbidden, colors_out, this->t_xadj, this->t_adj, numVerts));
            conflictTime += timer.seconds();
            break;
          case 8:
            timer.reset();
            Kokkos::parallel_for("NB D2 Coloring", range_policy_type(0, workBatches),
                                 NB_Coloring<8>(worklist, worklen, colorBase, forbidden, colors_out, this->xadj,
                                                this->adj, vertsPerThread, numCols));
            colorTime += timer.seconds();
            timer.reset();
            Kokkos::parallel_for("NB D2 Conflict Resolution", range_policy_type(0, numCols),
                                 NB_Conflict<8>(colorBase, forbidden, colors_out, this->t_xadj, this->t_adj, numVerts));
            conflictTime += timer.seconds();
            break;
          default:
            throw std::logic_error(
                "D2 symmetric color batch size is not a power-of-two, or is "
                "too big");
        }
        timer.reset();
        // Then build the next worklist
        Kokkos::parallel_scan("NB D2 worklist", range_policy_type(0, numVerts),
                              NB_Worklist(colors_out, worklist, worklen, numVerts), currentWork);
        worklistTime += timer.seconds();
        timer.reset();
        iter++;
      }
      // Will need to run with a different color base, so rebuild the work list
      Kokkos::parallel_scan("NB D2 Worklist Rebuild", range_policy_type(0, numVerts),
                            NB_UpdateBatch(colors_out, worklist, worklen, numVerts));
      Kokkos::deep_copy(currentWork, worklen);
      worklistTime += timer.seconds();
      timer.reset();
      if (currentWork == 0) {
        // Still have no work to do, meaning every vertex is colored
        break;
      }
      // Clear forbidden before continuing
      Kokkos::deep_copy(forbidden, 0U);
    }
    execution_space().fence();
    if (this->_ticToc) {
      std::cout << "~~ D2 timings ~~\n";
      std::cout << "Coloring: " << colorTime << '\n';
      std::cout << "Conflict: " << conflictTime << '\n';
      std::cout << "Forbidden: " << forbiddenTime << '\n';
      std::cout << "Worklist: " << worklistTime << '\n';
      std::cout << "** Total: " << colorTime + conflictTime + forbiddenTime + worklistTime << "\n\n";
    }
    if (this->_ticToc) {
      gc_handle->add_to_overall_coloring_time_phase1(timer.seconds());
      timer.reset();
    }

    // Save the number of phases and vertex colors to the graph coloring handle
    this->gc_handle->set_vertex_colors(colors_out);
    this->gc_handle->set_num_phases(iter);
  }

  void compute_d2_coloring_serial(const color_view_type& colors_out) {
    // Member data used:
    // gc_handle    = graph coloring handle
    // nr           = #vertices
    // nc           = #columns
    // xadj/adj     = graph where rows are vertices, and adjacent columns are
    // listed t_xadj/t_adj = graph where rows are columns, and adjacent vertices
    // are listed.
    //                Allowed to alias xadj/adj if same.
    if (this->_ticToc) {
      std::cout << "\tcolor_symmetric_graph_d2 params:\n"
                << "\t\t#vertices : " << this->nr << '\n'
                << "\t\t#edges: " << this->ne << '\n';
    }
    Kokkos::View<unsigned*, Kokkos::HostSpace> forbidden("Forbidden", this->nc);
    auto colors = Kokkos::create_mirror_view(colors_out);
    // Get the graph(s) in host space, if not already
    Kokkos::View<const size_type*, Kokkos::HostSpace> Vrowmap =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), this->xadj);
    Kokkos::View<const lno_t*, Kokkos::HostSpace> Vcolinds =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), this->adj);
    // Create worklist
    Kokkos::View<lno_t*, Kokkos::HostSpace> worklist(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Worklist"),
                                                     this->nr);
    int iter = 0;
    Kokkos::Timer timer;
    lno_t currentWork = this->nr;
    lno_t numCols     = this->nc;
    for (color_type colorBase = 1; currentWork > 0; colorBase += 32) {
      // Rebuilding the worklist in-place.
      lno_t worklistOutput = 0;
      for (lno_t i = 0; i < currentWork; i++) {
        lno_t v = i;
        if (iter > 0) v = worklist(i);
        // Compute v's forbidden for this batch
        unsigned forbid    = 0;
        size_type rowBegin = Vrowmap(v);
        size_type rowEnd   = Vrowmap(v + 1);
        // always include the diagonal (self-edge), to avoid distance-1
        // conflicts.
        if (!doing_bipartite) forbid |= forbidden(v);
        for (size_type j = rowBegin; j < rowEnd; j++) {
          lno_t nei = Vcolinds(j);
          if (nei < numCols) forbid |= forbidden(nei);
        }
        if (~forbid) {
          int bitOffset = KokkosKernels::Impl::least_set_bit(~forbid) - 1;
          colors(v)     = colorBase + bitOffset;
          // Together with including diagonal, setting forbidden(v)
          // with v's color will prevent all distance-1 conflicts
          if (v < numCols) forbidden(v) |= (1U << bitOffset);
          for (size_type j = rowBegin; j < rowEnd; j++) {
            lno_t nei = Vcolinds(j);
            // marking forbidden on out-neighbors of v prevents distance-2
            // conflicts
            if (nei < numCols) forbidden(nei) |= (1U << bitOffset);
          }
        } else {
          // Can't color in this batch, so add to worklist
          worklist(worklistOutput++) = v;
        }
      }
      currentWork = worklistOutput;
      // Clear all forbidden bits
      Kokkos::deep_copy(forbidden, 0U);
      iter++;
    }
    if (this->_ticToc) {
      gc_handle->add_to_overall_coloring_time_phase1(timer.seconds());
    }
    // Save the number of phases and vertex colors to the graph coloring handle
    Kokkos::deep_copy(colors_out, colors);
    this->gc_handle->set_vertex_colors(colors_out);
    this->gc_handle->set_num_phases(iter);
  }

 private:
  // -----------------------------------------------------------------
  //
  // GraphColorDistance2::colorGreedy()
  //
  // -----------------------------------------------------------------
  void colorGreedy(rowmap_t xadj_, entries_t adj_, rowmap_t t_xadj_, entries_t t_adj_, color_view_type vertex_colors_,
                   lno_view_t current_vertexList_, lno_t current_vertexListLength_) {
    lno_t chunkSize_ = this->_chunkSize;

    if (current_vertexListLength_ < 100 * chunkSize_) {
      chunkSize_ = 1;
    }

    // Pick the right coloring algorithm to use based on which algorithm we're
    // using
    switch (this->gc_handle->get_coloring_algo_type()) {
      // Single level parallelism on chunks
      // 1. [P] loop over vertices
      // 2. [S] loop over color offset blocks
      // 3. [S] loop over vertex neighbors
      // 4. [S] loop over vertex neighbors of neighbors
      case COLORING_D2_VB: {
        functorGreedyColorVB gc(this->nr, this->nc, xadj_, adj_, t_xadj_, t_adj_, vertex_colors_, current_vertexList_,
                                current_vertexListLength_);
        Kokkos::parallel_for("LoopOverChunks", range_policy_type(0, this->nr), gc);
      } break;

      // One level Perallelism, BIT Array for coloring
      // 1. [P] loop over vertices
      // 2. [S] loop over color offset blocks
      // 3. [S] loop over vertex neighbors
      // 4. [S] loop over vertex neighbors of neighbors
      case COLORING_D2_VB_BIT: {
        functorGreedyColorVB_BIT gc(this->nr, this->nc, xadj_, adj_, t_xadj_, t_adj_, vertex_colors_,
                                    current_vertexList_, current_vertexListLength_);
        Kokkos::parallel_for("LoopOverChunks", range_policy_type(0, this->nr), gc);
      } break;

      default:
        throw std::invalid_argument(
            "Unknown Distance-2 Algorithm Type or invalid for non Edge "
            "Filtering mode.");
    }

  }  // colorGreedy (end)

  // -----------------------------------------------------------------
  //
  // GraphColorDistance2::colorGreedyEF()
  //
  // -----------------------------------------------------------------
  void colorGreedyEF(rowmap_t xadj_, lno_view_t adj_copy_, rowmap_t t_xadj_, entries_t t_adj_copy_,
                     color_view_type vertex_colors_) {
    // Pick the right coloring algorithm to use based on which algorithm we're
    // using
    switch (this->gc_handle->get_coloring_algo_type()) {
      // One level parallelism, BIT Array for coloring + edge filtering
      // 1. [P] loop over vertices
      // 2. [S] loop over color offset blocks
      // 3. [S] loop over vertex neighbors
      // 4. [S] loop over vertex neighbors of neighbors
      case COLORING_D2_VB_BIT_EF: {
        functorGreedyColorVB_BIT_EF gc(this->nr, this->nc, xadj_, adj_copy_, t_xadj_, t_adj_copy_, vertex_colors_);
        Kokkos::parallel_for("LoopOverChunks", range_policy_type(0, this->nr), gc);
        // prettyPrint1DView(vertex_colors_, "COLORS_GC_VB_BIT",500);
      } break;

      default:
        throw std::invalid_argument(
            "Unknown Distance-2 Algorithm Type or algorithm does not use Edge "
            "Filtering.");
    }
  }  // colorGreedyEF (end)

  // -----------------------------------------------------------------
  //
  // GraphColorDistance2::findConflicts()
  //
  // -----------------------------------------------------------------
  lno_t findConflicts(bool& swap_work_arrays, rowmap_t xadj_, entries_t adj_, rowmap_t t_xadj_, entries_t t_adj_,
                      color_view_type vertex_colors_, lno_view_t current_vertexList_, lno_t current_vertexListLength_,
                      lno_view_t next_iteration_recolorList_, single_lno_view_t next_iteration_recolorListLength_) {
    swap_work_arrays          = true;
    lno_t output_numUncolored = 0;

    functorFindConflicts_Atomic conf(this->nr, this->nc, xadj_, adj_, t_xadj_, t_adj_, vertex_colors_,
                                     current_vertexList_, next_iteration_recolorList_,
                                     next_iteration_recolorListLength_);
    Kokkos::parallel_reduce("FindConflicts", range_policy_type(0, current_vertexListLength_), conf,
                            output_numUncolored);
    return output_numUncolored;
  }  // findConflicts (end)

  // -----------------------------------------------------------------
  //
  // GraphColorDistance2::resolveConflictsSerial()
  //
  // -----------------------------------------------------------------
  void resolveConflictsSerial(rowmap_t xadj_, entries_t adj_, rowmap_t t_xadj_, entries_t t_adj_,
                              color_view_type vertex_colors_, lno_view_t current_vertexList_,
                              size_type current_vertexListLength_) {
    color_type* forbidden = new color_type[nr];
    for (lno_t i = 0; i < nr; i++) forbidden[i] = nr;
    lno_t vid = 0;
    lno_t end = nr;

    typename lno_view_t::HostMirror h_recolor_list;

    end            = current_vertexListLength_;
    h_recolor_list = Kokkos::create_mirror_view(current_vertexList_);
    Kokkos::deep_copy(h_recolor_list, current_vertexList_);

    auto h_colors = Kokkos::create_mirror_view(vertex_colors_);

    auto h_idx = Kokkos::create_mirror_view(xadj_);
    auto h_adj = Kokkos::create_mirror_view(adj_);

    auto h_t_idx = Kokkos::create_mirror_view(t_xadj_);
    auto h_t_adj = Kokkos::create_mirror_view(t_adj_);

    Kokkos::deep_copy(h_colors, vertex_colors_);

    Kokkos::deep_copy(h_idx, xadj_);
    Kokkos::deep_copy(h_adj, adj_);

    Kokkos::deep_copy(h_t_idx, t_xadj_);
    Kokkos::deep_copy(h_t_adj, t_adj_);

    for (lno_t k = 0; k < end; k++) {
      vid = h_recolor_list(k);

      if (h_colors(vid) > 0) continue;

      // loop over distance-1 neighbors of vid
      for (size_type vid_d1_adj = h_idx(vid); vid_d1_adj < h_idx(vid + 1); vid_d1_adj++) {
        lno_t vid_d1 = h_adj(vid_d1_adj);
        if (vid_d1 < nc) {
          if (!doing_bipartite && vid_d1 != vid) {
            forbidden[h_colors(vid_d1)] = vid;
          }
          // loop over neighbors of vid_d1 (distance-2 from vid)
          for (size_type vid_d2_adj = h_t_idx(vid_d1); vid_d2_adj < h_t_idx(vid_d1 + 1); vid_d2_adj++) {
            lno_t vid_d2 = h_t_adj(vid_d2_adj);

            // skip over loops vid -- x -- vid, and filter out-of-bounds
            if (vid_d2 != vid && vid_d2 < nr) forbidden[h_colors(vid_d2)] = vid;
          }
        }
      }

      // color vertex vid with smallest available color
      int c = 1;
      while (forbidden[c] == vid) c++;
      h_colors(vid) = c;
    }
    Kokkos::deep_copy(vertex_colors_, h_colors);
    delete[] forbidden;
  }  // resolveConflictsSerial (end)

  // ------------------------------------------------------
  // Functions: Helpers
  // ------------------------------------------------------

 public:
  // pretty-print a 1D View with label
  template <typename kokkos_view_t>
  void prettyPrint1DView(kokkos_view_t& view, const char* label, const size_t max_entries = 500) const {
    int max_per_line = 20;
    int line_count   = 1;
    std::cout << label << " = [ \n\t";
    for (size_t i = 0; i < view.extent(0); i++) {
      std::cout << std::setw(5) << view(i) << " ";
      if (line_count >= max_per_line) {
        std::cout << std::endl << "\t";
        line_count = 0;
      }
      line_count++;
      if (i >= max_entries - 1) {
        std::cout << "<snip>";
        break;
      }
    }
    if (line_count > 1) std::cout << std::endl;
    std::cout << "\t ]" << std::endl;
  }  // prettyPrint1DView (end)

  // ------------------------------------------------------
  // Functors: Distance-2 Graph Coloring
  // ------------------------------------------------------

  /**
   * Functor to init a list sequentialy, that is list[i] = i
   */
  template <typename view_type>
  struct functorInitList {
    view_type _vertexList;
    functorInitList(view_type vertexList) : _vertexList(vertexList) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const lno_t i) const {
      // Natural order
      _vertexList(i) = i;
    }
  };  // struct functorInitList (end)

  /**
   * Functor for VB algorithm speculative coloring without edge filtering.
   * Single level parallelism
   */
  struct functorGreedyColorVB {
    lno_t nr;                 // num vertices
    lno_t nc;                 // num columns
    rowmap_t _idx;            // vertex degree list
    entries_t _adj;           // vertex adjacency list
    rowmap_t _t_idx;          // transpose vertex degree list
    entries_t _t_adj;         // transpose vertex adjacency list
    color_view_type _colors;  // vertex colors
    lno_view_t _vertexList;   //
    lno_t _vertexListLength;  //
    lno_t _chunkSize;         //

    functorGreedyColorVB(lno_t nr_, lno_t nc_, rowmap_t xadj_, entries_t adj_, rowmap_t t_xadj_, entries_t t_adj_,
                         color_view_type colors, lno_view_t vertexList, lno_t vertexListLength)
        : nr(nr_),
          nc(nc_),
          _idx(xadj_),
          _adj(adj_),
          _t_idx(t_xadj_),
          _t_adj(t_adj_),
          _colors(colors),
          _vertexList(vertexList),
          _vertexListLength(vertexListLength) {}

    // Color vertex i with smallest available color.
    //
    // Each thread colors a chunk of vertices to prevent all vertices getting
    // the same color.
    //
    // This version uses a bool array of size FORBIDDEN_SIZE.
    //
    // param: ii = vertex id
    //
    KOKKOS_INLINE_FUNCTION
    void operator()(const lno_t vid) const {
      // If vertex is not already colored...
      if (_colors(vid) <= 0) {
        const size_type vid_adj_begin = _idx(vid);
        const size_type vid_adj_end   = _idx(vid + 1);

        // Use forbidden array to find available color.
        // - should be small enough to fit into fast memory (use Kokkos
        // memoryspace?)
        bool forbidden[VB_D2_COLORING_FORBIDDEN_SIZE];  // Forbidden Colors

        // Do multiple passes if the forbidden array is too small.
        // * The Distance-1 code used the knowledge of the degree of the vertex
        // to cap the number of iterations
        //   but in distance-2 we'd need the total vertices at distance-2 which
        //   we don't easily have aprioi. This could be as big as all the
        //   vertices in the graph if diameter(G)=2...
        for (color_type offset = 1; offset <= nr; offset += VB_D2_COLORING_FORBIDDEN_SIZE) {
          // initialize
          for (int i = 0; i < VB_D2_COLORING_FORBIDDEN_SIZE; i++) {
            forbidden[i] = false;
          }
          // Check neighbors, fill forbidden array.
          for (size_type vid_adj = vid_adj_begin; vid_adj < vid_adj_end; vid_adj++) {
            const lno_t vid_d1 = _adj(vid_adj);
            if (vid_d1 < nc) {
              if (!doing_bipartite)  // note: compile-time branch (template
                                     // param)
              {
                if (vid_d1 != vid) {
                  const color_type c = _colors(vid_d1);
                  if ((c >= offset) && (c - offset < VB_D2_COLORING_FORBIDDEN_SIZE)) {
                    forbidden[c - offset] = true;
                  }
                }
              }
              const size_type vid_d1_adj_begin = _t_idx(vid_d1);
              const size_type vid_d1_adj_end   = _t_idx(vid_d1 + 1);
              for (size_type vid_d1_adj = vid_d1_adj_begin; vid_d1_adj < vid_d1_adj_end; vid_d1_adj++) {
                const lno_t vid_d2 = _t_adj(vid_d1_adj);

                // Skip distance-2-self-loops
                if (vid_d2 != vid && vid_d2 < nr) {
                  const color_type c = _colors(vid_d2);
                  if ((c >= offset) && (c - offset < VB_D2_COLORING_FORBIDDEN_SIZE)) {
                    forbidden[c - offset] = true;
                  }
                }
              }  // for vid_d1_adj...
            }
          }  // for vid_adj...

          // color vertex i with smallest available color (firstFit)
          for (int c = 0; c < VB_D2_COLORING_FORBIDDEN_SIZE; c++) {
            if (!forbidden[c]) {
              _colors(vid) = offset + c;
              return;
            }
          }  // for c...
        }    // for offset < nr
      }      // if _colors(vid) <= 0...
    }        // operator() (end)
  };         // struct functorGreedyColorVB (end)

  /**
   * Functor for VB_BIT algorithm coloring without edge filtering.
   * Single level parallelism
   */
  struct functorGreedyColorVB_BIT {
    lno_t nr;                 // num vertices
    lno_t nc;                 // num columns
    rowmap_t _idx;            // vertex degree list
    entries_t _adj;           // vertex adjacency list
    rowmap_t _t_idx;          // transpose vertex degree list
    entries_t _t_adj;         // transpose vertex adjacency list
    color_view_type _colors;  // vertex colors
    lno_view_t _vertexList;   //
    lno_t _vertexListLength;  //

    functorGreedyColorVB_BIT(lno_t nr_, lno_t nc_, rowmap_t xadj_, entries_t adj_, rowmap_t t_xadj_, entries_t t_adj_,
                             color_view_type colors, lno_view_t vertexList, lno_t vertexListLength)
        : nr(nr_),
          nc(nc_),
          _idx(xadj_),
          _adj(adj_),
          _t_idx(t_xadj_),
          _t_adj(t_adj_),
          _colors(colors),
          _vertexList(vertexList),
          _vertexListLength(vertexListLength) {}

    // Color vertex i with smallest available color.
    //
    // Each thread colors a chunk of vertices to prevent all vertices getting
    // the same color.
    //
    // This version uses a bool array of size FORBIDDEN_SIZE.
    //
    // param: ii = vertex id
    //
    KOKKOS_INLINE_FUNCTION
    void operator()(const lno_t vid) const {
      // If vertex is not colored yet...
      if (_colors(vid) == 0) {
        const size_type vid_adj_begin = _idx(vid);
        const size_type vid_adj_end   = _idx(vid + 1);

        for (color_type offset = 1; offset <= (nr + VBBIT_D2_COLORING_FORBIDDEN_SIZE);
             offset += VBBIT_D2_COLORING_FORBIDDEN_SIZE) {
          // Forbidden colors
          // - single long int for forbidden colors
          bit_64_forbidden_type forbidden = 0;

          // If all available colors for this range are unavailable we can break
          // out of the nested loops
          bool break_out = false;

          // Loop over distance-1 neighbors of vid
          for (size_type vid_adj = vid_adj_begin; !break_out && vid_adj < vid_adj_end; ++vid_adj) {
            const lno_t vid_d1 = _adj(vid_adj);
            if (vid_d1 < nc) {
              if (!doing_bipartite)  // note: compile-time branch (template
                                     // param)
              {
                // Check for dist-1 conflicts
                if (vid_d1 != vid) {
                  const color_type color        = _colors(vid_d1);
                  const color_type color_offset = color - offset;
                  if (color && color_offset <= VBBIT_D2_COLORING_FORBIDDEN_SIZE) {
                    // if it is in the current range, then add the color to the
                    // banned colors
                    if (color > offset) {
                      // convert color to bit representation
                      bit_64_forbidden_type ban_color_bit = 1;

                      ban_color_bit = ban_color_bit << color_offset;

                      // add it to forbidden colors
                      forbidden |= (bit_64_forbidden_type(1) << color_offset);
                    }
                  }
                }
              }
              const size_type vid_d1_adj_begin = _t_idx(vid_d1);
              const size_type vid_d1_adj_end   = _t_idx(vid_d1 + 1);

              // Loop over distance-2 neighbors of vid
              for (size_type vid_d1_adj = vid_d1_adj_begin; !break_out && vid_d1_adj < vid_d1_adj_end; ++vid_d1_adj) {
                const lno_t vid_d2 = _t_adj(vid_d1_adj);

                // Ignore Distance-2 Self Loops
                if (vid_d2 != vid && vid_d2 < nr) {
                  const color_type color        = _colors(vid_d2);
                  const color_type color_offset = color - offset;

                  // if color is within the current range, or if its color is in
                  // a previously traversed range
                  if (offset <= color && color_offset < VBBIT_D2_COLORING_FORBIDDEN_SIZE) {
                    // if it is in the current range, then add the color to the
                    // banned colors
                    forbidden |= (bit_64_forbidden_type(1) << color_offset);

                    // if there are no available colors in this range then exit
                    // early, no need to traverse the rest.
                    if (0 == ~forbidden) {
                      break_out = true;
                    }
                  }  // if color in current range ...
                }    // if vid_d2 ...
              }      // for vid_d1_adj ...
            }
          }  // for vid_adj ...

          // check if an available color exists.
          if (~forbidden) {
            bit_64_forbidden_type color_offset = KokkosKernels::Impl::least_set_bit(~forbidden) - 1;
            _colors(vid)                       = offset + color_offset;
            return;
          }
        }  // for offset <= (nr + VBBIT_D2_COLORING_FORBIDDEN_SIZE)
      }    // if _colors(vid)==0
    }      // operator() (end)
  };       // struct functorGreedyColorVB_BIT (end)

  /**
   * Functor for VB_BIT_EF algorithm coloring without edge filtering.
   * Single level parallelism
   */
  struct functorGreedyColorVB_BIT_EF {
    lno_t _nr;                // num vertices
    lno_t _nc;                // num vertices
    rowmap_t _idx;            // rowmap
    lno_view_t _adj;          // vertex adjacency list  (mutable)
    rowmap_t _t_idx;          // transpose vertex degree list
    entries_t _t_adj;         // transpose vertex adjacency list (NOT modified)
    color_view_type _colors;  // vertex colors

    functorGreedyColorVB_BIT_EF(lno_t nr_, lno_t nc_, rowmap_t xadj_, lno_view_t adj_, rowmap_t t_xadj_,
                                entries_t t_adj_, color_view_type colors)
        : _nr(nr_), _nc(nc_), _idx(xadj_), _adj(adj_), _t_idx(t_xadj_), _t_adj(t_adj_), _colors(colors) {}

    // Color vertex i with smallest available color.
    //
    // Each thread colors a chunk of vertices to prevent all vertices getting
    // the same color.
    //
    // This version uses a bool array of size FORBIDDEN_SIZE.
    //
    // param: ii = vertex id
    //
    KOKKOS_INLINE_FUNCTION
    void operator()(const lno_t vid) const {
      // If vertex is not colored yet..
      if (_colors(vid) == 0) {
        size_type vid_adj_begin = _idx(vid);
        size_type vid_adj_end   = _idx(vid + 1);

        bool foundColor = false;
        for (color_type offset = 0; !foundColor && offset <= (_nr + VBBIT_D2_COLORING_FORBIDDEN_SIZE);
             offset += VBBIT_D2_COLORING_FORBIDDEN_SIZE) {
          // Forbidden colors
          // - single long int for forbidden colors
          bit_64_forbidden_type forbidden = 0;

          // If all available colors for this range are unavailable we can break
          // out of the nested loops
          bool offset_colors_full = false;

          // Loop over distance-1 neighbors of vid
          for (size_type vid_adj = vid_adj_begin; !offset_colors_full && vid_adj < vid_adj_end; ++vid_adj) {
            const lno_t vid_d1 = _adj(vid_adj);
            if (vid_d1 < _nc) {
              if (!doing_bipartite)  // note: compile-time branch (template
                                     // param)
              {
                // Check for dist-1 conflict
                if (vid_d1 != vid) {
                  color_type color        = _colors(vid_d1);
                  color_type color_offset = color - offset;
                  // if color is within the current range, or if its color is in
                  // a previously traversed range
                  if (color && offset < color && color_offset <= VBBIT_D2_COLORING_FORBIDDEN_SIZE) {
                    // if it is in the current range, then add the color to the
                    // banned colors convert color to bit representation
                    bit_64_forbidden_type ban_color_bit = 1;
                    ban_color_bit                       = ban_color_bit << (color_offset - 1);
                    // add it to forbidden colors
                    forbidden = forbidden | ban_color_bit;
                  }
                }
              }

              size_type vid_d1_adj_begin            = _t_idx(vid_d1);
              const size_type vid_d1_adj_end        = _t_idx(vid_d1 + 1);
              const size_type degree_vid_d1         = vid_d1_adj_end - vid_d1_adj_begin;
              size_type num_vid_d2_colored_in_range = 0;

              // Store the maximum color value found in the vertices adjacent to
              // vid_d1
              color_type max_color_adj_to_d1 = 0;

              // Loop over distance-2 neighbors of vid
              for (size_type vid_d1_adj = vid_d1_adj_begin; !offset_colors_full && vid_d1_adj < vid_d1_adj_end;
                   ++vid_d1_adj) {
                const lno_t vid_d2 = _t_adj(vid_d1_adj);

                // Ignore Distance-2 Self Loops
                if (vid_d2 != vid && vid_d2 < _nr) {
                  color_type color        = _colors(vid_d2);
                  color_type color_offset = color - offset;  // color_offset < 0 means color is from a
                                                             // previous offset.

                  // Update maximum color adjacent to vid_d1 found so far.
                  max_color_adj_to_d1 = color > max_color_adj_to_d1 ? color : max_color_adj_to_d1;

                  // if color is within the current range, or if its color is in
                  // a previously traversed range
                  if (color && color_offset <= VBBIT_D2_COLORING_FORBIDDEN_SIZE) {
                    num_vid_d2_colored_in_range++;

                    // if it is in the current range, then add the color to the
                    // banned colors
                    if (color > offset) {
                      // convert color to bit representation
                      bit_64_forbidden_type ban_color_bit = 1;

                      ban_color_bit = ban_color_bit << (color_offset - 1);

                      // add it to forbidden colors
                      forbidden = forbidden | ban_color_bit;

                      // if there are no available colors in this range then
                      // exit early, no need to traverse the rest b/c they
                      // contribute no new information at this offset.
                      if (0 == ~forbidden) {
                        offset_colors_full = true;
                        // Note: with edge-filtering, this can short-circuit the
                        // loop over all
                        //       neighbors of VID and will reduce the number of
                        //       filtered edges.
                      }
                    }  // if color > offset
                  }    // if color && color_offset
                }      // if vid_d2 != vid ...
                else {
                  // If there's a self-loop then we should increment our
                  // 'colored in range' so we don't block filtering since we
                  // know there must be a (v2,v1) edge
                  num_vid_d2_colored_in_range++;
                }
              }  // for vid_d1_adj ...

              // Edge filtering on the neighbors of vid.  We can only do this if
              // ALL neighbors of vid_d1 have been visited and if all were
              // colored in current offset range or lower.
              if (degree_vid_d1 == num_vid_d2_colored_in_range) {
                if (vid_adj_begin > vid_adj) {
                  _adj(vid_adj)       = _adj(vid_adj_begin);
                  _adj(vid_adj_begin) = vid_d1;
                }
                vid_adj_begin++;
              }
            }

          }  // for vid_adj
          forbidden = ~(forbidden);

          // check if an available color exists.
          if (forbidden) {
            // if there is an available color, choose the first color, using 2s
            // complement.
            bit_64_forbidden_type new_color = forbidden & (-forbidden);
            color_type val                  = 1;

            // convert it back to decimal color.
            while ((new_color & 1) == 0) {
              ++val;
              new_color = new_color >> 1;
            }
            _colors(vid) = val + offset;
            foundColor   = true;
            break;
          }
        }  // for offset=0...
      }    // if _colors(vid)==0
    }      // operator() (end)
  };       // struct functorGreedyColorVB_BIT_EF (end)

  struct functorFindConflicts_Atomic {
    lno_t nr;  // num verts
    lno_t nc;  // num columns
    rowmap_t _idx;
    entries_t _adj;
    rowmap_t _t_idx;
    entries_t _t_adj;
    color_view_type _colors;
    lno_view_t _vertexList;
    lno_view_t _recolorList;
    single_lno_view_t _recolorListLength;

    functorFindConflicts_Atomic(lno_t nr_, lno_t nc_, rowmap_t xadj_, entries_t adj_, rowmap_t t_xadj_,
                                entries_t t_adj_, color_view_type colors, lno_view_t vertexList, lno_view_t recolorList,
                                single_lno_view_t recolorListLength)
        : nr(nr_),
          nc(nc_),
          _idx(xadj_),
          _adj(adj_),
          _t_idx(t_xadj_),
          _t_adj(t_adj_),
          _colors(colors),
          _vertexList(vertexList),
          _recolorList(recolorList),
          _recolorListLength(recolorListLength) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const lno_t i, lno_t& numConflicts) const {
      const lno_t vid                  = _vertexList(i);
      const color_type my_color        = _colors(vid);
      const size_type vid_d1_adj_begin = _idx(vid);
      const size_type vid_d1_adj_end   = _idx(vid + 1);
      // If vid is a valid column (vid < nc), check for column->vert conflicts
      for (size_type vid_d1_adj = vid_d1_adj_begin; vid_d1_adj < vid_d1_adj_end; vid_d1_adj++) {
        lno_t vid_d1 = _adj(vid_d1_adj);
        if (vid_d1 < nc) {
          if (!doing_bipartite)  // note: compile-time branch (template param)
          {
            // check for dist-1 conflict
            if (vid_d1 != vid && _colors(vid_d1) == my_color) {
              _colors(vid) = 0;  // uncolor vertex
              // Atomically add vertex to recolorList
              const lno_t k   = Kokkos::atomic_fetch_add(&_recolorListLength(), lno_t(1));
              _recolorList(k) = vid;
              numConflicts++;
              return;
            }
          }
          const size_type d2_adj_begin = _t_idx(vid_d1);
          const size_type d2_adj_end   = _t_idx(vid_d1 + 1);
          for (size_type vid_d2_adj = d2_adj_begin; vid_d2_adj < d2_adj_end; vid_d2_adj++) {
            const lno_t vid_d2 = _t_adj(vid_d2_adj);

            if (vid != vid_d2 && vid_d2 < nr) {
              if (_colors(vid_d2) == my_color) {
                _colors(vid) = 0;  // uncolor vertex
                // Atomically add vertex to recolorList
                const lno_t k   = Kokkos::atomic_fetch_add(&_recolorListLength(), lno_t(1));
                _recolorList(k) = vid;
                numConflicts++;
                return;
              }
            }  // if vid != vid_d2 ...
          }    // for vid_d2_adj ...
        }
      }  // for vid_d1_adj ...
    }    // operator() (end)
  };     // struct functorFindConflicts_Atomic (end)
};       // end class GraphColorDistance2

/**
 * Prints out a histogram of graph colors for Distance-2 Graph Coloring
 *
 * If the graph is symmetric, give the same value for col_map and row_map,
 * and for row_entries and col_entries.
 *
 * @param[in]  handle           The kernel handle
 * @param[in]  num_rows         Number of rows in the matrix (number of
 * vertices)
 * @param[in]  num_cols         Number of columns in the matrix
 * @param[in]  row_map          The row map
 * @param[in]  row_entries      The row entries
 * @param[in]  col_map          The column map
 * @param[in]  col_entries      The column entries
 * @param[out] validation_flags An array of 4 booleans.
 *                              validation_flags[0] : True IF the distance-2
 * coloring is invalid. validation_flags[1] : True IF the coloring is bad
 * because vertices are left uncolored. validation_flags[2] : True IF the
 * coloring is bad because at least one pair of vertices at distance=2 from each
 * other has the same color. validation_flags[3] : True IF a vertex has a color
 * greater than number of vertices in the graph. May not be an INVALID coloring,
 * but can indicate poor quality in coloring.
 * @param[in] csv               Output in CSV format? Default: false
 *
 * @return nothing
 */
template <class KernelHandle>
void graph_print_distance2_color_histogram(KernelHandle* handle, bool csv = false) {
  using lno_view_t      = typename KernelHandle::nnz_lno_temp_work_view_t;
  using lno_t           = typename KernelHandle::nnz_lno_t;
  using execution_space = typename KernelHandle::HandleExecSpace;
  using D2Handle        = typename KernelHandle::GraphColorDistance2HandleType;
  using color_view_t    = typename D2Handle::color_view_type;
  // Get handle
  D2Handle* gch_d2 = handle->get_distance2_graph_coloring_handle();
  // Get the coloring
  color_view_t colors = gch_d2->get_vertex_colors();
  lno_t num_colors    = gch_d2->get_num_colors();
  lno_view_t histogram("histogram", num_colors + 1);
  KokkosKernels::Impl::kk_get_histogram<color_view_t, lno_view_t, execution_space>(colors.extent(0), colors, histogram);
  auto h_histogram = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), histogram);
  // note: both modes ignore color 0 in output, since we assume the coloring is
  // valid
  if (csv) {
    size_t i = 1;
    for (; i < h_histogram.extent(0) - 1; i++) {
      std::cout << h_histogram(i) << ",";
    }
    std::cout << h_histogram(i);
  } else {
    auto histogram_slice = Kokkos::subview(histogram, std::make_pair((size_t)1, histogram.extent(0)));
    std::cout << "Distance-2 Color Histogram (1..N): " << std::endl;
    KokkosKernels::Impl::kk_print_1Dview(histogram_slice);
    std::cout << std::endl;
  }
}

}  // namespace Impl
}  // namespace KokkosGraph

#endif  // _KOKKOSGRAPH_DISTANCE2COLOR_IMPL_HPP
