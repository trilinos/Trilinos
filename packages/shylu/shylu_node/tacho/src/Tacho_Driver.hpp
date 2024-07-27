// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_DRIVER_HPP__
#define __TACHO_DRIVER_HPP__

/// \file Tacho_Driver.hpp
/// \brief temporary solver interface for refactoring
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>

namespace Tacho {

/// forward decl
class Graph;
#if defined(TACHO_HAVE_METIS)
class GraphTools_Metis;
#else
class GraphTools;
#endif

class SymbolicTools;
template <typename ValueType, typename DeviceType> class CrsMatrixBase;
template <typename ValueType, typename DeviceType> class NumericToolsBase;
template <typename ValueType, typename DeviceType> class NumericToolsSerial;
template <typename ValueType, typename DeviceType, int Var> class NumericToolsLevelSet;

///
/// Tacho Solver interface
///
template <typename ValueType, typename DeviceType> struct Driver {
public:
  using value_type = ValueType;
  using device_type = DeviceType;
  using exec_space = typename device_type::execution_space;
  using exec_memory_space = typename device_type::memory_space;

  using host_device_type = typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type;
  using host_space = typename host_device_type::execution_space;
  using host_memory_space = typename host_device_type::memory_space;

  using size_type_array = Kokkos::View<size_type *, device_type>;
  using ordinal_type_array = Kokkos::View<ordinal_type *, device_type>;
  using value_type_array = Kokkos::View<value_type *, device_type>;
  using value_type_matrix = Kokkos::View<value_type **, Kokkos::LayoutLeft, device_type>;

  using size_type_array_host = Kokkos::View<size_type *, host_device_type>;
  using ordinal_type_array_host = Kokkos::View<ordinal_type *, host_device_type>;
  using value_type_array_host = Kokkos::View<value_type *, host_device_type>;
  using value_type_matrix_host = Kokkos::View<value_type **, Kokkos::LayoutLeft, host_device_type>;

  using crs_matrix_type = CrsMatrixBase<value_type, device_type>;
  using crs_matrix_type_host = CrsMatrixBase<value_type, host_device_type>;

#if defined(TACHO_HAVE_METIS)
  using graph_tools_type = GraphTools_Metis;
#else
  using graph_tools_type = GraphTools;
#endif

  using symbolic_tools_type = SymbolicTools;
  using numeric_tools_base_type = NumericToolsBase<value_type, device_type>;
  using numeric_tools_serial_type = NumericToolsSerial<value_type, device_type>;
  using numeric_tools_levelset_var0_type = NumericToolsLevelSet<value_type, device_type, 0>;
  using numeric_tools_levelset_var1_type = NumericToolsLevelSet<value_type, device_type, 1>;
  using numeric_tools_levelset_var2_type = NumericToolsLevelSet<value_type, device_type, 2>;

private:
  enum : int { Cholesky = 1, LDL = 2, SymLU = 3, LU = 4 };

  // ** solver mode
  ordinal_type _method;

  // ** ordering options
  ordinal_type _order_connected_graph_separately;

  // ** problem
  ordinal_type _m;
  size_type _nnz;

  size_type_array _ap;
  size_type_array_host _h_ap;
  ordinal_type_array _aj;
  ordinal_type_array_host _h_aj;

  ordinal_type_array _perm;
  ordinal_type_array_host _h_perm;
  ordinal_type_array _peri;
  ordinal_type_array_host _h_peri;

  // ** condensed graph
  ordinal_type _m_graph;
  size_type _nnz_graph;

  size_type_array_host _h_ap_graph;
  ordinal_type_array_host _h_aj_graph;
  ordinal_type_array_host _h_aw_graph;

  ordinal_type_array_host _h_perm_graph;
  ordinal_type_array_host _h_peri_graph;

  // ** symbolic factorization output
  // supernodes output
  ordinal_type _nsupernodes;
  ordinal_type_array _supernodes;

  // dof mapping to sparse matrix
  size_type_array _gid_super_panel_ptr;
  ordinal_type_array _gid_super_panel_colidx;

  // supernode map and panel size configuration
  size_type_array _sid_super_panel_ptr;
  ordinal_type_array _sid_super_panel_colidx, _blk_super_panel_colidx;

  // supernode elimination tree (parent - children)
  size_type_array _stree_ptr;
  ordinal_type_array _stree_children;

  // supernode elimination tree (child - parent)
  ordinal_type_array _stree_parent;

  // roots of supernodes
  ordinal_type_array_host _stree_level, _stree_roots;

  // ** numeric factorization output
  numeric_tools_base_type *_N;

  // small dense matrix
  // - chol A is used
  // - ldl A D P are used
  value_type_matrix_host _A, _D;
  ordinal_type_array_host _P;

  // ** options
  ordinal_type _verbose;             // print
  ordinal_type _small_problem_thres; // smaller than this, use lapack

  // // ** tasking options
  ordinal_type _serial_thres_size; // serialization threshold size
  ordinal_type _mb;                // block size for byblocks algorithms
  ordinal_type _nb;                // panel size for panel algorithms
  ordinal_type _front_update_mode; // front update mode 0 - lock, 1 - atomic

  // ** levelset options
  bool _levelset;                    // use level set code instead of tasking
  ordinal_type _device_level_cut;    // above this level, matrices are computed on device
  ordinal_type _device_factor_thres; // bigger than this threshold, device function is used
  ordinal_type _device_solve_thres;  // bigger than this threshold, device function is used
  ordinal_type _variant;             // algorithmic variant in levelset 0: naive, 1: invert diagonals
  ordinal_type _nstreams;            // on cuda, multi streams are used

  // parallelism and memory constraint is made via this parameter
  ordinal_type _max_num_superblocks; // # of superblocks in the memoyrpool

public:
  Driver();
  /// delete copy constructor and assignment operator
  /// sharing numeric tools for different inputs does not make sense
  Driver(const Driver &) = default;
  Driver &operator=(const Driver &) = default;

  /// duplicate the solver with sharing symbolic factorization
  Driver duplicate();

  ///
  /// common options
  ///
  void setVerbose(const ordinal_type verbose = 1);
  void setSmallProblemThresholdsize(const ordinal_type small_problem_thres = 1024);
  void setMatrixType(const int symmetric, // 0 - unsymmetric, 1 - structure sym, 2 - symmetric
                     const bool is_positive_definite);
  void setSolutionMethod(const int method); /// 1 - cholesky, 2 - LDL, 3 - LU

  ///
  /// Graph options
  ///
  void setOrderConnectedGraphSeparately(const ordinal_type order_connected_graph_separately = 1);

  ///
  /// tasking options
  ///
  void setSerialThresholdsize(const ordinal_type serial_thres_size = -1);
  void setBlocksize(const ordinal_type mb = -1);
  void setPanelsize(const ordinal_type nb = -1);
  void setFrontUpdateMode(const ordinal_type front_update_mode = 1);
  void setMaxNumberOfSuperblocks(const ordinal_type max_num_superblocks = -1);

  ///
  /// Level set tools options
  ///
  void setLevelSetScheduling(const bool levelset);
  void setLevelSetOptionDeviceLevelCut(const ordinal_type device_level_cut);
  void setLevelSetOptionDeviceFunctionThreshold(const ordinal_type device_factor_thres,
                                                const ordinal_type device_solve_thres);
  void setLevelSetOptionNumStreams(const ordinal_type nstreams);
  void setLevelSetOptionAlgorithmVariant(const ordinal_type variant);

  ///
  /// get interface
  ///
  ordinal_type getNumSupernodes() const;
  ordinal_type_array getSupernodes() const;
  ordinal_type_array getPermutationVector() const;
  ordinal_type_array getInversePermutationVector() const;

  // internal only
  int analyze();
  int analyze_linear_system();
  int analyze_condensed_graph();

  template <typename arg_size_type_array, typename arg_ordinal_type_array>
  int analyze(const ordinal_type m, const arg_size_type_array &ap, const arg_ordinal_type_array &aj,
              const bool duplicate = false) {
    _m = m;

    if (duplicate) {
      /// for most cases, ap and aj are from host; so construct ap and aj and mirror to device
      _h_ap = size_type_array_host(Kokkos::ViewAllocateWithoutInitializing("h_ap"), ap.extent(0));
      Kokkos::deep_copy(_h_ap, ap);
      _h_aj = ordinal_type_array_host(Kokkos::ViewAllocateWithoutInitializing("h_aj"), aj.extent(0));
      Kokkos::deep_copy(_h_aj, aj);

      _ap = Kokkos::create_mirror_view(exec_memory_space(), _h_ap);
      Kokkos::deep_copy(_ap, _h_ap);
      _aj = Kokkos::create_mirror_view(exec_memory_space(), _h_aj);
      Kokkos::deep_copy(_aj, _h_aj);
    } else {
      /// this does not make any extra deep copy; users should hold the graph data
      _ap = Kokkos::create_mirror_view(exec_memory_space(), ap);
      Kokkos::deep_copy(_ap, ap);
      _aj = Kokkos::create_mirror_view(exec_memory_space(), aj);
      Kokkos::deep_copy(_aj, aj);

      _h_ap = Kokkos::create_mirror_view(host_memory_space(), ap);
      Kokkos::deep_copy(_h_ap, ap);
      _h_aj = Kokkos::create_mirror_view(host_memory_space(), aj);
      Kokkos::deep_copy(_h_aj, aj);
    }

    _h_perm = ordinal_type_array_host();
    _h_peri = ordinal_type_array_host();

    _nnz = _h_ap(m);

    _m_graph = 0;
    _nnz_graph = 0;

    _h_ap_graph = size_type_array_host();
    _h_aj_graph = ordinal_type_array_host();

    _h_perm_graph = ordinal_type_array_host();
    _h_peri_graph = ordinal_type_array_host();

    return analyze();
  }

  template <typename arg_size_type_array, typename arg_ordinal_type_array, typename arg_perm_type_array>
  int analyze(const ordinal_type m, const arg_size_type_array &ap, const arg_ordinal_type_array &aj,
              const arg_perm_type_array &perm, const arg_perm_type_array &peri, const bool duplicate = false) {
    _m = m;

    if (duplicate) {
      /// for most cases, ap and aj are from host; so construct ap and aj and mirror to device
      _h_ap = size_type_array_host(Kokkos::ViewAllocateWithoutInitializing("h_ap"), ap.extent(0));
      Kokkos::deep_copy(_h_ap, ap);
      _h_aj = ordinal_type_array_host(Kokkos::ViewAllocateWithoutInitializing("h_aj"), aj.extent(0));
      Kokkos::deep_copy(_h_aj, aj);

      _ap = Kokkos::create_mirror_view(exec_memory_space(), _h_ap);
      Kokkos::deep_copy(_ap, _h_ap);
      _aj = Kokkos::create_mirror_view(exec_memory_space(), _h_aj);
      Kokkos::deep_copy(_aj, _h_aj);

      _h_perm = ordinal_type_array_host(Kokkos::ViewAllocateWithoutInitializing("h_perm"), perm.extent(0));
      _h_peri = ordinal_type_array_host(Kokkos::ViewAllocateWithoutInitializing("h_peri"), peri.extent(0));
    } else {
      /// this does not make any extra deep copy; users should hold the graph data
      _ap = Kokkos::create_mirror_view(exec_memory_space(), ap);
      Kokkos::deep_copy(_ap, ap);
      _aj = Kokkos::create_mirror_view(exec_memory_space(), aj);
      Kokkos::deep_copy(_aj, aj);

      _h_ap = Kokkos::create_mirror_view(host_memory_space(), ap);
      Kokkos::deep_copy(_h_ap, ap);
      _h_aj = Kokkos::create_mirror_view(host_memory_space(), aj);
      Kokkos::deep_copy(_h_aj, aj);

      _h_perm = Kokkos::create_mirror_view(host_memory_space(), perm);
      _h_peri = Kokkos::create_mirror_view(host_memory_space(), peri);
    }

    Kokkos::deep_copy(_h_perm, perm);
    Kokkos::deep_copy(_h_peri, peri);

    _nnz = _h_ap(m);

    _m_graph = 0;
    _nnz_graph = 0;

    _h_ap_graph = size_type_array_host();
    _h_aj_graph = ordinal_type_array_host();

    _h_perm_graph = ordinal_type_array_host();
    _h_peri_graph = ordinal_type_array_host();

    return analyze();
  }

  template <typename arg_size_type_array, typename arg_ordinal_type_array>
  int analyze(const ordinal_type m, const arg_size_type_array &ap, const arg_ordinal_type_array &aj,
              const ordinal_type m_graph, const arg_size_type_array &ap_graph, const arg_ordinal_type_array &aj_graph,
              const arg_ordinal_type_array &aw_graph, const bool duplicate = false) {
    _m = m;

    if (duplicate) {
      /// for most cases, ap and aj are from host; so construct ap and aj and mirror to device
      _h_ap = size_type_array_host(Kokkos::ViewAllocateWithoutInitializing("h_ap"), ap.extent(0));
      Kokkos::deep_copy(_h_ap, ap);
      _h_aj = ordinal_type_array_host(Kokkos::ViewAllocateWithoutInitializing("h_aj"), aj.extent(0));
      Kokkos::deep_copy(_h_aj, aj);

      _ap = Kokkos::create_mirror_view(exec_memory_space(), _h_ap);
      Kokkos::deep_copy(_ap, _h_ap);
      _aj = Kokkos::create_mirror_view(exec_memory_space(), _h_aj);
      Kokkos::deep_copy(_aj, _h_aj);
    } else {
      /// this does not make any extra deep copy; users should hold the graph data
      _ap = Kokkos::create_mirror_view(exec_memory_space(), ap);
      Kokkos::deep_copy(_ap, ap);
      _aj = Kokkos::create_mirror_view(exec_memory_space(), aj);
      Kokkos::deep_copy(_aj, aj);

      _h_ap = Kokkos::create_mirror_view(host_memory_space(), ap);
      Kokkos::deep_copy(_h_ap, ap);
      _h_aj = Kokkos::create_mirror_view(host_memory_space(), aj);
      Kokkos::deep_copy(_h_aj, aj);
    }

    _h_perm = ordinal_type_array_host();
    _h_peri = ordinal_type_array_host();

    _nnz = _h_ap(m);

    _m_graph = m_graph;
    if (duplicate) {
      _h_ap_graph = size_type_array_host(Kokkos::ViewAllocateWithoutInitializing("h_ap_graph"), ap_graph.extent(0));
      _h_aj_graph = ordinal_type_array_host(Kokkos::ViewAllocateWithoutInitializing("h_aj_graph"), aj_graph.extent(0));
      _h_aw_graph = ordinal_type_array_host(Kokkos::ViewAllocateWithoutInitializing("h_aw_graph"), aw_graph.extent(0));
    } else {
      _h_ap_graph = Kokkos::create_mirror_view(host_memory_space(), ap_graph);
      _h_aj_graph = Kokkos::create_mirror_view(host_memory_space(), aj_graph);
      _h_aw_graph = Kokkos::create_mirror_view(host_memory_space(), aw_graph);
    }

    Kokkos::deep_copy(_h_ap_graph, ap_graph);
    Kokkos::deep_copy(_h_aj_graph, aj_graph);
    Kokkos::deep_copy(_h_aw_graph, aw_graph);

    _h_perm_graph = ordinal_type_array_host();
    _h_peri_graph = ordinal_type_array_host();

    _nnz_graph = _h_ap_graph(m_graph);

    return analyze();
  }

  int initialize();

  int factorize(const value_type_array &ax);
  int factorize_small_host(const value_type_array &ax);

  int solve(const value_type_matrix &x, const value_type_matrix &b, const value_type_matrix &t);
  int solve_small_host(const value_type_matrix &x, const value_type_matrix &b, const value_type_matrix &t);

  double computeRelativeResidual(const value_type_array &ax, const value_type_matrix &x, const value_type_matrix &b);
  void   computeSpMV(const value_type_array &ax, const value_type_matrix &x, value_type_matrix &b);

  int exportFactorsToCrsMatrix(crs_matrix_type &A);
  int release();
};

} // namespace Tacho

//#include "Tacho_Driver_Impl.hpp"

#endif
