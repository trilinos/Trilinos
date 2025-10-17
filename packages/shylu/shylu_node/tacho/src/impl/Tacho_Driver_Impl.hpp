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
#ifndef __TACHO_DRIVER_IMPL_HPP__
#define __TACHO_DRIVER_IMPL_HPP__

/// \file Tacho_Driver_Impl.hpp
/// \brief temporary solver interface for refactoring
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Driver.hpp"
#include "Tacho_Internal.hpp"

namespace Tacho {

template <typename VT, typename DT>
Driver<VT, DT>::Driver()
    : _method_setup(1), _method(1), _order_connected_graph_separately(true), _graph_algo_type(-1), _m(0), _nnz(0),
      _ap(), _h_ap(), _aj(), _h_aj(), _perm(), _h_perm(), _peri(), _h_peri(),
      _m_graph(0), _nnz_graph(0), _h_ap_graph(), _h_aj_graph(), _h_perm_graph(),
      _h_peri_graph(), _nnz_u(0), _nsupernodes(0), _N(nullptr), _verbose(0), _small_problem_thres(1024),
#ifdef TACHO_DEPRECATED_PARAMETERS
      _serial_thres_size(-1), _mb(-1), _nb(-1), _front_update_mode(-1), _levelset(0),
#endif
      _device_level_cut(0), _device_factor_thres(128), _device_solve_thres(128),
      #if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL)
      _variant(2),
      #else
      _variant(-1), // sequential by default
      #endif
      _nstreams(16), _pivot_tol(0.0),
#if defined(KOKKOS_ENABLE_HIP)
      _store_transpose(true)
#else
      _store_transpose(false)
#endif
#ifdef TACHO_DEPRECATED_PARAMETERS
      , _max_num_superblocks(-1)
#endif
     {}

///
/// duplicate the object
///
template <typename VT, typename DT> Driver<VT, DT> Driver<VT, DT>::duplicate() {
  /// input matrix should be given (m and nnz) and analysis is done (nsupernodes is greater than zero)
  const bool is_analysis_done = (_m > 0) && (_nnz > 0) && (_nsupernodes > 0);
  TACHO_TEST_FOR_EXCEPTION(!is_analysis_done, std::logic_error, "Analysis is not done yet");

  /// copy constructor of this
  Driver<VT, DT> r_val(*this);

  /// make sure numeric tool is null pointer
  r_val._N = nullptr;

  return r_val;
}

///
/// common options
///
template <typename VT, typename DT> void Driver<VT, DT>::setVerbose(const ordinal_type verbose) { _verbose = verbose; }

template <typename VT, typename DT>
void Driver<VT, DT>::setSmallProblemThresholdsize(const ordinal_type small_problem_thres) {
  _small_problem_thres = small_problem_thres;
}

template <typename VT, typename DT>
void Driver<VT, DT>::setMatrixType(const int symmetric, // 0 - unsymmetric, 1 - structure sym, 2 - symmetric
                                   const bool is_positive_definite) {
  switch (symmetric) {
  case 0: {
    _method_setup = LU;
    _method = LU;
    break;
  }
  case 1: {
    _method_setup = SymLU;
    _method = SymLU;
    break;
  }
  case 2: {
    if (is_positive_definite) {
      if (std::is_same<value_type, double>::value || std::is_same<value_type, float>::value ||
          std::is_same<value_type, Kokkos::complex<float>>::value ||
          std::is_same<value_type, Kokkos::complex<double>>::value) {
        // real symmetric posdef
        _method_setup = Cholesky;
        _method = Cholesky;
      }
    } else { // real or complex symmetric indef
      _method_setup = LDL;
      _method = LDL;
    }
    break;
  }
  default: {
    TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "symmetric argument is wrong");
  }
  }
}

template <typename VT, typename DT>
void Driver<VT, DT>::setSolutionMethod(const int method) { // 0 - LDL nopivot, 1 - Chol, 2 - LDL, 3 - LU
  {
    std::stringstream ss;
    ss << "Error: the given method (" << method << ") is not supported, 0 - LDL nopivot, 1 - Chol, 2 - LDL, 3 - SymLU";
    TACHO_TEST_FOR_EXCEPTION(method != LDL_nopiv && method != Cholesky && method != LDL && method != SymLU, std::logic_error,
                             ss.str().c_str());
  }
  _method_setup = method;
  _method = method;
}

template <typename VT, typename DT>
void Driver<VT, DT>::setFactorizationMethod(const int method) { // 0 - LDL nopivot, 1 - Chol, 2 - LDL, 3 - LU
  {
    std::stringstream ss;
    ss << "Error: the given method (" << method << ") is not supported, 0 - LDL nopivot, 1 - Chol, 2 - LDL, 3 - SymLU";
    TACHO_TEST_FOR_EXCEPTION(method != LDL_nopiv && method != Cholesky && method != LDL && method != SymLU, std::logic_error,
                             ss.str().c_str());
  }
  if (_method_setup == method || method == 1) {
    // switch only if it is the same as method_setup or chol
    _method = method;
    _N->setSolutionMethod(_method);
  }
}

template <typename VT, typename DT>
void Driver<VT, DT>::setOrderConnectedGraphSeparately(const bool order_connected_graph_separately) {
  _order_connected_graph_separately = order_connected_graph_separately;
}

template <typename VT, typename DT>
void Driver<VT, DT>::setGraphAlgorithmType(const int graph_algo_type) {
  _graph_algo_type = graph_algo_type;
}

#ifdef TACHO_DEPRECATED_PARAMETERS
///
/// tasking options
///
template <typename VT, typename DT> void Driver<VT, DT>::setSerialThresholdsize(const ordinal_type serial_thres_size) {
  _serial_thres_size = serial_thres_size;
}

template <typename VT, typename DT> void Driver<VT, DT>::setBlocksize(const ordinal_type mb) { _mb = mb; }

template <typename VT, typename DT> void Driver<VT, DT>::setPanelsize(const ordinal_type nb) { _nb = nb; }

template <typename VT, typename DT> void Driver<VT, DT>::setFrontUpdateMode(const ordinal_type front_update_mode) {
  _front_update_mode = front_update_mode;
}

template <typename VT, typename DT>
void Driver<VT, DT>::setMaxNumberOfSuperblocks(const ordinal_type max_num_superblocks) {
  _max_num_superblocks = max_num_superblocks;
}
#endif

///
/// Level set tools options
///
#ifdef TACHO_DEPRECATED_PARAMETERS
template <typename VT, typename DT> void Driver<VT, DT>::setLevelSetScheduling(const bool levelset) {
  _levelset = levelset;
}
#endif

template <typename VT, typename DT>
void Driver<VT, DT>::setLevelSetOptionDeviceLevelCut(const ordinal_type device_level_cut) {
  _device_level_cut = device_level_cut;
}

template <typename VT, typename DT>
void Driver<VT, DT>::setLevelSetOptionDeviceFunctionThreshold(const ordinal_type device_factor_thres,
                                                              const ordinal_type device_solve_thres) {
  _device_factor_thres = device_factor_thres;
  _device_solve_thres = device_solve_thres;
}

template <typename VT, typename DT> void Driver<VT, DT>::setLevelSetOptionAlgorithmVariant(const ordinal_type variant) {
  #if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL)
  if (variant > 3 || variant < 0) {
   TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "levelset algorithm variants range from 0 to 3");
  }
  #else
  if (variant > 3 || variant < -1) {
   TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "levelset algorithm variants range from -1 to 3 (-1 for serial)");
  }
  #endif
  _variant = variant;
}

template <typename VT, typename DT> void Driver<VT, DT>::setLevelSetOptionNumStreams(const ordinal_type nstreams) {
  _nstreams = nstreams;
}

template <typename VT, typename DT> void Driver<VT, DT>::setPivotTolerance(const mag_type pivot_tol) {
  _pivot_tol = pivot_tol;
}

template <typename VT, typename DT> void Driver<VT, DT>::useNoPivotTolerance() {
  _pivot_tol = 0.0;
}

template <typename VT, typename DT> void Driver<VT, DT>::useDefaultPivotTolerance() {
  using arith_traits = ArithTraits<value_type>;
  _pivot_tol = sqrt(arith_traits::epsilon());
}

template <typename VT, typename DT> void Driver<VT, DT>::storeExplicitTranspose(bool flag) {
  _store_transpose = flag;
}

///
/// get interface
///
template <typename VT, typename DT> ordinal_type Driver<VT, DT>::getNumNonZerosU() const { return _nnz_u; }
template <typename VT, typename DT> ordinal_type Driver<VT, DT>::getNumSupernodes() const { return _nsupernodes; }

template <typename VT, typename DT> typename Driver<VT, DT>::ordinal_type_array Driver<VT, DT>::getSupernodes() const {
  return _supernodes;
}

template <typename VT, typename DT>
typename Driver<VT, DT>::ordinal_type_array Driver<VT, DT>::getPermutationVector() const {
  return _perm;
}

template <typename VT, typename DT>
typename Driver<VT, DT>::ordinal_type_array Driver<VT, DT>::getInversePermutationVector() const {
  return _peri;
}

// internal only
template <typename VT, typename DT> int Driver<VT, DT>::analyze() {
  int r_val(0);
  if (_m <= _small_problem_thres) {
    /// do nothing
    if (_verbose) {
      printf("TachoSolver: Analyze (Small Problem)\n");
      printf("====================================\n");
      printf("  Linear system A\n");
      printf("             number of equations:                             %10d\n", _m);
      printf("\n");
      printf("  A is a small problem ( < %d ) and LAPACK is used\n", _small_problem_thres);
      printf("\n");
    }
  } else {
    const bool use_condensed_graph = (_m_graph > 0 && _m_graph < _m);
    if (use_condensed_graph) {
      Graph graph(_m_graph, _nnz_graph, _h_ap_graph, _h_aj_graph);
      graph_tools_type G(graph);
#if defined(TACHO_HAVE_METIS)
      if (_order_connected_graph_separately) {
        idx_t one_i = 1;
        G.setOption(METIS_OPTION_CCORDER, one_i);
      }
#endif
      if (_graph_algo_type >= 0) {
        G.setAlgorithm(_graph_algo_type);
      }
      G.reorder(_verbose);

      _h_perm_graph = G.PermVector();
      _h_peri_graph = G.InvPermVector();

      r_val = analyze_condensed_graph();
    } else {
      const bool use_graph_partitioner = (_h_perm.extent(0) == 0 && _h_peri.extent(0) == 0);
      if (use_graph_partitioner) {
        Graph graph(_m, _nnz, _h_ap, _h_aj);
        graph_tools_type G(graph);
#if defined(TACHO_HAVE_METIS)
        if (_order_connected_graph_separately) {
          idx_t one_i = 1;
          G.setOption(METIS_OPTION_CCORDER, one_i);
        }
#endif
        if (_graph_algo_type >= 0) {
          G.setAlgorithm(_graph_algo_type);
        }
        G.reorder(_verbose);

        _h_perm = G.PermVector();
        _h_peri = G.InvPermVector();

        r_val = analyze_linear_system();
      } else {
        r_val = analyze_linear_system();
      }
    }
  }
  return r_val;
}

template <typename VT, typename DT> int Driver<VT, DT>::analyze_linear_system() {
  if (_verbose) {
    printf("TachoSolver: Analyze Linear System\n");
    printf("==================================\n");
  }

  {
    symbolic_tools_type S(_m, _h_ap, _h_aj, _h_perm, _h_peri);
    S.symbolicFactorize(_verbose);

    _nnz_u = S.NumNonzerosU();
    _nsupernodes = S.NumSupernodes();
    _stree_level = S.SupernodesTreeLevel();
    _stree_roots = S.SupernodesTreeRoots();

    _supernodes = Kokkos::create_mirror_view(exec_memory_space(), S.Supernodes());
    _gid_super_panel_ptr = Kokkos::create_mirror_view(exec_memory_space(), S.gidSuperPanelPtr());
    _gid_super_panel_colidx = Kokkos::create_mirror_view(exec_memory_space(), S.gidSuperPanelColIdx());
    _sid_super_panel_ptr = Kokkos::create_mirror_view(exec_memory_space(), S.sidSuperPanelPtr());
    _sid_super_panel_colidx = Kokkos::create_mirror_view(exec_memory_space(), S.sidSuperPanelColIdx());
    _blk_super_panel_colidx = Kokkos::create_mirror_view(exec_memory_space(), S.blkSuperPanelColIdx());
    _stree_parent = Kokkos::create_mirror_view(exec_memory_space(), S.SupernodesTreeParent());
    _stree_ptr = Kokkos::create_mirror_view(exec_memory_space(), S.SupernodesTreePtr());
    _stree_children = Kokkos::create_mirror_view(exec_memory_space(), S.SupernodesTreeChildren());

    Kokkos::deep_copy(_supernodes, S.Supernodes());
    Kokkos::deep_copy(_gid_super_panel_ptr, S.gidSuperPanelPtr());
    Kokkos::deep_copy(_gid_super_panel_colidx, S.gidSuperPanelColIdx());
    Kokkos::deep_copy(_sid_super_panel_ptr, S.sidSuperPanelPtr());
    Kokkos::deep_copy(_sid_super_panel_colidx, S.sidSuperPanelColIdx());
    Kokkos::deep_copy(_blk_super_panel_colidx, S.blkSuperPanelColIdx());
    Kokkos::deep_copy(_stree_parent, S.SupernodesTreeParent());
    Kokkos::deep_copy(_stree_ptr, S.SupernodesTreePtr());
    Kokkos::deep_copy(_stree_children, S.SupernodesTreeChildren());

    // perm and peri is updated during symbolic factorization
    _perm = Kokkos::create_mirror_view(exec_memory_space(), _h_perm);
    _peri = Kokkos::create_mirror_view(exec_memory_space(), _h_peri);

    Kokkos::deep_copy(_perm, _h_perm);
    Kokkos::deep_copy(_peri, _h_peri);
  }
  return 0;
}

template <typename VT, typename DT> int Driver<VT, DT>::analyze_condensed_graph() {
  if (_verbose) {
    printf("TachoSolver: Analyze Condensed Graph and Evaporate the Graph\n");
    printf("============================================================\n");
  }

  {
    symbolic_tools_type S(_m_graph, _h_ap_graph, _h_aj_graph, _h_perm_graph, _h_peri_graph);
    S.symbolicFactorize(_verbose);
    S.evaporateSymbolicFactors(_h_aw_graph, _verbose);

    _nnz_u = S.NumNonzerosU();
    _nsupernodes = S.NumSupernodes();
    _stree_level = S.SupernodesTreeLevel();
    _stree_roots = S.SupernodesTreeRoots();

    _supernodes = Kokkos::create_mirror_view(exec_memory_space(), S.Supernodes());
    _gid_super_panel_ptr = Kokkos::create_mirror_view(exec_memory_space(), S.gidSuperPanelPtr());
    _gid_super_panel_colidx = Kokkos::create_mirror_view(exec_memory_space(), S.gidSuperPanelColIdx());
    _sid_super_panel_ptr = Kokkos::create_mirror_view(exec_memory_space(), S.sidSuperPanelPtr());
    _sid_super_panel_colidx = Kokkos::create_mirror_view(exec_memory_space(), S.sidSuperPanelColIdx());
    _blk_super_panel_colidx = Kokkos::create_mirror_view(exec_memory_space(), S.blkSuperPanelColIdx());
    _stree_parent = Kokkos::create_mirror_view(exec_memory_space(), S.SupernodesTreeParent());
    _stree_ptr = Kokkos::create_mirror_view(exec_memory_space(), S.SupernodesTreePtr());
    _stree_children = Kokkos::create_mirror_view(exec_memory_space(), S.SupernodesTreeChildren());
    _perm = Kokkos::create_mirror_view(exec_memory_space(), S.PermVector());
    _peri = Kokkos::create_mirror_view(exec_memory_space(), S.InvPermVector());

    Kokkos::deep_copy(_supernodes, S.Supernodes());
    Kokkos::deep_copy(_gid_super_panel_ptr, S.gidSuperPanelPtr());
    Kokkos::deep_copy(_gid_super_panel_colidx, S.gidSuperPanelColIdx());
    Kokkos::deep_copy(_sid_super_panel_ptr, S.sidSuperPanelPtr());
    Kokkos::deep_copy(_sid_super_panel_colidx, S.sidSuperPanelColIdx());
    Kokkos::deep_copy(_blk_super_panel_colidx, S.blkSuperPanelColIdx());
    Kokkos::deep_copy(_stree_parent, S.SupernodesTreeParent());
    Kokkos::deep_copy(_stree_ptr, S.SupernodesTreePtr());
    Kokkos::deep_copy(_stree_children, S.SupernodesTreeChildren());
    Kokkos::deep_copy(_perm, S.PermVector());
    Kokkos::deep_copy(_peri, S.InvPermVector());

    _h_perm = S.PermVector();
    _h_peri = S.InvPermVector();
  }
  return 0;
}

template <typename VT, typename DT> int Driver<VT, DT>::initialize() {
  if (_verbose) {
    printf("TachoSolver: Initialize(method = %d)\n",_method_setup);
    printf("====================================\n");
  }
  if (_method != _method_setup) {
    // Reset method
    setFactorizationMethod(_method_setup);
  }

  ///
  /// initialize numeric tools
  ///
  if (_m <= _small_problem_thres) {
    /// do nothing
  } else {
    ///
    /// create numeric tools serial for host space
    ///
    NumericToolsFactory<VT, DT> factory;
    factory.setBaseMember(_method, _m, _ap, _aj, _perm, _peri, _nsupernodes, _supernodes, _gid_super_panel_ptr,
                          _gid_super_panel_colidx, _sid_super_panel_ptr, _sid_super_panel_colidx,
                          _blk_super_panel_colidx, _stree_parent, _stree_ptr, _stree_children, _stree_level, _stree_roots,
                          _verbose);

    factory.setLevelSetMember(_variant, _device_level_cut, _device_factor_thres, _device_solve_thres, _nstreams);

    factory.createObject(_N);
  }
  return 0;
}

template <typename VT, typename DT> int Driver<VT, DT>::factorize(const value_type_array &ax) {
  return factorize(ax, _method_setup);
}

template <typename VT, typename DT> int Driver<VT, DT>::factorize(const value_type_array &ax, ordinal_type method) {
  if (method != _method) {
    setFactorizationMethod(method);
  }
  if (_verbose) {
    switch (_method) {
    case LDL_nopiv: {
      printf("TachoSolver: Factorize LDL (no pivot)\n");
      printf("=====================================\n");
      break;
    }
    case Cholesky: {
      printf("TachoSolver: Factorize Cholesky\n");
      printf("===============================\n");
      break;
    }
    case LDL: {
      printf("TachoSolver: Factorize LDL\n");
      printf("==========================\n");
      break;
    }
    case SymLU: {
      printf("TachoSolver: Factorize SymLU\n");
      printf("============================\n");
      break;
    }
    }
    if (_m <= _small_problem_thres) {
      printf( " Small matrix\n" );
    }
  }

  if (_m <= _small_problem_thres) {
    factorize_small_host(ax);
  } else {
    _N->factorize(ax, _store_transpose, _pivot_tol, _verbose);
  }
  return 0;
}

template <typename VT, typename DT> int Driver<VT, DT>::factorize_small_host(const value_type_array &ax) {
  double t_copy(0), t_factor(0);
  {
    Kokkos::Timer timer;

    timer.reset();
    _A = value_type_matrix_host("A", _m, _m);
    auto h_ax = Kokkos::create_mirror_view_and_copy(host_memory_space(), ax);
    for (ordinal_type i = 0; i < _m; ++i) {
      const size_type jbeg = _h_ap(i), jend = _h_ap(i + 1);
      for (size_type j = jbeg; j < jend; ++j) {
        const ordinal_type col = _h_aj(j);
        const bool flag = ((_method == LDL_nopiv && i <= col) || /// upper
                           (_method == Cholesky  && i <= col) || /// upper
                           (_method == LDL && i >= col) ||       /// lower
                           (_method == SymLU));                  /// full matrix
        if (flag)
          _A(i, col) = h_ax(j);
      }
    }
    t_copy = timer.seconds();

    timer.reset();
    switch (_method) {
    case LDL_nopiv:
    case Cholesky : {
      Tacho::Chol<Uplo::Upper, Algo::External>::invoke(_A);
      break;
    }
    case LDL: {
      _P = ordinal_type_array_host("P", 4 * _m);
      _D = value_type_matrix_host("D", _m, 2);
      auto W = value_type_array_host("W", 32 * _m);
      Tacho::LDL<Uplo::Lower, Algo::External>::invoke(_A, _P, W);
      Tacho::LDL<Uplo::Lower, Algo::External>::modify(_A, _P, _D);
      break;
    }
    case SymLU: {
      _P = ordinal_type_array_host("P", 4 * _m);
      Tacho::LU<Algo::External>::invoke(_A, _P);
      Tacho::LU<Algo::External>::modify(_m, _P);
      break;
    }
    default: {
      std::stringstream ss;
      ss << "Error: the solution method (" << _method << ") is not supported, 0 -  LDL no-pivot, 1 - Chol, 2 - LDL, 3 - SymLU";
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, ss.str().c_str());
    }
    }
    t_factor = timer.seconds();
  }

  if (_verbose) {
    switch (_method) {
    case Cholesky: {
      printf("TachoSolver: Factorize Cholesky (SmallDenseFactorization)\n");
      printf("=========================================================\n");
      break;
    }
    case LDL: {
      printf("TachoSolver: Factorize LDL (SmallDenseFactorization)\n");
      printf("====================================================\n");
      break;
    }
    case SymLU: {
      printf("TachoSolver: Factorize SymLU (SmallDenseFactorization)\n");
      printf("======================================================\n");
      break;
    }
    }
    printf("  Time\n");
    printf("             time for copying A into supernodes:              %10.6f s\n", t_copy);
    printf("             time for numeric factorization:                  %10.6f s\n", t_factor);
    printf("             total time spent:                                %10.6f s\n", (t_copy + t_factor));
    printf("\n");
  }

  return 0;
}

template <typename VT, typename DT>
int Driver<VT, DT>::solve(const value_type_matrix &x, const value_type_matrix &b, const value_type_matrix &t) {
  if (_verbose) {
    switch (_method) {
    case LDL_nopiv: {
      printf("TachoSolver: Solve LDL (no pivot)\n");
      printf("=================================\n");
      break;
    }
    case Cholesky: {
      printf("TachoSolver: Solve Cholesky\n");
      printf("===========================\n");
      break;
    }
    case LDL: {
      printf("TachoSolver: Solve LDL\n");
      printf("======================\n");
      break;
    }
    case SymLU: {
      printf("TachoSolver: Solve SymLU\n");
      printf("========================\n");
      break;
    }
    }
  }

  if (_m <= _small_problem_thres) {
    solve_small_host(x, b, t);
  } else {
    TACHO_TEST_FOR_EXCEPTION(t.extent(0) < x.extent(0) || t.extent(1) < x.extent(1), std::logic_error,
                             "Temporary rhs vector t is smaller than x");
    auto tt = Kokkos::subview(t, Kokkos::pair<ordinal_type, ordinal_type>(0, x.extent(0)),
                              Kokkos::pair<ordinal_type, ordinal_type>(0, x.extent(1)));
    _N->solve(x, b, tt, _verbose);
  }
  return 0;
}

template <typename VT, typename DT>
int Driver<VT, DT>::solve_small_host(const value_type_matrix &x, const value_type_matrix &b,
                                     const value_type_matrix &t) {
  Kokkos::Timer timer;
  double t_copy(0), t_solve(0);
  {
    timer.reset();
    Kokkos::deep_copy(x, b);
    t_copy = timer.seconds();

    timer.reset();
    switch (_method) {
    case LDL_nopiv:
    case Cholesky : {
      auto h_x = Kokkos::create_mirror_view_and_copy(host_memory_space(), x);
      Trsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, Algo::External>::invoke(Diag::NonUnit(), 1.0, _A, h_x);
      Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Algo::External>::invoke(Diag::NonUnit(), 1.0, _A, h_x);
      Kokkos::deep_copy(x, h_x);
      break;
    }
    case LDL: {
      auto perm = ordinal_type_array_host(_P.data() + 2 * _m, _m);
      auto peri = ordinal_type_array_host(_P.data() + 3 * _m, _m);
      auto h_x = Kokkos::create_mirror_view_and_copy(host_memory_space(), x);
      auto h_t = Kokkos::create_mirror_view(host_memory_space(), t);

      ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::Internal>::invoke(h_x, perm, h_t);
      Trsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Algo::External>::invoke(Diag::Unit(), 1.0, _A, h_t);
      Scale2x2_BlockInverseDiagonals<Side::Left, Algo::Internal>::invoke(_P, _D, h_t);
      Trsm<Side::Left, Uplo::Lower, Trans::Transpose, Algo::External>::invoke(Diag::Unit(), 1.0, _A, h_t);
      ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::Internal>::invoke(h_t, peri, h_x);
      Kokkos::deep_copy(x, h_x);
      break;
    }
    case SymLU: {
      auto perm = ordinal_type_array_host(_P.data() + 2 * _m, _m);
      auto h_x = Kokkos::create_mirror_view(host_memory_space(), x);
      auto h_t = Kokkos::create_mirror_view(host_memory_space(), t);
      Kokkos::deep_copy(h_t, x);
      ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::Internal>::invoke(h_t, perm, h_x);
      Trsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Algo::External>::invoke(Diag::Unit(), 1.0, _A, h_x);
      Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Algo::External>::invoke(Diag::NonUnit(), 1.0, _A, h_x);
      Kokkos::deep_copy(x, h_x);
      break;
    }
    }
    t_solve = timer.seconds();
  }

  if (_verbose) {
    printf("Summary: NumericTools (SmallDenseSolve)\n");
    printf("=======================================\n");
    printf("  Time\n");
    printf("             time for extra work e.g.,copy rhs:               %10.6f s\n", t_copy);
    printf("             time for numeric solve:                          %10.6f s\n", t_solve);
    printf("             total time spent:                                %10.6f s\n", (t_solve + t_copy));
    printf("\n");
  }
  return 0;
}

template <typename VT, typename DT>
double Driver<VT, DT>::computeRelativeResidual(const value_type_array &ax, const value_type_matrix &x,
                                               const value_type_matrix &b) {
  CrsMatrixBase<value_type, device_type> A;
  A.setExternalMatrix(_m, _m, _nnz, _ap, _aj, ax);

  return Tacho::computeRelativeResidual(A, x, b, _verbose);
}

template <typename VT, typename DT>
void Driver<VT, DT>::computeSpMV(const value_type_array &ax, const value_type_matrix &x, value_type_matrix &b) {
  CrsMatrixBase<value_type, device_type> A;
  A.setExternalMatrix(_m, _m, _nnz, _ap, _aj, ax);

  return Tacho::computeSpMV(A, x, b);
}

template <typename VT, typename DT> int Driver<VT, DT>::exportFactorsToCrsMatrix(crs_matrix_type &A) {
  if (_m <= _small_problem_thres) {
    typedef ArithTraits<value_type> ats;
    const typename ats::mag_type zero(0);

    /// count nonzero elements in dense U
    const ordinal_type m = _m;
    size_type_array_host h_ap("h_ap", m + 1);
    for (ordinal_type i = 0; i < m; ++i)
      for (ordinal_type j = 0; j < m; ++j)
        h_ap(i + 1) += (ats::abs(_A(i, j)) > zero);

    /// serial scan; this is a small problem
    h_ap(0) = 0;
    for (ordinal_type i = 0; i < m; ++i)
      h_ap(i + 1) += h_ap(i);

    /// create a host crs matrix
    const ordinal_type nnz = h_ap(m);
    ordinal_type_array_host h_aj(do_not_initialize_tag("h_aj"), nnz);
    value_type_array_host h_ax(do_not_initialize_tag("h_ax"), nnz);

    for (ordinal_type i = 0, k = 0; i < m; ++i)
      for (ordinal_type j = i; j < m; ++j)
        if (ats::abs(_A(i, j)) > zero) {
          h_aj(k) = j;
          h_ax(k) = _A(i, j);
          ++k;
        }

    crs_matrix_type_host h_A;
    h_A.setExternalMatrix(m, m, nnz, h_ap, h_aj, h_ax);
    /// h_A.showMe(std::cout, true);
    A.clear();
    A.createConfTo(h_A);
    A.copy(h_A);
  } else {
    _N->exportFactorsToCrsMatrix(A, false);
  }
  return 0;
}

template <typename VT, typename DT> int Driver<VT, DT>::release() {
  if (_verbose) {
    printf("TachoSolver: Release\n");
    printf("====================\n");
  }

  {
    if (_N != nullptr)
      _N->release(_verbose);
    delete _N;
    _N = nullptr;
  }
  {
    _method = 0;

    _m = 0;
    _nnz = 0;

    _ap = size_type_array();
    _h_ap = size_type_array_host();
    _aj = ordinal_type_array();
    _h_aj = ordinal_type_array_host();

    _perm = ordinal_type_array();
    _h_perm = ordinal_type_array_host();
    _peri = ordinal_type_array();
    _h_peri = ordinal_type_array_host();

    _m_graph = 0;
    _nnz_graph = 0;

    _h_ap_graph = size_type_array_host();
    _h_aj_graph = ordinal_type_array_host();

    _h_perm_graph = ordinal_type_array_host();
    _h_peri_graph = ordinal_type_array_host();

    _nnz_u = 0;
    _nsupernodes = 0;
    _supernodes = ordinal_type_array();

    _gid_super_panel_ptr = size_type_array();
    _gid_super_panel_colidx = ordinal_type_array();

    _sid_super_panel_ptr = size_type_array();

    _sid_super_panel_colidx = ordinal_type_array();
    _blk_super_panel_colidx = ordinal_type_array();

    _stree_ptr = size_type_array();
    _stree_children = ordinal_type_array();

    _stree_parent = ordinal_type_array();
    _stree_roots = ordinal_type_array_host();

    _A = value_type_matrix_host();

    _verbose = 0;
    _small_problem_thres = 1024;
  }
  return 0;
}

template <typename VT, typename DT> void Driver<VT, DT>::printParameters() {
    printf("\n");
    printf("TachoSolver: Parameters\n");
    printf("=======================\n");
    switch (_method) {
    case LDL_nopiv: {
      printf("Factorize LDL (no pivot, variant = %d)\n",_variant);
      break;
    }
    case Cholesky: {
      printf("Factorize Cholesky (variant = %d)\n",_variant);
      break;
    }
    case LDL: {
      printf("Factorize LDL (variant = %d)\n",_variant);
      break;
    }
    case SymLU: {
      printf("Factorize SymLU (variant = %d)\n",_variant);
      break;
    }
  }
  // ** options
  printf( " verbose             = %d\n", _verbose );             // print
  printf( " store_transpose     = %s\t (store transpose explicitly)\n", (_store_transpose ? "true" : "false"));
  printf( " small_problem_thres = %d\t (smaller than this, use lapack)\n\n", _small_problem_thres);

#ifdef TACHO_DEPRECATED_PARAMETERS
  // // ** tasking options
  printf( " **DEPRECATED** serial_thres_size = %d (serialization threshold size)\n", _serial_thres_size);
  printf( " **DEPRECATED** mb = %d (block size for byblocks algorithms)\n", _mb);
  printf( " **DEPRECATED** nb = %d (panel size for panel algorithms)\n", _nb);
  printf( " **DEPRECATED** max_num_superblocks = %d (superblocks in the memoyrpool)\n\n", _max_num_superblocks);
  printf( " **DEPRECATED** front_update_mode = %d (front update mode 0 - lock, 1 - atomic)\n", _front_update_mode);
#endif

  // ** levelset options
#ifdef TACHO_DEPRECATED_PARAMETERS
  printf( " **DEPRECATED** levelset = %s (use level set code instead of tasking)\n", (_levelset ? "true" : "false"));
#endif
  printf( " device_level_cut    = %d\t (above this level, matrices are computed on device)\n", _device_level_cut);
  printf( " device_factor_thres = %d\t (bigger than this threshold, device function is used)\n",  _device_factor_thres);
  printf( " device_solve_thres  = %d\t (bigger than this threshold, device function is used)\n", _device_solve_thres);
  printf( " nstreams            = %d\t (on device, multi streams are used)\n\n", _nstreams);

  printf( " pivot_tol           = %e\t (tolerance for tiny pivot perturbation)\n", _pivot_tol);


}

} // namespace Tacho

#endif
