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
#ifndef __TACHO_NUMERIC_TOOLS_BASE_HPP__
#define __TACHO_NUMERIC_TOOLS_BASE_HPP__

/// \file Tacho_NumericToolsBase.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_SupernodeInfo.hpp"

namespace Tacho {

template <typename ValueType, typename DeviceType> class NumericToolsBase {
public:
  using value_type = ValueType;
  using device_type = DeviceType;
  using exec_space = typename device_type::execution_space;
  using exec_memory_space = typename device_type::memory_space;

  using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

  using supernode_info_type = SupernodeInfo<value_type, device_type>;
  using crs_matrix_type = typename supernode_info_type::crs_matrix_type;

  using ordinal_type_array = typename supernode_info_type::ordinal_type_array;
  using size_type_array = typename supernode_info_type::size_type_array;
  using value_type_array = typename supernode_info_type::value_type_array;
  using int_type_array = typename supernode_info_type::int_type_array;

  using ordinal_pair_type_array = typename supernode_info_type::ordinal_pair_type_array;
  using value_type_matrix = typename supernode_info_type::value_type_matrix;
  using supernode_type_array = typename supernode_info_type::supernode_type_array;

  using host_device_type = typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type;
  using host_space = typename host_device_type::execution_space;
  using host_memory_space = typename host_device_type::memory_space;

  using ordinal_type_array_host = Kokkos::View<ordinal_type *, host_device_type>;
  using size_type_array_host = Kokkos::View<size_type *, host_device_type>;
  using supernode_type_array_host = Kokkos::View<typename supernode_info_type::supernode_type *, host_device_type>;

protected:
  ///
  /// supernode data structure memory "managed"
  /// this holds all necessary connectivity data
  ///

  // solution method
  ordinal_type _method; // 1 - cholesky, 2 - LDL, 3 - LU

  // matrix input
  ordinal_type _m;
  size_type_array _ap;
  ordinal_type_array _aj;
  value_type_array _ax;

  // graph ordering input
  ordinal_type_array _perm, _peri;

  // supernodes
  ordinal_type _nsupernodes;
  supernode_type_array _supernodes;

  // dof mapping to sparse matrix
  ordinal_type_array _gid_colidx;

  // supernode map and panel size configuration (sid and column blksize)
  ordinal_pair_type_array _sid_block_colidx;

  // supernode tree
  ordinal_type_array_host _stree_level, _stree_roots;

  // output : factors, pivot, diagonal blocks
  value_type_array _superpanel_buf;
  ordinal_type_array _piv;
  value_type_array _diag;

  ///
  /// supernode info: supernode data structure with "unamanged" view
  /// this is passed into computation algorithm without reference counting
  ///
  supernode_info_type _info;

  ///
  /// statistics
  ///
  struct {
    double t_init, t_mode_classification, t_copy, t_factor, t_solve, t_extra;
    double m_used, m_peak;
  } stat;

  inline void track_alloc(const double in) {
    stat.m_used += in;
    stat.m_peak = std::max(stat.m_peak, stat.m_used);
  }

  inline void track_free(const double out) { stat.m_used -= out; }

  inline void reset_stat() {
    stat.t_init = 0;
    stat.t_mode_classification = 0;
    stat.t_factor = 0;
    stat.t_solve = 0;
    stat.t_copy = 0;
    stat.t_extra = 0;
    stat.m_used = 0;
    stat.m_peak = 0;
  }

  virtual void print_stat_init() {
    /// nothing
  }

  virtual void print_stat_factor() {
    const double kilo(1024);
    printf("  Time\n");
    printf("             time for copying A into supernodes:              %10.6f s\n", stat.t_copy);
    printf("             time for numeric factorization:                  %10.6f s\n", stat.t_factor);
    printf("             total time spent:                                %10.6f s\n", (stat.t_copy + stat.t_factor));
    printf("\n");
    printf("  Memory\n");
    printf("             memory used in factorization:                    %10.3f MB\n", stat.m_used / kilo / kilo);
    printf("             peak memory used in factorization:               %10.3f MB\n", stat.m_peak / kilo / kilo);
    printf("\n");
  }

  virtual void print_stat_solve() {
    const double kilo(1024);
    printf("  Time\n");
    printf("             time for extra work e.g.,copy rhs:               %10.6f s\n", stat.t_extra);
    printf("             time for numeric solve:                          %10.6f s\n", stat.t_solve);
    printf("             total time spent:                                %10.6f s\n", (stat.t_solve + stat.t_extra));
    printf("  Memory\n");
    printf("             memory used in solve:                            %10.3f MB\n", stat.m_used / kilo / kilo);
    printf("\n");
  }

  inline void print_stat_memory() {
    const double kilo(1024);
    printf("  Memory\n");
    printf("             memory used now:                                 %10.3f MB\n", stat.m_used / kilo / kilo);
    printf("             peak memory used:                                %10.3f MB\n", stat.m_peak / kilo / kilo);
    printf("\n");
  }

public:
  NumericToolsBase() : _method(0), _m(0), stat() {}

  NumericToolsBase(const NumericToolsBase &b) = default;

  ///
  /// construction (assume input matrix and symbolic are from host)
  ///
  NumericToolsBase(const ordinal_type method,
                   // input matrix A
                   const ordinal_type m, const size_type_array &ap, const ordinal_type_array &aj,
                   // input permutation
                   const ordinal_type_array &perm, const ordinal_type_array &peri,
                   // supernodes
                   const ordinal_type nsupernodes, const ordinal_type_array &supernodes, const size_type_array &gid_ptr,
                   const ordinal_type_array &gid_colidx, const size_type_array &sid_ptr,
                   const ordinal_type_array &sid_colidx, const ordinal_type_array &blk_colidx,
                   const ordinal_type_array &stree_parent, const size_type_array &stree_ptr,
                   const ordinal_type_array &stree_children, const ordinal_type_array_host &stree_level,
                   const ordinal_type_array_host &stree_roots)
      : _method(method), _m(m), _ap(ap), _aj(aj), _perm(perm), _peri(peri), _nsupernodes(nsupernodes),
        _gid_colidx(gid_colidx), _stree_level(stree_level), _stree_roots(stree_roots) {

    reset_stat();

    ///
    /// symbolic input
    ///
    const bool allocate_l_buf = (_method == 3); /// for LU
    _info.initialize(_supernodes, _sid_block_colidx, _superpanel_buf, allocate_l_buf, supernodes, gid_ptr, gid_colidx,
                     sid_ptr, sid_colidx, blk_colidx, stree_parent, stree_ptr, stree_children);
    track_alloc(_superpanel_buf.span() * sizeof(value_type));

    _piv = ordinal_type_array("piv", 4 * _m);
    track_alloc(_piv.span() * sizeof(ordinal_type));

    _diag = value_type_array("diag", 2 * _m);
    track_alloc(_diag.span() * sizeof(value_type));
  }

  virtual ~NumericToolsBase() {}

  inline ordinal_type getSolutionMethod() const { return _method; }

  inline ordinal_type getNumRows() const { return _m; }

  inline ordinal_type getNumCols() const { return _m; }

  inline size_type_array getRowPtr() const { return _ap; }

  inline ordinal_type_array getCols() const { return _aj; }

  inline ordinal_type_array getPermutationVector() const { return _perm; }

  inline ordinal_type_array getInversePermutationVector() const { return _peri; }

  inline ordinal_type_array_host getSupernodesTreeLevel() const { return _stree_level; }

  inline supernode_info_type getSupernodesInfo() const { return _info; }

  inline virtual void release(const ordinal_type verbose = 0) {
    // release diagonal blocks and pivots
    track_free(_piv.span() * sizeof(ordinal_type));
    track_free(_diag.span() * sizeof(value_type));

    // release supernode buffer
    track_free(_superpanel_buf.span() * sizeof(value_type));
    _superpanel_buf = value_type_array();

    if (verbose) {
      printf("Summary: NumericTools (Release)\n");
      printf("===============================\n");
      print_stat_memory(); /// should report zero leak
    }
  }

  inline ordinal_type getMaxSupernodeSize() const { return _info.max_supernode_size; }

  inline ordinal_type getMaxSchurSize() const { return _info.max_schur_size; }

  inline void printMemoryStat(const ordinal_type verbose = 0) {
    if (verbose) {
      printf("Summary: NumericTools (Memory)\n");
      printf("==============================\n");
      print_stat_memory();
    }
  }

  inline virtual void factorize(const value_type_array &ax, const ordinal_type verbose = 0) {
    TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "The function should be overriden by derived classes");
  }

  inline virtual void solve(const value_type_matrix &x, // solution
                            const value_type_matrix &b, // right hand side
                            const value_type_matrix &t, // temporary workspace (store permuted vectors)
                            const ordinal_type verbose = 0) {
    TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "The function should be overriden by derived classes");
  }

  ///
  /// Utility on device
  ///
  inline void exportFactorsToCrsMatrix(crs_matrix_type &A, const bool replace_value_with_one = false) {
    _info.createCrsMatrix(A, replace_value_with_one);
  }

  inline double computeRelativeResidual(const value_type_matrix &x, const value_type_matrix &b) {
    crs_matrix_type A;
    auto d_last = Kokkos::subview(_ap, _m);
    auto h_last = Kokkos::create_mirror_view(host_memory_space(), d_last);
    Kokkos::deep_copy(h_last, d_last);
    A.setExternalMatrix(_m, _m, h_last(), //_ap(_m),
                        _ap, _aj, _ax);

    return Tacho::computeRelativeResidual(A, x, b);
  }
};

} // namespace Tacho
#endif
