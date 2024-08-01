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
#ifndef __TACHO_SYMBOLIC_TOOLS_HPP__
#define __TACHO_SYMBOLIC_TOOLS_HPP__

/// \file Tacho_SymbolicTools.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

class SymbolicTools {
public:
  typedef typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type host_device_type;
  typedef typename host_device_type::execution_space host_space;
  typedef typename host_device_type::memory_space host_memory_space;

  typedef Kokkos::View<ordinal_type *, host_device_type> ordinal_type_array;
  typedef Kokkos::View<size_type *, host_device_type> size_type_array;

  typedef Kokkos::pair<ordinal_type, ordinal_type> range_type;

  ///
  /// supernode tools
  ///

  // Tim Davis, Direct Methods for Sparse Linear Systems, Siam, p 42.
  static void computeEliminationTree(const ordinal_type m, const size_type_array &ap, const ordinal_type_array &aj,
                                     const ordinal_type_array &perm, const ordinal_type_array &peri,
                                     const ordinal_type_array &parent, const ordinal_type_array &ancestor);

  // Tim Davis, Direct Methods for Sparse Linear Systems, Siam, p 45.
  static ordinal_type TreeDepthFirstSearch(const ordinal_type j, const ordinal_type c, const ordinal_type_array &head,
                                           const ordinal_type_array &next, const ordinal_type_array &post,
                                           const ordinal_type_array &stack);

  // Tim Davis, Direct Methods for Sparse Linear Systems, Siam, p 45.
  static void computePostOrdering(const ordinal_type m, const ordinal_type_array &parent,
                                  const ordinal_type_array &post, const ordinal_type_array &work);

  // Tim Davis, Algorithm 849: A Concise Sparse Cholesky Factorization Package
  // ACM TOMS Vol 31 No. 4 pp 587--591.
  static void computeFillPatternUpper(const ordinal_type m, const size_type_array &ap, const ordinal_type_array &aj,
                                      const ordinal_type_array &perm, const ordinal_type_array &peri,
                                      /* */ size_type_array &up,
                                      /* */ ordinal_type_array &uj, const ordinal_type_array &work);

  // Joseph, W. H. Liu, Esmond G. Ng, and Barry W. Peyton,
  // "On Finding Supernodes for Sparse Matrix Computations,"
  // SIAM J. Matrix Anal. Appl., Vol. 14, No. 1, pp. 242-252.
  static void computeSupernodes(const ordinal_type m, const size_type_array &ap, const ordinal_type_array &aj,
                                const ordinal_type_array &perm, const ordinal_type_array &peri,
                                const ordinal_type_array &parent,
                                /* */ ordinal_type_array &supernodes, const ordinal_type_array &work);

  /// Based on the symbolic factors, allocate pannels
  static void allocateSupernodes(const ordinal_type m, const size_type_array &ap, const ordinal_type_array &aj,
                                 const ordinal_type_array &supernodes, const ordinal_type_array &work,
                                 /* */ size_type_array &gid_super_panel_ptr,
                                 /* */ ordinal_type_array &gid_super_panel_colidx,
                                 /* */ size_type_array &sid_super_panel_ptr,
                                 /* */ ordinal_type_array &sid_super_panel_colidx,
                                 /* */ ordinal_type_array &blk_super_panel_colidx);

  /// construct tree explicitly
  static void computeSupernodesAssemblyTree(const ordinal_type_array &parent, const ordinal_type_array &supernodes,
                                            /* */ ordinal_type_array &stree_level,
                                            /* */ ordinal_type_array &stree_parent,
                                            /* */ size_type_array &stree_ptr,
                                            /* */ ordinal_type_array &stree_children,
                                            /* */ ordinal_type_array &stree_roots, const ordinal_type_array &work);

  ///
  /// evaporation tools
  ///
  static void scanWeights(const ordinal_type m, const ordinal_type_array &aw, const ordinal_type_array &perm,
                          /* */ size_type_array &as,
                          /* */ size_type_array &aq);

  static void evaporateGraph(const ordinal_type m, const size_type_array &ap, const ordinal_type_array &aj,
                             const size_type_array &as,
                             /* */ size_type_array &ap_eva,
                             /* */ ordinal_type_array &aj_eva);

  static void evaporatePermutationVectors(const ordinal_type m, const ordinal_type_array &perm,
                                          const ordinal_type m_eva, const ordinal_type_array &aw,
                                          const size_type_array &aq,
                                          /* */ ordinal_type_array &perm_eva,
                                          /* */ ordinal_type_array &peri_eva);

  static void evaporateSupernodes(const ordinal_type_array &supernodes, const size_type_array &sid_super_panel_ptr,
                                  const size_type_array &gid_super_panel_ptr,
                                  const ordinal_type_array &gid_super_panel_colidx,
                                  const ordinal_type_array &blk_super_panel_colidx, const ordinal_type_array &perm,
                                  const size_type_array &as, const ordinal_type_array &peri_eva,
                                  /* */ ordinal_type_array &supernodes_eva,
                                  /* */ size_type_array &gid_super_panel_ptr_eva,
                                  /* */ ordinal_type_array &gid_super_panel_colidx_eva,
                                  /* */ ordinal_type_array &blk_super_panel_colidx_eva);

private:
  // matrix input
  ordinal_type _m;
  size_type_array _ap;
  ordinal_type_array _aj;

  // graph ordering input
  ordinal_type_array _perm, _peri;

  // supernodes output
  ordinal_type_array _supernodes;

  // dof mapping to sparse matrix
  size_type_array _gid_super_panel_ptr;
  ordinal_type_array _gid_super_panel_colidx;

  // supernode map and panel size configuration
  size_type_array _sid_super_panel_ptr;
  ordinal_type_array _sid_super_panel_colidx, _blk_super_panel_colidx;

  // supernode elimination tree (parent - children)
  size_type_array _stree_ptr;
  ordinal_type_array _stree_children, _stree_roots;

  // supernode elimination tree (child - parent)
  ordinal_type_array _stree_parent;

  // level information of supernodes
  ordinal_type_array _stree_level;

  // stat
  struct {
    ordinal_type nrows, nroots;
    size_type nnz_a, nnz_u;
    ordinal_type nsupernodes, max_nchildren, largest_supernode, largest_schur;
    ordinal_type nleaves, height; // tree
  } stat;

public:
  SymbolicTools();
  SymbolicTools(const SymbolicTools &b);

  ///
  /// construction
  ///
  SymbolicTools(const ordinal_type m, const size_type_array &ap, const ordinal_type_array &aj,
                const ordinal_type_array &perm, const ordinal_type_array &peri);

  template <typename CrsMatBaseType, typename GraphToolType> SymbolicTools(CrsMatBaseType &A, GraphToolType &G) {
    _m = A.NumRows();

    _ap = Kokkos::create_mirror_view(host_memory_space(), A.RowPtr());
    _aj = Kokkos::create_mirror_view(host_memory_space(), A.Cols());
    _perm = Kokkos::create_mirror_view(host_memory_space(), G.PermVector());
    _peri = Kokkos::create_mirror_view(host_memory_space(), G.InvPermVector());

    Kokkos::deep_copy(_ap, A.RowPtr());
    Kokkos::deep_copy(_aj, A.Cols());
    Kokkos::deep_copy(_perm, G.PermVector());
    Kokkos::deep_copy(_peri, G.InvPermVector());
  }

  ordinal_type NumSupernodes() const;
  ordinal_type_array Supernodes() const;
  size_type_array gidSuperPanelPtr() const;
  ordinal_type_array gidSuperPanelColIdx() const;
  size_type_array sidSuperPanelPtr() const;
  ordinal_type_array sidSuperPanelColIdx() const;
  ordinal_type_array blkSuperPanelColIdx() const;
  ordinal_type_array SupernodesTreeParent() const;
  size_type_array SupernodesTreePtr() const;
  ordinal_type_array SupernodesTreeChildren() const;
  ordinal_type_array SupernodesTreeRoots() const;
  ordinal_type_array SupernodesTreeLevel() const;
  ordinal_type_array PermVector() const;
  ordinal_type_array InvPermVector() const;

  void symbolicFactorize(const ordinal_type verbose = 0);
  void evaporateSymbolicFactors(const ordinal_type_array &aw, const ordinal_type verbose = 0);
};

} // namespace Tacho
#endif
