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
#ifndef __TACHO_TEST_SYMBOLIC_TOOLS_HPP__
#define __TACHO_TEST_SYMBOLIC_TOOLS_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_MatrixMarket.hpp"

#include "Tacho_Graph.hpp"
#include "Tacho_SymbolicTools.hpp"

#include "Tacho_GraphTools.hpp"

#if defined(TACHO_HAVE_METIS)
#include "Tacho_GraphTools_Metis.hpp"
#endif

namespace Test {

  TEST( Symbolic, constructor ) {
    const ordinal_type
      m = 4,
      n = 4,
      nnz = 16;

    crs_matrix_base_type_host A(m, n, nnz);

    ordinal_type cnt = 0;
    for (ordinal_type i=0;i<m;++i) {
      A.RowPtrBegin(i) = cnt;
      for (ordinal_type j=0;j<n;++j,++cnt) {
        A.Col(cnt) = j;
        A.Value(cnt) = i*n+j;
      }
      A.RowPtrEnd(i) = cnt;
    }

    ordinal_type_array_type_host idx("idx", m);
    for (ordinal_type i=0;i<m;++i) idx(i) = i;

    ///
    /// construction of symbolic tools
    ///
    SymbolicTools S(m, A.RowPtr(), A.Cols(), idx, idx);
    //std::cout << "num supernodes = " << S.NumSupernodes() << std::endl;
  }

  TEST( Symbolic, functions ) {
    std::string filename = MM_TEST_FILE;
    crs_matrix_base_type_host A;
    MatrixMarket<value_type>::read(filename, A);

    Graph G(A);

#if   defined(TACHO_HAVE_METIS)
    GraphTools_Metis T(G);
#else
    GraphTools T(G);
#endif
    T.reorder();

    ordinal_type m = A.NumRows();
    size_type_array_type_host ap = A.RowPtr();
    ordinal_type_array_type_host
      aj = A.Cols(),
      perm = T.PermVector(),
      peri = T.InvPermVector(),
      parent("parent", m),
      ancestor("ancestor", m);

    SymbolicTools::computeEliminationTree(m, ap, aj, perm, peri, parent, ancestor);

    ordinal_type_array_type_host work("work", m*4);

    using range_type = Kokkos::pair<ordinal_type,ordinal_type>;
    auto post = Kokkos::subview(work, range_type(0*m, 1*m));
    auto w    = Kokkos::subview(work, range_type(1*m, 4*m));
    SymbolicTools::computePostOrdering(m, parent, post, w);

    size_type_array_type_host up;
    ordinal_type_array_type_host uj;
    SymbolicTools::computeFillPatternUpper(m, ap, aj, perm, peri, up, uj, work);

    ordinal_type_array_type_host supernodes;
    SymbolicTools::computeSupernodes(m, ap, aj, perm, peri, parent, supernodes, work);

    // allocate supernodes
    size_type_array_type_host gid_super_panel_ptr, sid_super_panel_ptr;
    ordinal_type_array_type_host gid_super_panel_colidx, sid_super_panel_colidx, blk_super_panel_colidx;
    SymbolicTools::allocateSupernodes(m, up, uj, supernodes, work,
                                      gid_super_panel_ptr,
                                      gid_super_panel_colidx,
                                      sid_super_panel_ptr,
                                      sid_super_panel_colidx,
                                      blk_super_panel_colidx);

    size_type_array_type_host stree_ptr;
    ordinal_type_array_type_host stree_level, stree_parent, stree_children, stree_roots;
    SymbolicTools::computeSupernodesAssemblyTree(parent,
                                                 supernodes,
                                                 stree_level,
                                                 stree_parent,
                                                 stree_ptr,
                                                 stree_children,
                                                 stree_roots,
                                                 work);

    // const size_type numSupernodes = supernodes.extent(0) - 1;
    // printf("supernodes = \n");
    // for (size_type i=0;i<numSupernodes;++i) {
    //   printf("sid=  %d\n", supernodes(i));
    //   printf("-- gid = \n");
    //   for (size_type j=gid_super_panel_ptr(i);j<gid_super_panel_ptr(i+1);++j)
    //     printf("  %d", gid_super_panel_colidx(j));
    //   printf("\n");

    //   printf("-- supernodes id connected = \n");
    //   for (size_type j=sid_super_panel_ptr(i);j<(sid_super_panel_ptr(i+1)-1);++j)
    //     printf("  %d", sid_super_panel_colidx(j));
    //   printf("\n");
    //   printf("-- supernodes blocks connected = \n");
    //   for (size_type j=sid_super_panel_ptr(i);j<sid_super_panel_ptr(i+1);++j)
    //     printf("  %d", blk_super_panel_colidx(j));
    //   printf("\n");
    // }
  }

  TEST( Symbolic, interface ) {
    std::string filename = MM_TEST_FILE;

    crs_matrix_base_type_host A;
    MatrixMarket<value_type>::read(filename, A);

    Graph G(A);

#if   defined(TACHO_HAVE_METIS)
    GraphTools_Metis T(G);
#else
    GraphTools T(G);  
#endif
    T.reorder();

    {
      SymbolicTools S(A.NumRows(),
                      A.RowPtr(),
                      A.Cols(),
                      T.PermVector(),
                      T.InvPermVector());

      S.symbolicFactorize();
    }
    {
      SymbolicTools S(A, T);
      S.symbolicFactorize();
    }
  }

}
#endif
