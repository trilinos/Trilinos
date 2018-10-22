#ifndef __TACHO_TEST_SYMBOLIC_HPP__
#define __TACHO_TEST_SYMBOLIC_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_CrsMatrixBase.hpp"
#include "TachoExp_MatrixMarket.hpp"

#include "TachoExp_Graph.hpp"
#include "TachoExp_SymbolicTools.hpp"

#if defined(HAVE_SHYLU_NODETACHO_SCOTCH)
#include "TachoExp_GraphTools_Scotch.hpp"
#endif

#if defined(HAVE_SHYLU_NODETACHO_METIS)
#include "TachoExp_GraphTools_Metis.hpp"
#endif

#include "TachoExp_GraphTools_CAMD.hpp"

using namespace Tacho::Experimental;

typedef CrsMatrixBase<ValueType,HostSpaceType> CrsMatrixBaseHostType;
typedef CrsMatrixBase<ValueType,DeviceSpaceType> CrsMatrixBaseDeviceType;

TEST( Symbolic, constructor ) {
  const ordinal_type
    m = 4,
    n = 4,
    nnz = 16;

  CrsMatrixBaseHostType A("A", m, n, nnz);

  ordinal_type cnt = 0;
  for (ordinal_type i=0;i<m;++i) {
    A.RowPtrBegin(i) = cnt;
    for (ordinal_type j=0;j<n;++j,++cnt) {
      A.Col(cnt) = j;
      A.Value(cnt) = i*n+j;
    }
    A.RowPtrEnd(i) = cnt;
  }

  typedef Kokkos::View<ordinal_type*,HostSpaceType> ordinal_type_array;

  ordinal_type_array idx("idx", m);
  for (ordinal_type i=0;i<m;++i) idx(i) = i;

  ///
  /// construction of symbolic tools
  ///
  SymbolicTools S(m, A.RowPtr(), A.Cols(), idx, idx);
}

TEST( Symbolic, functions ) {
  CrsMatrixBaseHostType A("A");
  A = MatrixMarket<ValueType>::read("test.mtx");

  Graph G(A);

#if   defined(HAVE_SHYLU_NODETACHO_METIS)
  GraphTools_Metis T(G);
#elif defined(HAVE_SHYLU_NODETACHO_SCOTCH)
  GraphTools_Scotch T(G);
#else
  GraphTools_CAMD T(G);
#endif
  T.reorder();

  typedef Kokkos::View<ordinal_type*,HostSpaceType> ordinal_type_array;
  typedef Kokkos::View<size_type*,HostSpaceType> size_type_array;

  ordinal_type m = A.NumRows();
  size_type_array ap = A.RowPtr();
  ordinal_type_array
    aj = A.Cols(),
    perm = T.PermVector(),
    peri = T.InvPermVector(),
    parent("parent", m),
    ancestor("ancestor", m);

  SymbolicTools::computeEliminationTree(m, ap, aj, perm, peri, parent, ancestor);

  ordinal_type_array work("work", m*4);

  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
  auto post = Kokkos::subview(work, range_type(0*m, 1*m));
  auto w    = Kokkos::subview(work, range_type(1*m, 4*m));
  SymbolicTools::computePostOrdering(m, parent, post, w);

  size_type_array up;
  ordinal_type_array uj;
  SymbolicTools::computeFillPatternUpper(m, ap, aj, perm, peri, up, uj, work);

  ordinal_type_array supernodes;
  SymbolicTools::computeSuperNodes(m, ap, aj, perm, peri, parent, supernodes, work);

  // allocate supernodes
  size_type_array gid_super_panel_ptr, sid_super_panel_ptr;
  ordinal_type_array gid_super_panel_colidx, sid_super_panel_colidx, blk_super_panel_colidx;
  SymbolicTools::allocateSuperNodes(m, up, uj, supernodes, work,
                                    gid_super_panel_ptr,
                                    gid_super_panel_colidx,
                                    sid_super_panel_ptr,
                                    sid_super_panel_colidx,
                                    blk_super_panel_colidx);

  size_type_array stree_ptr;
  ordinal_type_array stree_children, stree_roots;
  SymbolicTools::computeSuperNodesAssemblyTree(parent,
                                               supernodes,
                                               stree_ptr,
                                               stree_children,
                                               stree_roots,
                                               work);

  // const size_type numSuperNodes = supernodes.extent(0) - 1;
  // printf("supernodes = \n");
  // for (size_type i=0;i<numSuperNodes;++i) {
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
  CrsMatrixBaseHostType A("A");
  A = MatrixMarket<ValueType>::read("test.mtx");

  Graph G(A);

#if   defined(HAVE_SHYLU_NODETACHO_METIS)
  GraphTools_Metis T(G);
#elif defined(HAVE_SHYLU_NODETACHO_SCOTCH)
  GraphTools_Scotch T(G);
#else
  GraphTools_CAMD T(G);
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


#endif
