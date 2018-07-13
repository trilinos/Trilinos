#ifndef __TACHO_TEST_NUMERIC_HPP__
#define __TACHO_TEST_NUMERIC_HPP__

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

#include "TachoExp_NumericTools.hpp"

using namespace Tacho::Experimental;

typedef CrsMatrixBase<ValueType,HostSpaceType> CrsMatrixBaseHostType;
typedef CrsMatrixBase<ValueType,DeviceSpaceType> CrsMatrixBaseDeviceType;

TEST( Numeric, constructor ) {
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

  SymbolicTools S(m, A.RowPtr(), A.Cols(), idx, idx);
  S.symbolicFactorize();

  NumericTools<ValueType,DeviceSpaceType> N(m, A.RowPtr(), A.Cols(), A.Values(),
                                            idx, idx,
                                            S.NumSuperNodes(), S.SuperNodes(),
                                            S.gidSuperPanelPtr(), S.gidSuperPanelColIdx(),
                                            S.sidSuperPanelPtr(), S.sidSuperPanelColIdx(), S.blkSuperPanelColIdx(),
                                            S.SuperNodesTreePtr(), S.SuperNodesTreeChildren(), S.SuperNodesTreeRoots());
}

TEST( Numeric, factorizeCholesky_Serial ) {
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

  SymbolicTools S(A, T);
  S.symbolicFactorize();

  NumericTools<ValueType,DeviceSpaceType> N(A.NumRows(), A.RowPtr(), A.Cols(), A.Values(),
                                            T.PermVector(), T.InvPermVector(),
                                            S.NumSuperNodes(), S.SuperNodes(),
                                            S.gidSuperPanelPtr(), S.gidSuperPanelColIdx(),
                                            S.sidSuperPanelPtr(), S.sidSuperPanelColIdx(), S.blkSuperPanelColIdx(),
                                            S.SuperNodesTreePtr(), S.SuperNodesTreeChildren(), S.SuperNodesTreeRoots());

  N.factorizeCholesky_Serial();
}

TEST( Numeric, factorizeCholesky_Parallel ) {
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

  SymbolicTools S(A, T);
  S.symbolicFactorize();

  NumericTools<ValueType,DeviceSpaceType> N(A.NumRows(), A.RowPtr(), A.Cols(), A.Values(),
                                            T.PermVector(), T.InvPermVector(),
                                            S.NumSuperNodes(), S.SuperNodes(),
                                            S.gidSuperPanelPtr(), S.gidSuperPanelColIdx(),
                                            S.sidSuperPanelPtr(), S.sidSuperPanelColIdx(), S.blkSuperPanelColIdx(),
                                            S.SuperNodesTreePtr(), S.SuperNodesTreeChildren(), S.SuperNodesTreeRoots());

  N.factorizeCholesky_Parallel();
}

#endif



// debug print

  // const auto m = A.NumRows();
  // const auto ap = A.RowPtr();
  // const auto aj = A.Cols();
  // const auto ax = A.Values();
  
  // const auto perm = T.PermVector();
  // const auto peri = T.InvPermVector();
  
  // const size_type ns = S.NumSuperNodes();
  // const auto supernodes = S.SuperNodes();
  // const auto gid_super_panel_ptr =  S.gidSuperPanelPtr();
  // const auto gid_super_panel_colidx = S.gidSuperPanelColIdx();
  // const auto sid_super_panel_ptr = S.sidSuperPanelPtr();
  // const auto sid_super_panel_colidx = S.sidSuperPanelColIdx();
  // const auto blk_super_panel_colidx = S.blkSuperPanelColIdx();

  // printf("supernodes = %d\n", ns);
  // for (size_type i=0;i<ns;++i) {
  //   printf("sid=  %d\n", supernodes(i));
  //   printf("-- gid = ");
  //   for (size_type j=gid_super_panel_ptr(i);j<gid_super_panel_ptr(i+1);++j)
  //     printf("  %d", gid_super_panel_colidx(j));
  //   printf("\n");

  //   printf("-- supernodes id connected = ");
  //   for (size_type j=sid_super_panel_ptr(i);j<(sid_super_panel_ptr(i+1)-1);++j)
  //     printf("  %d", sid_super_panel_colidx(j));
  //   printf("\n");
  //   printf("-- supernodes blocks connected = ");
  //   for (size_type j=sid_super_panel_ptr(i);j<sid_super_panel_ptr(i+1);++j)
  //     printf("  %d", blk_super_panel_colidx(j));
  //   printf("\n");
  // }

