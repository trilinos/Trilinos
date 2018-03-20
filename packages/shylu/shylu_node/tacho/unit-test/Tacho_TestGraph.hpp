#ifndef __TACHO_TEST_GRAPH_HPP__
#define __TACHO_TEST_GRAPH_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_MatrixMarket.hpp"

#include "Tacho_Graph.hpp"

#if defined(TACHO_HAVE_SCOTCH)
#include "Tacho_GraphTools_Scotch.hpp"
#endif

#if defined(TACHO_HAVE_METIS)
#include "Tacho_GraphTools_Metis.hpp"
#endif

#include "Tacho_GraphTools_CAMD.hpp"

using namespace Tacho;

typedef CrsMatrixBase<ValueType,HostSpaceType>   CrsMatrixBaseHostType;
//typedef CrsMatrixBase<ValueType,DeviceSpaceType> CrsMatrixBaseDeviceType;

TEST( Graph, constructor ) {  
  TEST_BEGIN;
  ///
  /// host space crs matrix base
  ///
  const ordinal_type 
    m = 4,
    n = 4,
    nnz = 16;

  CrsMatrixBaseHostType Ah(m, n, nnz);

  ordinal_type cnt = 0;
  for (ordinal_type i=0;i<m;++i) {
    Ah.RowPtrBegin(i) = cnt;
    for (ordinal_type j=0;j<n;++j,++cnt) {
      Ah.Col(cnt) = j;
      Ah.Value(cnt) = i*n+j;
    }
    Ah.RowPtrEnd(i) = cnt;
  }
  
  ///
  /// construction of graph
  ///
  Graph G(Ah);
  
  auto rptr = G.RowPtr();
  auto cidx = G.ColIdx();

  for (ordinal_type i=0;i<G.NumRows();++i) {
    const size_type jbeg = rptr(i), jend = rptr(i+1);
    for (size_type j=jbeg;j<jend;++j) {
      EXPECT_TRUE(cidx(j) != i);
    }
  }
  TEST_END;
}

#if defined(TACHO_HAVE_SCOTCH)
TEST( Graph, scotch ) {
  TEST_BEGIN;
  std::string inputfilename = MM_TEST_FILE + ".mtx";
  CrsMatrixBaseHostType Ah;
  MatrixMarket<ValueType>::read(inputfilename, Ah);

  Graph G(Ah);
  GraphTools_Scotch S(G);

  const ordinal_type m = G.NumRows();
  S.setTreeLevel(log2(m));
  S.setStrategy( SCOTCH_STRATSPEED
                 | SCOTCH_STRATSPEED
                 | SCOTCH_STRATLEVELMAX
                 | SCOTCH_STRATLEVELMIN
                 | SCOTCH_STRATLEAFSIMPLE
                 | SCOTCH_STRATSEPASIMPLE
                 );
  S.reorder();

  ///
  /// perm and invperm should be properly setup 
  ///
  const auto perm = S.PermVector();
  const auto peri = S.InvPermVector();

  for (ordinal_type i=0;i<m;++i) {
    EXPECT_EQ(i, peri(perm(i)));
  }
  TEST_END;
}
#endif

#if defined(TACHO_HAVE_METIS)
TEST( Graph, metis ) {
  TEST_BEGIN;
  std::string inputfilename = MM_TEST_FILE + ".mtx";
  CrsMatrixBaseHostType Ah;
  MatrixMarket<ValueType>::read(inputfilename, Ah);

  Graph G(Ah);
  GraphTools_Metis M(G);

  M.reorder();

  ///
  /// perm and invperm should be properly setup 
  ///
  const auto perm = M.PermVector();
  const auto peri = M.InvPermVector();

  const ordinal_type m = G.NumRows();
  for (ordinal_type i=0;i<m;++i) {
    EXPECT_EQ(i, peri(perm(i)));
  }
  TEST_END;  
}
#endif

#if defined(TACHO_HAVE_SCOTCH)
TEST( Graph, camd ) {
  TEST_BEGIN;
  std::string inputfilename = MM_TEST_FILE + ".mtx";

  CrsMatrixBaseHostType Ah;
  MatrixMarket<ValueType>::read(inputfilename, Ah);

  Graph G(Ah);
  GraphTools_Scotch S(G);

  const ordinal_type m = G.NumRows();
  S.setTreeLevel(log2(m));
  S.setStrategy( SCOTCH_STRATSPEED
                 | SCOTCH_STRATSPEED
                 | SCOTCH_STRATLEVELMAX
                 | SCOTCH_STRATLEVELMIN
                 | SCOTCH_STRATLEAFSIMPLE
                 | SCOTCH_STRATSEPASIMPLE
                 );
  S.reorder();

  GraphTools_CAMD C(G);
  C.setConstraint(S.NumBlocks(), 
                  S.RangeVector(), 
                  S.InvPermVector());
  C.reorder();

  const auto perm = C.PermVector();
  const auto peri = C.InvPermVector();
  
  ///
  /// perm and invperm should be properly setup 
  ///
  for (ordinal_type i=0;i<m;++i) {
    EXPECT_EQ(i, peri(perm(i)));
  }
  TEST_END;
}
#endif

#endif
