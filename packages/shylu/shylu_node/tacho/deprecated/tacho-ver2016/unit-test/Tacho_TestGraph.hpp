#ifndef __TACHO_TEST_GRAPH_HPP__
#define __TACHO_TEST_GRAPH_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_CrsMatrixBase.hpp"
#include "TachoExp_MatrixMarket.hpp"

#include "TachoExp_Graph.hpp"

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

TEST( Graph, constructor ) {  
  ///
  /// host space crs matrix base
  ///
  const ordinal_type 
    m = 4,
    n = 4,
    nnz = 16;

  CrsMatrixBaseHostType Ah("A host", m, n, nnz);

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
  /// device space crs matrix base and mirroring
  ///
  CrsMatrixBaseDeviceType Ad("A device");
  Ad.createMirror(Ah);
  Ad.copy(Ah);

  ///
  /// construction of graph
  ///
  Graph G(Ad);

  auto rptr = G.RowPtr();
  auto cidx = G.ColIdx();

  for (ordinal_type i=0;i<G.NumRows();++i) {
    const size_type jbeg = rptr(i), jend = rptr(i+1);
    for (size_type j=jbeg;j<jend;++j) {
      EXPECT_TRUE(cidx(j) != i);
    }
  }
}

#if defined(HAVE_SHYLU_NODETACHO_SCOTCH)
TEST( Graph, scotch ) {
  CrsMatrixBaseHostType Ah("A host");
  Ah = MatrixMarket<ValueType>::read("test.mtx");

  CrsMatrixBaseDeviceType Ad("A device");
  Ad.createMirror(Ah);
  Ad.copy(Ah);

  Graph G(Ad);
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
  
}
#endif

#if defined(HAVE_SHYLU_NODETACHO_METIS)
TEST( Graph, metis ) {
  CrsMatrixBaseHostType Ah("A host");
  Ah = MatrixMarket<ValueType>::read("test.mtx");

  CrsMatrixBaseDeviceType Ad("A device");
  Ad.createMirror(Ah);
  Ad.copy(Ah);

  Graph G(Ad);
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
  
}
#endif

#if defined(HAVE_SHYLU_NODETACHO_SCOTCH)
TEST( Graph, camd ) {

  CrsMatrixBaseHostType Ah("A host");
  Ah = MatrixMarket<ValueType>::read("test.mtx");

  CrsMatrixBaseDeviceType Ad("A device");
  Ad.createMirror(Ah);
  Ad.copy(Ah);

  Graph G(Ad);
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

}
#endif

#endif
