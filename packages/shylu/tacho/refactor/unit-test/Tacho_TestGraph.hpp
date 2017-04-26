#ifndef __TACHO_TEST_GRAPH_HPP__
#define __TACHO_TEST_GRAPH_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_CrsMatrixBase.hpp"
#include "TachoExp_MatrixMarket.hpp"

#include "TachoExp_Graph.hpp"
#if defined(HAVE_SHYLUTACHO_SCOTCH)
#include "TachoExp_GraphTools_Scotch.hpp"
#endif

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

#if defined(HAVE_SHYLUTACHO_SCOTCH)
TEST( Graph, scotch ) {
  CrsMatrixBaseHostType Ah("A host");
  Ah = MatrixMarket<ValueType>::read("test.mtx");

  CrsMatrixBaseDeviceType Ad("A device");
  Ad.createMirror(Ah);
  Ad.copy(Ah);

  Graph G(Ad);
  // GraphTools_Scotch S(G);
  // S.reorder();
  // std::ofstream file("test_graph_scotch.txt");
  // //S.showMe(file);
}
#endif

#endif
