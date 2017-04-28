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

#if defined(HAVE_SHYLUTACHO_SCOTCH)
#include "TachoExp_GraphTools_Scotch.hpp"
#endif

#if defined(HAVE_SHYLUTACHO_METIS)
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

TEST( Symbolic, eliminationtree ) {  
  CrsMatrixBaseHostType A("A");
  A = MatrixMarket<ValueType>::read("test.mtx");

  Graph G(A);

#if   defined(HAVE_SHYLUTACHO_METIS)
  GraphTools_Metis T(G);
#elif defined(HAVE_SHYLUTACHO_SCOTCH)
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
  
  ordinal_type_array 
    post("post", m),
    work("work", m*3);
  SymbolicTools::computePostOrdering(m, parent, post, work);

  size_type_array up;
  ordinal_type_array uj;
  SymbolicTools::computeFillPatternUpper(m, ap, aj, perm, peri, up, uj, work);
}


#endif
