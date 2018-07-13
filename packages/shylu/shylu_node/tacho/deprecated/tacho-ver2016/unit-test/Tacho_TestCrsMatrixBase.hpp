#ifndef __TACHO_TEST_CRS_MATRIX_BASE_HPP__
#define __TACHO_TEST_CRS_MATRIX_BASE_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_CrsMatrixBase.hpp"
#include "TachoExp_MatrixMarket.hpp"

using namespace Tacho::Experimental;

typedef CrsMatrixBase<ValueType,HostSpaceType> CrsMatrixBaseHostType;
typedef CrsMatrixBase<ValueType,DeviceSpaceType> CrsMatrixBaseDeviceType;

TEST( CrsMatrixBase, constructor ) {  
  ///
  /// host space crs matrix base
  ///
  const ordinal_type 
    m = 4,
    n = 4,
    nnz = 16;

  CrsMatrixBaseHostType Ah("A host", m, n, nnz);
  EXPECT_EQ(Ah.NumRows(), m);
  EXPECT_EQ(Ah.NumCols(), n);
  EXPECT_EQ(Ah.NumNonZeros(), nnz);

  ordinal_type cnt = 0;
  for (ordinal_type i=0;i<m;++i) {
    Ah.RowPtrBegin(i) = cnt;
    for (ordinal_type j=0;j<n;++j,++cnt) {
      Ah.Col(cnt) = j;
      Ah.Value(cnt) = i*n+j;
    }
    Ah.RowPtrEnd(i) = cnt;
  }

  for (ordinal_type i=0;i<m;++i) {
    EXPECT_EQ(Ah.RowPtrBegin(i), i*n);
    EXPECT_EQ(Ah.RowPtrEnd(i), i*n + n);
  }

  for (ordinal_type k=0;k<nnz;++k) {
    EXPECT_EQ(Ah.Col(k), k%n);
    EXPECT_EQ(Ah.Value(k), k);
  }
  
  ///
  /// device space crs matrix base and mirroring
  ///
  CrsMatrixBaseDeviceType Ad("A device");
  Ad.createMirror(Ah);
  Ad.copy(Ah);

  ///
  /// for checking copy back to host
  ///
  CrsMatrixBaseHostType Bh("B host");
  Bh.createConfTo(Ad);
  Bh.copy(Ad);
  
  for (ordinal_type i=0;i<m;++i) {
    EXPECT_EQ(Ah.RowPtrBegin(i), Bh.RowPtrBegin(i)); 
    EXPECT_EQ(Ah.RowPtrEnd(i), Bh.RowPtrEnd(i));
  }

  for (ordinal_type k=0;k<nnz;++k) {
    EXPECT_EQ(Ah.Col(k), Bh.Col(k));
    EXPECT_EQ(Ah.Value(k), Bh.Value(k));
  }
  
  ///
  /// clear
  ///
  Bh.clear();
  EXPECT_EQ(Bh.NumRows(), 0);
  EXPECT_EQ(Bh.NumCols(), 0);
  EXPECT_EQ(Bh.NumNonZeros(), 0);

  EXPECT_EQ(Bh.RowPtr().extent(0), 0);
  EXPECT_EQ(Bh.Cols().extent(0), 0);
  EXPECT_EQ(Bh.Values().extent(0), 0);

  EXPECT_TRUE(Bh.RowPtr().data() == NULL);
  EXPECT_TRUE(Bh.Cols().data() == NULL);
  EXPECT_TRUE(Bh.Values().data() == NULL);
}

TEST( CrsMatrixBase, matrixmarket ) {
  CrsMatrixBaseHostType Ah("A host");
  Ah = MatrixMarket<ValueType>::read("test.mtx");

  std::ofstream out("test_mm_read_output.mtx");
  MatrixMarket<ValueType>::write(out, Ah);

  CrsMatrixBaseHostType Bh("B host");
  Bh = MatrixMarket<ValueType>::read("test_mm_read_output.mtx");
  
  ///
  /// read and write the matrix and read again, 
  /// then check if they are same
  ///
  EXPECT_EQ(Ah.NumRows(), Bh.NumRows());
  EXPECT_EQ(Ah.NumCols(), Bh.NumCols());
  EXPECT_EQ(Ah.NumNonZeros(), Bh.NumNonZeros());

  const ordinal_type m = Ah.NumRows();
  for (ordinal_type i=0;i<m;++i) {
    EXPECT_EQ(Ah.RowPtrBegin(i), Bh.RowPtrBegin(i));
    const ordinal_type jbeg = Ah.RowPtrBegin(i), jend = Ah.RowPtrEnd(i);
    for (ordinal_type j=jbeg;j<jend;++j) {
      EXPECT_EQ(Ah.Col(j), Bh.Col(j));      
      EXPECT_EQ(Ah.Value(j), Bh.Value(j));      
    }
  }
}

#endif
