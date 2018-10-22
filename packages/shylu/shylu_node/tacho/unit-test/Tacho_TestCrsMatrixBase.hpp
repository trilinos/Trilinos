#ifndef __TACHO_TEST_CRS_MATRIX_BASE_HPP__
#define __TACHO_TEST_CRS_MATRIX_BASE_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_MatrixMarket.hpp"

using namespace Tacho;

typedef CrsMatrixBase<ValueType,HostSpaceType> CrsMatrixBaseHostType;
typedef CrsMatrixBase<ValueType,DeviceSpaceType> CrsMatrixBaseDeviceType;

TEST( CrsMatrixBase, constructor ) {  
  TEST_BEGIN;
  ///
  /// host space crs matrix base
  ///
  const ordinal_type 
    m = 4,
    n = 4,
    nnz = 16;

  CrsMatrixBaseHostType Ah(m, n, nnz);
  EXPECT_EQ(Ah.NumRows(), m);
  EXPECT_EQ(Ah.NumCols(), n);
  EXPECT_EQ(size_t(Ah.NumNonZeros()), size_t(nnz));

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
    EXPECT_EQ(size_t(Ah.RowPtrBegin(i)), size_t(i*n));
    EXPECT_EQ(size_t(Ah.RowPtrEnd(i)), size_t(i*n + n));
  }

  for (ordinal_type k=0;k<nnz;++k) {
    EXPECT_EQ(Ah.Col(k), k%n);
    EXPECT_EQ(Ah.Value(k), k);
  }
  
  ///
  /// device space crs matrix base and mirroring
  ///
  CrsMatrixBaseDeviceType Ad;
  Ad.createMirror(Ah);
  Ad.copy(Ah);

  ///
  /// for checking copy back to host
  ///
  CrsMatrixBaseHostType Bh;
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
  EXPECT_EQ(size_t(Bh.NumNonZeros()), size_t(0));

  EXPECT_EQ(size_t(Bh.RowPtr().extent(0)), size_t(0));
  EXPECT_EQ(size_t(Bh.Cols().extent(0)), size_t(0));
  EXPECT_EQ(size_t(Bh.Values().extent(0)), size_t(0));

  EXPECT_TRUE(Bh.RowPtr().data() == NULL);
  EXPECT_TRUE(Bh.Cols().data() == NULL);
  EXPECT_TRUE(Bh.Values().data() == NULL);
  TEST_END;
}

TEST( CrsMatrixBase, matrixmarket ) {
  TEST_BEGIN;
  std::string inputfilename  = MM_TEST_FILE + ".mtx";
  std::string outputfilename = MM_TEST_FILE + "_read_output.mtx";      

  CrsMatrixBaseHostType Ah;
  MatrixMarket<ValueType>::read(inputfilename, Ah);

  std::ofstream out(outputfilename);
  MatrixMarket<ValueType>::write(out, Ah);

  CrsMatrixBaseHostType Bh;
  MatrixMarket<ValueType>::read(outputfilename, Bh);

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
  TEST_END;
}

#endif
