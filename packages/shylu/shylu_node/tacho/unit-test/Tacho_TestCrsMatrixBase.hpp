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
#ifndef __TACHO_TEST_CRS_MATRIX_BASE_HPP__
#define __TACHO_TEST_CRS_MATRIX_BASE_HPP__

#include "Tacho_Util.hpp"
#include "Tacho_MatrixMarket.hpp"

namespace Test {

  TEST( CrsMatrixBase, coo ) {
    {
      auto a = Coo<double>();
      EXPECT_EQ(a.i, 0);
      EXPECT_EQ(a.j, 0);
      EXPECT_EQ(a.j, 0.0);
    }
    {
      auto a = Coo<double>(1,3, 3.0);
      auto b = Coo<double>(1,3,10.0);
      auto c = Coo<double>(2,3, 3.0);
      auto d = Coo<double>(1,1, 3.0);
    
      EXPECT_TRUE(a == b);
      EXPECT_TRUE(a != c);
      EXPECT_TRUE(a < c);
      EXPECT_FALSE(a < d);
    }
  }

  TEST( CrsMatrixBase, constructor ) {  
    ///
    /// host space crs matrix base
    ///
    const int 
      m = 4,
      n = 4,
      nnz = 16;

    crs_matrix_base_type_host Ah(m, n, nnz);
    EXPECT_EQ(Ah.NumRows(), m);
    EXPECT_EQ(Ah.NumCols(), n);
    EXPECT_EQ(size_t(Ah.NumNonZeros()), size_t(nnz));

    int cnt = 0;
    for (int i=0;i<m;++i) {
      Ah.RowPtrBegin(i) = cnt;
      for (int j=0;j<n;++j,++cnt) {
        Ah.Col(cnt) = j;
        Ah.Value(cnt) = i*n+j;
      }
      Ah.RowPtrEnd(i) = cnt;
    }

    for (int i=0;i<m;++i) {
      EXPECT_EQ(size_t(Ah.RowPtrBegin(i)), size_t(i*n));
      EXPECT_EQ(size_t(Ah.RowPtrEnd(i)), size_t(i*n + n));
    }

    for (int k=0;k<nnz;++k) {
      EXPECT_EQ(Ah.Col(k), k%n);
      EXPECT_EQ(Ah.Value(k), k);
    }
  
    ///
    /// device space crs matrix base and mirroring
    ///
    crs_matrix_base_type Ad;
    Ad.createMirror(Ah);
    Ad.copy(Ah);

    ///
    /// for checking copy back to host
    ///
    crs_matrix_base_type_host Bh;
    Bh.createConfTo(Ad);
    Bh.copy(Ad);
  
    for (int i=0;i<m;++i) {
      EXPECT_EQ(Ah.RowPtrBegin(i), Bh.RowPtrBegin(i)); 
      EXPECT_EQ(Ah.RowPtrEnd(i), Bh.RowPtrEnd(i));
    }

    for (int k=0;k<nnz;++k) {
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
  
  }

  TEST( CrsMatrixBase, permute ) {
    const std::string filename  = 
      std::is_same<value_type,magnitude_type>::value ? "test_double.mtx" : "test_dcomplex.mtx";

    ///
    /// host crs matrix read from matrix market
    ///
    crs_matrix_base_type_host Ah, Bh;
    MatrixMarket<value_type>::read(filename, Bh);
    Ah.createConfTo(Bh);

    ///
    /// device crs matrix
    ///
    crs_matrix_base_type Ad, Bd;
    Ad.createMirror(Ah);
    Bd.createMirror(Bh);
    Bd.copy(Bh);

    ///
    /// random permutation vector
    ///
    const int m = Ad.NumRows();
    ordinal_type_array_type_host perm("perm", m), peri("peri", m);

    for (int i=0;i<m;++i)
      perm(i) = i;

    {
      std::random_device rd;
      std::mt19937 g(rd());
      std::shuffle(perm.data(), perm.data()+m, g);
    }

    for (int i=0;i<m;++i)
      peri(perm(i)) = i;

    ///
    /// A = P B P^{-1}
    ///
    applyPermutationToCrsMatrix(Ad, Bd, perm, peri);

    ///
    /// check on host
    ///
    Ah.copy(Ad);
    for (int i=0;i<m;++i) {
      const int row = perm(i);
      const int jcnt = Bh.RowPtrEnd(row) - Bh.RowPtrBegin(row);
      const int boff = Bh.RowPtrBegin(row);
      const int aoff = Ah.RowPtrBegin(i);
      for (int j=0;j<jcnt;++j) {
        EXPECT_EQ(Ah.Col(aoff+j), peri(Bh.Col(boff+j)));
        EXPECT_EQ(Ah.Value(aoff+j), Bh.Value(boff+j));
      }
    }
  
  }
}


#endif
