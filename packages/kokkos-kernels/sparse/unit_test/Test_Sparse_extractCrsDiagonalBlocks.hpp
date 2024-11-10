//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#include "KokkosSparse_Utils.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosBlas1_nrm2.hpp"
#include "KokkosKernels_TestUtils.hpp"

namespace Test {
template <typename scalar_t, typename lno_t, typename size_type, typename device>
void run_test_extract_diagonal_blocks(int nrows, int nblocks) {
  using RowMapType     = Kokkos::View<size_type *, device>;
  using EntriesType    = Kokkos::View<lno_t *, device>;
  using ValuesType     = Kokkos::View<scalar_t *, device>;
  using RowMapType_hm  = typename RowMapType::HostMirror;
  using EntriesType_hm = typename EntriesType::HostMirror;
  using ValuesType_hm  = typename ValuesType::HostMirror;
  using crsMat_t       = CrsMatrix<scalar_t, lno_t, device, void, size_type>;

  crsMat_t A;
  std::vector<crsMat_t> DiagBlks(nblocks);
  std::vector<crsMat_t> DiagBlks_rcm(nblocks);

  if (nrows != 0) {
    // Generate test matrix
    const size_type nnz = 2 + (nrows - 2) * 3 + 2;
    RowMapType_hm hrow_map("hrow_map", nrows + 1);
    EntriesType_hm hentries("hentries", nnz);
    ValuesType_hm hvalues("hvalues", nnz);

    // first row
    hrow_map(0) = 0;
    hentries(0) = 0;
    hentries(1) = 1;
    hvalues(0)  = 0;
    hvalues(1)  = 1;
    // rows in between
    int cnt = 2;
    for (int i = 1; i <= (nrows - 2); i++) {
      hrow_map(i)       = cnt;
      hentries(cnt)     = -1 + i;
      hentries(cnt + 1) = 0 + i;
      hentries(cnt + 2) = 1 + i;
      hvalues(cnt)      = -1 + i;
      hvalues(cnt + 1)  = 0 + i;
      hvalues(cnt + 2)  = 1 + i;
      cnt += 3;
    }
    // last row
    hrow_map(nrows - 1) = cnt;
    hentries(nnz - 2)   = nrows - 2;
    hentries(nnz - 1)   = nrows - 1;
    hvalues(nnz - 2)    = nrows - 2;
    hvalues(nnz - 1)    = nrows - 1;
    // last element of row_map
    hrow_map(nrows) = nnz;

    // Allocate A on device memory
    RowMapType row_map("row_map", nrows + 1);
    EntriesType entries("entries", nnz);
    ValuesType values("values", nnz);

    // Copy from host to device
    Kokkos::deep_copy(row_map, hrow_map);
    Kokkos::deep_copy(entries, hentries);
    Kokkos::deep_copy(values, hvalues);

    // Construct a CRS matrix
    A = crsMat_t("CrsMatrix", nrows, nrows, nnz, values, row_map, entries);
  }

  // Extract
  KokkosSparse::Impl::kk_extract_diagonal_blocks_crsmatrix_sequential(A, DiagBlks);

  auto perm = KokkosSparse::Impl::kk_extract_diagonal_blocks_crsmatrix_sequential(A, DiagBlks_rcm, true);

  // Checking
  lno_t numRows = 0;
  lno_t numCols = 0;
  for (int i = 0; i < nblocks; i++) {
    numRows += DiagBlks[i].numRows();
    numCols += DiagBlks[i].numCols();
  }

  EXPECT_TRUE(numRows == static_cast<lno_t>(nrows));
  EXPECT_TRUE(numCols == static_cast<lno_t>(nrows));

  if (nrows > 0) {
    bool flag       = true;
    lno_t col_start = 0;
    for (int i = 0; i < nblocks; i++) {
      RowMapType_hm hrow_map_diagblk("hrow_map_diagblk", DiagBlks[i].numRows() + 1);
      EntriesType_hm hentries_diagblk("hentries_diagblk", DiagBlks[i].nnz());
      ValuesType_hm hvalues_diagblk("hvalues_diagblk", DiagBlks[i].nnz());

      Kokkos::deep_copy(hrow_map_diagblk, DiagBlks[i].graph.row_map);
      Kokkos::deep_copy(hentries_diagblk, DiagBlks[i].graph.entries);
      Kokkos::deep_copy(hvalues_diagblk, DiagBlks[i].values);

      for (int j = 0; j < static_cast<int>(DiagBlks[i].numRows()); j++) {
        size_type k1 = hrow_map_diagblk(j);
        size_type k2 = hrow_map_diagblk(j + 1);
        for (size_type k = k1; k < k2; k++) {
          scalar_t col = static_cast<scalar_t>(hentries_diagblk(k) + col_start);
          scalar_t val = hvalues_diagblk(k);
          if (Kokkos::abs(col - val) != 0) {
            flag = false;
            break;
          }
        }
        if (flag == false) break;
      }
      if (flag == false) break;
      col_start += DiagBlks[i].numCols();
    }
    EXPECT_TRUE(flag);

    // Checking RCM
    if (!perm.empty()) {
      scalar_t one  = scalar_t(1.0);
      scalar_t zero = scalar_t(0.0);
      scalar_t mone = scalar_t(-1.0);
      for (int i = 0; i < nblocks; i++) {
        ValuesType In("In", DiagBlks[i].numRows());
        ValuesType Out("Out", DiagBlks[i].numRows());

        ValuesType_hm h_Out     = Kokkos::create_mirror_view(Out);
        ValuesType_hm h_Out_tmp = Kokkos::create_mirror(Out);

        Kokkos::deep_copy(In, one);

        auto h_perm = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), perm[i]);

        KokkosSparse::spmv("N", one, DiagBlks_rcm[i], In, zero, Out);

        Kokkos::deep_copy(h_Out_tmp, Out);
        for (lno_t ii = 0; ii < static_cast<lno_t>(DiagBlks[i].numRows()); ii++) {
          lno_t rcm_ii = h_perm(ii);
          h_Out(ii)    = h_Out_tmp(rcm_ii);
        }
        Kokkos::deep_copy(Out, h_Out);

        KokkosSparse::spmv("N", one, DiagBlks[i], In, mone, Out);

        double nrm_val = KokkosBlas::nrm2(Out);
        EXPECT_LE(nrm_val, 1e-9);
      }
    }
  }
}
}  // namespace Test

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_extract_diagonal_blocks() {
  for (int s = 1; s <= 8; s++) {
    Test::run_test_extract_diagonal_blocks<scalar_t, lno_t, size_type, device>(0, s);
    Test::run_test_extract_diagonal_blocks<scalar_t, lno_t, size_type, device>(153, s);
    Test::run_test_extract_diagonal_blocks<scalar_t, lno_t, size_type, device>(1553, s);
  }
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                                       \
  TEST_F(TestCategory, sparse##_##extract_diagonal_blocks##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_extract_diagonal_blocks<SCALAR, ORDINAL, OFFSET, DEVICE>();                                      \
  }

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
