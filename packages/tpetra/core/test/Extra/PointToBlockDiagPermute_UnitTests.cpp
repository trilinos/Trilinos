// @HEADER
// *****************************************************************************
//      TpetraExt / Tpetra unit test
//
// Copyright 2025 NTESS
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Array.hpp>

#include <iostream>
#include <vector>
#include <algorithm>

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
#include "TpetraExt_PointToBlockDiagPermute_decl.hpp"
#include "TpetraCore_ETIHelperMacros.h"

using Teuchos::RCP;
using Teuchos::rcp;

namespace {

template <class ST, class LO, class GO, class NT>
RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> >
buildBlockDiagonalMatrix(const RCP<const Teuchos::Comm<int> >& comm,
                         int blockSize,
                         int numLocalBlocks) {
  using map_t = Tpetra::Map<LO, GO, NT>;
  using crs_t = Tpetra::CrsMatrix<ST, LO, GO, NT>;

  const GO numGlobal = static_cast<GO>(blockSize * numLocalBlocks * comm->getSize());
  auto map           = rcp(new map_t(numGlobal, 0, comm));

  auto A = rcp(new crs_t(map, 3));

  for (LO lrow = 0; lrow < static_cast<LO>(map->getLocalNumElements()); ++lrow) {
    GO gid = map->getGlobalElement(lrow);

    std::vector<GO> cols;
    std::vector<ST> vals;

    GO gidm1 = gid - 1;
    GO gidp1 = gid + 1;
    ST diag  = static_cast<ST>(10.0 + static_cast<double>(gid));
    ST left  = static_cast<ST>(-1.0);
    ST right = static_cast<ST>(-2.0);

    if ((lrow % blockSize) > 0) {
      cols.push_back(gidm1);
      vals.push_back(left);
    }

    cols.push_back(gid);
    vals.push_back(diag);

    if ((lrow % blockSize) < blockSize - 1) {
      cols.push_back(gidp1);
      vals.push_back(right);
    }

    A->insertGlobalValues(gid, Teuchos::ArrayView<const GO>(cols),
                          Teuchos::ArrayView<const ST>(vals));
  }

  A->fillComplete();
  return A;
}

template <class ST, class LO, class GO, class NT>
RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> >
buildExpectedReversedBlockDiagonal(const RCP<const Teuchos::Comm<int> >& comm,
                                   int blockSize,
                                   int numLocalBlocks) {
  using map_t = Tpetra::Map<LO, GO, NT>;
  using crs_t = Tpetra::CrsMatrix<ST, LO, GO, NT>;

  const GO numGlobal = static_cast<GO>(blockSize * numLocalBlocks * comm->getSize());
  auto map           = rcp(new map_t(numGlobal, 0, comm));

  auto A = rcp(new crs_t(map, blockSize));

  for (int blk = 0; blk < numLocalBlocks; ++blk) {
    std::vector<GO> gids(blockSize);
    for (int j = 0; j < blockSize; ++j) {
      gids[j] = map->getGlobalElement(static_cast<LO>(blk * blockSize + (blockSize - j - 1)));
    }

    for (int r = 0; r < blockSize; ++r) {
      GO rowGid = gids[r];
      std::vector<GO> cols;
      std::vector<ST> vals;

      if (r > 0) {
        cols.push_back(gids[r - 1]);
        vals.push_back(static_cast<ST>(-2.0));
      }

      cols.push_back(gids[r]);
      vals.push_back(static_cast<ST>(10.0 + static_cast<double>(rowGid)));

      if (r < blockSize - 1) {
        cols.push_back(gids[r + 1]);
        vals.push_back(static_cast<ST>(-1.0));
      }

      A->insertGlobalValues(rowGid, Teuchos::ArrayView<const GO>(cols),
                            Teuchos::ArrayView<const ST>(vals));
    }
  }

  A->fillComplete();
  return A;
}

template <class ST, class LO, class GO, class NT>
double diffNorm(const Tpetra::CrsMatrix<ST, LO, GO, NT>& A,
                const Tpetra::CrsMatrix<ST, LO, GO, NT>& B) {
  using crs_t = Tpetra::CrsMatrix<ST, LO, GO, NT>;

  auto map  = A.getRowMap();
  auto Diff = rcp(new crs_t(map, A.getGlobalMaxNumRowEntries() + B.getGlobalMaxNumRowEntries()));

  Tpetra::MatrixMatrix::Add(A, false, static_cast<ST>(1.0),
                            B, false, static_cast<ST>(-1.0), Diff);
  Diff->fillComplete(A.getDomainMap(), A.getRangeMap());

  return Diff->getFrobeniusNorm();
}

}  // namespace

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TpetraExt_PointToBlockDiagPermute, contiguous_blocks, ST, LO, GO, NT) {
  auto comm = Tpetra::getDefaultComm();

  const int blockSize      = 3;
  const int numLocalBlocks = 2;

  auto A = buildBlockDiagonalMatrix<ST, LO, GO, NT>(comm, blockSize, numLocalBlocks);

  Teuchos::ParameterList pl;
  pl.set("contiguous block size", blockSize);
  pl.set("number of local blocks", numLocalBlocks);

  Tpetra::Ext::PointToBlockDiagPermute<ST, LO, GO, NT> perm(*A);
  TEST_EQUALITY(perm.setParameters(pl), 0);
  TEST_EQUALITY(perm.compute(), 0);

  auto P = perm.createCrsMatrix();
  TEST_ASSERT(!P.is_null());

  const double norm = diffNorm<ST, LO, GO, NT>(*A, *P);
  TEST_COMPARE(norm, <, 1.0e-14);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TpetraExt_PointToBlockDiagPermute, explicit_reversed_gids, ST, LO, GO, NT) {
  auto comm = Tpetra::getDefaultComm();

  const int blockSize      = 3;
  const int numLocalBlocks = 2;

  auto A   = buildBlockDiagonalMatrix<ST, LO, GO, NT>(comm, blockSize, numLocalBlocks);
  auto map = A->getRowMap();

  const LO nLocal = static_cast<LO>(map->getLocalNumElements());

  Teuchos::Array<int> blockStarts(numLocalBlocks + 1);
  Teuchos::Array<GO> blockGids(nLocal);

  for (int blk = 0; blk < numLocalBlocks; ++blk) {
    blockStarts[blk] = blk * blockSize;
    for (int j = 0; j < blockSize; ++j) {
      blockGids[blk * blockSize + j] =
          map->getGlobalElement(static_cast<LO>(blk * blockSize + (blockSize - j - 1)));
    }
  }
  blockStarts[numLocalBlocks] = nLocal;

  Teuchos::ParameterList pl;
  pl.set("number of local blocks", numLocalBlocks);
  pl.set("block start index", blockStarts);
  pl.set("block entry gids", blockGids);

  Tpetra::Ext::PointToBlockDiagPermute<ST, LO, GO, NT> perm(*A);
  TEST_EQUALITY(perm.setParameters(pl), 0);
  TEST_EQUALITY(perm.compute(), 0);

  auto P = perm.createCrsMatrix();
  TEST_ASSERT(!P.is_null());

  auto E = buildExpectedReversedBlockDiagonal<ST, LO, GO, NT>(comm, blockSize, numLocalBlocks);

  const double norm = diffNorm<ST, LO, GO, NT>(*E, *P);
  TEST_COMPARE(norm, <, 1.0e-14);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TpetraExt_PointToBlockDiagPermute, explicit_identity_order_gids, ST, LO, GO, NT) {
  auto comm = Tpetra::getDefaultComm();

  const int blockSize      = 3;
  const int numLocalBlocks = 2;

  auto A   = buildBlockDiagonalMatrix<ST, LO, GO, NT>(comm, blockSize, numLocalBlocks);
  auto map = A->getRowMap();

  const LO nLocal = static_cast<LO>(map->getLocalNumElements());

  Teuchos::Array<int> blockStarts(numLocalBlocks + 1);
  Teuchos::Array<GO> blockGids(nLocal);

  for (int blk = 0; blk < numLocalBlocks; ++blk) {
    blockStarts[blk] = blk * blockSize;
    for (int j = 0; j < blockSize; ++j) {
      blockGids[blk * blockSize + j] =
          map->getGlobalElement(static_cast<LO>(blk * blockSize + j));
    }
  }
  blockStarts[numLocalBlocks] = nLocal;

  Teuchos::ParameterList pl;
  pl.set("number of local blocks", numLocalBlocks);
  pl.set("block start index", blockStarts);
  pl.set("block entry gids", blockGids);

  Tpetra::Ext::PointToBlockDiagPermute<ST, LO, GO, NT> perm(*A);
  TEST_EQUALITY(perm.setParameters(pl), 0);
  TEST_EQUALITY(perm.compute(), 0);

  auto P = perm.createCrsMatrix();
  TEST_ASSERT(!P.is_null());

  const double norm = diffNorm<ST, LO, GO, NT>(*A, *P);
  TEST_COMPARE(norm, <, 1.0e-14);
}

#define UNIT_TEST_GROUP(ST, LO, GO, NT)                                                                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TpetraExt_PointToBlockDiagPermute, contiguous_blocks, ST, LO, GO, NT)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TpetraExt_PointToBlockDiagPermute, explicit_reversed_gids, ST, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TpetraExt_PointToBlockDiagPermute, explicit_identity_order_gids, ST, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP)
