// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Thyra testing tools
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"

// Thyra includes
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "tBlockedTpetraOperator.hpp"

#include "Trilinos_Util_CrsMatrixGallery.h"

#include "Teko_BlockedTpetraOperator.hpp"
#include "Teko_TpetraHelpers.hpp"

#include <random>

namespace Teko {
namespace Test {

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Thyra::createMember;
using Thyra::LinearOpBase;
using Thyra::LinearOpTester;
using Thyra::VectorBase;

void tBlockedTpetraOperator::buildBlockGIDs(std::vector<std::vector<GO>>& gids,
                                            const Tpetra::Map<LO, GO, NT>& map,
                                            bool singleBlock) const {
  LO numLocal = map.getLocalNumElements();
  LO numHalf  = numLocal / 2;
  numHalf += ((numHalf % 2 == 0) ? 0 : 1);

  gids.clear();
  gids.resize(singleBlock ? 1 : 3);

  std::vector<GO>& blk0 = gids[0];
  std::vector<GO>& blk1 = singleBlock ? gids[0] : gids[1];
  std::vector<GO>& blk2 = singleBlock ? gids[0] : gids[2];

  // loop over global IDs: treat first block as strided
  GO gid = -1;
  for (LO i = 0; i < numHalf; i += 2) {
    gid = map.getGlobalElement(i);
    blk0.push_back(gid);

    gid = map.getGlobalElement(i + 1);
    blk1.push_back(gid);
  }

  // loop over global IDs: treat remainder as contiguous
  for (LO i = numHalf; i < numLocal; i++) {
    gid = map.getGlobalElement(i);
    blk2.push_back(gid);
  }

  TEUCHOS_ASSERT(LO(singleBlock ? blk0.size() : blk0.size() + blk1.size() + blk2.size()) ==
                 numLocal);
}

void tBlockedTpetraOperator::initializeTest() { tolerance_ = 1e-14; }

int tBlockedTpetraOperator::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                                    int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tBlockedTpetraOperator";

  status = test_vector_constr(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"vector_constr\" ... PASSED", "   \"vector_constr\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_single_block(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"test_single_block\" ... PASSED",
                "   \"test_single_block\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_noncontig(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"test_noncontig\" ... PASSED", "   \"test_noncontig\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_reorder(verbosity, failstrm, 0);
  Teko_TEST_MSG(stdstrm, 1, "   \"reorder(flat reorder)\" ... PASSED",
                "   \"reorder(flat reorder)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_reorder(verbosity, failstrm, 1);
  Teko_TEST_MSG(stdstrm, 1, "   \"reorder(composite reorder = " << 1 << ")\" ... PASSED",
                "   \"reorder(composite reorder)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_reorder(verbosity, failstrm, 2);
  Teko_TEST_MSG(stdstrm, 1, "   \"reorder(composite reorder = " << 2 << ")\" ... PASSED",
                "   \"reorder(composite reorder)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG(failstrm, 0, "tBlockedTpetraOperator...PASSED",
                  "tBlockedTpetraOperator...FAILED");
  } else {  // Normal Operating Procedures (NOP)
    Teko_TEST_MSG(failstrm, 0, "...PASSED", "tBlockedTpetraOperator...FAILED");
  }

  return failcount;
}

bool tBlockedTpetraOperator::test_vector_constr(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  const Epetra_Comm& comm_epetra            = *GetComm();
  RCP<const Teuchos::Comm<int>> comm_tpetra = GetComm_tpetra();

  TEST_MSG("\n   tBlockedTpetraOperator::test_vector_constr: "
           << "Running on " << comm_epetra.NumProc() << " processors");

  // pick
  int nx = 5;  // * comm_epetra.NumProc();
  int ny = 5;  // * comm_epetra.NumProc();

  // create a big matrix to play with
  // note: this matrix is not really strided
  //       however, I just need a nontrivial
  //       matrix to play with
  Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d", comm_epetra, false);
  FGallery.Set("nx", nx);
  FGallery.Set("ny", ny);
  RCP<Epetra_CrsMatrix> epetraA = rcp(FGallery.GetMatrix(), false);
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> A =
      Teko::TpetraHelpers::nonConstEpetraCrsMatrixToTpetra(epetraA, comm_tpetra);
  ST beforeNorm = A->getFrobeniusNorm();

  int width = 3;
  Tpetra::MultiVector<ST, LO, GO, NT> x(A->getDomainMap(), width);
  Tpetra::MultiVector<ST, LO, GO, NT> ys(A->getRangeMap(), width);
  Tpetra::MultiVector<ST, LO, GO, NT> y(A->getRangeMap(), width);

  std::vector<std::vector<GO>> vars;
  buildBlockGIDs(vars, *A->getRowMap(), false);

  Teko::TpetraHelpers::BlockedTpetraOperator shell(vars, A);

  // test the operator against a lot of random vectors
  int numtests = 50;
  ST max       = 0.0;
  ST min       = 1.0;
  for (int i = 0; i < numtests; i++) {
    std::vector<ST> norm(width);
    std::vector<ST> rel(width);
    x.randomize();

    shell.apply(x, y);
    A->apply(x, ys);

    Tpetra::MultiVector<ST, LO, GO, NT> e(y, Teuchos::Copy);
    e.update(-1.0, ys, 1.0);
    e.norm2(Teuchos::ArrayView<ST>(norm));

    // compute relative error
    ys.norm2(Teuchos::ArrayView<ST>(rel));
    for (int j = 0; j < width; j++) {
      max = max > norm[j] / rel[j] ? max : norm[j] / rel[j];
      min = min < norm[j] / rel[j] ? min : norm[j] / rel[j];
    }
  }
  TEST_ASSERT(max >= min, "\n   tBlockedTpetraOperator::test_vector_constr: "
                              << toString(status) << ": "
                              << "sanity checked - " << max << " >= " << min);
  TEST_ASSERT(max <= tolerance_, "\n   tBlockedTpetraOperator::test_vector_constr: "
                                     << toString(status) << ": "
                                     << "testing tolerance over many matrix vector multiplies ( "
                                     << max << " <= "
                                     << "testing tolerance over many matrix vector multiplies ( "
                                     << max << " <= " << tolerance_ << " )");

  A->scale(2.0);  // double everything

  ST afterNorm = A->getFrobeniusNorm();
  TEST_ASSERT(beforeNorm != afterNorm, "\n   tBlockedTpetraOperator::test_vector_constr "
                                           << toString(status) << ": "
                                           << "verify matrix has been modified");

  shell.RebuildOps();

  // test the operator against a lot of random vectors
  numtests = 50;
  max      = 0.0;
  min      = 1.0;
  for (int i = 0; i < numtests; i++) {
    std::vector<ST> norm(width);
    std::vector<ST> rel(width);
    x.randomize();

    shell.apply(x, y);
    A->apply(x, ys);

    Tpetra::MultiVector<ST, LO, GO, NT> e(y, Teuchos::Copy);
    e.update(-1.0, ys, 1.0);
    e.norm2(Teuchos::ArrayView<ST>(norm));

    // compute relative error
    ys.norm2(Teuchos::ArrayView<ST>(rel));
    for (int j = 0; j < width; j++) {
      max = max > norm[j] / rel[j] ? max : norm[j] / rel[j];
      min = min < norm[j] / rel[j] ? min : norm[j] / rel[j];
    }
  }
  TEST_ASSERT(max >= min, "\n   tBlockedTpetraOperator::test_vector_constr (rebuild): "
                              << toString(status) << ": "
                              << "sanity checked - " << max << " >= " << min);
  TEST_ASSERT(max <= tolerance_, "\n   tBlockedTpetraOperator::test_vector_constr (rebuild): "
                                     << toString(status) << ": "
                                     << "testing tolerance over many matrix vector multiplies ( "
                                     << max << " <= " << tolerance_ << " )");

  return allPassed;
}

std::vector<GO> random_gids(const RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>>& contigMat) {
  std::random_device rd;
  std::mt19937 g(rd());

  const auto contigRowMap          = contigMat->getRowMap();
  const auto numGlobalElements     = contigRowMap->getGlobalNumElements();
  const auto numLargestPossibleGid = 10 * numGlobalElements;
  std::vector<GO> randomGids(numLargestPossibleGid);
  std::iota(randomGids.begin(), randomGids.end(), 0);
  std::shuffle(randomGids.begin(), randomGids.end(), g);

  const auto numLocalGids = contigRowMap->getLocalNumElements();
  const auto gidStart     = contigRowMap->getMinGlobalIndex();
  const auto gidEnd       = contigRowMap->getMaxGlobalIndex();

  size_t numGids = gidEnd - gidStart + 1;
  TEUCHOS_ASSERT(numGids == numLocalGids);

  std::vector<GO> localRandomGids(numLocalGids);
  std::copy(randomGids.begin() + gidStart, randomGids.begin() + gidStart + numGids,
            localRandomGids.begin());

  return localRandomGids;
}

RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> assemble_noncontig_matrix(
    const RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>>& contigMat) {
  const auto contigRowMap        = contigMat->getRowMap();
  const auto comm                = contigRowMap->getComm();
  const auto contigGlobalIndices = contigRowMap->getMyGlobalIndices();

  auto randomGlobalGids = random_gids(contigMat);

  decltype(contigGlobalIndices)::non_const_type noncontigGlobalIndices(
      "noncontig", contigGlobalIndices.extent(0));
  for (auto i = 0U; i < contigGlobalIndices.extent(0); ++i) {
    noncontigGlobalIndices(i) = randomGlobalGids[i];
  }

  const auto indexBase       = 0;
  const auto invalid         = Teuchos::OrdinalTraits<GO>::invalid();
  const auto noncontigRowMap = Teuchos::make_rcp<Tpetra::Map<LO, GO, NT>>(
      invalid,
      Teuchos::ArrayView<const GO>(noncontigGlobalIndices.data(), noncontigGlobalIndices.extent(0)),
      indexBase, comm);

  auto nrowsLocal = contigRowMap->getLocalNumElements();
  std::vector<size_t> numEntPerRow(nrowsLocal);
  for (size_t row = 0; row < nrowsLocal; ++row) {
    numEntPerRow[row] = contigMat->getNumEntriesInLocalRow(row);
  }

  auto noncontigMat = Teuchos::make_rcp<Tpetra::CrsMatrix<ST, LO, GO, NT>>(
      noncontigRowMap, Teuchos::ArrayView<const size_t>(numEntPerRow));

  noncontigMat->resumeFill();

  const auto globalMaxNumRowEntries = contigMat->getGlobalMaxNumRowEntries();
  Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_global_inds_host_view_type contigColumnIndices(
      "contigColumnIndices", globalMaxNumRowEntries);
  Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_global_inds_host_view_type noncontigColumnIndices(
      "noncontigColumnIndices", globalMaxNumRowEntries);
  Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_values_host_view_type columnValues(
      "columnValues", globalMaxNumRowEntries);

  for (size_t row = 0; row < nrowsLocal; ++row) {
    const auto contigGlobalRow = contigRowMap->getGlobalElement(row);
    size_t numEntries          = Teuchos::OrdinalTraits<size_t>::invalid();
    contigMat->getGlobalRowCopy(contigGlobalRow, contigColumnIndices, columnValues, numEntries);

    for (auto index = 0U; index < numEntries; ++index) {
      auto localCol                 = contigRowMap->getLocalElement(contigGlobalIndices(index));
      noncontigColumnIndices(index) = noncontigRowMap->getGlobalElement(localCol);
    }

    const auto noncontigGlobalRow = noncontigRowMap->getGlobalElement(row);
    noncontigMat->insertGlobalValues(
        noncontigGlobalRow, Teuchos::ArrayView<const GO>(noncontigColumnIndices.data(), numEntries),
        Teuchos::ArrayView<ST>(columnValues.data(), numEntries));
  }

  noncontigMat->fillComplete(noncontigRowMap, noncontigRowMap);

  return noncontigMat;
}

bool tBlockedTpetraOperator::test_single_block(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  const Epetra_Comm& comm_epetra            = *GetComm();
  RCP<const Teuchos::Comm<int>> comm_tpetra = GetComm_tpetra();

  TEST_MSG("\n   tBlockedTpetraOperator::test_single_block: "
           << "Running on " << comm_epetra.NumProc() << " processors");

  // pick
  int nx = 5;  // * comm_epetra.NumProc();
  int ny = 5;  // * comm_epetra.NumProc();

  // create a big matrix to play with
  // note: this matrix is not really strided
  //       however, I just need a nontrivial
  //       matrix to play with
  Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d", comm_epetra, false);
  FGallery.Set("nx", nx);
  FGallery.Set("ny", ny);
  RCP<Epetra_CrsMatrix> epetraA = rcp(FGallery.GetMatrix(), false);
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> contigA =
      Teko::TpetraHelpers::nonConstEpetraCrsMatrixToTpetra(epetraA, comm_tpetra);

  auto A = assemble_noncontig_matrix(contigA);

  ST beforeNorm = A->getFrobeniusNorm();

  int width = 3;
  Tpetra::MultiVector<ST, LO, GO, NT> x(A->getDomainMap(), width);
  Tpetra::MultiVector<ST, LO, GO, NT> ys(A->getRangeMap(), width);
  Tpetra::MultiVector<ST, LO, GO, NT> y(A->getRangeMap(), width);

  std::vector<std::vector<GO>> vars;
  buildBlockGIDs(vars, *A->getRowMap(), true);

  Teko::TpetraHelpers::BlockedTpetraOperator shell(vars, A);

  // test the operator against a lot of random vectors
  int numtests = 50;
  ST max       = 0.0;
  ST min       = 1.0;
  for (int i = 0; i < numtests; i++) {
    std::vector<ST> norm(width);
    std::vector<ST> rel(width);
    x.randomize();

    shell.apply(x, y);
    A->apply(x, ys);

    Tpetra::MultiVector<ST, LO, GO, NT> e(y, Teuchos::Copy);
    e.update(-1.0, ys, 1.0);
    e.norm2(Teuchos::ArrayView<ST>(norm));

    // compute relative error
    ys.norm2(Teuchos::ArrayView<ST>(rel));
    for (int j = 0; j < width; j++) {
      max = max > norm[j] / rel[j] ? max : norm[j] / rel[j];
      min = min < norm[j] / rel[j] ? min : norm[j] / rel[j];
    }
  }
  TEST_ASSERT(max >= min, "\n   tBlockedTpetraOperator::test_single_block: "
                              << toString(status) << ": "
                              << "sanity checked - " << max << " >= " << min);
  TEST_ASSERT(max <= tolerance_, "\n   tBlockedTpetraOperator::test_single_block: "
                                     << toString(status) << ": "
                                     << "testing tolerance over many matrix vector multiplies ( "
                                     << max << " <= "
                                     << "testing tolerance over many matrix vector multiplies ( "
                                     << max << " <= " << tolerance_ << " )");

  A->scale(2.0);  // double everything

  ST afterNorm = A->getFrobeniusNorm();
  TEST_ASSERT(beforeNorm != afterNorm, "\n   tBlockedTpetraOperator::test_single_block "
                                           << toString(status) << ": "
                                           << "verify matrix has been modified");

  shell.RebuildOps();

  // test the operator against a lot of random vectors
  numtests = 50;
  max      = 0.0;
  min      = 1.0;
  for (int i = 0; i < numtests; i++) {
    std::vector<ST> norm(width);
    std::vector<ST> rel(width);
    x.randomize();

    shell.apply(x, y);
    A->apply(x, ys);

    Tpetra::MultiVector<ST, LO, GO, NT> e(y, Teuchos::Copy);
    e.update(-1.0, ys, 1.0);
    e.norm2(Teuchos::ArrayView<ST>(norm));

    // compute relative error
    ys.norm2(Teuchos::ArrayView<ST>(rel));
    for (int j = 0; j < width; j++) {
      max = max > norm[j] / rel[j] ? max : norm[j] / rel[j];
      min = min < norm[j] / rel[j] ? min : norm[j] / rel[j];
    }
  }
  TEST_ASSERT(max >= min, "\n   tBlockedTpetraOperator::test_single_block (rebuild): "
                              << toString(status) << ": "
                              << "sanity checked - " << max << " >= " << min);
  TEST_ASSERT(max <= tolerance_, "\n   tBlockedTpetraOperator::test_single_block (rebuild): "
                                     << toString(status) << ": "
                                     << "testing tolerance over many matrix vector multiplies ( "
                                     << max << " <= " << tolerance_ << " )");

  return allPassed;
}

bool tBlockedTpetraOperator::test_noncontig(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  const Epetra_Comm& comm_epetra            = *GetComm();
  RCP<const Teuchos::Comm<int>> comm_tpetra = GetComm_tpetra();

  TEST_MSG("\n   tBlockedTpetraOperator::test_noncontig: "
           << "Running on " << comm_epetra.NumProc() << " processors");

  // pick
  int nx = 15;  // * comm_epetra.NumProc();
  int ny = 15;  // * comm_epetra.NumProc();

  // create a big matrix to play with
  // note: this matrix is not really strided
  //       however, I just need a nontrivial
  //       matrix to play with
  Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d", comm_epetra, false);
  FGallery.Set("nx", nx);
  FGallery.Set("ny", ny);
  RCP<Epetra_CrsMatrix> epetraA = rcp(FGallery.GetMatrix(), false);
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> contigA =
      Teko::TpetraHelpers::nonConstEpetraCrsMatrixToTpetra(epetraA, comm_tpetra);

  auto A = assemble_noncontig_matrix(contigA);

  ST beforeNorm = A->getFrobeniusNorm();

  int width = 3;
  Tpetra::MultiVector<ST, LO, GO, NT> x(A->getDomainMap(), width);
  Tpetra::MultiVector<ST, LO, GO, NT> ys(A->getRangeMap(), width);
  Tpetra::MultiVector<ST, LO, GO, NT> y(A->getRangeMap(), width);

  std::vector<std::vector<GO>> vars;
  buildBlockGIDs(vars, *A->getRowMap(), false);

  Teko::TpetraHelpers::BlockedTpetraOperator shell(vars, A);

  // test the operator against a lot of random vectors
  int numtests = 50;
  ST max       = 0.0;
  ST min       = 1.0;
  for (int i = 0; i < numtests; i++) {
    std::vector<ST> norm(width);
    std::vector<ST> rel(width);
    x.randomize();

    shell.apply(x, y);
    A->apply(x, ys);

    Tpetra::MultiVector<ST, LO, GO, NT> e(y, Teuchos::Copy);
    e.update(-1.0, ys, 1.0);
    e.norm2(Teuchos::ArrayView<ST>(norm));

    // compute relative error
    ys.norm2(Teuchos::ArrayView<ST>(rel));
    for (int j = 0; j < width; j++) {
      max = max > norm[j] / rel[j] ? max : norm[j] / rel[j];
      min = min < norm[j] / rel[j] ? min : norm[j] / rel[j];
    }
  }
  TEST_ASSERT(max >= min, "\n   tBlockedTpetraOperator::test_noncontig: "
                              << toString(status) << ": "
                              << "sanity checked - " << max << " >= " << min);
  TEST_ASSERT(max <= tolerance_, "\n   tBlockedTpetraOperator::test_noncontig: "
                                     << toString(status) << ": "
                                     << "testing tolerance over many matrix vector multiplies ( "
                                     << max << " <= "
                                     << "testing tolerance over many matrix vector multiplies ( "
                                     << max << " <= " << tolerance_ << " )");

  A->scale(2.0);  // double everything

  ST afterNorm = A->getFrobeniusNorm();
  TEST_ASSERT(beforeNorm != afterNorm, "\n   tBlockedTpetraOperator::test_noncontig "
                                           << toString(status) << ": "
                                           << "verify matrix has been modified");

  shell.RebuildOps();

  // test the operator against a lot of random vectors
  numtests = 50;
  max      = 0.0;
  min      = 1.0;
  for (int i = 0; i < numtests; i++) {
    std::vector<ST> norm(width);
    std::vector<ST> rel(width);
    x.randomize();

    shell.apply(x, y);
    A->apply(x, ys);

    Tpetra::MultiVector<ST, LO, GO, NT> e(y, Teuchos::Copy);
    e.update(-1.0, ys, 1.0);
    e.norm2(Teuchos::ArrayView<ST>(norm));

    // compute relative error
    ys.norm2(Teuchos::ArrayView<ST>(rel));
    for (int j = 0; j < width; j++) {
      max = max > norm[j] / rel[j] ? max : norm[j] / rel[j];
      min = min < norm[j] / rel[j] ? min : norm[j] / rel[j];
    }
  }
  TEST_ASSERT(max >= min, "\n   tBlockedTpetraOperator::test_noncontig (rebuild): "
                              << toString(status) << ": "
                              << "sanity checked - " << max << " >= " << min);
  TEST_ASSERT(max <= tolerance_, "\n   tBlockedTpetraOperator::test_noncontig (rebuild): "
                                     << toString(status) << ": "
                                     << "testing tolerance over many matrix vector multiplies ( "
                                     << max << " <= " << tolerance_ << " )");

  return allPassed;
}

bool tBlockedTpetraOperator::test_reorder(int verbosity, std::ostream& os, int total) {
  bool status    = false;
  bool allPassed = true;

  const Epetra_Comm& comm_epetra            = *GetComm();
  RCP<const Teuchos::Comm<int>> comm_tpetra = GetComm_tpetra();

  std::string tstr = total ? "(composite reorder)" : "(flat reorder)";

  TEST_MSG("\n   tBlockedTpetraOperator::test_reorder" << tstr << ": "
                                                       << "Running on " << comm_epetra.NumProc()
                                                       << " processors");

  // pick
  int nx = 5;  // 3 * 25 * comm_epetra.NumProc();
  int ny = 5;  // 3 * 50 * comm_epetra.NumProc();

  // create a big matrix to play with
  // note: this matrix is not really strided
  //       however, I just need a nontrivial
  //       matrix to play with
  Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d", comm_epetra, false);
  FGallery.Set("nx", nx);
  FGallery.Set("ny", ny);
  RCP<Epetra_CrsMatrix> epetraA = rcp(FGallery.GetMatrix(), false);
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> A =
      Teko::TpetraHelpers::nonConstEpetraCrsMatrixToTpetra(epetraA, comm_tpetra);

  int width = 3;
  Tpetra::MultiVector<ST, LO, GO, NT> x(A->getDomainMap(), width);
  Tpetra::MultiVector<ST, LO, GO, NT> yf(A->getRangeMap(), width);
  Tpetra::MultiVector<ST, LO, GO, NT> yr(A->getRangeMap(), width);

  std::vector<std::vector<GO>> vars;
  buildBlockGIDs(vars, *A->getRowMap(), false);

  Teko::TpetraHelpers::BlockedTpetraOperator flatShell(vars, A, "Af");
  Teko::TpetraHelpers::BlockedTpetraOperator reorderShell(vars, A, "Ar");

  Teko::BlockReorderManager brm;
  switch (total) {
    case 0:
      brm.SetNumBlocks(3);
      brm.SetBlock(0, 1);
      brm.SetBlock(1, 0);
      brm.SetBlock(2, 2);
      break;
    case 1:
      brm.SetNumBlocks(2);
      brm.SetBlock(0, 1);
      brm.GetBlock(1)->SetNumBlocks(2);
      brm.GetBlock(1)->SetBlock(0, 0);
      brm.GetBlock(1)->SetBlock(1, 2);
      break;
    case 2:
      brm.SetNumBlocks(2);
      brm.GetBlock(0)->SetNumBlocks(2);
      brm.GetBlock(0)->SetBlock(0, 0);
      brm.GetBlock(0)->SetBlock(1, 2);
      brm.SetBlock(1, 1);
      break;
  }
  reorderShell.Reorder(brm);
  TEST_MSG("\n   tBlockedTpetraOperator::test_reorder" << tstr << ": patern = " << brm.toString());

  TEST_MSG("\n   tBlockedTpetraOperator::test_reorder" << tstr << ":\n");
  TEST_MSG("\n      " << Teuchos::describe(*reorderShell.getThyraOp(), Teuchos::VERB_HIGH)
                      << std::endl);

  // test the operator against a lot of random vectors
  int numtests = 10;
  ST max       = 0.0;
  ST min       = 1.0;
  for (int i = 0; i < numtests; i++) {
    std::vector<ST> norm(width);
    std::vector<ST> rel(width);
    x.randomize();

    flatShell.apply(x, yf);
    reorderShell.apply(x, yr);

    Tpetra::MultiVector<ST, LO, GO, NT> e(yf, Teuchos::Copy);
    e.update(-1.0, yr, 1.0);
    e.norm2(Teuchos::ArrayView<ST>(norm));

    // compute relative error
    yf.norm2(Teuchos::ArrayView<ST>(rel));
    for (int j = 0; j < width; j++) {
      max = max > norm[j] / rel[j] ? max : norm[j] / rel[j];
      min = min < norm[j] / rel[j] ? min : norm[j] / rel[j];
    }
  }
  TEST_ASSERT(max >= min, "   tBlockedTpetraOperator::test_reorder"
                              << tstr << ": " << toString(status) << ": "
                              << "sanity checked - " << max << " >= " << min);
  TEST_ASSERT(max <= tolerance_, "   tBlockedTpetraOperator::test_reorder"
                                     << tstr << ": " << toString(status) << ": "
                                     << "testing tolerance over many matrix vector multiplies ( "
                                     << max << " <= " << tolerance_ << " )");

  return allPassed;
}

}  // namespace Test
}  // namespace Teko
