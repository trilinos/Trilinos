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
#include "Teuchos_CommHelpers.hpp"

#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "tBlockedTpetraOperator.hpp"

#include "Galeri_XpetraMaps.hpp"
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraParameters.hpp"

#include "Teko_BlockedTpetraOperator.hpp"
#include "Teko_TpetraHelpers.hpp"

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

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

namespace {

using ST = Teko::ST;
using LO = Teko::LO;
using GO = Teko::GO;
using NT = Teko::NT;

using map_t = Tpetra::Map<LO, GO, NT>;
using crs_t = Tpetra::CrsMatrix<ST, LO, GO, NT>;
using vec_t = Tpetra::Vector<ST, LO, GO, NT>;
using mv_t  = Tpetra::MultiVector<ST, LO, GO, NT>;

RCP<crs_t> buildRecirc2DMatrix(const RCP<const Teuchos::Comm<int>>& comm, GO nx, GO ny) {
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  galeriList.set("mx", comm->getSize());
  galeriList.set("my", 1);

  auto tMap = Galeri::Xpetra::CreateMap<LO, GO, map_t>("Cartesian2D", comm, galeriList);

  auto problem =
      Galeri::Xpetra::BuildProblem<ST, LO, GO, map_t, crs_t, mv_t>("Recirc2D", tMap, galeriList);

  return problem->BuildMatrix();
}

GO gcd_go(GO a, GO b) {
  if (a < 0) a = -a;
  if (b < 0) b = -b;

  while (b != 0) {
    const GO t = a % b;
    a          = b;
    b          = t;
  }

  return a;
}

GO chooseCoprimeMultiplier(GO n) {
  if (n <= 2) return 1;

  for (GO a = 3; a < n; a += 2) {
    if (gcd_go(a, n) == 1) return a;
  }

  return 1;
}

// Deterministic global permutation.
//
// This is used to create a globally noncontiguous map without assuming that
// each MPI rank owns a contiguous range of GIDs.  The permutation is bijective
// as long as the multiplier is coprime with numGlobalElements.
GO permute_gid(GO gid, GO numGlobalElements, GO indexBase) {
  if (numGlobalElements <= 1) return gid;

  const GO zeroBasedGid = gid - indexBase;

  const GO multiplier = chooseCoprimeMultiplier(numGlobalElements);
  const GO shift      = 1 % numGlobalElements;

  const GO permutedZeroBasedGid = (multiplier * zeroBasedGid + shift) % numGlobalElements;

  return indexBase + permutedZeroBasedGid;
}

std::vector<GO> permuted_local_gids(const RCP<const map_t>& contigRowMap) {
  const auto numLocalElements  = contigRowMap->getLocalNumElements();
  const auto numGlobalElements = contigRowMap->getGlobalNumElements();
  const auto indexBase         = contigRowMap->getIndexBase();

  std::vector<GO> localPermutedGids(numLocalElements);

  for (size_t localRow = 0; localRow < numLocalElements; ++localRow) {
    const GO oldGid             = contigRowMap->getGlobalElement(static_cast<LO>(localRow));
    localPermutedGids[localRow] = permute_gid(oldGid, numGlobalElements, indexBase);
  }

  return localPermutedGids;
}

RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> assemble_noncontig_matrix(
    const RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>>& contigMat) {
  const auto contigRowMap      = contigMat->getRowMap();
  const auto comm              = contigRowMap->getComm();
  const auto numGlobalElements = contigRowMap->getGlobalNumElements();
  const auto indexBase         = contigRowMap->getIndexBase();
  const auto invalid           = Teuchos::OrdinalTraits<GO>::invalid();

  // Build the new local GID list by applying a deterministic global
  // permutation to exactly the GIDs owned by this rank.  Do not assume that
  // minGID:maxGID is a contiguous local interval.
  auto localPermutedGids = permuted_local_gids(contigRowMap);

  const auto noncontigRowMap = Teuchos::make_rcp<Tpetra::Map<LO, GO, NT>>(
      invalid, Teuchos::ArrayView<const GO>(localPermutedGids.data(), localPermutedGids.size()),
      indexBase, comm);

  const auto nrowsLocal = contigRowMap->getLocalNumElements();

  std::vector<size_t> numEntPerRow(nrowsLocal);
  for (size_t localRow = 0; localRow < nrowsLocal; ++localRow) {
    numEntPerRow[localRow] = contigMat->getNumEntriesInLocalRow(localRow);
  }

  auto noncontigMat = Teuchos::make_rcp<Tpetra::CrsMatrix<ST, LO, GO, NT>>(
      noncontigRowMap, Teuchos::ArrayView<const size_t>(numEntPerRow));

  noncontigMat->resumeFill();

  const auto globalMaxNumRowEntries = contigMat->getGlobalMaxNumRowEntries();

  crs_t::nonconst_global_inds_host_view_type contigColumnIndices("contigColumnIndices",
                                                                 globalMaxNumRowEntries);
  crs_t::nonconst_global_inds_host_view_type noncontigColumnIndices("noncontigColumnIndices",
                                                                    globalMaxNumRowEntries);
  crs_t::nonconst_values_host_view_type columnValues("columnValues", globalMaxNumRowEntries);

  for (size_t localRow = 0; localRow < nrowsLocal; ++localRow) {
    const GO contigGlobalRow = contigRowMap->getGlobalElement(static_cast<LO>(localRow));

    size_t numEntries = Teuchos::OrdinalTraits<size_t>::invalid();

    contigMat->getGlobalRowCopy(contigGlobalRow, contigColumnIndices, columnValues, numEntries);

    for (size_t entry = 0; entry < numEntries; ++entry) {
      const GO contigGlobalCol      = contigColumnIndices(entry);
      noncontigColumnIndices(entry) = permute_gid(contigGlobalCol, numGlobalElements, indexBase);
    }

    const GO noncontigGlobalRow = permute_gid(contigGlobalRow, numGlobalElements, indexBase);

    noncontigMat->insertGlobalValues(
        noncontigGlobalRow, Teuchos::ArrayView<const GO>(noncontigColumnIndices.data(), numEntries),
        Teuchos::ArrayView<const ST>(columnValues.data(), numEntries));
  }

  noncontigMat->fillComplete(noncontigRowMap, noncontigRowMap);

  return noncontigMat;
}

// Running summary of the relative error ||test - ref|| / ||ref|| (per column)
// over many trials, along with enough detail about the worst column to
// diagnose a failure.
struct ApplyErrorStats {
  ST maxRelErr   = 0.0;
  bool allFinite = true;
  int numTrials  = 0;

  // details of the column with the largest relative error
  int worstTrial  = -1;
  int worstColumn = -1;
  ST worstErrNorm = 0.0;
  ST worstRefNorm = 0.0;

  // details of the entry with the largest absolute error in the worst column
  GO worstGid     = -1;
  int worstOwner  = -1;
  ST worstTestVal = 0.0;
  ST worstRefVal  = 0.0;
};

// Find the entry of column j with the largest absolute error.  Collective.
void locateWorstEntry(const mv_t& e, const mv_t& test, const mv_t& ref, int j,
                      ApplyErrorStats& stats) {
  const auto comm   = ref.getMap()->getComm();
  const int myRank  = comm->getRank();
  const auto eH     = e.getLocalViewHost(Tpetra::Access::ReadOnly);
  const auto testH  = test.getLocalViewHost(Tpetra::Access::ReadOnly);
  const auto refH   = ref.getLocalViewHost(Tpetra::Access::ReadOnly);
  const auto inf    = std::numeric_limits<ST>::infinity();
  const LO numLocal = static_cast<LO>(ref.getLocalLength());

  ST localMax = -1.0;
  LO localIdx = -1;
  for (LO i = 0; i < numLocal; i++) {
    const ST absErr = std::isfinite(eH(i, j)) ? std::abs(eH(i, j)) : inf;
    if (absErr > localMax) {
      localMax = absErr;
      localIdx = i;
    }
  }

  ST globalMax = -1.0;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, localMax, Teuchos::outArg(globalMax));

  const int candidate = (localIdx >= 0 && localMax == globalMax) ? myRank : comm->getSize();
  int owner           = comm->getSize();
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, candidate, Teuchos::outArg(owner));
  if (owner >= comm->getSize()) return;  // no rows anywhere

  GO gid    = -1;
  ST vals[] = {0.0, 0.0};
  if (myRank == owner) {
    gid     = ref.getMap()->getGlobalElement(localIdx);
    vals[0] = testH(localIdx, j);
    vals[1] = refH(localIdx, j);
  }
  Teuchos::broadcast(*comm, owner, Teuchos::outArg(gid));
  Teuchos::broadcast(*comm, owner, 2, vals);

  stats.worstGid     = gid;
  stats.worstOwner   = owner;
  stats.worstTestVal = vals[0];
  stats.worstRefVal  = vals[1];
}

// Accumulate the per-column relative error of test against ref.  Collective.
void accumulateApplyError(const mv_t& test, const mv_t& ref, int trial, ApplyErrorStats& stats) {
  const int width = static_cast<int>(ref.getNumVectors());
  std::vector<ST> errNorm(width);
  std::vector<ST> refNorm(width);

  mv_t e(test, Teuchos::Copy);
  e.update(-1.0, ref, 1.0);
  e.norm2(Teuchos::ArrayView<ST>(errNorm));
  ref.norm2(Teuchos::ArrayView<ST>(refNorm));

  stats.numTrials++;
  for (int j = 0; j < width; j++) {
    // fall back to the absolute error if the reference is identically zero
    const ST relErr   = refNorm[j] > 0.0 ? errNorm[j] / refNorm[j] : errNorm[j];
    const bool finite = std::isfinite(relErr);

    // norms are global, so every rank takes the same branch here
    const bool firstNonFinite = !finite && stats.allFinite;
    if (firstNonFinite || (finite && stats.allFinite && relErr > stats.maxRelErr)) {
      stats.maxRelErr    = relErr;
      stats.worstTrial   = trial;
      stats.worstColumn  = j;
      stats.worstErrNorm = errNorm[j];
      stats.worstRefNorm = refNorm[j];
      locateWorstEntry(e, test, ref, j, stats);
    }
    stats.allFinite &= finite;
  }
}

bool passed(const ApplyErrorStats& stats, ST tolerance) {
  return stats.allFinite && stats.maxRelErr <= tolerance;
}

std::string describe(const ApplyErrorStats& stats, ST tolerance) {
  std::ostringstream ss;
  ss << std::scientific << std::setprecision(16);
  ss << "max relative error = " << stats.maxRelErr << " (tolerance = " << tolerance << ") over "
     << stats.numTrials << " trials";
  if (!stats.allFinite) ss << " [NON-FINITE VALUES DETECTED]";
  if (stats.worstTrial >= 0) {
    ss << "\n      worst column: trial " << stats.worstTrial << ", column " << stats.worstColumn
       << ", ||test - ref|| = " << stats.worstErrNorm << ", ||ref|| = " << stats.worstRefNorm;
  }
  if (stats.worstOwner >= 0) {
    ss << "\n      worst entry:  GID " << stats.worstGid << " (rank " << stats.worstOwner
       << "), test = " << stats.worstTestVal << ", ref = " << stats.worstRefVal
       << ", test - ref = " << (stats.worstTestVal - stats.worstRefVal);
  }
  return ss.str();
}

}  // namespace

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

// The blocked and reordered operators sum the block contributions in a different
// order than a monolithic SpMV, and the SpMV kernel itself (e.g. cuSPARSE) need not be
// bitwise reproducible, so only agreement to within a modest multiple of machine
// epsilon can be expected.
void tBlockedTpetraOperator::initializeTest() { tolerance_ = 1e-12; }

int tBlockedTpetraOperator::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                                    int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tBlockedTpetraOperator";

  status = test_vector_constr(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"vector_constr\" ... PASSED",
                       "   \"vector_constr\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_single_block(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"test_single_block\" ... PASSED",
                       "   \"test_single_block\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_noncontig(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"test_noncontig\" ... PASSED",
                       "   \"test_noncontig\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_reorder(verbosity, failstrm, 0);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"reorder(flat reorder)\" ... PASSED",
                       "   \"reorder(flat reorder)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_reorder(verbosity, failstrm, 1);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"reorder(composite reorder = " << 1 << ")\" ... PASSED",
                       "   \"reorder(composite reorder = " << 1 << ")\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_reorder(verbosity, failstrm, 2);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"reorder(composite reorder = " << 2 << ")\" ... PASSED",
                       "   \"reorder(composite reorder = " << 2 << ")\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG_tpetra(failstrm, 0, "tBlockedTpetraOperator...PASSED",
                         "tBlockedTpetraOperator...FAILED");
  } else {  // Normal Operating Procedures (NOP)
    Teko_TEST_MSG_tpetra(failstrm, 0, "...PASSED", "tBlockedTpetraOperator...FAILED");
  }

  return failcount;
}

bool tBlockedTpetraOperator::test_vector_constr(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  RCP<const Teuchos::Comm<int>> comm_tpetra = GetComm_tpetra();

  TEST_MSG("\n   tBlockedTpetraOperator::test_vector_constr: "
           << "Running on " << comm_tpetra->getSize() << " processors");

  // pick
  GO nx = 5;  // * comm_epetra.NumProc();
  GO ny = 5;  // * comm_epetra.NumProc();

  // create a big matrix to play with
  // note: this matrix is not really strided
  //       however, I just need a nontrivial
  //       matrix to play with
  auto A        = buildRecirc2DMatrix(comm_tpetra, nx, ny);
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
  ApplyErrorStats stats;
  for (int i = 0; i < numtests; i++) {
    x.randomize();

    shell.apply(x, y);
    A->apply(x, ys);

    accumulateApplyError(y, ys, i, stats);
  }
  TEST_ASSERT(passed(stats, tolerance_), "\n   tBlockedTpetraOperator::test_vector_constr: "
                                             << "testing tolerance over many matrix vector "
                                                "multiplies: "
                                             << describe(stats, tolerance_));

  A->scale(2.0);  // double everything

  ST afterNorm = A->getFrobeniusNorm();
  TEST_ASSERT(beforeNorm != afterNorm, "\n   tBlockedTpetraOperator::test_vector_constr "
                                           << toString(status) << ": "
                                           << "verify matrix has been modified");

  shell.RebuildOps();

  // test the operator against a lot of random vectors
  numtests = 50;
  stats    = ApplyErrorStats();
  for (int i = 0; i < numtests; i++) {
    x.randomize();

    shell.apply(x, y);
    A->apply(x, ys);

    accumulateApplyError(y, ys, i, stats);
  }
  TEST_ASSERT(passed(stats, tolerance_),
              "\n   tBlockedTpetraOperator::test_vector_constr (rebuild): "
                  << "testing tolerance over many matrix vector "
                     "multiplies: "
                  << describe(stats, tolerance_));

  return allPassed;
}

bool tBlockedTpetraOperator::test_single_block(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  RCP<const Teuchos::Comm<int>> comm_tpetra = GetComm_tpetra();

  TEST_MSG("\n   tBlockedTpetraOperator::test_single_block: "
           << "Running on " << comm_tpetra->getSize() << " processors");

  // pick
  GO nx = 5;  // * comm_epetra.NumProc();
  GO ny = 5;  // * comm_epetra.NumProc();

  // create a big matrix to play with
  // note: this matrix is not really strided
  //       however, I just need a nontrivial
  //       matrix to play with
  auto contigA = buildRecirc2DMatrix(comm_tpetra, nx, ny);

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
  ApplyErrorStats stats;
  for (int i = 0; i < numtests; i++) {
    x.randomize();

    shell.apply(x, y);
    A->apply(x, ys);

    accumulateApplyError(y, ys, i, stats);
  }
  TEST_ASSERT(passed(stats, tolerance_), "\n   tBlockedTpetraOperator::test_single_block: "
                                             << "testing tolerance over many matrix vector "
                                                "multiplies: "
                                             << describe(stats, tolerance_));

  A->scale(2.0);  // double everything

  ST afterNorm = A->getFrobeniusNorm();
  TEST_ASSERT(beforeNorm != afterNorm, "\n   tBlockedTpetraOperator::test_single_block "
                                           << toString(status) << ": "
                                           << "verify matrix has been modified");

  shell.RebuildOps();

  // test the operator against a lot of random vectors
  numtests = 50;
  stats    = ApplyErrorStats();
  for (int i = 0; i < numtests; i++) {
    x.randomize();

    shell.apply(x, y);
    A->apply(x, ys);

    accumulateApplyError(y, ys, i, stats);
  }
  TEST_ASSERT(passed(stats, tolerance_),
              "\n   tBlockedTpetraOperator::test_single_block (rebuild): "
                  << "testing tolerance over many matrix vector "
                     "multiplies: "
                  << describe(stats, tolerance_));

  return allPassed;
}

bool tBlockedTpetraOperator::test_noncontig(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  RCP<const Teuchos::Comm<int>> comm_tpetra = GetComm_tpetra();

  TEST_MSG("\n   tBlockedTpetraOperator::test_noncontig: "
           << "Running on " << comm_tpetra->getSize() << " processors");

  // pick
  GO nx = 15;  // * comm_epetra.NumProc();
  GO ny = 15;  // * comm_epetra.NumProc();

  // create a big matrix to play with
  // note: this matrix is not really strided
  //       however, I just need a nontrivial
  //       matrix to play with
  auto contigA = buildRecirc2DMatrix(comm_tpetra, nx, ny);

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
  ApplyErrorStats stats;
  for (int i = 0; i < numtests; i++) {
    x.randomize();

    shell.apply(x, y);
    A->apply(x, ys);

    accumulateApplyError(y, ys, i, stats);
  }
  TEST_ASSERT(passed(stats, tolerance_), "\n   tBlockedTpetraOperator::test_noncontig: "
                                             << "testing tolerance over many matrix vector "
                                                "multiplies: "
                                             << describe(stats, tolerance_));

  A->scale(2.0);  // double everything

  ST afterNorm = A->getFrobeniusNorm();
  TEST_ASSERT(beforeNorm != afterNorm, "\n   tBlockedTpetraOperator::test_noncontig "
                                           << toString(status) << ": "
                                           << "verify matrix has been modified");

  shell.RebuildOps();

  // test the operator against a lot of random vectors
  numtests = 50;
  stats    = ApplyErrorStats();
  for (int i = 0; i < numtests; i++) {
    x.randomize();

    shell.apply(x, y);
    A->apply(x, ys);

    accumulateApplyError(y, ys, i, stats);
  }
  TEST_ASSERT(passed(stats, tolerance_), "\n   tBlockedTpetraOperator::test_noncontig (rebuild): "
                                             << "testing tolerance over many matrix vector "
                                                "multiplies: "
                                             << describe(stats, tolerance_));

  return allPassed;
}

bool tBlockedTpetraOperator::test_reorder(int verbosity, std::ostream& os, int total) {
  bool status    = false;
  bool allPassed = true;

  RCP<const Teuchos::Comm<int>> comm_tpetra = GetComm_tpetra();

  std::string tstr =
      total ? "(composite reorder = " + std::to_string(total) + ")" : "(flat reorder)";

  TEST_MSG("\n   tBlockedTpetraOperator::test_reorder" << tstr << ": "
                                                       << "Running on " << comm_tpetra->getSize()
                                                       << " processors");

  // pick
  GO nx = 5;  // 3 * 25 * comm_epetra.NumProc();
  GO ny = 5;  // 3 * 50 * comm_epetra.NumProc();

  // create a big matrix to play with
  // note: this matrix is not really strided
  //       however, I just need a nontrivial
  //       matrix to play with
  auto A = buildRecirc2DMatrix(comm_tpetra, nx, ny);

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
  Tpetra::MultiVector<ST, LO, GO, NT> ys(A->getRangeMap(), width);

  // second applies of the same input, used to detect run-to-run nondeterminism
  Tpetra::MultiVector<ST, LO, GO, NT> yf2(A->getRangeMap(), width);
  Tpetra::MultiVector<ST, LO, GO, NT> ys2(A->getRangeMap(), width);

  int numtests = 10;
  ApplyErrorStats reorderVsFlat;
  ApplyErrorStats flatVsMatrix;
  ApplyErrorStats reorderVsMatrix;
  ApplyErrorStats flatRepeat;
  ApplyErrorStats matrixRepeat;
  for (int i = 0; i < numtests; i++) {
    x.randomize();

    flatShell.apply(x, yf);
    reorderShell.apply(x, yr);
    A->apply(x, ys);

    flatShell.apply(x, yf2);
    A->apply(x, ys2);

    accumulateApplyError(yr, yf, i, reorderVsFlat);
    accumulateApplyError(yf, ys, i, flatVsMatrix);
    accumulateApplyError(yr, ys, i, reorderVsMatrix);
    accumulateApplyError(yf2, yf, i, flatRepeat);
    accumulateApplyError(ys2, ys, i, matrixRepeat);
  }

  // On failure, report all comparisons so it is clear which operator deviates.  The repeat
  // comparisons are informational only (a nonzero value means the SpMV, or the blocked apply
  // built on it, is not bitwise reproducible from run to run); they are not asserted.
  const bool allOk = passed(reorderVsFlat, tolerance_) && passed(flatVsMatrix, tolerance_) &&
                     passed(reorderVsMatrix, tolerance_);
  TEST_ASSERT(allOk, "   tBlockedTpetraOperator::test_reorder"
                         << tstr << " (pattern = " << brm.toString() << "): "
                         << "testing tolerance over many matrix vector multiplies"
                         << "\n    reordered vs. flat:   " << describe(reorderVsFlat, tolerance_)
                         << "\n    flat vs. matrix:      " << describe(flatVsMatrix, tolerance_)
                         << "\n    reordered vs. matrix: " << describe(reorderVsMatrix, tolerance_)
                         << "\n    (informational) flat vs. flat, repeated apply: "
                         << describe(flatRepeat, tolerance_)
                         << "\n    (informational) matrix vs. matrix, repeated apply: "
                         << describe(matrixRepeat, tolerance_));

  return allPassed;
}

}  // namespace Test
}  // namespace Teko
