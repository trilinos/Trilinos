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
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
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
#include "Tpetra_Details_Behavior.hpp"

#include "KokkosKernels_config.h"

#include "tBlockedTpetraOperator.hpp"

#include "Galeri_XpetraMaps.hpp"
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraParameters.hpp"

#include "Teko_BlockedTpetraOperator.hpp"
#include "Teko_TpetraHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
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

  // relative error of every column of every trial, in order
  std::vector<ST> history;
  int width = 0;
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
  stats.width = width;
  for (int j = 0; j < width; j++) {
    // fall back to the absolute error if the reference is identically zero
    const ST relErr   = refNorm[j] > 0.0 ? errNorm[j] / refNorm[j] : errNorm[j];
    const bool finite = std::isfinite(relErr);
    stats.history.push_back(relErr);

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
  if (!passed(stats, tolerance) && stats.width > 0) {
    ss << std::setprecision(3) << "\n      relative error per trial (one value per column):";
    for (size_t k = 0; k < stats.history.size(); k++) {
      if (k % stats.width == 0) ss << "\n        trial " << k / stats.width << ":";
      ss << " " << stats.history[k];
    }
  }
  return ss.str();
}

// The test driver sends failure messages to a log file (failure.log by default), which CI
// does not capture.  Also echo failures to stdout on rank 0 so they show up in the CI log.
void echoFailure(const Teuchos::Comm<int>& comm, bool ok, const std::string& msg) {
  if (!ok && comm.getRank() == 0) std::cout << msg << std::endl;
}

// Gather onto rank 0 the partitioning (local rows and block sizes per rank) and every entry
// where any two of yf (flat), yr (reordered), ys (monolithic) differ by more than
// tolerance * ||ys(:, j)||.  Returns the formatted report on rank 0 and an empty string
// elsewhere.  Collective.
std::string describeLargeEntries(const mv_t& yf, const mv_t& yr, const mv_t& ys,
                                 const std::vector<std::vector<GO>>& vars, ST tolerance,
                                 int maxEntries = 100) {
  const auto map     = ys.getMap();
  const auto comm    = map->getComm();
  const int myRank   = comm->getRank();
  const int numRanks = comm->getSize();
  const int width    = static_cast<int>(ys.getNumVectors());
  const int numBlks  = static_cast<int>(vars.size());
  const LO numLocal  = static_cast<LO>(ys.getLocalLength());

  std::vector<ST> refNorm(width);
  ys.norm2(Teuchos::ArrayView<ST>(refNorm));

  std::map<GO, int> blockOf;
  for (int b = 0; b < numBlks; b++)
    for (const GO gid : vars[b]) blockOf[gid] = b;

  const auto yfH = yf.getLocalViewHost(Tpetra::Access::ReadOnly);
  const auto yrH = yr.getLocalViewHost(Tpetra::Access::ReadOnly);
  const auto ysH = ys.getLocalViewHost(Tpetra::Access::ReadOnly);

  // each record: gid, block, column, flat, reordered, monolithic
  const int recLen = 6;
  std::vector<double> records;
  for (LO i = 0; i < numLocal; i++) {
    const GO gid  = map->getGlobalElement(i);
    const auto it = blockOf.find(gid);
    for (int j = 0; j < width; j++) {
      const ST scale = tolerance * (refNorm[j] > 0.0 ? refNorm[j] : 1.0);
      const ST f = yfH(i, j), r = yrH(i, j), s = ysH(i, j);
      const ST d = std::max({std::abs(r - f), std::abs(f - s), std::abs(r - s)});
      if (!std::isfinite(d) || d > scale) {
        records.insert(records.end(), {static_cast<double>(gid),
                                       static_cast<double>(it == blockOf.end() ? -1 : it->second),
                                       static_cast<double>(j), f, r, s});
      }
    }
  }

  // per rank: number of local rows, size of each block, number of doubles in records
  std::vector<int> layout(numBlks + 2);
  layout[0] = numLocal;
  for (int b = 0; b < numBlks; b++) layout[b + 1] = static_cast<int>(vars[b].size());
  layout[numBlks + 1] = static_cast<int>(records.size());

  const int layoutLen = numBlks + 2;
  std::vector<int> allLayout(myRank == 0 ? layoutLen * numRanks : 0);
  Teuchos::gather<int, int>(layout.data(), layoutLen, allLayout.data(), layoutLen, 0, *comm);

  std::vector<int> counts(numRanks, 0), displs(numRanks, 0);
  if (myRank == 0) {
    for (int p = 0; p < numRanks; p++) {
      counts[p] = allLayout[p * layoutLen + numBlks + 1];
      displs[p] = (p == 0) ? 0 : displs[p - 1] + counts[p - 1];
    }
  }
  std::vector<double> allRecords(myRank == 0 ? displs.back() + counts.back() : 0);
  Teuchos::gatherv<int, double>(records.data(), static_cast<int>(records.size()), allRecords.data(),
                                counts.data(), displs.data(), 0, *comm);

  if (myRank != 0) return "";

  std::ostringstream ss;
  ss << std::scientific << std::setprecision(16);
  ss << "\n      partitioning:";
  for (int p = 0; p < numRanks; p++) {
    ss << "\n        rank " << p << ": " << allLayout[p * layoutLen]
       << " local rows, block sizes [";
    for (int b = 0; b < numBlks; b++) ss << " " << allLayout[p * layoutLen + b + 1];
    ss << " ], " << counts[p] / recLen << " entries out of tolerance";
  }

  const int numRecords = static_cast<int>(allRecords.size()) / recLen;
  ss << "\n      entries where flat/reordered/matrix disagree by more than tolerance * ||matrix||"
     << " (" << numRecords << " total, showing up to " << maxEntries << "):";
  int shown = 0;
  for (int p = 0; p < numRanks && shown < maxEntries; p++) {
    for (int k = displs[p]; k < displs[p] + counts[p] && shown < maxEntries; k += recLen, shown++) {
      const double* rec = &allRecords[k];
      ss << "\n        rank " << p << " GID " << static_cast<GO>(rec[0]) << " block "
         << static_cast<int>(rec[1]) << " col " << static_cast<int>(rec[2]) << ": flat = " << rec[3]
         << ", reordered = " << rec[4] << ", matrix = " << rec[5]
         << ", reordered - flat = " << (rec[4] - rec[3])
         << ", flat - matrix = " << (rec[3] - rec[5]);
    }
  }
  return ss.str();
}

// ---------------------------------------------------------------------------------------
// Diagnostics run by test_reorder on failure.  Each one isolates one stage of the blocked
// apply (block assembly, vector plumbing, per-block SpMV, synchronization).
// ---------------------------------------------------------------------------------------

std::string describeConfiguration(const Teuchos::Comm<int>& comm) {
  using exec_space = NT::execution_space;
  using mem_space  = NT::memory_space;

  const char* launchBlocking = std::getenv("CUDA_LAUNCH_BLOCKING");

  std::ostringstream ss;
  ss << "\n      configuration: " << comm.getSize() << " ranks, Kokkos " << KOKKOS_VERSION
     << ", Kokkos Kernels " << KOKKOSKERNELS_VERSION << ", execution space " << exec_space::name()
     << ", memory space " << mem_space::name();
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  ss << ", cuSPARSE TPL enabled";
#else
  ss << ", cuSPARSE TPL disabled";
#endif
  ss << ", Tpetra assumes GPU-aware MPI = "
     << (Tpetra::Details::Behavior::assumeMpiIsGPUAware() ? "yes" : "no")
     << ", CUDA_LAUNCH_BLOCKING = " << (launchBlocking ? launchBlocking : "(unset)");
  return ss.str();
}

// Report whether a value computed from global reductions agrees on every rank.  Collective.
std::string checkRankConsistency(const Teuchos::Comm<int>& comm, ST value,
                                 const std::string& name) {
  ST minVal = 0.0, maxVal = 0.0;
  Teuchos::reduceAll(comm, Teuchos::REDUCE_MIN, value, Teuchos::outArg(minVal));
  Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, value, Teuchos::outArg(maxVal));

  std::ostringstream ss;
  ss << std::scientific << std::setprecision(16);
  ss << "\n      cross-rank consistency of " << name << ": ";
  if (minVal == maxVal)
    ss << "identical on all ranks";
  else
    ss << "DIFFERS ACROSS RANKS (min = " << minVal << ", max = " << maxVal << ")";
  return ss.str();
}

RCP<const crs_t> getBlockMatrix(const Teko::TpetraHelpers::BlockedTpetraOperator& shell, int i,
                                int j) {
  return rcp_dynamic_cast<const crs_t>(shell.GetBlock(i, j), true);
}

std::vector<RCP<const crs_t>> getShellBlocks(
    const Teko::TpetraHelpers::BlockedTpetraOperator& shell, int numBlocks) {
  std::vector<RCP<const crs_t>> blocks;
  for (int i = 0; i < numBlocks; i++)
    for (int j = 0; j < numBlocks; j++) blocks.push_back(getBlockMatrix(shell, i, j));
  return blocks;
}

// Compare the sub-blocks built by two shells from identical input: they must be bitwise
// identical, otherwise sub-block assembly is nondeterministic.  Must be called before either
// shell is reordered (GetBlock indexes the current, possibly reordered, operator).
// Collective.
bool compareShellBlocks(const Teko::TpetraHelpers::BlockedTpetraOperator& a,
                        const Teko::TpetraHelpers::BlockedTpetraOperator& b, int numBlocks,
                        std::string& report) {
  const auto comm = a.getRangeMap()->getComm();

  std::ostringstream ss;
  ss << "\n      sub-block assembly determinism (flat shell vs. reordered shell, before reorder):";
  bool allSame = true;
  for (int i = 0; i < numBlocks; i++) {
    for (int j = 0; j < numBlocks; j++) {
      const auto Ma = getBlockMatrix(a, i, j);
      const auto Mb = getBlockMatrix(b, i, j);

      const bool sameMaps = Ma->getRowMap()->isSameAs(*Mb->getRowMap()) &&
                            Ma->getColMap()->isSameAs(*Mb->getColMap()) &&
                            Ma->getDomainMap()->isSameAs(*Mb->getDomainMap()) &&
                            Ma->getRangeMap()->isSameAs(*Mb->getRangeMap());

      const auto la = Ma->getLocalMatrixHost();
      const auto lb = Mb->getLocalMatrixHost();

      int localStructDiffs = 0;
      int localValueDiffs  = 0;
      if (la.numRows() != lb.numRows() || la.nnz() != lb.nnz()) {
        localStructDiffs++;
      } else {
        for (LO r = 0; r <= la.numRows(); r++)
          if (la.graph.row_map(r) != lb.graph.row_map(r)) localStructDiffs++;
        for (size_t k = 0; localStructDiffs == 0 && k < la.nnz(); k++) {
          const GO ca = Ma->getColMap()->getGlobalElement(la.graph.entries(k));
          const GO cb = Mb->getColMap()->getGlobalElement(lb.graph.entries(k));
          if (ca != cb) localStructDiffs++;
          if (la.values(k) != lb.values(k)) localValueDiffs++;
        }
      }

      int structDiffs = 0, valueDiffs = 0;
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, localStructDiffs,
                         Teuchos::outArg(structDiffs));
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, localValueDiffs, Teuchos::outArg(valueDiffs));

      if (!sameMaps || structDiffs != 0 || valueDiffs != 0) {
        allSame = false;
        ss << "\n        block (" << i << "," << j << "): DIFFERS - maps "
           << (sameMaps ? "same" : "DIFFERENT") << ", " << structDiffs
           << " structural differences, " << valueDiffs << " value differences";
      }
    }
  }
  if (allSame) ss << " all " << numBlocks * numBlocks << " blocks bitwise identical";
  report = ss.str();
  return allSame;
}

// Copy x into the shell's Thyra product layout and back; the result must equal x bitwise.
// Isolates the vector plumbing (one2many/many2one, Import/Export, reorder flattening).
// Collective.
std::string checkMappingRoundTrip(const Teko::TpetraHelpers::BlockedTpetraOperator& shell,
                                  const mv_t& x, const std::string& name) {
  const auto strategy = shell.getMapStrategy();
  const auto op       = shell.getThyraOp();

  ApplyErrorStats stats;
  for (int trial = 0; trial < 2; trial++) {
    RCP<Thyra::MultiVectorBase<ST>> tX = Thyra::createMembers(op->domain(), x.getNumVectors());
    Thyra::assign(tX.ptr(), 0.0);
    strategy->copyTpetraIntoThyra(x, tX.ptr());

    mv_t x2(x.getMap(), x.getNumVectors());
    strategy->copyThyraIntoTpetra(tX, x2);
    accumulateApplyError(x2, x, trial, stats);
  }

  std::ostringstream ss;
  ss << "\n      Tpetra -> Thyra -> Tpetra round trip (" << name
     << " shell, expected exactly 0): " << describe(stats, 0.0);
  return ss.str();
}

// For every sub-block, compare the device SpMV against a host CSR reference with beta = 0
// and beta = 1 (the blocked apply accumulates with beta = 1), and check that repeated device
// applies are bitwise reproducible.  Isolates the SpMV from the Teko/Thyra layers.
// Collective.
std::string checkBlockSpmv(const std::vector<RCP<const crs_t>>& blocks, int numBlocks, ST tolerance,
                           const std::string& name) {
  std::ostringstream ss;
  ss << "\n      per-block device SpMV vs. host CSR reference (" << name << " shell):";

  bool anyBad = false;
  for (int i = 0; i < numBlocks; i++) {
    for (int j = 0; j < numBlocks; j++) {
      const auto M = blocks[i * numBlocks + j];
      const int w  = 3;

      mv_t xj(M->getDomainMap(), w);
      xj.randomize();
      mv_t y0(M->getRangeMap(), w);
      mv_t y0repeat(M->getRangeMap(), w);
      mv_t r(M->getRangeMap(), w);
      r.randomize();
      mv_t y1(r, Teuchos::Copy);

      M->apply(xj, y0);
      M->apply(xj, y0repeat);
      M->apply(xj, y1, Teuchos::NO_TRANS, 1.0, 1.0);

      // host reference: import x into the column map and do a serial CSR product
      mv_t xcol(M->getColMap(), w);
      Tpetra::Import<LO, GO, NT> importer(M->getDomainMap(), M->getColMap());
      xcol.doImport(xj, importer, Tpetra::INSERT);

      mv_t ref0(M->getRangeMap(), w);
      mv_t ref1(M->getRangeMap(), w);
      {
        const auto lcl   = M->getLocalMatrixHost();
        const auto xh    = xcol.getLocalViewHost(Tpetra::Access::ReadOnly);
        const auto rh    = r.getLocalViewHost(Tpetra::Access::ReadOnly);
        const auto ref0h = ref0.getLocalViewHost(Tpetra::Access::OverwriteAll);
        const auto ref1h = ref1.getLocalViewHost(Tpetra::Access::OverwriteAll);
        for (LO row = 0; row < lcl.numRows(); row++) {
          for (int c = 0; c < w; c++) {
            ST sum = 0.0;
            for (auto k = lcl.graph.row_map(row); k < lcl.graph.row_map(row + 1); k++)
              sum += lcl.values(k) * xh(lcl.graph.entries(k), c);
            ref0h(row, c) = sum;
            ref1h(row, c) = rh(row, c) + sum;
          }
        }
      }

      ApplyErrorStats beta0, beta1, repeat;
      accumulateApplyError(y0, ref0, 0, beta0);
      accumulateApplyError(y1, ref1, 0, beta1);
      accumulateApplyError(y0repeat, y0, 0, repeat);

      const bool bad =
          !passed(beta0, tolerance) || !passed(beta1, tolerance) || !passed(repeat, 0.0);
      anyBad |= bad;
      ss << "\n        block (" << i << "," << j << ") " << M->getGlobalNumRows() << "x"
         << M->getGlobalNumCols() << ", nnz " << M->getGlobalNumEntries() << ": "
         << (bad ? "MISMATCH" : "ok") << std::scientific << std::setprecision(3)
         << " (beta=0 rel err " << beta0.maxRelErr << ", beta=1 rel err " << beta1.maxRelErr
         << ", repeat rel diff " << repeat.maxRelErr << ")";
      if (bad) {
        ss << "\n          beta=0: " << describe(beta0, tolerance)
           << "\n          beta=1: " << describe(beta1, tolerance)
           << "\n          repeat: " << describe(repeat, 0.0);
      }
    }
  }
  ss << "\n      per-block SpMV summary (" << name
     << " shell): " << (anyBad ? "MISMATCHES FOUND" : "all blocks ok");
  return ss.str();
}

// Number of non-finite entries in v, summed over all ranks.  Collective.
int countNonFinite(const mv_t& v) {
  const auto vH  = v.getLocalViewHost(Tpetra::Access::ReadOnly);
  int localCount = 0;
  for (size_t i = 0; i < vH.extent(0); i++)
    for (size_t j = 0; j < vH.extent(1); j++)
      if (!std::isfinite(vH(i, j))) localCount++;
  int count = 0;
  Teuchos::reduceAll(*v.getMap()->getComm(), Teuchos::REDUCE_SUM, localCount,
                     Teuchos::outArg(count));
  return count;
}

// Rank-local layout of v (reported from rank 0 only).
std::string describeLayout(const mv_t& v) {
  const auto vD = v.getLocalViewDevice(Tpetra::Access::ReadOnly);
  std::ostringstream ss;
  ss << "local " << v.getLocalLength() << "x" << v.getNumVectors() << ", stride "
     << (v.isConstantStride() ? std::to_string(v.getStride()) : std::string("non-constant"))
     << ", view extents " << vD.extent(0) << "x" << vD.extent(1) << ", device ptr "
     << static_cast<const void*>(vD.data());
  return ss.str();
}

// Gather (rank, GID, column) of every non-finite entry onto rank 0.  Collective.
std::string listNonFinite(const mv_t& v, int maxEntries = 20) {
  const auto map  = v.getMap();
  const auto comm = map->getComm();
  const auto vH   = v.getLocalViewHost(Tpetra::Access::ReadOnly);

  std::vector<double> records;  // rank, gid, column
  for (size_t i = 0; i < vH.extent(0); i++)
    for (size_t j = 0; j < vH.extent(1); j++)
      if (!std::isfinite(vH(i, j)))
        records.insert(records.end(),
                       {static_cast<double>(comm->getRank()),
                        static_cast<double>(map->getGlobalElement(static_cast<LO>(i))),
                        static_cast<double>(j)});

  const int numRanks = comm->getSize();
  int myCount        = static_cast<int>(records.size());
  std::vector<int> counts(numRanks, 0), displs(numRanks, 0);
  Teuchos::gather<int, int>(&myCount, 1, counts.data(), 1, 0, *comm);
  for (int p = 1; p < numRanks; p++) displs[p] = displs[p - 1] + counts[p - 1];
  std::vector<double> all(comm->getRank() == 0 ? displs.back() + counts.back() : 0);
  Teuchos::gatherv<int, double>(records.data(), myCount, all.data(), counts.data(), displs.data(),
                                0, *comm);

  std::ostringstream ss;
  for (size_t k = 0; k < all.size() && static_cast<int>(k / 3) < maxEntries; k += 3)
    ss << " (rank " << static_cast<int>(all[k]) << ", GID " << static_cast<GO>(all[k + 1])
       << ", col " << static_cast<int>(all[k + 2]) << ")";
  return ss.str();
}

// y = alpha * M * x + beta * y0 computed on the host from an explicit column-map import of x.
// Collective.
mv_t hostReferenceApply(const crs_t& M, const mv_t& x, const mv_t& y0, ST alpha, ST beta) {
  const int w = static_cast<int>(x.getNumVectors());
  mv_t xcol(M.getColMap(), w);
  Tpetra::Import<LO, GO, NT> importer(M.getDomainMap(), M.getColMap());
  xcol.doImport(x, importer, Tpetra::INSERT);

  mv_t y(M.getRangeMap(), w);
  {
    const auto lcl = M.getLocalMatrixHost();
    const auto xh  = xcol.getLocalViewHost(Tpetra::Access::ReadOnly);
    const auto y0h = y0.getLocalViewHost(Tpetra::Access::ReadOnly);
    const auto yh  = y.getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (LO row = 0; row < lcl.numRows(); row++) {
      for (int c = 0; c < w; c++) {
        ST sum = 0.0;
        for (auto k = lcl.graph.row_map(row); k < lcl.graph.row_map(row + 1); k++)
          sum += lcl.values(k) * xh(lcl.graph.entries(k), c);
        yh(row, c) = alpha * sum + (beta == 0.0 ? 0.0 : beta * y0h(row, c));
      }
    }
  }
  return y;
}

// Re-run a leaf apply that turned finite inputs into non-finite output in several ways, to
// separate the possible causes.  Collective.
void diagnoseLeafApply(const crs_t& M, const mv_t& X, const mv_t& Ybefore, const mv_t& Yafter,
                       ST alpha, ST beta, std::ostringstream& log) {
  const int w = static_cast<int>(X.getNumVectors());
  log << "\n          non-finite entries in output:" << listNonFinite(Yafter);

  // (a) same matrix, fresh copies of the inputs
  {
    mv_t Xa(X, Teuchos::Copy);
    mv_t Ya(Ybefore, Teuchos::Copy);
    M.apply(Xa, Ya, Teuchos::NO_TRANS, alpha, beta);
    log << "\n          (a) same matrix, fresh input copies:       " << countNonFinite(Ya)
        << " non-finite";
  }
  // (b) deep copy of the matrix, i.e. a fresh SpMV handle and column-map import cache
  {
    crs_t Mb(M, Teuchos::Copy);
    mv_t Yb(Ybefore, Teuchos::Copy);
    Mb.apply(X, Yb, Teuchos::NO_TRANS, alpha, beta);
    log << "\n          (b) deep-copied matrix (fresh handle):     " << countNonFinite(Yb)
        << " non-finite";
  }
  // (c) one column at a time (rank-1 SpMV path instead of the multi-vector path)
  {
    int nonFinite = 0;
    for (int c = 0; c < w; c++) {
      vec_t yc(*Ybefore.getVector(c), Teuchos::Copy);
      M.apply(*X.getVector(c), yc, Teuchos::NO_TRANS, alpha, beta);
      nonFinite += countNonFinite(yc);
    }
    log << "\n          (c) same matrix, column by column:         " << nonFinite << " non-finite";
  }
  // (d) same matrix, multi-vector apply with beta = 1 on a zeroed output
  {
    mv_t Yd(Ybefore.getMap(), w);
    if (beta != 0.0) Yd.update(beta, Ybefore, 0.0);
    M.apply(X, Yd, Teuchos::NO_TRANS, alpha, 1.0);
    log << "\n          (d) same matrix, beta = 1 on prescaled y:  " << countNonFinite(Yd)
        << " non-finite";
  }
  // (e) host reference
  {
    const mv_t Yref = hostReferenceApply(M, X, Ybefore, alpha, beta);
    ApplyErrorStats stats;
    accumulateApplyError(Yafter, Yref, 0, stats);
    log << "\n          (e) host reference: " << countNonFinite(Yref)
        << " non-finite; device result vs. host reference: " << describe(stats, 1e-12);
  }
}

// Replay a (possibly nested) Thyra blocked apply leaf by leaf, exactly as
// Thyra::DefaultBlockedLinearOp::applyImpl does, checking every leaf CrsMatrix::apply for
// non-finite inputs and outputs.  Collective.
void instrumentedApply(const Thyra::LinearOpBase<ST>& op, const Thyra::MultiVectorBase<ST>& X,
                       Thyra::MultiVectorBase<ST>& Y, ST alpha, ST beta, const std::string& path,
                       std::ostringstream& log, bool& introduced,
                       const std::vector<RCP<const crs_t>>& blocks, int numBlocks) {
  using TOE = Thyra::TpetraOperatorVectorExtraction<ST, LO, GO, NT>;

  if (const auto blk = dynamic_cast<const Thyra::BlockedLinearOpBase<ST>*>(&op)) {
    const auto range =
        rcp_dynamic_cast<const Thyra::DefaultProductVectorSpace<ST>>(blk->productRange(), true);
    const auto domain =
        rcp_dynamic_cast<const Thyra::DefaultProductVectorSpace<ST>>(blk->productDomain(), true);
    const auto Xp =
        Thyra::castOrCreateSingleBlockProductMultiVector<ST>(domain, Teuchos::rcpFromRef(X));
    const auto Yp =
        Thyra::nonconstCastOrCreateSingleBlockProductMultiVector<ST>(range, Teuchos::rcpFromRef(Y));

    for (int i = 0; i < range->numBlocks(); i++) {
      const auto Yi = Yp->getNonconstMultiVectorBlock(i);
      for (int j = 0; j < domain->numBlocks(); j++) {
        const auto opij      = blk->getBlock(i, j);
        const std::string ij = path + "(" + std::to_string(i) + "," + std::to_string(j) + ")";
        const ST b           = (j == 0) ? beta : 1.0;
        if (nonnull(opij))
          instrumentedApply(*opij, *Xp->getMultiVectorBlock(j), *Yi, alpha, b, ij, log, introduced,
                            blocks, numBlocks);
        else if (j == 0)
          Thyra::scale(beta, Yi.ptr());
      }
    }
    return;
  }

  const auto tpetraOp = dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT>*>(&op);
  TEUCHOS_ASSERT(tpetraOp != nullptr);
  const auto M  = rcp_dynamic_cast<const crs_t>(tpetraOp->getConstTpetraOperator(), true);
  const auto tX = TOE::getConstTpetraMultiVector(Teuchos::rcpFromRef(X));
  const auto tY = TOE::getTpetraMultiVector(Teuchos::rcpFromRef(Y));

  const mv_t Xcopy(*tX, Teuchos::Copy);
  const mv_t Ybefore(*tY, Teuchos::Copy);
  const int nfX      = countNonFinite(*tX);
  const int nfBefore = countNonFinite(*tY);

  M->apply(*tX, *tY, Teuchos::NO_TRANS, alpha, beta);
  const int nfAfter = countNonFinite(*tY);

  const bool aliased = tX->getLocalViewDevice(Tpetra::Access::ReadOnly).data() ==
                       tY->getLocalViewDevice(Tpetra::Access::ReadOnly).data();

  // identify the leaf by its original (unreordered) block indices
  std::string original = "?";
  for (size_t b = 0; b < blocks.size(); b++)
    if (blocks[b].get() == M.get())
      original = "(" + std::to_string(b / numBlocks) + "," + std::to_string(b % numBlocks) + ")";

  log << "\n        leaf " << path << " = original block " << original << ": "
      << M->getGlobalNumRows() << "x" << M->getGlobalNumCols() << " nnz "
      << M->getGlobalNumEntries() << ", beta " << beta
      << ", non-finite X/Y-before/Y-after = " << nfX << "/" << nfBefore << "/" << nfAfter
      << (aliased ? ", X AND Y ALIAS" : "") << "\n          X: " << describeLayout(*tX)
      << "\n          Y: " << describeLayout(*tY);

  const bool inputsFinite = nfX == 0 && (beta == 0.0 || nfBefore == 0);
  if (inputsFinite && nfAfter > 0 && !introduced) {
    introduced = true;
    log << "\n          ^^^ NON-FINITE VALUES INTRODUCED BY THIS CrsMatrix::apply ^^^";
    diagnoseLeafApply(*M, Xcopy, Ybefore, *tY, alpha, beta, log);
  }
}

// Replay shell.apply(x, y) with y initially zero, exactly as TpetraOperatorWrapper::apply
// does, but with the Thyra apply instrumented leaf by leaf.  Collective.
std::string instrumentedShellApply(const Teko::TpetraHelpers::BlockedTpetraOperator& shell,
                                   const mv_t& x, const std::string& name,
                                   const std::vector<RCP<const crs_t>>& blocks, int numBlocks) {
  const auto strategy = shell.getMapStrategy();
  const auto op       = shell.getThyraOp();
  const int w         = static_cast<int>(x.getNumVectors());

  RCP<Thyra::MultiVectorBase<ST>> tX = Thyra::createMembers(op->domain(), w);
  RCP<Thyra::MultiVectorBase<ST>> tY = Thyra::createMembers(op->range(), w);
  Thyra::assign(tX.ptr(), 0.0);
  Thyra::assign(tY.ptr(), 0.0);

  mv_t y(shell.getRangeMap(), w);
  strategy->copyTpetraIntoThyra(x, tX.ptr());
  strategy->copyTpetraIntoThyra(y, tY.ptr());

  std::ostringstream log;
  log << "\n      instrumented leaf-by-leaf replay of the " << name
      << " shell apply (leaf paths are in the shell's own, possibly reordered, block indices):";
  bool introduced = false;
  instrumentedApply(*op, *tX, *tY, 1.0, 0.0, "", log, introduced, blocks, numBlocks);

  strategy->copyThyraIntoTpetra(tY, y);
  const int nfReplay = countNonFinite(y);

  mv_t yActual(shell.getRangeMap(), w);
  shell.apply(x, yActual);
  const int nfActual = countNonFinite(yActual);

  log << "\n        " << name << " replay result: " << nfReplay
      << " non-finite; actual shell.apply on same input: " << nfActual << " non-finite"
      << (introduced ? "" : " (no leaf apply introduced non-finite values)");
  return log.str();
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
  std::string failMsg =
      std::string("\n   tBlockedTpetraOperator::test_vector_constr: ") +
      "testing tolerance over many matrix vector multiplies: " + describe(stats, tolerance_);
  TEST_ASSERT(passed(stats, tolerance_), failMsg);
  echoFailure(*comm_tpetra, status, failMsg);

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
  failMsg = std::string("\n   tBlockedTpetraOperator::test_vector_constr (rebuild): ") +
            "testing tolerance over many matrix vector multiplies: " + describe(stats, tolerance_);
  TEST_ASSERT(passed(stats, tolerance_), failMsg);
  echoFailure(*comm_tpetra, status, failMsg);

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
  std::string failMsg =
      std::string("\n   tBlockedTpetraOperator::test_single_block: ") +
      "testing tolerance over many matrix vector multiplies: " + describe(stats, tolerance_);
  TEST_ASSERT(passed(stats, tolerance_), failMsg);
  echoFailure(*comm_tpetra, status, failMsg);

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
  failMsg = std::string("\n   tBlockedTpetraOperator::test_single_block (rebuild): ") +
            "testing tolerance over many matrix vector multiplies: " + describe(stats, tolerance_);
  TEST_ASSERT(passed(stats, tolerance_), failMsg);
  echoFailure(*comm_tpetra, status, failMsg);

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
  std::string failMsg =
      std::string("\n   tBlockedTpetraOperator::test_noncontig: ") +
      "testing tolerance over many matrix vector multiplies: " + describe(stats, tolerance_);
  TEST_ASSERT(passed(stats, tolerance_), failMsg);
  echoFailure(*comm_tpetra, status, failMsg);

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
  failMsg = std::string("\n   tBlockedTpetraOperator::test_noncontig (rebuild): ") +
            "testing tolerance over many matrix vector multiplies: " + describe(stats, tolerance_);
  TEST_ASSERT(passed(stats, tolerance_), failMsg);
  echoFailure(*comm_tpetra, status, failMsg);

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

  // the two shells are built from identical input, so their blocks must be bitwise identical
  const int numBlocks = static_cast<int>(vars.size());
  std::string blockReport;
  const bool blocksSame = compareShellBlocks(flatShell, reorderShell, numBlocks, blockReport);

  // capture the sub-block matrices now: after reordering, GetBlock indexes the reordered op
  const auto flatBlocks    = getShellBlocks(flatShell, numBlocks);
  const auto reorderBlocks = getShellBlocks(reorderShell, numBlocks);

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
  Tpetra::MultiVector<ST, LO, GO, NT> yr2(A->getRangeMap(), width);
  Tpetra::MultiVector<ST, LO, GO, NT> ys2(A->getRangeMap(), width);

  int numtests = 10;
  ApplyErrorStats reorderVsFlat;
  ApplyErrorStats flatVsMatrix;
  ApplyErrorStats reorderVsMatrix;
  ApplyErrorStats flatRepeat;
  ApplyErrorStats reorderRepeat;
  ApplyErrorStats matrixRepeat;

  // input and entry-level report for the first trial that exceeds the tolerance
  RCP<mv_t> xFirstFail;
  int firstFailTrial = -1;
  std::string largeEntries;

  for (int i = 0; i < numtests; i++) {
    x.randomize();

    flatShell.apply(x, yf);
    reorderShell.apply(x, yr);
    A->apply(x, ys);

    flatShell.apply(x, yf2);
    reorderShell.apply(x, yr2);
    A->apply(x, ys2);

    accumulateApplyError(yr, yf, i, reorderVsFlat);
    accumulateApplyError(yf, ys, i, flatVsMatrix);
    accumulateApplyError(yr, ys, i, reorderVsMatrix);
    accumulateApplyError(yf2, yf, i, flatRepeat);
    accumulateApplyError(yr2, yr, i, reorderRepeat);
    accumulateApplyError(ys2, ys, i, matrixRepeat);

    // the stats are built from global norms, so every rank agrees on this branch
    const bool okSoFar = passed(reorderVsFlat, tolerance_) && passed(flatVsMatrix, tolerance_) &&
                         passed(reorderVsMatrix, tolerance_);
    if (!okSoFar && firstFailTrial < 0) {
      firstFailTrial = i;
      xFirstFail     = rcp(new mv_t(x, Teuchos::Copy));
      largeEntries   = describeLargeEntries(yf, yr, ys, vars, tolerance_);
    }
  }

  const bool allOk = passed(reorderVsFlat, tolerance_) && passed(flatVsMatrix, tolerance_) &&
                     passed(reorderVsMatrix, tolerance_);

  // On failure, run every diagnostic and report to stdout as well as the failure log.  The
  // repeat comparisons and the diagnostics are informational; only allOk and blocksSame are
  // asserted.
  std::string failMsg;
  if (!allOk || !blocksSame) {
    std::ostringstream diag;
    diag << "\n   tBlockedTpetraOperator::test_reorder" << tstr << " (pattern = " << brm.toString()
         << "): FAILED; diagnostics follow" << describeConfiguration(*comm_tpetra)
         << "\n    reordered vs. flat:   " << describe(reorderVsFlat, tolerance_)
         << "\n    flat vs. matrix:      " << describe(flatVsMatrix, tolerance_)
         << "\n    reordered vs. matrix: " << describe(reorderVsMatrix, tolerance_)
         << "\n    flat vs. flat, repeated apply (expected exactly 0): "
         << describe(flatRepeat, 0.0)
         << "\n    reordered vs. reordered, repeated apply (expected exactly 0): "
         << describe(reorderRepeat, 0.0)
         << "\n    matrix vs. matrix, repeated apply (expected exactly 0): "
         << describe(matrixRepeat, 0.0)
         << checkRankConsistency(*comm_tpetra, reorderVsFlat.maxRelErr, "reordered vs. flat error")
         << checkRankConsistency(*comm_tpetra, flatVsMatrix.maxRelErr, "flat vs. matrix error")
         << checkRankConsistency(*comm_tpetra, reorderVsMatrix.maxRelErr,
                                 "reordered vs. matrix error")
         << blockReport;

    if (firstFailTrial >= 0) {
      diag << "\n      first failing trial: " << firstFailTrial << largeEntries;

      // re-run the failing input, once as-is and once with a full fence around every apply,
      // to see whether the failure is reproducible and whether it is synchronization related
      for (const bool fenced : {false, true}) {
        mv_t yf3(A->getRangeMap(), width);
        mv_t yr3(A->getRangeMap(), width);
        mv_t ys3(A->getRangeMap(), width);
        if (fenced) Kokkos::fence();
        flatShell.apply(*xFirstFail, yf3);
        if (fenced) Kokkos::fence();
        reorderShell.apply(*xFirstFail, yr3);
        if (fenced) Kokkos::fence();
        A->apply(*xFirstFail, ys3);
        if (fenced) Kokkos::fence();

        ApplyErrorStats rf, fm, rm;
        accumulateApplyError(yr3, yf3, 0, rf);
        accumulateApplyError(yf3, ys3, 0, fm);
        accumulateApplyError(yr3, ys3, 0, rm);
        diag << "\n      re-run of failing input" << (fenced ? " with Kokkos::fence()" : "")
             << ":\n        reordered vs. flat:   " << describe(rf, tolerance_)
             << "\n        flat vs. matrix:      " << describe(fm, tolerance_)
             << "\n        reordered vs. matrix: " << describe(rm, tolerance_);
      }

      diag << checkMappingRoundTrip(flatShell, *xFirstFail, "flat")
           << checkMappingRoundTrip(reorderShell, *xFirstFail, "reordered")
           << instrumentedShellApply(reorderShell, *xFirstFail, "reordered", reorderBlocks,
                                     numBlocks)
           << instrumentedShellApply(flatShell, *xFirstFail, "flat", flatBlocks, numBlocks);
    }

    diag << checkBlockSpmv(reorderBlocks, numBlocks, tolerance_, "reordered")
         << checkBlockSpmv(flatBlocks, numBlocks, tolerance_, "flat");
    failMsg = diag.str();
  }
  TEST_ASSERT(allOk && blocksSame, failMsg);
  echoFailure(*comm_tpetra, status, failMsg);

  return allPassed;
}

}  // namespace Test
}  // namespace Teko
