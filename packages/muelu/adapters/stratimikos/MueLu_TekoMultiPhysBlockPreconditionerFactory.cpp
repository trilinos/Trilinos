// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_TekoMultiPhysBlockPreconditionerFactory.hpp"

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA) && defined(HAVE_MUELU_THYRATPETRAADAPTERS) && defined(HAVE_MUELU_TEKO)

#include "MueLu_Exceptions.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_TpetraOperator_decl.hpp"
#include "MueLu_Utilities_decl.hpp"

#include "Teko_PreconditionerFactory.hpp"
#include "Teko_TpetraOperatorWrapper.hpp"
#include "Teko_TpetraThyraConverter.hpp"

#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

#include "Teuchos_Comm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Xpetra_BlockedCrsMatrix_decl.hpp"
#include "Xpetra_BlockedMap_decl.hpp"
#include "Xpetra_BlockedMultiVector_decl.hpp"
#include "Xpetra_MapExtractor_decl.hpp"
#include "Xpetra_MultiVectorFactory_decl.hpp"
#include "Xpetra_ThyraUtils.hpp"
#include "Xpetra_TpetraCrsMatrix_decl.hpp"
#include "Xpetra_TpetraMultiVector.hpp"

#include "Tpetra_Details_makeColMap_decl.hpp"
#include "KokkosSparse_SortCrs.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <tuple>
#include <utility>

namespace MueLu {

namespace {
using ST                       = Teko::ST;
using LO                       = Teko::LO;
using GO                       = Teko::GO;
using NT                       = Teko::NT;
using TpetraMatrix             = Tpetra::CrsMatrix<ST, LO, GO, NT>;
using TpetraMultiVector        = Tpetra::MultiVector<ST, LO, GO, NT>;
using XpetraMatrix             = Xpetra::Matrix<ST, LO, GO, NT>;
using XpetraMultiVector        = Xpetra::MultiVector<ST, LO, GO, NT>;
using XpetraBlockedMatrix      = Xpetra::BlockedCrsMatrix<ST, LO, GO, NT>;
using XpetraBlockedMultiVector = Xpetra::BlockedMultiVector<ST, LO, GO, NT>;
using XpetraMapExtractor       = Xpetra::MapExtractor<ST, LO, GO, NT>;
using XpetraMap                = Xpetra::Map<LO, GO, NT>;
using XpetraBlockedMap         = Xpetra::BlockedMap<LO, GO, NT>;
using MueLuHierarchy           = MueLu::Hierarchy<ST, LO, GO, NT>;
using MultiPhysOperator        = MueLu::MultiPhys<ST, LO, GO, NT>;
using Teuchos::RCP;

struct AMGComplexityStats {
  int levels{0};
  double operatorComplexity{-1.0};
  double smootherComplexity{-1.0};
  double cycleComplexity{-1.0};
  double gridComplexity{-1.0};
};

AMGComplexityStats computeAMGComplexity(const RCP<MueLuHierarchy>& hierarchy) {
  AMGComplexityStats stats;
  if (Teuchos::is_null(hierarchy)) return stats;

  stats.levels             = hierarchy->GetNumLevels();
  stats.operatorComplexity = hierarchy->GetOperatorComplexity();
  stats.smootherComplexity = hierarchy->GetSmootherComplexity();
  stats.cycleComplexity    = stats.operatorComplexity;
  if (stats.smootherComplexity >= 0.0) stats.cycleComplexity += stats.smootherComplexity;

  double fineRows        = 0.0;
  double rowsAcrossLevel = 0.0;
  for (int levelId = 0; levelId < stats.levels; ++levelId) {
    auto level = hierarchy->GetLevel(levelId);
    if (level->IsAvailable("A")) {
      auto A          = level->template Get<RCP<XpetraMatrix> >("A");
      const double nr = static_cast<double>(A->getGlobalNumRows());
      if (levelId == 0) fineRows = nr;
      rowsAcrossLevel += nr;
    }
  }

  if (fineRows > 0.0) stats.gridComplexity = rowsAcrossLevel / fineRows;
  return stats;
}

void printAMGComplexity(const std::string& kind,
                        int index,
                        const AMGComplexityStats& stats,
                        const RCP<const Teuchos::Comm<int> >& comm,
                        bool enabled) {
  if (!enabled) return;
  if (!Teuchos::is_null(comm) && comm->getRank() != 0) return;

  std::cout << "AMG_COMPLEXITY kind=" << kind
            << " index=" << index
            << " levels=" << stats.levels
            << " operator=" << stats.operatorComplexity
            << " smoother=" << stats.smootherComplexity
            << " cycle=" << stats.cycleComplexity
            << " grid=" << stats.gridComplexity << std::endl;
}

RCP<const TpetraMatrix> getCrsMatrix(int i, int j, const Teko::BlockedLinearOp& blocked) {
  bool transposed = false;
  ST scalar       = 0.0;
  auto matrix     = Teko::TpetraHelpers::getTpetraCrsMatrix(Teko::getBlock(i, j, blocked), &scalar, &transposed);
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(matrix), MueLu::Exceptions::RuntimeError,
                             "MueLu::TekoMultiPhysBlockPreconditionerFactory: failed to extract Tpetra::CrsMatrix for block (" << i << "," << j << ").");
  TEUCHOS_TEST_FOR_EXCEPTION(transposed, MueLu::Exceptions::RuntimeError,
                             "MueLu::TekoMultiPhysBlockPreconditionerFactory: transposed Teko blocks are not supported.");
  TEUCHOS_TEST_FOR_EXCEPTION(scalar != Teuchos::ScalarTraits<ST>::one(), MueLu::Exceptions::RuntimeError,
                             "MueLu::TekoMultiPhysBlockPreconditionerFactory: scaled Teko blocks are not supported.");
  return matrix;
}

RCP<const Teuchos::Comm<int> > extractCommunicator(const Teko::BlockedLinearOp& blocked) {
  return getCrsMatrix(0, 0, blocked)->getRowMap()->getComm();
}

auto computeBlockOffsets(const Teko::BlockedLinearOp& blocked) {
  using device_type    = typename NT::device_type;
  const int nBlockRows = Teko::blockRowCount(blocked);

  LO numLocalRows  = 0;
  GO numGlobalRows = 0;
  std::vector<LO> scanLocalBlockOffset(nBlockRows + 1, 0);
  std::vector<GO> scanGlobalBlockOffset(nBlockRows + 1, 0);

  for (int row = 0; row < nBlockRows; ++row) {
    auto Aii    = getCrsMatrix(row, row, blocked);
    auto rowMap = Aii->getRowMap();
    numGlobalRows += rowMap->getGlobalNumElements();
    numLocalRows += static_cast<LO>(rowMap->getLocalNumElements());
    scanLocalBlockOffset[row + 1]  = numLocalRows;
    scanGlobalBlockOffset[row + 1] = numGlobalRows;
  }

  Kokkos::View<LO*, device_type> scanLocalBlockOffsetView(
      Kokkos::ViewAllocateWithoutInitializing("scanLocalBlockOffsetView"), nBlockRows + 1);
  {
    auto hostMirror = Kokkos::create_mirror_view(scanLocalBlockOffsetView);
    std::copy(scanLocalBlockOffset.begin(), scanLocalBlockOffset.end(), hostMirror.data());
    Kokkos::deep_copy(scanLocalBlockOffsetView, hostMirror);
  }

  Kokkos::View<GO*, device_type> scanGlobalBlockOffsetView(
      Kokkos::ViewAllocateWithoutInitializing("scanGlobalBlockOffsetView"), nBlockRows + 1);
  {
    auto hostMirror = Kokkos::create_mirror_view(scanGlobalBlockOffsetView);
    std::copy(scanGlobalBlockOffset.begin(), scanGlobalBlockOffset.end(), hostMirror.data());
    Kokkos::deep_copy(scanGlobalBlockOffsetView, hostMirror);
  }

  return std::make_tuple(numLocalRows, numGlobalRows, scanLocalBlockOffsetView, scanGlobalBlockOffsetView);
}

auto extractViewOfLocalSubblockMatricesAndMaps(const Teko::BlockedLinearOp& blocked) {
  using local_matrix_type = typename TpetraMatrix::local_matrix_device_type;
  using local_map_type    = typename Tpetra::Map<LO, GO, NT>::local_map_type;

  const int nBlockRows = Teko::blockRowCount(blocked);
  const int nBlockCols = Teko::blockColCount(blocked);
  Kokkos::View<local_matrix_type**, Kokkos::SharedSpace> subblockMatrices(
      Kokkos::view_alloc("subblockMatrices", Kokkos::SequentialHostInit), nBlockRows, nBlockCols);
  Kokkos::View<local_map_type**, Kokkos::SharedSpace> subblockColMaps(
      Kokkos::view_alloc("subblockColMaps", Kokkos::SequentialHostInit), nBlockRows, nBlockCols);

  for (int row = 0; row < nBlockRows; ++row) {
    for (int col = 0; col < nBlockCols; ++col) {
      auto Aij                   = getCrsMatrix(row, col, blocked);
      subblockMatrices(row, col) = local_matrix_type(Aij->getLocalMatrixDevice());
      subblockColMaps(row, col)  = local_map_type(Aij->getColMap()->getLocalMap());
    }
  }

  return std::make_tuple(subblockMatrices, subblockColMaps);
}

// localRow always lies in [0, scanOffsetView(nBlockRows)), so a block row is always found.
template <class OffsetViewType>
KOKKOS_INLINE_FUNCTION int computeBlockRow(const OffsetViewType& scanOffsetView, LO localRow) {
  int blockRow          = -1;
  const auto nBlockRows = scanOffsetView.extent(0) - 1;
  for (std::size_t row = 0; row < nBlockRows; ++row) {
    const auto lower = scanOffsetView(row);
    const auto upper = scanOffsetView(row + 1);
    if (lower <= localRow && localRow < upper) {
      blockRow = static_cast<int>(row);
      return blockRow;
    }
  }
  return blockRow;
}

RCP<TpetraMatrix> convertBlockMatrixToCrsMatrix(const Teko::BlockedLinearOp& blocked) {
  using local_matrix_type = typename TpetraMatrix::local_matrix_device_type;
  using local_map_type    = typename Tpetra::Map<LO, GO, NT>::local_map_type;
  using row_map_type      = typename local_matrix_type::row_map_type::non_const_type;
  using values_type       = typename local_matrix_type::values_type::non_const_type;
  using index_type        = typename local_matrix_type::index_type::non_const_type;
  using size_type         = typename row_map_type::non_const_value_type;
  using execution_space   = typename local_matrix_type::execution_space;
  using device_type       = typename NT::device_type;

  const int nBlockCols = Teko::blockColCount(blocked);

  LO numLocalRows{0};
  GO numGlobalRows{0};
  Kokkos::View<LO*, device_type> scanLocalBlockOffsetView;
  Kokkos::View<GO*, device_type> scanGlobalBlockOffsetView;
  std::tie(numLocalRows, numGlobalRows, scanLocalBlockOffsetView, scanGlobalBlockOffsetView) =
      computeBlockOffsets(blocked);

  const auto invalidLO = Teuchos::OrdinalTraits<LO>::invalid();
  const LO numRows     = numLocalRows;
  auto nEntriesPerRow  = Kokkos::View<LO*, device_type>(Kokkos::ViewAllocateWithoutInitializing("nEntriesPerRow"), numRows);
  Kokkos::deep_copy(nEntriesPerRow, invalidLO);
  auto prefixSumEntriesPerRow = row_map_type(
      Kokkos::ViewAllocateWithoutInitializing("prefixSumEntriesPerRow"), numRows + 1);

  Kokkos::View<local_matrix_type**, Kokkos::SharedSpace> subblockMatrices;
  Kokkos::View<local_map_type**, Kokkos::SharedSpace> subblockColMaps;
  std::tie(subblockMatrices, subblockColMaps) = extractViewOfLocalSubblockMatricesAndMaps(blocked);

  size_type totalNumOwnedCols = 0;
  Kokkos::parallel_scan(
      Kokkos::RangePolicy<Kokkos::Schedule<Kokkos::Dynamic>, execution_space>(0, numRows),
      KOKKOS_LAMBDA(const LO localRow, size_type& sumNumEntries, bool finalPass) {
        auto numOwnedCols = nEntriesPerRow(localRow);
        if (numOwnedCols == invalidLO) {
          const auto blockRow = computeBlockRow(scanLocalBlockOffsetView, localRow);
          numOwnedCols        = 0;
          const auto lid      = localRow - scanLocalBlockOffsetView(blockRow);
          for (int col = 0; col < nBlockCols; ++col) numOwnedCols += subblockMatrices(blockRow, col).row(lid).length;
          nEntriesPerRow(localRow) = numOwnedCols;
        }
        if (finalPass) {
          prefixSumEntriesPerRow(localRow) = sumNumEntries;
          if (localRow == (numRows - 1)) prefixSumEntriesPerRow(numRows) = prefixSumEntriesPerRow(localRow) + numOwnedCols;
        }
        sumNumEntries += numOwnedCols;
      },
      totalNumOwnedCols);

  auto columnIndices = Kokkos::View<GO*, device_type>(Kokkos::ViewAllocateWithoutInitializing("columnIndices"), totalNumOwnedCols);
  auto values        = values_type(Kokkos::ViewAllocateWithoutInitializing("values"), totalNumOwnedCols);

  Kokkos::parallel_for(
      Kokkos::RangePolicy<Kokkos::Schedule<Kokkos::Dynamic>, execution_space>(0, numRows),
      KOKKOS_LAMBDA(const LO localRow) {
        const auto blockRow        = computeBlockRow(scanLocalBlockOffsetView, localRow);
        const size_type colIdStart = prefixSumEntriesPerRow(localRow);
        LO colId                   = 0;
        const auto lid             = localRow - scanLocalBlockOffsetView(blockRow);
        for (int blockCol = 0; blockCol < nBlockCols; ++blockCol) {
          auto colMap              = subblockColMaps(blockRow, blockCol);
          const auto gidStart      = scanGlobalBlockOffsetView(blockCol);
          const auto sparseRowView = subblockMatrices(blockRow, blockCol).row(lid);
          for (LO col = 0; col < sparseRowView.length; ++col) {
            auto colidx                       = colMap.getGlobalElement(sparseRowView.colidx(col));
            values[colIdStart + colId]        = sparseRowView.value(col);
            columnIndices[colIdStart + colId] = colidx + gidStart;
            ++colId;
          }
        }
      });

  Teko::LinearOp blockedOpLo = blocked;
  auto blockedOpTpetra       = Teuchos::rcp(new Teko::TpetraHelpers::TpetraOperatorWrapper(blockedOpLo));
  auto rangeMap              = blockedOpTpetra->getRangeMap();
  auto domainMap             = blockedOpTpetra->getDomainMap();

  Teuchos::RCP<const Tpetra::Map<LO, GO, NT> > colMap;
  Tpetra::Details::makeColMap<LO, GO, NT>(colMap, domainMap, columnIndices);
  TEUCHOS_ASSERT(colMap);

  auto colMapDev          = colMap->getLocalMap();
  auto localColumnIndices = index_type(Kokkos::ViewAllocateWithoutInitializing("localColumnIndices"), totalNumOwnedCols);
  Kokkos::parallel_for(
      Kokkos::RangePolicy<execution_space, Kokkos::IndexType<size_type> >(0, totalNumOwnedCols),
      KOKKOS_LAMBDA(const size_type index) { localColumnIndices(index) = colMapDev.getLocalElement(columnIndices(index)); });

  KokkosSparse::sort_crs_matrix<execution_space, row_map_type, index_type, values_type>(
      prefixSumEntriesPerRow, localColumnIndices, values);

  const LO numCols = static_cast<LO>(colMap->getLocalNumElements());
  auto localMatrix = local_matrix_type(
      "localMat", numRows, numCols, totalNumOwnedCols, values, prefixSumEntriesPerRow, localColumnIndices);

  return Teuchos::rcp(new TpetraMatrix(localMatrix, rangeMap, colMap, domainMap, rangeMap));
}

// Update the values of a matrix previously built by convertBlockMatrixToCrsMatrix in place.
// Returns false if the row layout or sparsity pattern of the blocks changed, in which case
// the matrix must be rebuilt.
bool updateBlockMatrixToCrsMatrix(const Teko::BlockedLinearOp& blocked,
                                  Teuchos::RCP<TpetraMatrix>& matrix) {
  using local_matrix_type = typename TpetraMatrix::local_matrix_device_type;
  using local_map_type    = typename Tpetra::Map<LO, GO, NT>::local_map_type;
  using execution_space   = typename local_matrix_type::execution_space;
  using device_type       = typename NT::device_type;

  const int nBlockCols = Teko::blockColCount(blocked);

  LO numLocalRows{0};
  GO numGlobalRows{0};
  Kokkos::View<LO*, device_type> scanLocalBlockOffsetView;
  Kokkos::View<GO*, device_type> scanGlobalBlockOffsetView;
  std::tie(numLocalRows, numGlobalRows, scanLocalBlockOffsetView, scanGlobalBlockOffsetView) =
      computeBlockOffsets(blocked);

  // A nonzero count means the row layout or sparsity pattern changed on this rank.
  const LO numRows     = static_cast<LO>(matrix->getRowMap()->getLocalNumElements());
  LO numMissingEntries = (numRows == numLocalRows) ? 0 : 1;

  if (numMissingEntries == 0) {
    matrix->resumeFill();
    matrix->setAllToScalar(0.0);

    auto matrixDevice = matrix->getLocalMatrixDevice();
    auto colMapDev    = matrix->getColMap()->getLocalMap();

    Kokkos::View<local_matrix_type**, Kokkos::SharedSpace> subblockMatrices;
    Kokkos::View<local_map_type**, Kokkos::SharedSpace> subblockColMaps;
    std::tie(subblockMatrices, subblockColMaps) = extractViewOfLocalSubblockMatricesAndMaps(blocked);

    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<Kokkos::Schedule<Kokkos::Dynamic>, execution_space>(0, numRows),
        KOKKOS_LAMBDA(const LO localRow, LO& numMissing) {
          const auto blockRow = computeBlockRow(scanLocalBlockOffsetView, localRow);
          const auto lid      = localRow - scanLocalBlockOffsetView(blockRow);
          for (int blockCol = 0; blockCol < nBlockCols; ++blockCol) {
            const GO gidStart        = scanGlobalBlockOffsetView(blockCol);
            const auto sparseRowView = subblockMatrices(blockRow, blockCol).row(lid);
            const auto mapij         = subblockColMaps(blockRow, blockCol);
            for (LO col = 0; col < sparseRowView.length; ++col) {
              auto value       = sparseRowView.value(col);
              const auto gidij = mapij.getGlobalElement(sparseRowView.colidx(col));
              const auto gid   = gidij + gidStart;
              auto colId       = colMapDev.getLocalElement(gid);
              if (matrixDevice.sumIntoValues(localRow, &colId, 1, &value, true, false) != 1) ++numMissing;
            }
          }
        },
        numMissingEntries);
  }

  // The layout may only have changed on some ranks; agree globally before the collective fillComplete.
  LO globalNumMissingEntries = 0;
  Teuchos::reduceAll(*matrix->getRowMap()->getComm(), Teuchos::REDUCE_MAX, numMissingEntries, Teuchos::outArg(globalNumMissingEntries));
  if (globalNumMissingEntries > 0) return false;

  Teko::LinearOp blockedOpLo = blocked;
  auto blockedOpTpetra       = Teuchos::rcp(new Teko::TpetraHelpers::TpetraOperatorWrapper(blockedOpLo));
  matrix->fillComplete(blockedOpTpetra->getDomainMap(), blockedOpTpetra->getRangeMap());
  return true;
}

RCP<const Teuchos::Comm<int> > extractCommunicator(const Teko::BlockedMultiVector& blocked) {
  auto mv = Thyra::TpetraOperatorVectorExtraction<ST, LO, GO, NT>::getConstTpetraMultiVector(Teko::getBlock(0, blocked));
  return mv->getMap()->getComm();
}

Teko::MultiVector convertBlockVectorToVector(const Teko::BlockedMultiVector& blocked) {
  auto comm = extractCommunicator(blocked);
  auto map  = Teko::TpetraHelpers::thyraVSToTpetraMap(*blocked->productSpace(), Thyra::convertTpetraToThyraComm(comm));
  auto mv   = Teuchos::rcp(new TpetraMultiVector(map, Teuchos::as<size_t>(blocked->domain()->dim())));
  Teko::TpetraHelpers::blockThyraToTpetra(blocked, *mv);
  return Thyra::createMultiVector(mv);
}

std::string toLower(const std::string& input) {
  std::string lowered = input;
  std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return lowered;
}

bool parseSubblockIndex(const std::string& name, int& index) {
  const std::string prefix = "subblockList";
  if (name.rfind(prefix, 0) != 0) return false;

  const std::string suffix = name.substr(prefix.size());
  if (suffix.empty()) return false;
  for (char c : suffix)
    if (!std::isdigit(static_cast<unsigned char>(c))) return false;

  index = std::stoi(suffix);
  return true;
}

bool detectBlockSmootherRecursive(const Teuchos::ParameterList& pl) {
  for (Teuchos::ParameterList::ConstIterator it = pl.begin(); it != pl.end(); ++it) {
    const auto& name  = pl.name(it);
    const auto& entry = pl.entry(it);
    if (entry.isList()) {
      if (detectBlockSmootherRecursive(Teuchos::getValue<Teuchos::ParameterList>(entry))) return true;
      continue;
    }
    const auto loweredName = toLower(name);
    if (loweredName.find("smoother:") == std::string::npos || loweredName.find("type") == std::string::npos) continue;
    if (!entry.isType<std::string>()) continue;
    if (toLower(Teuchos::getValue<std::string>(entry)).find("teko") != std::string::npos) return true;
  }
  return false;
}

class WrapMultiPhysAsTekoBlockOperator : public Teko::BlockImplicitLinearOp {
 public:
  WrapMultiPhysAsTekoBlockOperator(Teko::BlockedLinearOp& A,
                                   const RCP<MultiPhysOperator>& Minv,
                                   bool useBlockSmoothers,
                                   const RCP<const Teuchos::Comm<int> >& comm,
                                   const RCP<const XpetraMapExtractor>& domainMapExtractor)
    : A_(A)
    , Minv_(Minv)
    , productRange_(A->productRange())
    , productDomain_(A->productDomain())
    , useBlockSmoothers_(useBlockSmoothers)
    , comm_(comm)
    , domainMapExtractor_(domainMapExtractor) {
    if (useBlockSmoothers_) {
      TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(comm_), MueLu::Exceptions::RuntimeError,
                                 "MueLu::TekoMultiPhysBlockPreconditionerFactory: comm is null for blocked apply.");
      TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(domainMapExtractor_), MueLu::Exceptions::RuntimeError,
                                 "MueLu::TekoMultiPhysBlockPreconditionerFactory: domain map extractor is null for blocked apply.");
      // The vector spaces are fixed for the lifetime of this operator, so build the Xpetra blocked maps once
      // instead of on every apply (which is collective and allocates).
      rBlockedMap_ = Teuchos::rcp_dynamic_cast<const XpetraBlockedMap>(Xpetra::ThyraUtils<ST, LO, GO, NT>::toXpetra(productDomain_, comm_), true);
      zBlockedMap_ = Teuchos::rcp_dynamic_cast<const XpetraBlockedMap>(Xpetra::ThyraUtils<ST, LO, GO, NT>::toXpetra(productRange_, comm_), true);
    }
  }

  Teko::VectorSpace range() const override { return productRange_; }
  Teko::VectorSpace domain() const override { return productDomain_; }

  void implicitApply(const Teko::BlockedMultiVector& r,
                     Teko::BlockedMultiVector& z,
                     const double alpha = 1.0,
                     const double beta  = 0.0) const override {
    // MultiPhys::apply ignores alpha and beta and always uses a zero initial guess,
    // so compute M r separately and form z = alpha * M r + beta * z here.
    const Teuchos::Ptr<Thyra::MultiVectorBase<ST> > zPtr = z.ptr();
    if (beta == 0.0) {
      applyMultiPhys(r, z);
      if (alpha != 1.0) Thyra::scale(alpha, zPtr);
    } else {
      auto Mr = Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<ST> >(
          Thyra::createMembers(z->range(), z->domain()->dim()), true);
      applyMultiPhys(r, Mr);
      Thyra::scale(beta, zPtr);
      Thyra::update(alpha, *Mr, zPtr);
    }
  }

  void implicitApply(const Thyra::EOpTransp M_trans,
                     const Teko::BlockedMultiVector& r,
                     Teko::BlockedMultiVector& z,
                     const double alpha = 1.0,
                     const double beta  = 0.0) const override {
    TEUCHOS_TEST_FOR_EXCEPTION(M_trans != Thyra::NOTRANS, std::logic_error,
                               "MueLu::TekoMultiPhysBlockPreconditionerFactory does not support transpose apply.");
    implicitApply(r, z, alpha, beta);
  }

 private:
  // Computes z = M r.
  void applyMultiPhys(const Teko::BlockedMultiVector& r, const Teko::BlockedMultiVector& z) const {
    using TpetraExtraction = Thyra::TpetraOperatorVectorExtraction<ST, LO, GO, NT>;

    if (useBlockSmoothers_) {
      // Wrap the Thyra blocks without copying them.
      const int numBlocks = r->productSpace()->numBlocks();
      std::vector<RCP<XpetraMultiVector> > rBlocks(numBlocks);
      std::vector<RCP<XpetraMultiVector> > zBlocks(numBlocks);
      for (int block = 0; block < numBlocks; ++block) {
        rBlocks[block] = Xpetra::toXpetra(Teuchos::rcp_const_cast<TpetraMultiVector>(
            TpetraExtraction::getConstTpetraMultiVector(r->getMultiVectorBlock(block))));
        zBlocks[block] = Xpetra::toXpetra(TpetraExtraction::getTpetraMultiVector(z->getNonconstMultiVectorBlock(block)));
      }
      auto rXpetra = Teuchos::rcp(new XpetraBlockedMultiVector(rBlockedMap_, rBlocks));
      auto zXpetra = Teuchos::rcp(new XpetraBlockedMultiVector(zBlockedMap_, zBlocks));

      Minv_->apply(*rXpetra, *zXpetra, Teuchos::NO_TRANS);
      Xpetra::ThyraUtils<ST, LO, GO, NT>::updateThyra(
          Teuchos::rcp_implicit_cast<const XpetraMultiVector>(zXpetra), domainMapExtractor_,
          Teuchos::rcp_implicit_cast<Thyra::MultiVectorBase<ST> >(z));
      return;
    }

    // Reuse the flat scratch vectors unless the number of columns changed. The previous contents of z
    // do not need to be copied since MultiPhys::apply uses a zero initial guess.
    if (flattenedR_.is_null() || flattenedR_->domain()->dim() != r->domain()->dim()) {
      flattenedR_ = convertBlockVectorToVector(r);
      flattenedZ_ = convertBlockVectorToVector(z);
    } else {
      Teko::TpetraHelpers::blockThyraToTpetra(r, *TpetraExtraction::getTpetraMultiVector(flattenedR_));
    }

    auto zMv = TpetraExtraction::getTpetraMultiVector(flattenedZ_);
    auto rMv = TpetraExtraction::getTpetraMultiVector(flattenedR_);
    auto zXp = Xpetra::toXpetra(zMv);
    auto rXp = Xpetra::toXpetra(rMv);
    Minv_->apply(*rXp, *zXp, Teuchos::NO_TRANS);
    Teko::TpetraHelpers::blockTpetraToThyra(*zMv, z.ptr());
  }

  Teko::BlockedLinearOp A_;
  RCP<MultiPhysOperator> Minv_;
  Teko::VectorSpace productRange_;
  Teko::VectorSpace productDomain_;
  bool useBlockSmoothers_;
  RCP<const Teuchos::Comm<int> > comm_;
  RCP<const XpetraMapExtractor> domainMapExtractor_;
  RCP<const XpetraBlockedMap> rBlockedMap_;
  RCP<const XpetraBlockedMap> zBlockedMap_;
  mutable Teko::MultiVector flattenedR_;
  mutable Teko::MultiVector flattenedZ_;
};

}  // namespace

TekoMultiPhysBlockPreconditionerState::TekoMultiPhysBlockPreconditionerState()
  : Teko::BlockPreconditionerState()
  , numBlocks(-1) {}

TekoMultiPhysBlockPreconditionerFactory::TekoMultiPhysBlockPreconditionerFactory()
  : Teko::BlockPreconditionerFactory()
  , amgSettings_(Teuchos::rcp(new Teuchos::ParameterList()))
  , useBlockSmoothers_(false)
  , printComplexity_(false) {}

Teuchos::RCP<Teko::PreconditionerState>
TekoMultiPhysBlockPreconditionerFactory::buildPreconditionerState() const {
  return Teuchos::rcp(new TekoMultiPhysBlockPreconditionerState());
}

Teuchos::RCP<TekoMultiPhysBlockPreconditionerState>
TekoMultiPhysBlockPreconditionerFactory::ensureConcreteState(Teko::BlockPreconditionerState& state) const {
  auto concrete = Teuchos::rcp_dynamic_cast<TekoMultiPhysBlockPreconditionerState>(Teuchos::rcpFromRef(state), false);
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(concrete), MueLu::Exceptions::RuntimeError,
                             "MueLu::TekoMultiPhysBlockPreconditionerFactory requires TekoMultiPhysBlockPreconditionerState.");
  return concrete;
}

Teuchos::RCP<XpetraMultiVector>
TekoMultiPhysBlockPreconditionerFactory::convertToXpetraMultiVector(const Teuchos::ParameterList& pl,
                                                                    const std::string& name) const {
  if (!pl.isParameter(name)) return Teuchos::null;
  if (pl.isType<RCP<XpetraMultiVector> >(name)) return pl.get<RCP<XpetraMultiVector> >(name);
  if (pl.isType<RCP<const XpetraMultiVector> >(name)) return Teuchos::rcp_const_cast<XpetraMultiVector>(pl.get<RCP<const XpetraMultiVector> >(name));
  if (pl.isType<RCP<TpetraMultiVector> >(name)) return Xpetra::toXpetra(pl.get<RCP<TpetraMultiVector> >(name));
  if (pl.isType<RCP<const TpetraMultiVector> >(name)) return Xpetra::toXpetra(Teuchos::rcp_const_cast<TpetraMultiVector>(pl.get<RCP<const TpetraMultiVector> >(name)));

  TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
                             "MueLu::TekoMultiPhysBlockPreconditionerFactory: parameter '" << name
                                                                                           << "' must be an Xpetra or Tpetra multivector.");
  return Teuchos::null;
}

void TekoMultiPhysBlockPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList& pl) {
  const auto& settings = pl.isSublist("AMG Settings") ? pl.sublist("AMG Settings") : pl;
  amgSettings_         = Teuchos::rcp(new Teuchos::ParameterList(settings));
  amgSettings_->remove("Type", false);

  std::set<int> subblockIndices;
  for (Teuchos::ParameterList::ConstIterator it = amgSettings_->begin(); it != amgSettings_->end(); ++it) {
    const auto& name = amgSettings_->name(it);
    if (!amgSettings_->isSublist(name)) continue;
    int index = -1;
    if (parseSubblockIndex(name, index)) subblockIndices.insert(index);
  }

  int numBlocks = 0;
  if (!subblockIndices.empty()) {
    numBlocks = *subblockIndices.rbegin() + 1;
    for (int block = 0; block < numBlocks; ++block)
      TEUCHOS_TEST_FOR_EXCEPTION(subblockIndices.count(block) == 0, MueLu::Exceptions::RuntimeError,
                                 "MueLu::TekoMultiPhysBlockPreconditionerFactory: non-contiguous subblockList entries; missing subblockList" << block << ".");
  }

  coordinates_.assign(numBlocks, Teuchos::null);
  material_.assign(numBlocks, Teuchos::null);
  nullspace_.assign(numBlocks, Teuchos::null);

  if (pl.isSublist("user data")) {
    const auto& userData = pl.sublist("user data");
    for (int block = 0; block < numBlocks; ++block) {
      coordinates_[block] = convertToXpetraMultiVector(userData, "Coordinates" + std::to_string(block));
      material_[block]    = convertToXpetraMultiVector(userData, "Material" + std::to_string(block));
      nullspace_[block]   = convertToXpetraMultiVector(userData, "Nullspace" + std::to_string(block));
    }
  }

  for (int block = 0; block < numBlocks; ++block) {
    const auto subblockName = "subblockList" + std::to_string(block);
    if (!amgSettings_->isSublist(subblockName)) continue;
    auto& subblockParams = amgSettings_->sublist(subblockName);
    if (!subblockParams.isSublist("user data")) continue;
    const auto& userData = subblockParams.sublist("user data");
    if (userData.isParameter("Coordinates")) coordinates_[block] = convertToXpetraMultiVector(userData, "Coordinates");
    if (userData.isParameter("Material")) material_[block] = convertToXpetraMultiVector(userData, "Material");
    if (userData.isParameter("Nullspace")) nullspace_[block] = convertToXpetraMultiVector(userData, "Nullspace");
  }

  useBlockSmoothers_ = detectBlockSmootherRecursive(*amgSettings_);
  printComplexity_   = false;
  if (amgSettings_->isParameter("multiphys: print complexity")) {
    printComplexity_ = amgSettings_->get<bool>("multiphys: print complexity");
    amgSettings_->remove("multiphys: print complexity");
  }

  if (useBlockSmoothers_) {
    const bool useSubcommunicators = false;
    amgSettings_->set("repartition: use subcommunicators", useSubcommunicators);
    for (int subblock = 0; subblock < numBlocks; ++subblock)
      amgSettings_->sublist("subblockList" + std::to_string(subblock), true)
          .set("repartition: use subcommunicators", useSubcommunicators);
  }
}

Teko::LinearOp
TekoMultiPhysBlockPreconditionerFactory::buildPreconditionerOperator(Teko::BlockedLinearOp& blo,
                                                                     Teko::BlockPreconditionerState& state) const {
  Teuchos::TimeMonitor timer(*Teuchos::TimeMonitor::getNewTimer("MueLu::TekoMultiPhysBlockPreconditionerFactory::buildPreconditionerOperator"));

  const int numBlockRows = Teko::blockRowCount(blo);
  const int numBlockCols = Teko::blockColCount(blo);
  TEUCHOS_TEST_FOR_EXCEPTION(numBlockRows != numBlockCols, MueLu::Exceptions::RuntimeError,
                             "MueLu::TekoMultiPhysBlockPreconditionerFactory requires a square blocked operator.");

  auto comm          = extractCommunicator(blo);
  auto concreteState = ensureConcreteState(state);

  if (!coordinates_.empty())
    TEUCHOS_TEST_FOR_EXCEPTION(static_cast<int>(coordinates_.size()) != numBlockRows, MueLu::Exceptions::RuntimeError,
                               "MueLu::TekoMultiPhysBlockPreconditionerFactory has " << coordinates_.size()
                                                                                     << " configured blocks, but the operator has " << numBlockRows << ".");

  auto localCoordinates = coordinates_;
  auto localMaterial    = material_;
  auto localNullspace   = nullspace_;
  if (static_cast<int>(localCoordinates.size()) < numBlockRows) {
    localCoordinates.resize(numBlockRows);
    localMaterial.resize(numBlockRows);
    localNullspace.resize(numBlockRows);
  }

  for (int block = 0; block < numBlockRows; ++block)
    TEUCHOS_TEST_FOR_EXCEPTION(!amgSettings_->isSublist("subblockList" + std::to_string(block)), MueLu::Exceptions::RuntimeError,
                               "MueLu::TekoMultiPhysBlockPreconditionerFactory: missing parameter sublist subblockList" << block << ".");

  const bool sameNumBlocks = concreteState->numBlocks == numBlockRows;
  bool reuseMultiPhys      = sameNumBlocks && !Teuchos::is_null(concreteState->multiphysPreconditioner);

  if (useBlockSmoothers_) {
    if (Teuchos::is_null(concreteState->blockedMatrix) || !sameNumBlocks) {
      concreteState->blockedMatrix = Teuchos::rcp(new XpetraBlockedMatrix(blo, comm));
    } else {
      for (int row = 0; row < numBlockRows; ++row)
        for (int col = 0; col < numBlockCols; ++col)
          concreteState->blockedMatrix->setMatrix(
              row, col,
              Xpetra::toXpetra(Teuchos::rcp_const_cast<TpetraMatrix>(getCrsMatrix(row, col, blo))));
    }
    concreteState->monolithicMatrix = Teuchos::null;
  } else {
    if (Teuchos::is_null(concreteState->monolithicMatrix) || !sameNumBlocks ||
        !updateBlockMatrixToCrsMatrix(blo, concreteState->monolithicMatrix)) {
      // First build, or the layout/sparsity pattern of the blocks changed: rebuild from scratch.
      concreteState->monolithicMatrix = convertBlockMatrixToCrsMatrix(blo);
      reuseMultiPhys                  = false;
    }
    concreteState->blockedMatrix = Teuchos::null;
  }

  Teuchos::ArrayRCP<RCP<XpetraMatrix> > arrayOfAuxMatrices(numBlockRows);
  Teuchos::ArrayRCP<RCP<XpetraMultiVector> > arrayOfCoords(numBlockRows);
  Teuchos::ArrayRCP<RCP<XpetraMultiVector> > arrayOfMaterials(numBlockRows);
  Teuchos::ArrayRCP<RCP<XpetraMultiVector> > arrayOfNullspaces(numBlockRows);

  for (int block = 0; block < numBlockRows; ++block) {
    arrayOfAuxMatrices[block] = Xpetra::toXpetra(Teuchos::rcp_const_cast<TpetraMatrix>(getCrsMatrix(block, block, blo)));
    arrayOfCoords[block]      = localCoordinates[block];
    arrayOfMaterials[block]   = localMaterial[block];
    arrayOfNullspaces[block]  = localNullspace[block];
  }

  RCP<XpetraMatrix> A = useBlockSmoothers_
                            ? Teuchos::rcp_implicit_cast<XpetraMatrix>(concreteState->blockedMatrix)
                            : Xpetra::toXpetra(concreteState->monolithicMatrix);

  if (reuseMultiPhys) {
    concreteState->multiphysPreconditioner->resetMatrix(A, true);
  } else {
    concreteState->multiphysPreconditioner = Teuchos::rcp(new MultiPhysOperator(
        A, arrayOfAuxMatrices, arrayOfNullspaces, arrayOfCoords, numBlockRows, *amgSettings_, true,
        arrayOfMaterials, true));
  }
  concreteState->numBlocks = numBlockRows;

  printAMGComplexity("multiphysics", -1,
                     computeAMGComplexity(concreteState->multiphysPreconditioner->multiphysicsHierarchy()),
                     comm, printComplexity_);
  const auto subblockHierarchies = concreteState->multiphysPreconditioner->subblockHierarchies();
  for (int block = 0; block < static_cast<int>(subblockHierarchies.size()); ++block)
    printAMGComplexity("subblock", block, computeAMGComplexity(subblockHierarchies[block]), comm,
                       printComplexity_);

  RCP<const XpetraMapExtractor> domainMapExtractor;
  if (!Teuchos::is_null(concreteState->blockedMatrix)) domainMapExtractor = concreteState->blockedMatrix->getDomainMapExtractor();

  return Teuchos::rcp(new WrapMultiPhysAsTekoBlockOperator(
      blo, concreteState->multiphysPreconditioner, useBlockSmoothers_, comm, domainMapExtractor));
}

void registerTekoMultiPhysBlockPreconditionerFactory(const std::string& name) {
  Teko::PreconditionerFactory::addPreconditionerFactory(
      name,
      Teuchos::rcp(new Teko::AutoClone<MueLu::TekoMultiPhysBlockPreconditionerFactory>()));
}

}  // namespace MueLu

#endif