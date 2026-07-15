// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_STRUCTUREDRAPFACTORY_DEF_HPP
#define MUELU_STRUCTUREDRAPFACTORY_DEF_HPP

#include <algorithm>
#include <sstream>
#include <vector>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <stdexcept>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_TripleMatrixMultiply.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_MapFactory.hpp>

#include "MueLu_StructuredRAPFactory_decl.hpp"

#include "MueLu_Utilities.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Behavior.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_CommHelpers.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::StructuredRAPFactory()
  : hasDeclaredInput_(false) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~StructuredRAPFactory() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  validParamList->set<std::string>(
      "rap: matrix type", "", "Galeri matrix type used to infer the structured RAP graph.");
  validParamList->set<int>("rap: processor grid x", -1, "Number of processor subdomains in x direction.");
  validParamList->set<int>("rap: processor grid y", -1, "Number of processor subdomains in y direction.");
  validParamList->set<int>("rap: processor grid z", -1, "Number of processor subdomains in z direction.");

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("rap: triple product");         // in the long term this has to be the only option for multiplication
  SET_VALID_ENTRY("rap: prebuild coarse graph");  // if true, the coarse graph is prebuilt and passed to the triple matrix product kernel. This can be used to optimize the symbolic phase of the triple matrix product.
  SET_VALID_ENTRY("rap: fix zero diagonals");
  SET_VALID_ENTRY("rap: fix zero diagonals threshold");
  SET_VALID_ENTRY("rap: fix zero diagonals replacement");
  SET_VALID_ENTRY("rap: relative diagonal floor");
#undef SET_VALID_ENTRY
  validParamList->set<bool>(
      "transpose: use implicit", true,
      "Use P^T as the restriction operator. StructuredRAPFactory requires this option to be true.");
  validParamList->set<RCP<const FactoryBase>>("A", null, "Generating factory of the matrix A used during the prolongator smoothing process");
  validParamList->set<RCP<const FactoryBase>>("P", null, "Prolongator factory");
  validParamList->set<RCP<const FactoryBase>>("lCoarseNodesPerDim", null, "Number of nodes per spatial dimension on the coarse grid.");
  validParamList->set<RCP<const FactoryBase>>("structuredInterpolationOrder", null, "Interpolation order used to construct the structured prolongator.");
  validParamList->set<RCP<const FactoryBase>>(
      "matrixType", null, "Galeri matrix type used to infer the structured RAP graph.");

  validParamList->set<bool>("CheckMainDiagonal", false, "Check main diagonal for zeros");
  validParamList->set<bool>("RepairMainDiagonal", false, "Repair zeros on main diagonal");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  const Teuchos::ParameterList& pL = GetParameterList();
  TEUCHOS_TEST_FOR_EXCEPTION(
      !pL.get<bool>("transpose: use implicit"), Exceptions::RuntimeError,
      "StructuredRAPFactory requires \"transpose: use implicit\" = true because "
      "the prebuilt coarse graph assumes the Galerkin product P^T A P.");

  Input(fineLevel, "A");
  Input(coarseLevel, "P");

  const bool prebuildCoarseGraph = pL.get<bool>("rap: prebuild coarse graph");

  if (prebuildCoarseGraph) {
    Input(fineLevel, "lCoarseNodesPerDim");
    Input(fineLevel, "structuredInterpolationOrder");

    if (pL.get<std::string>("rap: matrix type").empty())
      Input(fineLevel, "matrixType");
  }

  // call DeclareInput of all user-given transfer factories
  for (std::vector<RCP<const FactoryBase>>::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it)
    (*it)->CallDeclareInput(coarseLevel);

  hasDeclaredInput_ = true;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetStructured1D(RCP<Matrix>& Ac, RCP<Matrix> P, const Teuchos::Array<LocalOrdinal> lCoarseNodesPerDim, LocalOrdinal dofsPerNode, const std::string& matrixType) const {
  const GO localNx       = Teuchos::as<GO>(lCoarseNodesPerDim[0]);
  const GO dofsPerNodeGO = Teuchos::as<GO>(dofsPerNode);
  TEUCHOS_TEST_FOR_EXCEPTION(localNx <= 0, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured1D(" << matrixType
                                                                      << "): local coarse dimension must be positive.");
  TEUCHOS_TEST_FOR_EXCEPTION(dofsPerNode <= 0, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured1D(" << matrixType
                                                                      << "): dofsPerNode must be positive.");

  RCP<ParameterList> paramList = rcp(new ParameterList);
  paramList->set("No Nonlocal Changes", true);
  paramList->set("Optimize Storage", true);
  paramList->set("compute global constants", true);

  // P->getDomainMap() contains all properties of the coarse dofs owned by this MPI rank.
  auto rowMap                = P->getDomainMap();
  const size_t localNumRows  = rowMap->getLocalNumElements();
  const GO expectedLocalRows = localNx * dofsPerNodeGO;
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<GO>(localNumRows) != expectedLocalRows, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured1D(" << matrixType
                                                                      << "): local coarse dimension " << localNx << " with " << dofsPerNode
                                                                      << " dofs per node does not match local coarse row count " << localNumRows << ".");

  const GO numGlobalRows = Teuchos::as<GO>(rowMap->getGlobalNumElements());
  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalRows % dofsPerNodeGO != 0, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured1D(" << matrixType
                                                                      << "): global coarse row count " << numGlobalRows
                                                                      << " is not divisible by dofsPerNode " << dofsPerNode << ".");

  const GO numGlobalNodes = numGlobalRows / dofsPerNodeGO;
  const GO lowerBound     = rowMap->getMinAllGlobalIndex();
  const GO upperBound     = rowMap->getMaxAllGlobalIndex() + 1;
  TEUCHOS_TEST_FOR_EXCEPTION(upperBound - lowerBound != numGlobalRows, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured1D(" << matrixType
                                                                      << "): coarse row map GIDs must form a contiguous range.");

  const Teuchos::ArrayView<const GO> localRowGids = rowMap->getLocalElementList();
  RCP<const Teuchos::Comm<int>> comm              = rowMap->getComm();
  const int numRanks                              = comm->getSize();

  Teuchos::Array<GO> localRankData(2);
  localRankData[0] = localNx;
  localRankData[1] = Teuchos::as<GO>(localNumRows);
  Teuchos::Array<GO> rankData(2 * numRanks);
  Teuchos::gatherAll(*comm, 2, localRankData.getRawPtr(), 2 * numRanks, rankData.getRawPtr());
  GO globalRankBlockedNx   = 0;
  GO globalRankBlockedRows = 0;
  for (int rank = 0; rank < numRanks; ++rank) {
    globalRankBlockedNx += rankData[2 * rank];
    globalRankBlockedRows += rankData[2 * rank + 1];
  }
  TEUCHOS_TEST_FOR_EXCEPTION(globalRankBlockedNx != numGlobalNodes, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured1D(" << matrixType
                                                                      << "): rank-blocked coarse dimension " << globalRankBlockedNx
                                                                      << " does not match coarse node count " << numGlobalNodes << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(globalRankBlockedRows != numGlobalRows, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured1D(" << matrixType
                                                                      << "): rank-blocked coarse row count " << globalRankBlockedRows
                                                                      << " does not match coarse row count " << numGlobalRows << ".");

  const auto getNodeOrdinal = [&](const GO rowGid) -> GO {
    const GO relativeGid = rowGid - lowerBound;
    TEUCHOS_TEST_FOR_EXCEPTION(relativeGid < 0 || relativeGid >= numGlobalRows, Exceptions::RuntimeError,
                               "StructuredRAPFactory::GetStructured1D(" << matrixType
                                                                        << "): row GID " << rowGid << " is outside the contiguous coarse row range.");
    return relativeGid / dofsPerNodeGO;
  };

  const auto getDofGid = [&](const GO nodeOrdinal, const LO dof) -> GO {
    return lowerBound + nodeOrdinal * dofsPerNodeGO + Teuchos::as<GO>(dof);
  };

  std::vector<GO> remoteColGids;
  remoteColGids.reserve(Teuchos::as<size_t>(2 * dofsPerNode));
  for (int rowLid = 0; rowLid < localRowGids.size(); ++rowLid) {
    const GO rowNode = getNodeOrdinal(localRowGids[rowLid]);
    for (GO colNode = rowNode - 1; colNode <= rowNode + 1; ++colNode) {
      if (colNode < 0 || colNode >= numGlobalNodes)
        continue;
      for (LO colDof = 0; colDof < dofsPerNode; ++colDof) {
        const GO colGid = getDofGid(colNode, colDof);
        if (rowMap->getLocalElement(colGid) == Teuchos::OrdinalTraits<LO>::invalid())
          remoteColGids.push_back(colGid);
      }
    }
  }
  std::sort(remoteColGids.begin(), remoteColGids.end());
  remoteColGids.erase(std::unique(remoteColGids.begin(), remoteColGids.end()), remoteColGids.end());

  Array<GO> colMapGids;
  colMapGids.reserve(localRowGids.size() + Teuchos::as<int>(remoteColGids.size()));
  for (int i = 0; i < localRowGids.size(); ++i)
    colMapGids.push_back(localRowGids[i]);
  for (typename std::vector<GO>::const_iterator it = remoteColGids.begin(); it != remoteColGids.end(); ++it)
    colMapGids.push_back(*it);

  RCP<const Map> colMap = MapFactory::Build(rowMap->lib(),
                                            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                            colMapGids(),
                                            rowMap->getIndexBase(),
                                            comm);

  ArrayRCP<size_t> rowptr(localNumRows + 1);
  rowptr[0]       = 0;
  size_t localNnz = 0;
  for (int rowLid = 0; rowLid < localRowGids.size(); ++rowLid) {
    const GO rowNode = getNodeOrdinal(localRowGids[rowLid]);
    if (rowNode > 0)
      localNnz += Teuchos::as<size_t>(dofsPerNode);
    localNnz += Teuchos::as<size_t>(dofsPerNode);
    if (rowNode + 1 < numGlobalNodes)
      localNnz += Teuchos::as<size_t>(dofsPerNode);
    rowptr[rowLid + 1] = localNnz;
  }

  ArrayRCP<LO> colind(localNnz);
  const auto getColLid = [&](const GO colGid) -> LO {
    const LO colLid = colMap->getLocalElement(colGid);
    TEUCHOS_TEST_FOR_EXCEPTION(colLid == Teuchos::OrdinalTraits<LO>::invalid(), Exceptions::RuntimeError,
                               "StructuredRAPFactory::GetStructured1D(" << matrixType
                                                                        << "): column GID " << colGid
                                                                        << " was not found in the prebuilt coarse column map.");
    return colLid;
  };

  size_t k = 0;
  for (int rowLid = 0; rowLid < localRowGids.size(); ++rowLid) {
    const GO rowNode = getNodeOrdinal(localRowGids[rowLid]);
    for (GO colNode = rowNode - 1; colNode <= rowNode + 1; ++colNode) {
      if (colNode < 0 || colNode >= numGlobalNodes)
        continue;
      for (LO colDof = 0; colDof < dofsPerNode; ++colDof)
        colind[k++] = getColLid(getDofGid(colNode, colDof));
    }
  }
  TEUCHOS_ASSERT(k == localNnz);

  RCP<CrsGraph> myGraph = CrsGraphFactory::Build(rowMap, colMap, rowptr, colind, paramList);
  if (!myGraph->isFillComplete())
    myGraph->fillComplete(rowMap, rowMap, paramList);
  Ac = MatrixFactory::Build(myGraph, paramList);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetLaplace1D(RCP<Matrix>& Ac, RCP<Matrix> P, Teuchos::Array<LocalOrdinal> lCoarseNodesPerDim) const {
  GetStructured1D(Ac, P, lCoarseNodesPerDim, Teuchos::as<LO>(1), "Laplace1D");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetElasticity1D(RCP<Matrix>& Ac, RCP<Matrix> P, Teuchos::Array<LocalOrdinal> lCoarseNodesPerDim) const {
  // Galeri builds exactly the same scalar 1D operator for Laplace1D and Elasticity1D, so we can just call GetStructured1D with dofsPerNode=1
  GetStructured1D(Ac, P, lCoarseNodesPerDim, Teuchos::as<LO>(1), "Elasticity1D");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetStructured2D(RCP<Matrix>& Ac, RCP<Matrix> P, Teuchos::Array<LocalOrdinal> lCoarseNodesPerDim, LocalOrdinal dofsPerNode, bool includeDiagonalNeighbors, int procX, int procY, const std::string& matrixType) const {
  const GO localNx       = Teuchos::as<GO>(lCoarseNodesPerDim[0]);
  const GO localNy       = Teuchos::as<GO>(lCoarseNodesPerDim[1]);
  const GO dofsPerNodeGO = Teuchos::as<GO>(dofsPerNode);

  TEUCHOS_TEST_FOR_EXCEPTION(localNx <= 0 || localNy <= 0, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured2D(" << matrixType
                                                                      << "): local coarse dimensions must be positive.");
  TEUCHOS_TEST_FOR_EXCEPTION(dofsPerNode <= 0, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured2D(" << matrixType
                                                                      << "): dofsPerNode must be positive.");

  RCP<ParameterList> paramList = rcp(new ParameterList);
  paramList->set("No Nonlocal Changes", true);
  paramList->set("Optimize Storage", true);
  paramList->set("compute global constants", true);

  // P->getDomainMap() contains all properties of the coarse dofs owned by this MPI rank.
  auto rowMap                                     = P->getDomainMap();
  const size_t localNumRows                       = rowMap->getLocalNumElements();
  const Teuchos::ArrayView<const GO> localRowGids = rowMap->getLocalElementList();
  const GO numGlobalCoarseDofs                    = Teuchos::as<GO>(rowMap->getGlobalNumElements());
  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalCoarseDofs % dofsPerNodeGO != 0, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured2D(" << matrixType
                                                                      << "): coarse domain map size " << numGlobalCoarseDofs
                                                                      << " is not divisible by dofsPerNode " << dofsPerNode << ".");
  const GO numGlobalCoarseNodes = numGlobalCoarseDofs / dofsPerNodeGO;
  const GO expectedLocalRows    = localNx * localNy * dofsPerNodeGO;
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<GO>(localNumRows) != expectedLocalRows, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured2D(" << matrixType
                                                                      << "): local coarse dimensions " << localNx << "x" << localNy
                                                                      << " with " << dofsPerNode << " dofs per node do not match local coarse row count "
                                                                      << localNumRows << ".");

  const GO localMinGid       = rowMap->getMinGlobalIndex();
  const GO localMaxGid       = rowMap->getMaxGlobalIndex();
  bool localRowMapContiguous = true;
  for (int i = 0; i < localRowGids.size(); ++i) {
    if (localRowGids[i] != localMinGid + Teuchos::as<GO>(i)) {
      localRowMapContiguous = false;
      break;
    }
  }

  RCP<const Teuchos::Comm<int>> comm = rowMap->getComm();
  const int myRank                   = comm->getRank();
  const int numRanks                 = comm->getSize();
  if (procX <= 0) procX = 1;
  if (procY <= 0) procY = numRanks / procX;
  TEUCHOS_TEST_FOR_EXCEPTION(procX <= 0 || procY <= 0 || procX * procY != numRanks, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured2D(" << matrixType
                                                                      << "): invalid processor grid " << procX << "x" << procY
                                                                      << " for " << numRanks << " ranks.");

  const int myProcX        = myRank % procX;
  const int myProcY        = myRank / procX;
  const size_t rowsPerNode = Teuchos::as<size_t>(dofsPerNode);
  const GO localMinNodeGid = localMinGid / dofsPerNodeGO;
  Teuchos::Array<GO> localRankData(4);
  localRankData[0] = localMinNodeGid;
  localRankData[1] = localNx;
  localRankData[2] = localNy;
  localRankData[3] = Teuchos::as<GO>(localNumRows / rowsPerNode);
  Teuchos::Array<GO> rankData(4 * numRanks);
  Teuchos::gatherAll(*comm, 4, localRankData.getRawPtr(), 4 * numRanks, rankData.getRawPtr());

  GO globalRankBlockedNx = 0;
  for (int px = 0; px < procX; ++px)
    globalRankBlockedNx += rankData[4 * px + 1];
  GO globalRankBlockedNy = 0;
  for (int py = 0; py < procY; ++py)
    globalRankBlockedNy += rankData[4 * (py * procX) + 2];
  TEUCHOS_TEST_FOR_EXCEPTION(globalRankBlockedNx * globalRankBlockedNy != numGlobalCoarseNodes, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured2D(" << matrixType
                                                                      << "): processor-grid coarse dimensions " << globalRankBlockedNx << "x" << globalRankBlockedNy
                                                                      << " do not match coarse node count " << numGlobalCoarseNodes << ".");

  const bool pairedContiguousRows =
      localRowMapContiguous &&
      (localMinGid % dofsPerNodeGO == 0 && localNumRows % rowsPerNode == 0);

  const auto useStencilOffset = [&](const int dx, const int dy) -> bool {
    return includeDiagonalNeighbors || dx == 0 || dy == 0;
  };

  const auto resolveNeighbor = [&](const GO x, const GO y, const int dx, const int dy,
                                   int& neighborRank, GO& neighborX, GO& neighborY) -> bool {
    int neighborProcX = myProcX;
    int neighborProcY = myProcY;
    neighborX         = x + Teuchos::as<GO>(dx);
    neighborY         = y + Teuchos::as<GO>(dy);

    if (neighborX < 0)
      --neighborProcX;
    else if (neighborX >= localNx)
      ++neighborProcX;
    if (neighborY < 0)
      --neighborProcY;
    else if (neighborY >= localNy)
      ++neighborProcY;

    if (neighborProcX < 0 || neighborProcX >= procX || neighborProcY < 0 || neighborProcY >= procY)
      return false;

    neighborRank             = neighborProcY * procX + neighborProcX;
    const GO neighborLocalNx = rankData[4 * neighborRank + 1];
    const GO neighborLocalNy = rankData[4 * neighborRank + 2];
    if (neighborX < 0)
      neighborX = neighborLocalNx - 1;
    else if (neighborX >= localNx)
      neighborX = 0;
    if (neighborY < 0)
      neighborY = neighborLocalNy - 1;
    else if (neighborY >= localNy)
      neighborY = 0;

    return true;
  };

  const auto countStencilColumns = [&](const GO x, const GO y) -> size_t {
    size_t nnz = 0;
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dx = -1; dx <= 1; ++dx) {
        if (!useStencilOffset(dx, dy))
          continue;
        int neighborRank = myRank;
        GO neighborX     = x;
        GO neighborY     = y;
        if (resolveNeighbor(x, y, dx, dy, neighborRank, neighborX, neighborY))
          nnz += rowsPerNode;
      }
    }
    return nnz;
  };

  ArrayRCP<size_t> rowptr(localNumRows + 1);
  rowptr[0] = 0;

  size_t localNnz    = 0;
  GO previousNodeGid = Teuchos::OrdinalTraits<GO>::invalid();
  size_t previousNnz = 0;

  for (size_t rowLid = 0; rowLid < localNumRows; rowLid += (pairedContiguousRows ? rowsPerNode : 1)) {
    const GO rowGid  = localRowMapContiguous ? localMinGid + Teuchos::as<GO>(rowLid) : rowMap->getGlobalElement(Teuchos::as<LO>(rowLid));
    const GO nodeGid = rowGid / dofsPerNodeGO;

    if (!pairedContiguousRows && nodeGid == previousNodeGid) {
      localNnz += previousNnz;
      rowptr[rowLid + 1] = localNnz;
      continue;
    }

    const GO localNodeLid = nodeGid - localMinNodeGid;
    const GO x            = localNodeLid % localNx;
    const GO y            = localNodeLid / localNx;
    const size_t nnz      = countStencilColumns(x, y);

    previousNodeGid = nodeGid;
    previousNnz     = nnz;
    if (pairedContiguousRows) {
      for (size_t rowDof = 0; rowDof < rowsPerNode; ++rowDof) {
        localNnz += nnz;
        rowptr[rowLid + rowDof + 1] = localNnz;
      }
    } else {
      localNnz += nnz;
      rowptr[rowLid + 1] = localNnz;
    }
  }

  std::vector<GO> remoteColGids;
  const size_t maxNodalNeighbors = includeDiagonalNeighbors ? 9 : 5;
  remoteColGids.reserve(Teuchos::as<size_t>((2 * localNx + 2 * localNy) * dofsPerNodeGO * Teuchos::as<GO>(maxNodalNeighbors)));

  const auto addRemoteColumnsForNode = [&](const GO x, const GO y) {
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dx = -1; dx <= 1; ++dx) {
        if (!useStencilOffset(dx, dy))
          continue;

        int neighborRank = myRank;
        GO neighborX     = x;
        GO neighborY     = y;
        if (!resolveNeighbor(x, y, dx, dy, neighborRank, neighborX, neighborY))
          continue;
        if (neighborRank == myRank)
          continue;

        const GO neighborLocalNx = rankData[4 * neighborRank + 1];
        const GO colNodeGid      = rankData[4 * neighborRank] + neighborY * neighborLocalNx + neighborX;
        for (LO colDof = 0; colDof < dofsPerNode; ++colDof)
          remoteColGids.push_back(dofsPerNodeGO * colNodeGid + Teuchos::as<GO>(colDof));
      }
    }
  };

  for (GO x = 0; x < localNx; ++x)
    addRemoteColumnsForNode(x, 0);
  if (localNy > 1)
    for (GO x = 0; x < localNx; ++x)
      addRemoteColumnsForNode(x, localNy - 1);
  for (GO y = 1; y + 1 < localNy; ++y) {
    addRemoteColumnsForNode(0, y);
    if (localNx > 1)
      addRemoteColumnsForNode(localNx - 1, y);
  }

  std::sort(remoteColGids.begin(), remoteColGids.end());
  remoteColGids.erase(std::unique(remoteColGids.begin(), remoteColGids.end()), remoteColGids.end());

  Array<GO> colMapGids;
  colMapGids.reserve(localRowGids.size() + Teuchos::as<int>(remoteColGids.size()));
  for (int i = 0; i < localRowGids.size(); ++i)
    colMapGids.push_back(localRowGids[i]);
  for (typename std::vector<GO>::const_iterator it = remoteColGids.begin(); it != remoteColGids.end(); ++it)
    colMapGids.push_back(*it);

  RCP<const Map> colMap = MapFactory::Build(rowMap->lib(),
                                            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                            colMapGids(),
                                            rowMap->getIndexBase(),
                                            rowMap->getComm());

  ArrayRCP<LO> colind(localNnz);
  const size_t maxColumnsPerRow = rowsPerNode * maxNodalNeighbors;
  Array<LO> colLids(Teuchos::as<int>(maxColumnsPerRow));
  Array<LO> previousColLids(Teuchos::as<int>(maxColumnsPerRow));
  size_t k        = 0;
  previousNodeGid = Teuchos::OrdinalTraits<GO>::invalid();
  previousNnz     = 0;

  for (size_t rowLid = 0; rowLid < localNumRows; rowLid += (pairedContiguousRows ? rowsPerNode : 1)) {
    const GO rowGid  = localRowMapContiguous ? localMinGid + Teuchos::as<GO>(rowLid) : rowMap->getGlobalElement(Teuchos::as<LO>(rowLid));
    const GO nodeGid = rowGid / dofsPerNodeGO;

    if (!pairedContiguousRows && nodeGid == previousNodeGid) {
      for (size_t i = 0; i < previousNnz; ++i)
        colind[k++] = previousColLids[i];
      continue;
    }

    const GO localNodeLid = nodeGid - localMinNodeGid;
    const GO x            = localNodeLid % localNx;
    const GO y            = localNodeLid / localNx;

    if (pairedContiguousRows && x > 0 && x + 1 < localNx && y > 0 && y + 1 < localNy) {
      size_t nnz = 0;
      for (int dy = -1; dy <= 1; ++dy) {
        const GO yy = y + Teuchos::as<GO>(dy);
        for (int dx = -1; dx <= 1; ++dx) {
          if (!useStencilOffset(dx, dy))
            continue;
          const GO xx      = x + Teuchos::as<GO>(dx);
          const LO colBase = Teuchos::as<LO>(dofsPerNodeGO * (yy * localNx + xx));
          for (LO colDof = 0; colDof < dofsPerNode; ++colDof)
            colLids[nnz++] = colBase + colDof;
        }
      }

      previousNodeGid = nodeGid;
      previousNnz     = nnz;
      for (size_t rowDof = 0; rowDof < rowsPerNode; ++rowDof)
        for (size_t i = 0; i < nnz; ++i)
          colind[k++] = colLids[i];
      continue;
    }

    size_t nnz = 0;
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dx = -1; dx <= 1; ++dx) {
        if (!useStencilOffset(dx, dy))
          continue;

        int neighborRank = myRank;
        GO neighborX     = x;
        GO neighborY     = y;
        if (!resolveNeighbor(x, y, dx, dy, neighborRank, neighborX, neighborY))
          continue;

        const GO neighborLocalNx = rankData[4 * neighborRank + 1];
        const GO colNodeGid      = rankData[4 * neighborRank] + neighborY * neighborLocalNx + neighborX;
        for (LO colDof = 0; colDof < dofsPerNode; ++colDof) {
          const GO colGid          = dofsPerNodeGO * colNodeGid + Teuchos::as<GO>(colDof);
          const bool localFastPath = localRowMapContiguous && colGid >= localMinGid && colGid <= localMaxGid;
          const LO colLid          = localFastPath ? Teuchos::as<LO>(colGid - localMinGid) : colMap->getLocalElement(colGid);
          TEUCHOS_TEST_FOR_EXCEPTION(colLid == Teuchos::OrdinalTraits<LO>::invalid(), Exceptions::RuntimeError,
                                     "StructuredRAPFactory::GetStructured2D(" << matrixType
                                                                              << "): column GID " << colGid
                                                                              << " was not found in the prebuilt coarse column map.");
          colLids[nnz++] = colLid;
        }
      }
    }

    std::sort(colLids.getRawPtr(), colLids.getRawPtr() + nnz);

    previousNodeGid = nodeGid;
    previousNnz     = nnz;
    if (pairedContiguousRows) {
      for (size_t rowDof = 0; rowDof < rowsPerNode; ++rowDof)
        for (size_t i = 0; i < nnz; ++i)
          colind[k++] = colLids[i];
    } else {
      for (size_t i = 0; i < nnz; ++i) {
        previousColLids[i] = colLids[i];
        colind[k++]        = colLids[i];
      }
    }
  }
  TEUCHOS_ASSERT(k == localNnz);

  RCP<CrsGraph> myGraph = CrsGraphFactory::Build(rowMap, colMap, rowptr, colind, paramList);
  if (!myGraph->isFillComplete())
    myGraph->fillComplete(rowMap, rowMap, paramList);
  Ac = MatrixFactory::Build(myGraph, paramList);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetLaplace2D(RCP<Matrix>& Ac, RCP<Matrix> P, Teuchos::Array<LocalOrdinal> lCoarseNodesPerDim, int procX, int procY, int interpolationOrder) const {
  TEUCHOS_TEST_FOR_EXCEPTION(interpolationOrder < 0 || interpolationOrder > 1, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetLaplace2D: interpolation order "
                                 << interpolationOrder << " is not supported.");
  // Piecewise-linear interpolation can create diagonal coarse couplings for scalar Laplace2D
  // So we set includeDiagonalNeighbors to true for interpolationOrder=1, and false for interpolationOrder=0
  const bool includeDiagonalNeighbors = (interpolationOrder == 1);
  GetStructured2D(Ac, P, lCoarseNodesPerDim, Teuchos::as<LO>(1), includeDiagonalNeighbors, procX, procY, "Laplace2D");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetElasticity2D(RCP<Matrix>& Ac, RCP<Matrix> P, Teuchos::Array<LocalOrdinal> lCoarseNodesPerDim, int procX, int procY) const {
  // Call GetStructured2D with dofsPerNode=2 and includeDiagonalNeighbors=true for Elasticity2D (9 point stencil)
  GetStructured2D(Ac, P, lCoarseNodesPerDim, Teuchos::as<LO>(2), true, procX, procY, "Elasticity2D");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetStructured3D(RCP<Matrix>& Ac, RCP<Matrix> P, Teuchos::Array<LocalOrdinal> lCoarseNodesPerDim, LocalOrdinal dofsPerNode, bool includeDiagonalNeighbors, int procX, int procY, int procZ, const std::string& matrixType) const {
  const GO localNx       = Teuchos::as<GO>(lCoarseNodesPerDim[0]);
  const GO localNy       = Teuchos::as<GO>(lCoarseNodesPerDim[1]);
  const GO localNz       = Teuchos::as<GO>(lCoarseNodesPerDim[2]);
  const GO dofsPerNodeGO = Teuchos::as<GO>(dofsPerNode);

  TEUCHOS_TEST_FOR_EXCEPTION(localNx <= 0 || localNy <= 0 || localNz <= 0, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured3D(" << matrixType
                                                                      << "): local coarse dimensions must be positive.");
  TEUCHOS_TEST_FOR_EXCEPTION(dofsPerNode <= 0, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured3D(" << matrixType
                                                                      << "): dofsPerNode must be positive.");

  RCP<ParameterList> paramList = rcp(new ParameterList);
  paramList->set("No Nonlocal Changes", true);
  paramList->set("Optimize Storage", true);
  paramList->set("compute global constants", true);

  // P->getDomainMap() contains all properties of the coarse dofs owned by this MPI rank.
  auto rowMap                                     = P->getDomainMap();
  const size_t localNumRows                       = rowMap->getLocalNumElements();
  const Teuchos::ArrayView<const GO> localRowGids = rowMap->getLocalElementList();
  const GO numGlobalCoarseDofs                    = Teuchos::as<GO>(rowMap->getGlobalNumElements());
  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalCoarseDofs % dofsPerNodeGO != 0, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured3D(" << matrixType
                                                                      << "): coarse domain map size " << numGlobalCoarseDofs
                                                                      << " is not divisible by dofsPerNode " << dofsPerNode << ".");
  const GO numGlobalCoarseNodes = numGlobalCoarseDofs / dofsPerNodeGO;
  const GO expectedLocalRows    = localNx * localNy * localNz * dofsPerNodeGO;
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<GO>(localNumRows) != expectedLocalRows, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured3D(" << matrixType
                                                                      << "): local coarse dimensions " << localNx << "x" << localNy << "x" << localNz
                                                                      << " with " << dofsPerNode << " dofs per node do not match local coarse row count "
                                                                      << localNumRows << ".");

  const GO localMinGid       = rowMap->getMinGlobalIndex();
  const GO localMaxGid       = rowMap->getMaxGlobalIndex();
  bool localRowMapContiguous = true;
  for (int i = 0; i < localRowGids.size(); ++i) {
    if (localRowGids[i] != localMinGid + Teuchos::as<GO>(i)) {
      localRowMapContiguous = false;
      break;
    }
  }

  RCP<const Teuchos::Comm<int>> comm = rowMap->getComm();
  const int myRank                   = comm->getRank();
  const int numRanks                 = comm->getSize();
  if (procX <= 0) procX = 1;
  if (procY <= 0) procY = 1;
  if (procZ <= 0) procZ = numRanks / (procX * procY);
  TEUCHOS_TEST_FOR_EXCEPTION(procX <= 0 || procY <= 0 || procZ <= 0 || procX * procY * procZ != numRanks,
                             Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured3D(" << matrixType
                                                                      << "): invalid processor grid " << procX << "x" << procY << "x" << procZ
                                                                      << " for " << numRanks << " ranks.");

  const int procXY         = procX * procY;
  const int myProcX        = myRank % procX;
  const int myProcY        = (myRank % procXY) / procX;
  const int myProcZ        = myRank / procXY;
  const size_t rowsPerNode = Teuchos::as<size_t>(dofsPerNode);
  const GO localMinNodeGid = localMinGid / dofsPerNodeGO;
  Teuchos::Array<GO> localRankData(5);
  localRankData[0] = localMinNodeGid;
  localRankData[1] = localNx;
  localRankData[2] = localNy;
  localRankData[3] = localNz;
  localRankData[4] = Teuchos::as<GO>(localNumRows / rowsPerNode);
  Teuchos::Array<GO> rankData(5 * numRanks);
  Teuchos::gatherAll(*comm, 5, localRankData.getRawPtr(), 5 * numRanks, rankData.getRawPtr());

  GO globalRankBlockedNx = 0;
  for (int px = 0; px < procX; ++px)
    globalRankBlockedNx += rankData[5 * px + 1];
  GO globalRankBlockedNy = 0;
  for (int py = 0; py < procY; ++py)
    globalRankBlockedNy += rankData[5 * (py * procX) + 2];
  GO globalRankBlockedNz = 0;
  for (int pz = 0; pz < procZ; ++pz)
    globalRankBlockedNz += rankData[5 * (pz * procXY) + 3];
  TEUCHOS_TEST_FOR_EXCEPTION(globalRankBlockedNx * globalRankBlockedNy * globalRankBlockedNz != numGlobalCoarseNodes,
                             Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructured3D(" << matrixType
                                                                      << "): processor-grid coarse dimensions " << globalRankBlockedNx << "x"
                                                                      << globalRankBlockedNy << "x" << globalRankBlockedNz
                                                                      << " do not match coarse node count " << numGlobalCoarseNodes << ".");

  const bool groupedContiguousRows =
      localRowMapContiguous &&
      (localMinGid % dofsPerNodeGO == 0 && localNumRows % rowsPerNode == 0);

  const auto useStencilOffset = [&](const int dx, const int dy, const int dz) -> bool {
    return includeDiagonalNeighbors ||
           (dx == 0 && dy == 0) ||
           (dx == 0 && dz == 0) ||
           (dy == 0 && dz == 0);
  };

  const auto resolveNeighbor = [&](const GO x, const GO y, const GO z,
                                   const int dx, const int dy, const int dz,
                                   int& neighborRank, GO& neighborX, GO& neighborY, GO& neighborZ) -> bool {
    int neighborProcX = myProcX;
    int neighborProcY = myProcY;
    int neighborProcZ = myProcZ;
    neighborX         = x + Teuchos::as<GO>(dx);
    neighborY         = y + Teuchos::as<GO>(dy);
    neighborZ         = z + Teuchos::as<GO>(dz);

    if (neighborX < 0)
      --neighborProcX;
    else if (neighborX >= localNx)
      ++neighborProcX;
    if (neighborY < 0)
      --neighborProcY;
    else if (neighborY >= localNy)
      ++neighborProcY;
    if (neighborZ < 0)
      --neighborProcZ;
    else if (neighborZ >= localNz)
      ++neighborProcZ;

    if (neighborProcX < 0 || neighborProcX >= procX ||
        neighborProcY < 0 || neighborProcY >= procY ||
        neighborProcZ < 0 || neighborProcZ >= procZ)
      return false;

    neighborRank             = neighborProcZ * procXY + neighborProcY * procX + neighborProcX;
    const GO neighborLocalNx = rankData[5 * neighborRank + 1];
    const GO neighborLocalNy = rankData[5 * neighborRank + 2];
    const GO neighborLocalNz = rankData[5 * neighborRank + 3];
    if (neighborX < 0)
      neighborX = neighborLocalNx - 1;
    else if (neighborX >= localNx)
      neighborX = 0;
    if (neighborY < 0)
      neighborY = neighborLocalNy - 1;
    else if (neighborY >= localNy)
      neighborY = 0;
    if (neighborZ < 0)
      neighborZ = neighborLocalNz - 1;
    else if (neighborZ >= localNz)
      neighborZ = 0;

    return true;
  };

  const auto countStencilColumns = [&](const GO x, const GO y, const GO z) -> size_t {
    size_t nnz = 0;
    for (int dz = -1; dz <= 1; ++dz) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
          if (!useStencilOffset(dx, dy, dz))
            continue;
          int neighborRank = myRank;
          GO neighborX     = x;
          GO neighborY     = y;
          GO neighborZ     = z;
          if (resolveNeighbor(x, y, z, dx, dy, dz, neighborRank, neighborX, neighborY, neighborZ))
            nnz += rowsPerNode;
        }
      }
    }
    return nnz;
  };

  ArrayRCP<size_t> rowptr(localNumRows + 1);
  rowptr[0] = 0;

  size_t localNnz    = 0;
  GO previousNodeGid = Teuchos::OrdinalTraits<GO>::invalid();
  size_t previousNnz = 0;

  for (size_t rowLid = 0; rowLid < localNumRows; rowLid += (groupedContiguousRows ? rowsPerNode : 1)) {
    const GO rowGid  = localRowMapContiguous ? localMinGid + Teuchos::as<GO>(rowLid) : rowMap->getGlobalElement(Teuchos::as<LO>(rowLid));
    const GO nodeGid = rowGid / dofsPerNodeGO;

    if (!groupedContiguousRows && nodeGid == previousNodeGid) {
      localNnz += previousNnz;
      rowptr[rowLid + 1] = localNnz;
      continue;
    }

    const GO localNodeLid = nodeGid - localMinNodeGid;
    const GO x            = localNodeLid % localNx;
    const GO y            = (localNodeLid / localNx) % localNy;
    const GO z            = localNodeLid / (localNx * localNy);
    const size_t nnz      = countStencilColumns(x, y, z);

    previousNodeGid = nodeGid;
    previousNnz     = nnz;
    if (groupedContiguousRows) {
      for (size_t rowDof = 0; rowDof < rowsPerNode; ++rowDof) {
        localNnz += nnz;
        rowptr[rowLid + rowDof + 1] = localNnz;
      }
    } else {
      localNnz += nnz;
      rowptr[rowLid + 1] = localNnz;
    }
  }

  std::vector<GO> remoteColGids;
  const auto addRemoteColumnsForNode = [&](const GO x, const GO y, const GO z) {
    for (int dz = -1; dz <= 1; ++dz) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
          if (!useStencilOffset(dx, dy, dz))
            continue;

          int neighborRank = myRank;
          GO neighborX     = x;
          GO neighborY     = y;
          GO neighborZ     = z;
          if (!resolveNeighbor(x, y, z, dx, dy, dz, neighborRank, neighborX, neighborY, neighborZ))
            continue;
          if (neighborRank == myRank)
            continue;

          const GO neighborLocalNx = rankData[5 * neighborRank + 1];
          const GO neighborLocalNy = rankData[5 * neighborRank + 2];
          const GO colNodeGid      = rankData[5 * neighborRank] +
                                neighborZ * neighborLocalNx * neighborLocalNy +
                                neighborY * neighborLocalNx + neighborX;
          for (LO colDof = 0; colDof < dofsPerNode; ++colDof)
            remoteColGids.push_back(dofsPerNodeGO * colNodeGid + Teuchos::as<GO>(colDof));
        }
      }
    }
  };

  for (GO z = 0; z < localNz; ++z) {
    for (GO y = 0; y < localNy; ++y) {
      for (GO x = 0; x < localNx; ++x) {
        if (x == 0 || x + 1 == localNx ||
            y == 0 || y + 1 == localNy ||
            z == 0 || z + 1 == localNz)
          addRemoteColumnsForNode(x, y, z);
      }
    }
  }

  std::sort(remoteColGids.begin(), remoteColGids.end());
  remoteColGids.erase(std::unique(remoteColGids.begin(), remoteColGids.end()), remoteColGids.end());

  Array<GO> colMapGids;
  colMapGids.reserve(localRowGids.size() + Teuchos::as<int>(remoteColGids.size()));
  for (int i = 0; i < localRowGids.size(); ++i)
    colMapGids.push_back(localRowGids[i]);
  for (typename std::vector<GO>::const_iterator it = remoteColGids.begin(); it != remoteColGids.end(); ++it)
    colMapGids.push_back(*it);

  RCP<const Map> colMap = MapFactory::Build(rowMap->lib(),
                                            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                            colMapGids(),
                                            rowMap->getIndexBase(),
                                            rowMap->getComm());

  ArrayRCP<LO> colind(localNnz);
  const size_t maxNodalNeighbors = includeDiagonalNeighbors ? 27 : 7;
  const size_t maxColumnsPerRow  = rowsPerNode * maxNodalNeighbors;
  Array<LO> colLids(Teuchos::as<int>(maxColumnsPerRow));
  Array<LO> previousColLids(Teuchos::as<int>(maxColumnsPerRow));
  size_t k        = 0;
  previousNodeGid = Teuchos::OrdinalTraits<GO>::invalid();
  previousNnz     = 0;

  for (size_t rowLid = 0; rowLid < localNumRows; rowLid += (groupedContiguousRows ? rowsPerNode : 1)) {
    const GO rowGid  = localRowMapContiguous ? localMinGid + Teuchos::as<GO>(rowLid) : rowMap->getGlobalElement(Teuchos::as<LO>(rowLid));
    const GO nodeGid = rowGid / dofsPerNodeGO;

    if (!groupedContiguousRows && nodeGid == previousNodeGid) {
      for (size_t i = 0; i < previousNnz; ++i)
        colind[k++] = previousColLids[i];
      continue;
    }

    const GO localNodeLid = nodeGid - localMinNodeGid;
    const GO x            = localNodeLid % localNx;
    const GO y            = (localNodeLid / localNx) % localNy;
    const GO z            = localNodeLid / (localNx * localNy);

    if (groupedContiguousRows &&
        x > 0 && x + 1 < localNx &&
        y > 0 && y + 1 < localNy &&
        z > 0 && z + 1 < localNz) {
      size_t nnz = 0;
      for (int dz = -1; dz <= 1; ++dz) {
        const GO zz = z + Teuchos::as<GO>(dz);
        for (int dy = -1; dy <= 1; ++dy) {
          const GO yy = y + Teuchos::as<GO>(dy);
          for (int dx = -1; dx <= 1; ++dx) {
            if (!useStencilOffset(dx, dy, dz))
              continue;
            const GO xx      = x + Teuchos::as<GO>(dx);
            const LO colBase = Teuchos::as<LO>(dofsPerNodeGO *
                                               (zz * localNx * localNy + yy * localNx + xx));
            for (LO colDof = 0; colDof < dofsPerNode; ++colDof)
              colLids[nnz++] = colBase + colDof;
          }
        }
      }

      previousNodeGid = nodeGid;
      previousNnz     = nnz;
      for (size_t rowDof = 0; rowDof < rowsPerNode; ++rowDof)
        for (size_t i = 0; i < nnz; ++i)
          colind[k++] = colLids[i];
      continue;
    }

    size_t nnz = 0;
    for (int dz = -1; dz <= 1; ++dz) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
          if (!useStencilOffset(dx, dy, dz))
            continue;

          int neighborRank = myRank;
          GO neighborX     = x;
          GO neighborY     = y;
          GO neighborZ     = z;
          if (!resolveNeighbor(x, y, z, dx, dy, dz, neighborRank, neighborX, neighborY, neighborZ))
            continue;

          const GO neighborLocalNx = rankData[5 * neighborRank + 1];
          const GO neighborLocalNy = rankData[5 * neighborRank + 2];
          const GO colNodeGid      = rankData[5 * neighborRank] +
                                neighborZ * neighborLocalNx * neighborLocalNy +
                                neighborY * neighborLocalNx + neighborX;
          for (LO colDof = 0; colDof < dofsPerNode; ++colDof) {
            const GO colGid          = dofsPerNodeGO * colNodeGid + Teuchos::as<GO>(colDof);
            const bool localFastPath = localRowMapContiguous && colGid >= localMinGid && colGid <= localMaxGid;
            const LO colLid          = localFastPath ? Teuchos::as<LO>(colGid - localMinGid) : colMap->getLocalElement(colGid);
            TEUCHOS_TEST_FOR_EXCEPTION(colLid == Teuchos::OrdinalTraits<LO>::invalid(), Exceptions::RuntimeError,
                                       "StructuredRAPFactory::GetStructured3D(" << matrixType
                                                                                << "): column GID " << colGid
                                                                                << " was not found in the prebuilt coarse column map.");
            colLids[nnz++] = colLid;
          }
        }
      }
    }

    std::sort(colLids.getRawPtr(), colLids.getRawPtr() + nnz);

    previousNodeGid = nodeGid;
    previousNnz     = nnz;
    if (groupedContiguousRows) {
      for (size_t rowDof = 0; rowDof < rowsPerNode; ++rowDof)
        for (size_t i = 0; i < nnz; ++i)
          colind[k++] = colLids[i];
    } else {
      for (size_t i = 0; i < nnz; ++i) {
        previousColLids[i] = colLids[i];
        colind[k++]        = colLids[i];
      }
    }
  }
  TEUCHOS_ASSERT(k == localNnz);

  RCP<CrsGraph> myGraph = CrsGraphFactory::Build(rowMap, colMap, rowptr, colind, paramList);
  if (!myGraph->isFillComplete())
    myGraph->fillComplete(rowMap, rowMap, paramList);
  Ac = MatrixFactory::Build(myGraph, paramList);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetLaplace3D(RCP<Matrix>& Ac, RCP<Matrix> P, Teuchos::Array<LocalOrdinal> lCoarseNodesPerDim, int procX, int procY, int procZ, int interpolationOrder) const {
  TEUCHOS_TEST_FOR_EXCEPTION(interpolationOrder < 0 || interpolationOrder > 1, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetLaplace3D: interpolation order "
                                 << interpolationOrder << " is not supported.");
  // Piecewise-linear interpolation can create edge and corner coarse couplings for scalar Laplace3D.
  const bool includeDiagonalNeighbors = (interpolationOrder == 1);
  GetStructured3D(Ac, P, lCoarseNodesPerDim, Teuchos::as<LO>(1), includeDiagonalNeighbors,
                  procX, procY, procZ, "Laplace3D");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetElasticity3D(RCP<Matrix>& Ac, RCP<Matrix> P, Teuchos::Array<LocalOrdinal> lCoarseNodesPerDim, int procX, int procY, int procZ) const {
  // Elasticity3D has a full 27-point nodal stencil with three dofs per node.
  GetStructured3D(Ac, P, lCoarseNodesPerDim, Teuchos::as<LO>(3), true,
                  procX, procY, procZ, "Elasticity3D");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  const bool doTranspose           = true;
  const bool doFillComplete        = true;
  const bool doOptimizeStorage     = true;
  const Teuchos::ParameterList& pL = GetParameterList();
  const bool prebuildCoarseGraph   = pL.get<bool>("rap: prebuild coarse graph");

  TEUCHOS_TEST_FOR_EXCEPTION(
      !pL.get<bool>("transpose: use implicit"), Exceptions::RuntimeError,
      "StructuredRAPFactory requires \"transpose: use implicit\" = true because "
      "the prebuilt coarse graph assumes the Galerkin product P^T A P.");

  RCP<Matrix> Ac;

  {
    FactoryMonitor m(*this, "Computing Ac", coarseLevel);

    std::ostringstream levelstr;
    levelstr << coarseLevel.GetLevelID();
    std::string labelstr = FormattingHelper::getColonLabel(coarseLevel.getObjectLabel());

    TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_ == false, Exceptions::RuntimeError,
                               "MueLu::RAPFactory::Build(): CallDeclareInput has not been called before Build!");

    RCP<Matrix> A = Get<RCP<Matrix>>(fineLevel, "A");
    RCP<Matrix> P = Get<RCP<Matrix>>(coarseLevel, "P");
    // We don't have a valid P (e.g., # global aggregates = 0) so we bail.
    // This level will ultimately be removed in MueLu_Hierarchy_defs.h via a resize()
    if (P.is_null()) {
      Ac = Teuchos::null;
      Set(coarseLevel, "A", Ac);
      return;
    }

#ifdef KOKKOS_ENABLE_CUDA
    bool isCuda = typeid(Node).name() == typeid(Kokkos::Compat::KokkosCudaWrapperNode).name();
#else
    bool isCuda = false;
#endif

    TEUCHOS_TEST_FOR_EXCEPTION(
        pL.get<bool>("rap: triple product") == false, Exceptions::RuntimeError,
        "StructuredRAPFactory requires \"rap: triple product\" = true.");
    TEUCHOS_TEST_FOR_EXCEPTION(
        isCuda, Exceptions::RuntimeError,
        "StructuredRAPFactory does not currently support CUDA.");

    {
      RCP<ParameterList> RAPparams = rcp(new ParameterList);
      if (pL.isSublist("matrixmatrix: kernel params"))
        RAPparams->sublist("matrixmatrix: kernel params") = pL.sublist("matrixmatrix: kernel params");

      if (coarseLevel.IsAvailable("RAP reuse data", this)) {
        GetOStream(static_cast<MsgType>(Runtime0 | Test)) << "Reusing previous RAP data" << std::endl;

        RAPparams = coarseLevel.Get<RCP<ParameterList>>("RAP reuse data", this);

        if (RAPparams->isParameter("graph"))
          Ac = RAPparams->get<RCP<Matrix>>("graph");

        // Some eigenvalue may have been cached with the matrix in the previous run.
        // As the matrix values will be updated, we need to reset the eigenvalue.
        Ac->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());

        // If we want to prebuild the coarse graph, do that here. Otherwise, we will get it in the symbolic phase of the triple matrix product,
        // but that will be more expensive
      } else if (prebuildCoarseGraph) {
        // if reuse data not available, try to get sparse fill graph via the knowledge of the matrix structure
        std::string matrixType = pL.get<std::string>("rap: matrix type");
        if (matrixType.empty())
          matrixType = Get<std::string>(fineLevel, "matrixType");
        Teuchos::Array<LocalOrdinal> lCoarseNodesPerDim =
            Get<Teuchos::Array<LocalOrdinal>>(fineLevel, "lCoarseNodesPerDim");
        const int interpolationOrder = Get<int>(fineLevel, "structuredInterpolationOrder");
        const int procX              = pL.get<int>("rap: processor grid x");
        const int procY              = pL.get<int>("rap: processor grid y");
        const int procZ              = pL.get<int>("rap: processor grid z");

        if (matrixType == "Laplace1D") {
          GetOStream(Statistics1) << "StructuredRAP: Using Laplace1D pattern determination routine." << std::endl;
          GetLaplace1D(Ac, P, lCoarseNodesPerDim);
        } else if (matrixType == "Elasticity1D") {
          GetOStream(Statistics1) << "StructuredRAP: Using Elasticity1D pattern determination routine." << std::endl;
          GetElasticity1D(Ac, P, lCoarseNodesPerDim);
        } else if (matrixType == "Laplace2D") {
          GetOStream(Statistics1) << "StructuredRAP: Using Laplace2D pattern determination routine.\n";
          GetLaplace2D(Ac, P, lCoarseNodesPerDim, procX, procY, interpolationOrder);
        } else if (matrixType == "Elasticity2D") {
          GetOStream(Statistics1) << "StructuredRAP: Using Elasticity2D pattern determination routine.\n";
          GetElasticity2D(Ac, P, lCoarseNodesPerDim, procX, procY);
        } else if (matrixType == "Laplace3D") {
          GetOStream(Statistics1) << "StructuredRAP: Using Laplace3D pattern determination routine.\n";
          GetLaplace3D(Ac, P, lCoarseNodesPerDim, procX, procY, procZ, interpolationOrder);
        } else if (matrixType == "Elasticity3D") {
          GetOStream(Statistics1) << "StructuredRAP: Using Elasticity3D pattern determination routine.\n";
          GetElasticity3D(Ac, P, lCoarseNodesPerDim, procX, procY, procZ);
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(
              true, Exceptions::RuntimeError,
              "StructuredRAPFactory: matrixType \"" << matrixType
                                                    << "\" is not supported for prebuilt Ac graph.");
        }
      }

      // We *always* need global constants for the RAP, but not for the temps
      RAPparams->set("compute global constants: temporaries", RAPparams->get("compute global constants: temporaries", false));
      RAPparams->set("compute global constants", true);

      if (Ac.is_null())
        Ac = MatrixFactory::Build(P->getDomainMap(), Teuchos::as<LocalOrdinal>(0));

      SubFactoryMonitor m2(*this, "MxMxM: P^T x A x P (implicit)", coarseLevel);

      Xpetra::TripleMatrixMultiply<SC, LO, GO, NO>::
          MultiplyRAP(*P, doTranspose, *A, !doTranspose, *P, !doTranspose, *Ac, doFillComplete,
                      doOptimizeStorage, labelstr + std::string("MueLu::P^T*A*P-implicit-") + levelstr.str(),
                      RAPparams);

      GetOStream(Statistics1) << "StructuredRAP: Ac nnz (prebuild coarse graph = "
                              << (prebuildCoarseGraph ? "true" : "false")
                              << "): local = " << Ac->getLocalNumEntries()
                              << ", global = " << Ac->getGlobalNumEntries() << std::endl;

      Teuchos::ArrayView<const double> relativeFloor = pL.get<Teuchos::Array<double>>("rap: relative diagonal floor")();
      if (relativeFloor.size() > 0) {
        Xpetra::MatrixUtils<SC, LO, GO, NO>::RelativeDiagonalBoost(Ac, relativeFloor, GetOStream(Statistics2));
      }

      bool repairZeroDiagonals = pL.get<bool>("RepairMainDiagonal") || pL.get<bool>("rap: fix zero diagonals");
      bool checkAc             = pL.get<bool>("CheckMainDiagonal") || pL.get<bool>("rap: fix zero diagonals");
      if (checkAc || repairZeroDiagonals) {
        using magnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
        magnitudeType threshold;
        if (pL.isType<magnitudeType>("rap: fix zero diagonals threshold"))
          threshold = pL.get<magnitudeType>("rap: fix zero diagonals threshold");
        else
          threshold = Teuchos::as<magnitudeType>(pL.get<double>("rap: fix zero diagonals threshold"));
        Scalar replacement = Teuchos::as<Scalar>(pL.get<double>("rap: fix zero diagonals replacement"));
        Xpetra::MatrixUtils<SC, LO, GO, NO>::CheckRepairMainDiagonal(Ac, repairZeroDiagonals, GetOStream(Warnings1), threshold, replacement);
      }

      if (IsPrint(Statistics2)) {
        RCP<ParameterList> params = rcp(new ParameterList());
        params->set("printLoadBalancingInfo", true);
        params->set("printCommInfo", true);

        GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*Ac, "Ac", params);
      }

      if (!Ac.is_null()) {
        std::ostringstream oss;
        oss << "A_" << coarseLevel.GetLevelID();
        Ac->setObjectLabel(oss.str());
      }
      Set(coarseLevel, "A", Ac);

      RAPparams->set("graph", Ac);
      Set(coarseLevel, "RAP reuse data", RAPparams);
    }
  }

  if (Behavior::debug())
    MatrixUtils::checkLocalRowMapMatchesColMap(*Ac);

  if (transferFacts_.begin() != transferFacts_.end()) {
    SubFactoryMonitor m(*this, "Projections", coarseLevel);

    // call Build of all user-given transfer factories
    for (std::vector<RCP<const FactoryBase>>::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it) {
      RCP<const FactoryBase> fac = *it;
      GetOStream(Runtime0) << "RAPFactory: call transfer factory: " << fac->description() << std::endl;
      fac->CallBuild(coarseLevel);
      // Coordinates transfer is marginally different from all other operations
      // because it is *optional*, and not required. For instance, we may need
      // coordinates only on level 4 if we start repartitioning from that level,
      // but we don't need them on level 1,2,3. As our current Hierarchy setup
      // assumes propagation of dependencies only through three levels, this
      // means that we need to rely on other methods to propagate optional data.
      //
      // The method currently used is through RAP transfer factories, which are
      // simply factories which are called at the end of RAP with a single goal:
      // transfer some fine data to coarser level. Because these factories are
      // kind of outside of the mainline factories, they behave different. In
      // particular, we call their Build method explicitly, rather than through
      // Get calls. This difference is significant, as the Get call is smart
      // enough to know when to release all factory dependencies, and Build is
      // dumb. This led to the following CoordinatesTransferFactory sequence:
      // 1. Request level 0
      // 2. Request level 1
      // 3. Request level 0
      // 4. Release level 0
      // 5. Release level 1
      //
      // The problem is missing "6. Release level 0". Because it was missing,
      // we had outstanding request on "Coordinates", "Aggregates" and
      // "CoarseMap" on level 0.
      //
      // This was fixed by explicitly calling Release on transfer factories in
      // RAPFactory. I am still unsure how exactly it works, but now we have
      // clear data requests for all levels.
      coarseLevel.Release(*fac);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddTransferFactory(const RCP<const FactoryBase>& factory) {
  // check if it's a TwoLevelFactoryBase based transfer factory
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast,
                             "MueLu::StructuredRAPFactory::AddTransferFactory: Transfer factory is not derived from TwoLevelFactoryBase. "
                             "This is very strange. (Note: you can remove this exception if there's a good reason for)");
  TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_, Exceptions::RuntimeError, "MueLu::StructuredRAPFactory::AddTransferFactory: Factory is being added after we have already declared input");
  transferFacts_.push_back(factory);
}

}  // namespace MueLu

#define MUELU_STRUCTUREDRAPFACTORY_SHORT
#endif  // MUELU_STRUCTUREDRAPFACTORY_DEF_HPP
