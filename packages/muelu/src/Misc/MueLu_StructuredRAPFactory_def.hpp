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
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_TripleMatrixMultiply.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_MapFactory.hpp>

#include "MueLu_StructuredRAPFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Behavior.hpp"
#include "MueLu_RAPFactory_def.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_CommHelpers.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::StructuredRAPFactory()
  : hasDeclaredInput_(false)
  , rapFactoryDelegate_(rcp(new MueLu::RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>())) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~StructuredRAPFactory() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  validParamList->set<std::string>(
      "rap: matrix type", "", "Galeri matrix type used to infer the structured RAP graph.");

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
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ConfigureRAPFactoryDelegate() const {
  const Teuchos::ParameterList& pL = GetParameterList();

  ParameterList rapParams;
  RCP<const ParameterList> validRAPParams = rapFactoryDelegate_->GetValidParameterList();
  for (ParameterList::ConstIterator it = validRAPParams->begin(); it != validRAPParams->end(); ++it) {
    const std::string& paramName = validRAPParams->name(it);
    if (pL.isParameter(paramName))
      rapParams.setEntry(paramName, pL.getEntry(paramName));
    else if (pL.isSublist(paramName))
      rapParams.sublist(paramName) = pL.sublist(paramName);
  }

  rapFactoryDelegate_->SetParameterList(rapParams);
  rapFactoryDelegate_->SetFactory("A", GetFactory("A"));
  rapFactoryDelegate_->SetFactory("P", GetFactory("P"));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  const Teuchos::ParameterList& pL = GetParameterList();
  TEUCHOS_TEST_FOR_EXCEPTION(
      !pL.get<bool>("transpose: use implicit"), Exceptions::RuntimeError,
      "StructuredRAPFactory requires \"transpose: use implicit\" = true because "
      "the prebuilt coarse graph assumes the Galerkin product P^T A P.");

  const bool prebuildCoarseGraph = pL.get<bool>("rap: prebuild coarse graph");

  if (!prebuildCoarseGraph) {
    ConfigureRAPFactoryDelegate();
    coarseLevel.DeclareInput("A", rapFactoryDelegate_.get(), this);
    coarseLevel.DeclareInput("RAP reuse data", rapFactoryDelegate_.get(), this);
    hasDeclaredInput_ = true;
    return;
  }

  Input(fineLevel, "A");
  Input(coarseLevel, "P");

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
typename StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::StructuredGraphSpec
StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetStructuredGraphSpec(
    const std::string& matrixType, const int interpolationOrder) const {
  StructuredGraphSpec graphSpec;
  graphSpec.description = matrixType;

  bool useFullStencil = false;
  if (matrixType == "Laplace1D" || matrixType == "Elasticity1D") {
    graphSpec.numDimensions = 1;
    graphSpec.dofsPerNode   = Teuchos::as<LO>(1);
  } else if (matrixType == "Laplace2D") {
    TEUCHOS_TEST_FOR_EXCEPTION(interpolationOrder < 0 || interpolationOrder > 1, Exceptions::RuntimeError,
                               "StructuredRAPFactory::GetStructuredGraphSpec: interpolation order "
                                   << interpolationOrder << " is not supported for " << matrixType << ".");
    graphSpec.numDimensions = 2;
    graphSpec.dofsPerNode   = Teuchos::as<LO>(1);
    useFullStencil          = interpolationOrder == 1;
  } else if (matrixType == "Elasticity2D") {
    graphSpec.numDimensions = 2;
    graphSpec.dofsPerNode   = Teuchos::as<LO>(2);
    useFullStencil          = true;
  } else if (matrixType == "Laplace3D") {
    TEUCHOS_TEST_FOR_EXCEPTION(interpolationOrder < 0 || interpolationOrder > 1, Exceptions::RuntimeError,
                               "StructuredRAPFactory::GetStructuredGraphSpec: interpolation order "
                                   << interpolationOrder << " is not supported for " << matrixType << ".");
    graphSpec.numDimensions = 3;
    graphSpec.dofsPerNode   = Teuchos::as<LO>(1);
    useFullStencil          = interpolationOrder == 1;
  } else if (matrixType == "Elasticity3D") {
    graphSpec.numDimensions = 3;
    graphSpec.dofsPerNode   = Teuchos::as<LO>(3);
    useFullStencil          = true;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError,
                               "StructuredRAPFactory: matrixType \"" << matrixType
                                                                       << "\" is not supported for prebuilt Ac graph.");
  }

  const int minY = graphSpec.numDimensions > 1 ? -1 : 0;
  const int maxY = graphSpec.numDimensions > 1 ? 1 : 0;
  const int minZ = graphSpec.numDimensions > 2 ? -1 : 0;
  const int maxZ = graphSpec.numDimensions > 2 ? 1 : 0;
  for (int dz = minZ; dz <= maxZ; ++dz) {
    for (int dy = minY; dy <= maxY; ++dy) {
      for (int dx = -1; dx <= 1; ++dx) {
        const int numChangedDimensions = (dx != 0 ? 1 : 0) + (dy != 0 ? 1 : 0) + (dz != 0 ? 1 : 0);
        if (useFullStencil || numChangedDimensions <= 1)
          graphSpec.stencilOffsets.push_back(StencilOffset{dx, dy, dz});
      }
    }
  }

  return graphSpec;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetStructuredGraph(
    RCP<Matrix>& Ac, RCP<Matrix> P,
    const Teuchos::Array<LocalOrdinal>& lCoarseNodesPerDim,
    const StructuredGraphSpec& graphSpec) const {
  TEUCHOS_TEST_FOR_EXCEPTION(graphSpec.numDimensions < 1 || graphSpec.numDimensions > 3, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                         << "): number of dimensions must be between one and three.");
  TEUCHOS_TEST_FOR_EXCEPTION(lCoarseNodesPerDim.size() < graphSpec.numDimensions, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                         << "): insufficient local coarse-grid dimensions.");
  TEUCHOS_TEST_FOR_EXCEPTION(graphSpec.dofsPerNode <= 0, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                         << "): dofsPerNode must be positive.");
  TEUCHOS_TEST_FOR_EXCEPTION(graphSpec.stencilOffsets.empty(), Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                         << "): stencil must contain at least one offset.");

  Teuchos::Array<GO> localNodes(3, Teuchos::as<GO>(1));
  for (int dim = 0; dim < graphSpec.numDimensions; ++dim) {
    localNodes[dim] = Teuchos::as<GO>(lCoarseNodesPerDim[dim]);
    TEUCHOS_TEST_FOR_EXCEPTION(localNodes[dim] <= 0, Exceptions::RuntimeError,
                               "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                           << "): local coarse dimensions must be positive.");
  }
  for (typename std::vector<StencilOffset>::const_iterator offset = graphSpec.stencilOffsets.begin();
       offset != graphSpec.stencilOffsets.end(); ++offset) {
    TEUCHOS_TEST_FOR_EXCEPTION(offset->x < -1 || offset->x > 1 ||
                                   offset->y < -1 || offset->y > 1 ||
                                   offset->z < -1 || offset->z > 1,
                               Exceptions::RuntimeError,
                               "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                           << "): only radius-one stencil offsets are currently supported.");
    TEUCHOS_TEST_FOR_EXCEPTION((graphSpec.numDimensions < 2 && offset->y != 0) ||
                                   (graphSpec.numDimensions < 3 && offset->z != 0),
                               Exceptions::RuntimeError,
                               "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                           << "): stencil contains an offset in an inactive dimension.");
  }

  RCP<ParameterList> paramList = rcp(new ParameterList);
  paramList->set("No Nonlocal Changes", true);
  paramList->set("Optimize Storage", true);
  paramList->set("compute global constants", true);

  auto rowMap                                     = P->getDomainMap();
  const size_t localNumRows                       = rowMap->getLocalNumElements();
  const Teuchos::ArrayView<const GO> localRowGids = rowMap->getLocalElementList();
  const GO dofsPerNodeGO                          = Teuchos::as<GO>(graphSpec.dofsPerNode);
  const size_t rowsPerNode                        = Teuchos::as<size_t>(graphSpec.dofsPerNode);
  const GO localNumNodes                          = localNodes[0] * localNodes[1] * localNodes[2];
  const GO expectedLocalRows                      = localNumNodes * dofsPerNodeGO;
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<GO>(localNumRows) != expectedLocalRows, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                         << "): local coarse dimensions with " << graphSpec.dofsPerNode
                                                                         << " dofs per node do not match local coarse row count "
                                                                         << localNumRows << ".");

  const GO numGlobalRows = Teuchos::as<GO>(rowMap->getGlobalNumElements());
  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalRows % dofsPerNodeGO != 0, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                         << "): global coarse row count is not divisible by dofsPerNode.");
  const GO numGlobalNodes = numGlobalRows / dofsPerNodeGO;
  const GO globalMinGid   = rowMap->getMinAllGlobalIndex();
  const GO globalMaxGid   = rowMap->getMaxAllGlobalIndex();
  TEUCHOS_TEST_FOR_EXCEPTION(globalMaxGid - globalMinGid + 1 != numGlobalRows, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                         << "): coarse row-map GIDs must form a contiguous range.");

  RCP<const Teuchos::Comm<int>> comm = rowMap->getComm();
  const int myRank                   = comm->getRank();
  const int numRanks                 = comm->getSize();
  const GO localMinGid = rowMap->getMinGlobalIndex();
  const GO localMaxGid = rowMap->getMaxGlobalIndex();
  TEUCHOS_TEST_FOR_EXCEPTION((localMinGid - globalMinGid) % dofsPerNodeGO != 0, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                         << "): local row range does not begin at a nodal boundary.");
  const GO firstLocalNode = (localMinGid - globalMinGid) / dofsPerNodeGO;

  bool localRowMapContiguous = true;
  for (int rowLid = 0; rowLid < localRowGids.size(); ++rowLid) {
    if (localRowGids[rowLid] != localMinGid + Teuchos::as<GO>(rowLid)) {
      localRowMapContiguous = false;
      break;
    }
  }

  Teuchos::Array<GO> localRankData(5);
  localRankData[0] = firstLocalNode;
  localRankData[1] = localNodes[0];
  localRankData[2] = localNodes[1];
  localRankData[3] = localNodes[2];
  localRankData[4] = localNumNodes;
  Teuchos::Array<GO> rankData(5 * numRanks);
  Teuchos::gatherAll(*comm, 5, localRankData.getRawPtr(), 5 * numRanks, rankData.getRawPtr());

  Teuchos::Array<int> procGrid(3, 1);
  if (graphSpec.numDimensions == 1) {
    procGrid[0] = numRanks;
  } else if (graphSpec.numDimensions == 2) {
    procGrid[0] = 1;
    while ((procGrid[0] + 1) * (procGrid[0] + 1) <= numRanks)
      ++procGrid[0];
    procGrid[1] = numRanks / procGrid[0];
    while (procGrid[0] * procGrid[1] != numRanks)
      procGrid[1] = numRanks / (--procGrid[0]);
  } else {
    int cubeRootNumRanks = 1;
    while ((cubeRootNumRanks + 1) * (cubeRootNumRanks + 1) * (cubeRootNumRanks + 1) <= numRanks)
      ++cubeRootNumRanks;
    procGrid[0] = cubeRootNumRanks;
    procGrid[1] = cubeRootNumRanks;
    procGrid[2] = cubeRootNumRanks;

    if (procGrid[0] * procGrid[1] * procGrid[2] != numRanks) {
      procGrid[0] = 1;
      procGrid[1] = 1;
      procGrid[2] = 1;

      int procTemp = numRanks;
      int factors[50];
      for (int factor = 0; factor < 50; ++factor)
        factors[factor] = 0;
      for (int factor = 2; factor < 50; ++factor) {
        while (procTemp % factor == 0) {
          ++factors[factor];
          procTemp /= factor;
        }
      }
      procGrid[0] = procTemp;
      for (int factor = 49; factor > 0; --factor) {
        while (factors[factor] != 0) {
          if (procGrid[0] <= procGrid[1] && procGrid[0] <= procGrid[2])
            procGrid[0] *= factor;
          else if (procGrid[1] <= procGrid[0] && procGrid[1] <= procGrid[2])
            procGrid[1] *= factor;
          else
            procGrid[2] *= factor;
          --factors[factor];
        }
      }
    }
  }

  const int procXY = procGrid[0] * procGrid[1];
  TEUCHOS_TEST_FOR_EXCEPTION(procGrid[0] <= 0 || procGrid[1] <= 0 || procGrid[2] <= 0 ||
                                 procXY * procGrid[2] != numRanks,
                             Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                         << "): invalid inferred processor grid "
                                                                         << procGrid[0] << "x" << procGrid[1] << "x" << procGrid[2]
                                                                         << " for " << numRanks << " ranks.");

  const int myProcX = myRank % procGrid[0];
  const int myProcY = (myRank % procXY) / procGrid[0];
  const int myProcZ = myRank / procXY;

  Teuchos::Array<GO> globalNodes(3, Teuchos::as<GO>(0));
  for (int px = 0; px < procGrid[0]; ++px)
    globalNodes[0] += rankData[5 * px + 1];
  for (int py = 0; py < procGrid[1]; ++py)
    globalNodes[1] += rankData[5 * (py * procGrid[0]) + 2];
  for (int pz = 0; pz < procGrid[2]; ++pz)
    globalNodes[2] += rankData[5 * (pz * procXY) + 3];
  TEUCHOS_TEST_FOR_EXCEPTION(globalNodes[0] * globalNodes[1] * globalNodes[2] != numGlobalNodes,
                             Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                         << "): processor-grid coarse dimensions do not match coarse node count.");

  const auto getDofGid = [&](const GO nodeOrdinal, const LO dof) -> GO {
    return globalMinGid + nodeOrdinal * dofsPerNodeGO + Teuchos::as<GO>(dof);
  };

  const auto getLocalNodeCoordinates = [&](const GO localNode, GO& x, GO& y, GO& z) {
    x = localNode % localNodes[0];
    y = (localNode / localNodes[0]) % localNodes[1];
    z = localNode / (localNodes[0] * localNodes[1]);
  };

  const auto resolveNeighbor = [&](const GO x, const GO y, const GO z, const StencilOffset& offset,
                                   int& neighborRank, GO& neighborNode) -> bool {
    int neighborProcX = myProcX;
    int neighborProcY = myProcY;
    int neighborProcZ = myProcZ;
    GO neighborX      = x + Teuchos::as<GO>(offset.x);
    GO neighborY      = y + Teuchos::as<GO>(offset.y);
    GO neighborZ      = z + Teuchos::as<GO>(offset.z);

    if (neighborX < 0)
      --neighborProcX;
    else if (neighborX >= localNodes[0])
      ++neighborProcX;
    if (neighborY < 0)
      --neighborProcY;
    else if (neighborY >= localNodes[1])
      ++neighborProcY;
    if (neighborZ < 0)
      --neighborProcZ;
    else if (neighborZ >= localNodes[2])
      ++neighborProcZ;

    if (neighborProcX < 0 || neighborProcX >= procGrid[0] ||
        neighborProcY < 0 || neighborProcY >= procGrid[1] ||
        neighborProcZ < 0 || neighborProcZ >= procGrid[2])
      return false;

    neighborRank             = neighborProcZ * procXY + neighborProcY * procGrid[0] + neighborProcX;
    const GO neighborLocalNx = rankData[5 * neighborRank + 1];
    const GO neighborLocalNy = rankData[5 * neighborRank + 2];
    const GO neighborLocalNz = rankData[5 * neighborRank + 3];
    if (neighborX < 0)
      neighborX = neighborLocalNx - 1;
    else if (neighborX >= localNodes[0])
      neighborX = 0;
    if (neighborY < 0)
      neighborY = neighborLocalNy - 1;
    else if (neighborY >= localNodes[1])
      neighborY = 0;
    if (neighborZ < 0)
      neighborZ = neighborLocalNz - 1;
    else if (neighborZ >= localNodes[2])
      neighborZ = 0;

    neighborNode = rankData[5 * neighborRank] +
                   neighborZ * neighborLocalNx * neighborLocalNy +
                   neighborY * neighborLocalNx + neighborX;
    return true;
  };

  const size_t maxStencilSize = graphSpec.stencilOffsets.size();
  Teuchos::Array<int> neighborRanks(Teuchos::as<int>(maxStencilSize));
  Teuchos::Array<GO> neighborNodes(Teuchos::as<int>(maxStencilSize));
  const auto getNeighbors = [&](const GO x, const GO y, const GO z) -> size_t {
    size_t numNeighbors = 0;
    for (typename std::vector<StencilOffset>::const_iterator offset = graphSpec.stencilOffsets.begin();
         offset != graphSpec.stencilOffsets.end(); ++offset) {
      int neighborRank = myRank;
      GO neighborNode  = 0;
      if (resolveNeighbor(x, y, z, *offset, neighborRank, neighborNode)) {
        neighborRanks[numNeighbors] = neighborRank;
        neighborNodes[numNeighbors] = neighborNode;
        ++numNeighbors;
      }
    }
    return numNeighbors;
  };

  const bool groupedContiguousRows =
      localRowMapContiguous && localNumRows % rowsPerNode == 0;
  const auto isInteriorNode = [&](const GO x, const GO y, const GO z) -> bool {
    return x > 0 && x + 1 < localNodes[0] &&
           (graphSpec.numDimensions < 2 || (y > 0 && y + 1 < localNodes[1])) &&
           (graphSpec.numDimensions < 3 || (z > 0 && z + 1 < localNodes[2]));
  };

  Teuchos::Array<GO> interiorNodeOffsets(Teuchos::as<int>(maxStencilSize));
  for (size_t stencil = 0; stencil < maxStencilSize; ++stencil) {
    const StencilOffset& offset = graphSpec.stencilOffsets[stencil];
    interiorNodeOffsets[stencil] =
        Teuchos::as<GO>(offset.x) +
        localNodes[0] * (Teuchos::as<GO>(offset.y) +
                         localNodes[1] * Teuchos::as<GO>(offset.z));
  }

  ArrayRCP<size_t> rowptr(localNumRows + 1);
  rowptr[0]       = 0;
  size_t localNnz = 0;
  std::vector<GO> remoteColGids;

  if (groupedContiguousRows) {
    // Every DOF row at a node has the same columns, so resolve the stencil once per node.
    for (GO localNode = 0; localNode < localNumNodes; ++localNode) {
      GO x = 0, y = 0, z = 0;
      getLocalNodeCoordinates(localNode, x, y, z);
      const bool interior       = isInteriorNode(x, y, z);
      const size_t numNeighbors = interior ? maxStencilSize : getNeighbors(x, y, z);
      const size_t rowNnz       = numNeighbors * rowsPerNode;
      const size_t firstRow     = Teuchos::as<size_t>(localNode) * rowsPerNode;
      for (size_t rowDof = 0; rowDof < rowsPerNode; ++rowDof) {
        localNnz += rowNnz;
        rowptr[firstRow + rowDof + 1] = localNnz;
      }

      // Interior nodes cannot reference remote columns.
      if (!interior) {
        for (size_t neighbor = 0; neighbor < numNeighbors; ++neighbor) {
          if (neighborRanks[neighbor] == myRank)
            continue;
          for (LO colDof = 0; colDof < graphSpec.dofsPerNode; ++colDof)
            remoteColGids.push_back(getDofGid(neighborNodes[neighbor], colDof));
        }
      }
    }
  } else {
    // Preserve support for maps whose local rows are not grouped by node.
    for (size_t rowLid = 0; rowLid < localNumRows; ++rowLid) {
      const GO rowGid      = rowMap->getGlobalElement(Teuchos::as<LO>(rowLid));
      const GO nodeOrdinal = (rowGid - globalMinGid) / dofsPerNodeGO;
      const GO localNode   = nodeOrdinal - firstLocalNode;
      TEUCHOS_TEST_FOR_EXCEPTION(localNode < 0 || localNode >= localNumNodes, Exceptions::RuntimeError,
                                 "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                             << "): row GID does not belong to the rank-local structured block.");
      GO x = 0, y = 0, z = 0;
      getLocalNodeCoordinates(localNode, x, y, z);
      const size_t numNeighbors = getNeighbors(x, y, z);
      localNnz += numNeighbors * rowsPerNode;
      rowptr[rowLid + 1] = localNnz;
    }

    for (GO localNode = 0; localNode < localNumNodes; ++localNode) {
      GO x = 0, y = 0, z = 0;
      getLocalNodeCoordinates(localNode, x, y, z);
      const size_t numNeighbors = getNeighbors(x, y, z);
      for (size_t neighbor = 0; neighbor < numNeighbors; ++neighbor) {
        if (neighborRanks[neighbor] == myRank)
          continue;
        for (LO colDof = 0; colDof < graphSpec.dofsPerNode; ++colDof)
          remoteColGids.push_back(getDofGid(neighborNodes[neighbor], colDof));
      }
    }
  }
  std::sort(remoteColGids.begin(), remoteColGids.end());
  remoteColGids.erase(std::unique(remoteColGids.begin(), remoteColGids.end()), remoteColGids.end());

  Array<GO> colMapGids;
  colMapGids.reserve(localRowGids.size() + Teuchos::as<int>(remoteColGids.size()));
  for (int rowLid = 0; rowLid < localRowGids.size(); ++rowLid)
    colMapGids.push_back(localRowGids[rowLid]);
  for (typename std::vector<GO>::const_iterator gid = remoteColGids.begin(); gid != remoteColGids.end(); ++gid)
    colMapGids.push_back(*gid);

  RCP<const Map> colMap = MapFactory::Build(rowMap->lib(),
                                            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                            colMapGids(), rowMap->getIndexBase(), comm);

  const size_t maxColumnsPerRow = maxStencilSize * rowsPerNode;
  Array<LO> rowColLids(Teuchos::as<int>(maxColumnsPerRow));
  ArrayRCP<LO> colind(localNnz);

  if (groupedContiguousRows) {
    for (GO localNode = 0; localNode < localNumNodes; ++localNode) {
      GO x = 0, y = 0, z = 0;
      getLocalNodeCoordinates(localNode, x, y, z);
      const bool interior = isInteriorNode(x, y, z);
      size_t rowNnz       = 0;

      if (interior) {
        // Stencil offsets are in x-fastest order, so these local column LIDs are already sorted.
        for (size_t stencil = 0; stencil < maxStencilSize; ++stencil) {
          const GO colLocalNode = localNode + interiorNodeOffsets[stencil];
          const LO colBase      = Teuchos::as<LO>(colLocalNode * dofsPerNodeGO);
          for (LO colDof = 0; colDof < graphSpec.dofsPerNode; ++colDof)
            rowColLids[rowNnz++] = colBase + colDof;
        }
      } else {
        const size_t numNeighbors = getNeighbors(x, y, z);
        for (size_t neighbor = 0; neighbor < numNeighbors; ++neighbor) {
          for (LO colDof = 0; colDof < graphSpec.dofsPerNode; ++colDof) {
            const GO colGid = getDofGid(neighborNodes[neighbor], colDof);
            const bool localFastPath = colGid >= localMinGid && colGid <= localMaxGid;
            const LO colLid = localFastPath
                                  ? Teuchos::as<LO>(colGid - localMinGid)
                                  : colMap->getLocalElement(colGid);
            TEUCHOS_TEST_FOR_EXCEPTION(colLid == Teuchos::OrdinalTraits<LO>::invalid(), Exceptions::RuntimeError,
                                       "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                                   << "): column GID " << colGid
                                                                                   << " was not found in the coarse column map.");
            rowColLids[rowNnz++] = colLid;
          }
        }
        std::sort(rowColLids.getRawPtr(), rowColLids.getRawPtr() + rowNnz);
      }

      const size_t firstRow = Teuchos::as<size_t>(localNode) * rowsPerNode;
      for (size_t rowDof = 0; rowDof < rowsPerNode; ++rowDof) {
        const size_t rowStart = rowptr[firstRow + rowDof];
        TEUCHOS_ASSERT(rowptr[firstRow + rowDof + 1] - rowStart == rowNnz);
        for (size_t column = 0; column < rowNnz; ++column)
          colind[rowStart + column] = rowColLids[column];
      }
    }
  } else {
    size_t entry = 0;
    for (size_t rowLid = 0; rowLid < localNumRows; ++rowLid) {
      const GO rowGid      = rowMap->getGlobalElement(Teuchos::as<LO>(rowLid));
      const GO nodeOrdinal = (rowGid - globalMinGid) / dofsPerNodeGO;
      const GO localNode   = nodeOrdinal - firstLocalNode;
      GO x = 0, y = 0, z = 0;
      getLocalNodeCoordinates(localNode, x, y, z);
      const size_t numNeighbors = getNeighbors(x, y, z);

      size_t rowNnz = 0;
      for (size_t neighbor = 0; neighbor < numNeighbors; ++neighbor) {
        for (LO colDof = 0; colDof < graphSpec.dofsPerNode; ++colDof) {
          const GO colGid = getDofGid(neighborNodes[neighbor], colDof);
          const LO colLid = colMap->getLocalElement(colGid);
          TEUCHOS_TEST_FOR_EXCEPTION(colLid == Teuchos::OrdinalTraits<LO>::invalid(), Exceptions::RuntimeError,
                                     "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                                 << "): column GID " << colGid
                                                                                 << " was not found in the coarse column map.");
          rowColLids[rowNnz++] = colLid;
        }
      }
      std::sort(rowColLids.getRawPtr(), rowColLids.getRawPtr() + rowNnz);
      for (size_t column = 0; column < rowNnz; ++column)
        colind[entry++] = rowColLids[column];
    }
    TEUCHOS_ASSERT(entry == localNnz);
  }

  RCP<CrsGraph> graph = CrsGraphFactory::Build(rowMap, colMap, rowptr, colind, paramList);
  if (!graph->isFillComplete())
    graph->fillComplete(rowMap, rowMap, paramList);
  Ac = MatrixFactory::Build(graph, paramList);
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

  TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_ == false, Exceptions::RuntimeError,
                             "MueLu::RAPFactory::Build(): CallDeclareInput has not been called before Build!");

  if (!prebuildCoarseGraph) {
    if (coarseLevel.IsAvailable("RAP reuse data", this)) {
      RCP<ParameterList> RAPparams = coarseLevel.Get<RCP<ParameterList>>("RAP reuse data", this);
      coarseLevel.Set("RAP reuse data", RAPparams, rapFactoryDelegate_.get());
    }

    // If prebuildCoarseGraph is false, we delegate the whole RAP computation to RAPFactory.
    rapFactoryDelegate_->Build(fineLevel, coarseLevel);

    Ac = coarseLevel.Get<RCP<Matrix>>("A", rapFactoryDelegate_.get());
    Set(coarseLevel, "A", Ac);

    if (coarseLevel.IsAvailable("RAP reuse data", rapFactoryDelegate_.get())) {
      RCP<ParameterList> RAPparams = coarseLevel.Get<RCP<ParameterList>>("RAP reuse data", rapFactoryDelegate_.get());
      Set(coarseLevel, "RAP reuse data", RAPparams);
    }

    return;
  }

  {
    FactoryMonitor m(*this, "Computing Ac", coarseLevel);

    std::ostringstream levelstr;
    levelstr << coarseLevel.GetLevelID();
    std::string labelstr = FormattingHelper::getColonLabel(coarseLevel.getObjectLabel());

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

    RCP<Matrix> A = Get<RCP<Matrix>>(fineLevel, "A");
    RCP<Matrix> P = Get<RCP<Matrix>>(coarseLevel, "P");
    // We don't have a valid P (e.g., # global aggregates = 0) so we bail.
    // This level will ultimately be removed in MueLu_Hierarchy_defs.h via a resize()
    if (P.is_null()) {
      Ac = Teuchos::null;
      Set(coarseLevel, "A", Ac);
      return;
    }

    {
      RCP<ParameterList> RAPparams = rcp(new ParameterList);
      if (pL.isSublist("matrixmatrix: kernel params"))
        RAPparams->sublist("matrixmatrix: kernel params") = pL.sublist("matrixmatrix: kernel params");

      if (coarseLevel.IsAvailable("RAP reuse data", this)) {
        GetOStream(static_cast<MsgType>(Runtime0 | Test)) << "Reusing previous RAP data" << std::endl;

        RAPparams = coarseLevel.Get<RCP<ParameterList>>("RAP reuse data", this);

        TEUCHOS_TEST_FOR_EXCEPTION(!RAPparams->isParameter("graph"), Exceptions::RuntimeError,
                                   "StructuredRAPFactory::Build(): \"RAP reuse data\" does not contain the expected graph.");
        Ac = RAPparams->get<RCP<Matrix>>("graph");
        TEUCHOS_TEST_FOR_EXCEPTION(Ac.is_null(), Exceptions::RuntimeError,
                                   "StructuredRAPFactory::Build(): \"RAP reuse data\" graph is null.");

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
        const StructuredGraphSpec graphSpec = GetStructuredGraphSpec(matrixType, interpolationOrder);
        GetOStream(Statistics1) << "StructuredRAP: Using " << graphSpec.description
                                << " stencil with " << graphSpec.stencilOffsets.size()
                                << " nodal entries." << std::endl;
        GetStructuredGraph(Ac, P, lCoarseNodesPerDim, graphSpec);
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
  rapFactoryDelegate_->AddTransferFactory(factory);
}

}  // namespace MueLu

#define MUELU_STRUCTUREDRAPFACTORY_SHORT
#endif  // MUELU_STRUCTUREDRAPFACTORY_DEF_HPP
