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

#include <sstream>
#include <vector>

#include <Kokkos_Core.hpp>

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
#include "Teuchos_TimeMonitor.hpp"

namespace MueLu {

namespace StructuredRAPFactoryDetails {

template <class GO, class LocalNodes>
KOKKOS_INLINE_FUNCTION void getLocalNodeIndices(const GO localNode,
                                                const LocalNodes& localNodes,
                                                GO& x, GO& y, GO& z) {
  x = localNode % localNodes[0];
  y = (localNode / localNodes[0]) % localNodes[1];
  z = localNode / (localNodes[0] * localNodes[1]);
}

template <class GO, class Offset, class LocalNodes, class ProcGrid, class RankData>
KOKKOS_INLINE_FUNCTION bool resolveNeighbor(const GO x, const GO y, const GO z,
                                            const Offset& offset,
                                            const LocalNodes& localNodes,
                                            const ProcGrid& procGrid,
                                            const RankData& rankData,
                                            const int myProcX,
                                            const int myProcY,
                                            const int myProcZ,
                                            int& neighborRank,
                                            GO& neighborNode) {
  int neighborProcX = myProcX;
  int neighborProcY = myProcY;
  int neighborProcZ = myProcZ;
  GO neighborX      = x + static_cast<GO>(offset.x);
  GO neighborY      = y + static_cast<GO>(offset.y);
  GO neighborZ      = z + static_cast<GO>(offset.z);

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

  neighborRank             = neighborProcZ * procGrid[0] * procGrid[1] + neighborProcY * procGrid[0] + neighborProcX;
  const GO neighborLocalNx = rankData(4 * neighborRank + 1);
  const GO neighborLocalNy = rankData(4 * neighborRank + 2);
  const GO neighborLocalNz = rankData(4 * neighborRank + 3);
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

  neighborNode = rankData(4 * neighborRank) +
                 neighborZ * neighborLocalNx * neighborLocalNy +
                 neighborY * neighborLocalNx + neighborX;
  return true;
}

}  // namespace StructuredRAPFactoryDetails

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
      "matrixType", null, "Matrix type used to infer the structured RAP graph.");

  validParamList->set<bool>("CheckMainDiagonal", false, "Check main diagonal for zeros");
  validParamList->set<bool>("RepairMainDiagonal", false, "Repair zeros on main diagonal");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

// Configure RAPFactory to delegate to if coarse graph prebuilding is disabled
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
  const bool useRAPDelegate      = !prebuildCoarseGraph;

  if (useRAPDelegate) {
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

// Describe the expected coarse-matrix sparsity pattern based on the matrix type and interpolation order
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

// Prebuild sparsity structure of coarse matrix
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetStructuredGraph(
    RCP<Matrix>& Ac, RCP<Matrix> P,
    const Teuchos::Array<LocalOrdinal>& lCoarseNodesPerDim,
    const StructuredGraphSpec& graphSpec) const {
  using local_graph_type = typename CrsGraph::local_graph_type;
  using row_map_type     = typename local_graph_type::row_map_type::non_const_type;
  using entries_type     = typename local_graph_type::entries_type::non_const_type;
  using device_type      = typename Node::device_type;
  using execution_space  = typename device_type::execution_space;
  using range_policy     = Kokkos::RangePolicy<execution_space, Kokkos::IndexType<size_t>>;

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
                                                                         << "): coarse-grid stencil is empty; it must define at least one relative node position (dx, dy, dz).");
  constexpr size_t maxSupportedStencilSize = 27;
  constexpr size_t maxSupportedDofsPerNode = 3;
  TEUCHOS_TEST_FOR_EXCEPTION(graphSpec.stencilOffsets.size() > maxSupportedStencilSize, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                         << "): a radius-one stencil may contain at most 27 entries.");
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(graphSpec.dofsPerNode) > maxSupportedDofsPerNode, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                         << "): at most three DOFs per node are supported.");

  Kokkos::Array<GO, 3> localNodes{};
  for (int dim = 0; dim < 3; ++dim)
    localNodes[dim] = Teuchos::as<GO>(1);
  for (int dim = 0; dim < graphSpec.numDimensions; ++dim) {
    localNodes[dim] = Teuchos::as<GO>(lCoarseNodesPerDim[dim]);
    TEUCHOS_TEST_FOR_EXCEPTION(localNodes[dim] <= 0, Exceptions::RuntimeError,
                               "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                           << "): local coarse dimensions must be positive.");
  }

  Kokkos::Array<StencilOffset, maxSupportedStencilSize> stencilOffsets{};
  Kokkos::Array<GO, maxSupportedStencilSize> localStencilNodeOffsets{};
  for (size_t stencil = 0; stencil < graphSpec.stencilOffsets.size(); ++stencil) {
    const StencilOffset& offset = graphSpec.stencilOffsets[stencil];
    TEUCHOS_TEST_FOR_EXCEPTION(offset.x < -1 || offset.x > 1 ||
                                   offset.y < -1 || offset.y > 1 ||
                                   offset.z < -1 || offset.z > 1,
                               Exceptions::RuntimeError,
                               "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                           << "): only radius-one stencil offsets are currently supported.");
    TEUCHOS_TEST_FOR_EXCEPTION((graphSpec.numDimensions < 2 && offset.y != 0) ||
                                   (graphSpec.numDimensions < 3 && offset.z != 0),
                               Exceptions::RuntimeError,
                               "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                           << "): stencil contains an offset in an inactive dimension.");
    stencilOffsets[stencil] = offset;
    localStencilNodeOffsets[stencil] =
        static_cast<GO>(offset.x) + static_cast<GO>(offset.y) * localNodes[0] +
        static_cast<GO>(offset.z) * localNodes[0] * localNodes[1];
  }

  RCP<ParameterList> paramList = rcp(new ParameterList);
  paramList->set("No Nonlocal Changes", true);
  paramList->set("Optimize Storage", true);
  paramList->set("compute global constants", true);

  // Columns of P represent coarse-grid nodes
  auto rowMap                = P->getDomainMap();
  const size_t localNumRows  = rowMap->getLocalNumElements();
  const LO dofsPerNode       = graphSpec.dofsPerNode;
  const GO dofsPerNodeGO     = Teuchos::as<GO>(dofsPerNode);
  const size_t rowsPerNode   = Teuchos::as<size_t>(dofsPerNode);
  const GO localNumNodesGO   = localNodes[0] * localNodes[1] * localNodes[2];
  const size_t localNumNodes = Teuchos::as<size_t>(localNumNodesGO);
  const GO expectedLocalRows = localNumNodesGO * dofsPerNodeGO;
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<GO>(localNumRows) != expectedLocalRows, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                         << "): local coarse dimensions with " << graphSpec.dofsPerNode
                                                                         << " dofs per node do not match local coarse row count "
                                                                         << localNumRows << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(rowMap->lib() != Xpetra::UseTpetra, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph requires the Tpetra backend.");

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
  const GO localMinGid               = rowMap->getMinGlobalIndex();
  TEUCHOS_TEST_FOR_EXCEPTION((localMinGid - globalMinGid) % dofsPerNodeGO != 0, Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                         << "): local row range does not begin at a nodal boundary.");
  const GO firstLocalNode = (localMinGid - globalMinGid) / dofsPerNodeGO;

  Teuchos::Array<GO> localRankData(4);
  localRankData[0] = firstLocalNode;
  localRankData[1] = localNodes[0];
  localRankData[2] = localNodes[1];
  localRankData[3] = localNodes[2];
  Teuchos::Array<GO> rankData(4 * numRanks);
  {
    // For RankData on every rank
    Teuchos::gatherAll(*comm, 4, localRankData.getRawPtr(), 4 * numRanks, rankData.getRawPtr());
  }

  // Get mx, my and mz (number of ranks per dimension) like they are defined in Galeri (compare Galeri_XpetraMaps_def.hpp)
  Kokkos::Array<int, 3> procGrid{};
  for (int dim = 0; dim < 3; ++dim)
    procGrid[dim] = 1;
  if (graphSpec.numDimensions == 1) {
    procGrid[0] = numRanks;
  } else if (graphSpec.numDimensions == 2) {
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

  Kokkos::Array<GO, 3> globalNodes{};
  for (int dim = 0; dim < 3; ++dim)
    globalNodes[dim] = Teuchos::as<GO>(0);
  for (int px = 0; px < procGrid[0]; ++px)
    globalNodes[0] += rankData[4 * px + 1];
  for (int py = 0; py < procGrid[1]; ++py)
    globalNodes[1] += rankData[4 * (py * procGrid[0]) + 2];
  for (int pz = 0; pz < procGrid[2]; ++pz)
    globalNodes[2] += rankData[4 * (pz * procXY) + 3];
  TEUCHOS_TEST_FOR_EXCEPTION(globalNodes[0] * globalNodes[1] * globalNodes[2] != numGlobalNodes,
                             Exceptions::RuntimeError,
                             "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                         << "): processor-grid coarse dimensions do not match coarse node count.");

  execution_space executionSpace;
  Kokkos::View<GO*, device_type> rankDataDevice(
      Kokkos::ViewAllocateWithoutInitializing("StructuredRAP: rank metadata"), rankData.size());
  auto rankDataHost = Kokkos::create_mirror_view(rankDataDevice);
  for (int entry = 0; entry < rankData.size(); ++entry)
    rankDataHost(entry) = rankData[entry];
  Kokkos::deep_copy(executionSpace, rankDataDevice, rankDataHost);

  const size_t stencilSize = graphSpec.stencilOffsets.size();
  const int numDimensions  = graphSpec.numDimensions;
  const auto rowLocalMap   = rowMap->getLocalMap();

  const bool debug = Behavior::debug();
  if (debug) {
    size_t invalidLocalRows = 0;
    Kokkos::parallel_reduce(
        "StructuredRAP: validate local rows", range_policy(executionSpace, 0, localNumRows),
        KOKKOS_LAMBDA(const size_t rowLid, size_t& invalid) {
          const GO rowGid         = rowLocalMap.getGlobalElement(static_cast<LO>(rowLid));
          const GO expectedRowGid = localMinGid + static_cast<GO>(rowLid);
          if (rowGid != expectedRowGid)
            ++invalid;
        },
        invalidLocalRows);
    TEUCHOS_TEST_FOR_EXCEPTION(invalidLocalRows != 0, Exceptions::RuntimeError,
                               "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                           << "): coarse row-map GIDs are not locally contiguous and ordered.");
  }

  // A radius-one stencil reaches disjoint faces, edges, and corners of at most
  // 26 neighboring ranks, so generate each remote column exactly once.
  constexpr size_t maxNeighborRegions = 26;
  Kokkos::Array<int, maxNeighborRegions> remoteRegionRanks{};
  Kokkos::Array<StencilOffset, maxNeighborRegions> remoteRegionDirections{};
  Kokkos::Array<GO, maxNeighborRegions> remoteRegionFirstNodes{};
  Kokkos::Array<GO, maxNeighborRegions> remoteRegionNx{};
  Kokkos::Array<GO, maxNeighborRegions> remoteRegionNy{};
  Kokkos::Array<GO, maxNeighborRegions> remoteRegionNz{};
  Kokkos::Array<size_t, maxNeighborRegions + 1> remoteRegionOffsets{};
  size_t numRemoteRegions = 0;

  const int minProcY = graphSpec.numDimensions > 1 ? -1 : 0;
  const int maxProcY = graphSpec.numDimensions > 1 ? 1 : 0;
  const int minProcZ = graphSpec.numDimensions > 2 ? -1 : 0;
  const int maxProcZ = graphSpec.numDimensions > 2 ? 1 : 0;
  for (int procDz = minProcZ; procDz <= maxProcZ; ++procDz) {
    for (int procDy = minProcY; procDy <= maxProcY; ++procDy) {
      for (int procDx = -1; procDx <= 1; ++procDx) {
        if (procDx == 0 && procDy == 0 && procDz == 0)
          continue;

        bool stencilReachesRegion = false;
        for (size_t stencil = 0; stencil < stencilSize; ++stencil) {
          const StencilOffset& offset = stencilOffsets[stencil];
          const bool reachesX         = procDx == 0 || offset.x == procDx;
          const bool reachesY         = procDy == 0 || offset.y == procDy;
          const bool reachesZ         = procDz == 0 || offset.z == procDz;
          if (reachesX && reachesY && reachesZ) {
            stencilReachesRegion = true;
            break;
          }
        }
        if (!stencilReachesRegion)
          continue;

        const int neighborProcX = myProcX + procDx;
        const int neighborProcY = myProcY + procDy;
        const int neighborProcZ = myProcZ + procDz;
        if (neighborProcX < 0 || neighborProcX >= procGrid[0] ||
            neighborProcY < 0 || neighborProcY >= procGrid[1] ||
            neighborProcZ < 0 || neighborProcZ >= procGrid[2])
          continue;

        const int neighborRank   = neighborProcZ * procXY + neighborProcY * procGrid[0] + neighborProcX;
        const GO neighborNx      = rankData[4 * neighborRank + 1];
        const GO neighborNy      = rankData[4 * neighborRank + 2];
        const GO neighborNz      = rankData[4 * neighborRank + 3];
        const size_t regionNodes = Teuchos::as<size_t>((procDx == 0 ? neighborNx : 1) *
                                                       (procDy == 0 ? neighborNy : 1) *
                                                       (procDz == 0 ? neighborNz : 1));

        remoteRegionRanks[numRemoteRegions]      = neighborRank;
        remoteRegionDirections[numRemoteRegions] = StencilOffset{procDx, procDy, procDz};
        remoteRegionFirstNodes[numRemoteRegions] = rankData[4 * neighborRank];
        remoteRegionNx[numRemoteRegions]         = neighborNx;
        remoteRegionNy[numRemoteRegions]         = neighborNy;
        remoteRegionNz[numRemoteRegions]         = neighborNz;
        remoteRegionOffsets[numRemoteRegions + 1] =
            remoteRegionOffsets[numRemoteRegions] + regionNodes * rowsPerNode;
        ++numRemoteRegions;
      }
    }
  }
  const size_t numRemoteColumns = remoteRegionOffsets[numRemoteRegions];

  Kokkos::View<GO*, device_type> colMapGids(
      Kokkos::ViewAllocateWithoutInitializing("StructuredRAP: column map GIDs"),
      localNumRows + numRemoteColumns);
  const auto localRowGidsDevice = rowMap->getMyGlobalIndicesDevice();
  Kokkos::parallel_for(
      "StructuredRAP: fill column map", range_policy(executionSpace, 0, localNumRows + numRemoteColumns),
      KOKKOS_LAMBDA(const size_t entry) {
        if (entry < localNumRows) {
          colMapGids(entry) = localRowGidsDevice(entry);
          return;
        }

        const size_t remoteEntry = entry - localNumRows;
        size_t region            = 0;
        while (remoteEntry >= remoteRegionOffsets[region + 1])
          ++region;

        const StencilOffset direction = remoteRegionDirections[region];
        const GO neighborNx           = remoteRegionNx[region];
        const GO neighborNy           = remoteRegionNy[region];
        const GO neighborNz           = remoteRegionNz[region];
        size_t regionEntry            = remoteEntry - remoteRegionOffsets[region];
        const GO dof                  = static_cast<GO>(regionEntry % rowsPerNode);
        size_t regionNode             = regionEntry / rowsPerNode;

        GO neighborX = direction.x < 0 ? neighborNx - 1 : 0;
        GO neighborY = direction.y < 0 ? neighborNy - 1 : 0;
        GO neighborZ = direction.z < 0 ? neighborNz - 1 : 0;
        if (direction.x == 0) {
          neighborX = static_cast<GO>(regionNode % static_cast<size_t>(neighborNx));
          regionNode /= static_cast<size_t>(neighborNx);
        }
        if (direction.y == 0) {
          neighborY = static_cast<GO>(regionNode % static_cast<size_t>(neighborNy));
          regionNode /= static_cast<size_t>(neighborNy);
        }
        if (direction.z == 0)
          neighborZ = static_cast<GO>(regionNode);

        const GO neighborNode = remoteRegionFirstNodes[region] +
                                neighborZ * neighborNx * neighborNy + neighborY * neighborNx + neighborX;
        colMapGids(entry) = globalMinGid + neighborNode * dofsPerNodeGO + dof;
      });
  executionSpace.fence("StructuredRAP: column map GIDs ready");

  typename Map::global_indices_array_device_type constColMapGids = colMapGids;
  // colMap is filled here
  RCP<const Map> colMap;
  {
    colMap = MapFactory::Build(
        rowMap->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
        constColMapGids, rowMap->getIndexBase(), comm);
  }

  // Loop to fill rowptr
  row_map_type rowptr;
  size_t localNnz = 0;
  {
    rowptr = row_map_type(
        Kokkos::ViewAllocateWithoutInitializing("StructuredRAP: graph row pointers"), localNumRows + 1);
    Kokkos::parallel_scan(
        "StructuredRAP: compute graph row pointers", range_policy(executionSpace, 0, localNumNodes + 1),
        KOKKOS_LAMBDA(const size_t localNodeIndex, size_t& update, const bool final) {
          if (localNodeIndex == localNumNodes) {
            if (final)
              rowptr(localNumRows) = update;
            return;
          }

          const GO localNode = static_cast<GO>(localNodeIndex);
          GO x = 0, y = 0, z = 0;
          StructuredRAPFactoryDetails::getLocalNodeIndices(localNode, localNodes, x, y, z);
          const bool isInterior =
              x > 0 && x + 1 < localNodes[0] &&
              (numDimensions < 2 || (y > 0 && y + 1 < localNodes[1])) &&
              (numDimensions < 3 || (z > 0 && z + 1 < localNodes[2]));
          size_t rowLength = isInterior ? stencilSize * rowsPerNode : 0;
          if (!isInterior) {
            for (size_t stencil = 0; stencil < stencilSize; ++stencil) {
              int neighborRank = myRank;
              GO neighborNode  = 0;
              if (StructuredRAPFactoryDetails::resolveNeighbor(
                      x, y, z, stencilOffsets[stencil], localNodes, procGrid, rankDataDevice,
                      myProcX, myProcY, myProcZ, neighborRank, neighborNode))
                rowLength += rowsPerNode;
            }
          }

          if (final) {
            const size_t firstRow = localNodeIndex * rowsPerNode;
            for (size_t rowDof = 0; rowDof < rowsPerNode; ++rowDof)
              rowptr(firstRow + rowDof) = update + rowDof * rowLength;
          }
          update += rowLength * rowsPerNode;
        },
        localNnz);
  }

  entries_type colind;
  {
    colind = entries_type(
        Kokkos::ViewAllocateWithoutInitializing("StructuredRAP: graph column indices"), localNnz);
  }

  const LO invalidLocalOrdinal = Teuchos::OrdinalTraits<LO>::invalid();
  Kokkos::View<int, device_type> invalidColumn;
  if (debug) {
    invalidColumn = Kokkos::View<int, device_type>("StructuredRAP: invalid column");
    Kokkos::deep_copy(executionSpace, invalidColumn, 0);
  }

  // Loop to fill colind
  {
    Kokkos::parallel_for(
        "StructuredRAP: fill graph column indices", range_policy(executionSpace, 0, localNumNodes),
        KOKKOS_LAMBDA(const size_t localNodeIndex) {
          const GO localNode = static_cast<GO>(localNodeIndex);
          GO x = 0, y = 0, z = 0;
          StructuredRAPFactoryDetails::getLocalNodeIndices(localNode, localNodes, x, y, z);
          size_t columnOffset = 0;
          const bool isInterior =
              x > 0 && x + 1 < localNodes[0] &&
              (numDimensions < 2 || (y > 0 && y + 1 < localNodes[1])) &&
              (numDimensions < 3 || (z > 0 && z + 1 < localNodes[2]));

          if (isInterior) {
            const size_t firstRow = localNodeIndex * rowsPerNode;
            for (size_t stencil = 0; stencil < stencilSize; ++stencil) {
              const GO neighborLocalNode = localNode + localStencilNodeOffsets[stencil];
              for (size_t colDof = 0; colDof < rowsPerNode; ++colDof) {
                const LO colLid = static_cast<LO>(
                    neighborLocalNode * dofsPerNodeGO + static_cast<GO>(colDof));
                for (size_t rowDof = 0; rowDof < rowsPerNode; ++rowDof)
                  colind(rowptr(firstRow + rowDof) + columnOffset) = colLid;
                ++columnOffset;
              }
            }
            return;
          }

          Kokkos::Array<LO, maxSupportedStencilSize * maxSupportedDofsPerNode> nodeColumns;
          for (size_t stencil = 0; stencil < stencilSize; ++stencil) {
            int neighborRank = myRank;
            GO neighborNode  = 0;
            if (!StructuredRAPFactoryDetails::resolveNeighbor(
                    x, y, z, stencilOffsets[stencil], localNodes, procGrid, rankDataDevice,
                    myProcX, myProcY, myProcZ, neighborRank, neighborNode))
              continue;

            for (size_t colDof = 0; colDof < rowsPerNode; ++colDof) {
              const GO colGid = globalMinGid + neighborNode * dofsPerNodeGO + static_cast<GO>(colDof);
              LO colLid       = invalidLocalOrdinal;
              if (neighborRank == myRank) {
                colLid = static_cast<LO>(colGid - localMinGid);
              } else {
                size_t region = 0;
                while (region < numRemoteRegions && remoteRegionRanks[region] != neighborRank)
                  ++region;
                if (region < numRemoteRegions) {
                  const GO neighborNx           = remoteRegionNx[region];
                  const GO neighborNy           = remoteRegionNy[region];
                  const GO localNeighborNode    = neighborNode - remoteRegionFirstNodes[region];
                  const GO neighborX            = localNeighborNode % neighborNx;
                  const GO neighborY            = (localNeighborNode / neighborNx) % neighborNy;
                  const GO neighborZ            = localNeighborNode / (neighborNx * neighborNy);
                  const StencilOffset direction = remoteRegionDirections[region];
                  size_t regionNode             = 0;
                  size_t stride                 = 1;
                  if (direction.x == 0) {
                    regionNode += static_cast<size_t>(neighborX) * stride;
                    stride *= static_cast<size_t>(neighborNx);
                  }
                  if (direction.y == 0) {
                    regionNode += static_cast<size_t>(neighborY) * stride;
                    stride *= static_cast<size_t>(neighborNy);
                  }
                  if (direction.z == 0)
                    regionNode += static_cast<size_t>(neighborZ) * stride;
                  colLid = static_cast<LO>(localNumRows + remoteRegionOffsets[region] +
                                           regionNode * rowsPerNode + colDof);
                }
              }

              if (debug && colLid == invalidLocalOrdinal)
                Kokkos::atomic_exchange(&invalidColumn(), 1);

              nodeColumns[columnOffset++] = colLid;
            }
          }

          for (size_t entry = 1; entry < columnOffset; ++entry) {
            const LO value   = nodeColumns[entry];
            size_t insertion = entry;
            while (insertion > 0 && value < nodeColumns[insertion - 1]) {
              nodeColumns[insertion] = nodeColumns[insertion - 1];
              --insertion;
            }
            nodeColumns[insertion] = value;
          }

          const size_t firstRow = localNodeIndex * rowsPerNode;
          for (size_t rowDof = 0; rowDof < rowsPerNode; ++rowDof) {
            const size_t rowStart = rowptr(firstRow + rowDof);
            for (size_t entry = 0; entry < columnOffset; ++entry)
              colind(rowStart + entry) = nodeColumns[entry];
          }
        });
    executionSpace.fence("StructuredRAP: graph column indices ready");
  }

  if (debug) {
    int invalidColumnHost = 0;
    Kokkos::deep_copy(executionSpace, invalidColumnHost, invalidColumn);
    TEUCHOS_TEST_FOR_EXCEPTION(invalidColumnHost != 0, Exceptions::RuntimeError,
                               "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                           << "): a generated column GID was not found in the coarse column map.");

    size_t unsortedRows = 0;
    Kokkos::parallel_reduce(
        "StructuredRAP: validate sorted graph rows", range_policy(executionSpace, 0, localNumRows),
        KOKKOS_LAMBDA(const size_t row, size_t& invalid) {
          for (size_t entry = rowptr(row) + 1; entry < rowptr(row + 1); ++entry) {
            if (colind(entry) < colind(entry - 1)) {
              ++invalid;
              break;
            }
          }
        },
        unsortedRows);
    TEUCHOS_TEST_FOR_EXCEPTION(unsortedRows != 0, Exceptions::RuntimeError,
                               "StructuredRAPFactory::GetStructuredGraph(" << graphSpec.description
                                                                           << "): generated graph rows are not sorted.");
  }

  local_graph_type localGraph(colind, rowptr);
  RCP<CrsGraph> graph;
  {
    graph = CrsGraphFactory::Build(
        localGraph, rowMap, colMap, rowMap, rowMap, paramList);
  }
  {
    using values_type = typename Matrix::local_matrix_type::values_type;
    values_type values(
        Kokkos::ViewAllocateWithoutInitializing("StructuredRAP: matrix values"), localNnz);
    Ac = MatrixFactory::Build(graph, values, paramList);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  const bool doTranspose           = true;
  const bool doFillComplete        = true;
  const bool doOptimizeStorage     = true;
  const Teuchos::ParameterList& pL = GetParameterList();
  const bool prebuildCoarseGraph   = pL.get<bool>("rap: prebuild coarse graph");
  const bool useRAPDelegate        = !prebuildCoarseGraph;

  TEUCHOS_TEST_FOR_EXCEPTION(
      !pL.get<bool>("transpose: use implicit"), Exceptions::RuntimeError,
      "StructuredRAPFactory requires \"transpose: use implicit\" = true because "
      "the prebuilt coarse graph assumes the Galerkin product P^T A P.");

  RCP<Matrix> Ac;

  TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_ == false, Exceptions::RuntimeError,
                             "MueLu::RAPFactory::Build(): CallDeclareInput has not been called before Build!");

  if (useRAPDelegate) {
    if (coarseLevel.IsAvailable("RAP reuse data", this)) {
      RCP<ParameterList> RAPparams = coarseLevel.Get<RCP<ParameterList>>("RAP reuse data", this);
      coarseLevel.Set("RAP reuse data", RAPparams, rapFactoryDelegate_.get());
    }

    // RAPFactory does not accept a prebuilt Ac graph, so delegate only when
    // coarse-graph prebuilding is disabled.
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

    TEUCHOS_TEST_FOR_EXCEPTION(
        pL.get<bool>("rap: triple product") == false, Exceptions::RuntimeError,
        "StructuredRAPFactory requires \"rap: triple product\" = true.");

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
        const int interpolationOrder        = Get<int>(fineLevel, "structuredInterpolationOrder");
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
