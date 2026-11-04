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

#include <cstdint>
#include <sstream>
#include <vector>

#include <Kokkos_Core.hpp>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_TripleMatrixMultiply.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_TpetraCrsGraph.hpp>
#include <Xpetra_TpetraMap.hpp>

#include <Tpetra_CreateCrsGraphFromStencil.hpp>

#include "MueLu_StructuredRAPFactory_decl.hpp"
#include "MueLu_StructuredRAPKernel.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_NoFactory.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Behavior.hpp"
#include "MueLu_RAPFactory_def.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_CommHelpers.hpp"

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

}  // namespace StructuredRAPFactoryDetails

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::StructuredRAPFactory()
  : hasDeclaredInput_(false)
  , rapFactoryDelegate_(rcp(new MueLu::RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>())) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  // If no triple product implementation is prescribed, use the structured kernel by default.
  validParamList->set<std::string>(
      "rap: triple product implementation", "structured",
      "Implementation used with a prebuilt coarse graph: unstructured or structured.");

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
  validParamList->set<RCP<const FactoryBase>>("numDimensions", null, "Number of spatial dimensions.");
  validParamList->set<RCP<const FactoryBase>>("lNodesPerDim", null, "Local number of fine-grid nodes per spatial dimension.");
  validParamList->set<RCP<const FactoryBase>>("lCoarseNodesPerDim", null, "Local number of coarse-grid nodes per spatial dimension.");
  validParamList->set<RCP<const FactoryBase>>("structuredInterpolationOrder", null, "Interpolation order used to construct the structured prolongator.");

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
RCP<const FactoryBase>
StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetStructuredMetadataFactory(
    // Used to get lCoarseNodesPerDim and structuredInterpolationOrder from the prolongator factory
    const std::string& varName) const {
  const RCP<const Factory> pFactory =
      rcp_dynamic_cast<const Factory>(GetFactory("P"));
  if (!pFactory.is_null()) {
    const RCP<const FactoryBase> metadataFactory = pFactory->GetFactory(varName);
    if (!metadataFactory.is_null())
      return metadataFactory;
  }

  return GetFactory(varName);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  const Teuchos::ParameterList& pL = GetParameterList();
  TEUCHOS_TEST_FOR_EXCEPTION(
      !pL.get<bool>("transpose: use implicit"), Exceptions::RuntimeError,
      "StructuredRAPFactory requires \"transpose: use implicit\" = true because "
      "the prebuilt coarse graph assumes the Galerkin product P^T A P.");

  const bool useRAPDelegate =
      !pL.get<bool>("rap: prebuild coarse graph");

  if (useRAPDelegate) {
    ConfigureRAPFactoryDelegate();
    coarseLevel.DeclareInput("A", rapFactoryDelegate_.get(), this);
    coarseLevel.DeclareInput("RAP reuse data", rapFactoryDelegate_.get(), this);
    hasDeclaredInput_ = true;
    return;
  }

  Input(fineLevel, "A");
  Input(coarseLevel, "P");

  const RCP<const FactoryBase> coarseDimensionsFactory =
      GetStructuredMetadataFactory("lCoarseNodesPerDim");
  const RCP<const FactoryBase> interpolationOrderFactory =
      GetStructuredMetadataFactory("structuredInterpolationOrder");
  fineLevel.DeclareInput("lCoarseNodesPerDim", coarseDimensionsFactory.get(), this);
  fineLevel.DeclareInput("structuredInterpolationOrder", interpolationOrderFactory.get(), this);
  if (fineLevel.GetLevelID() == 0) {
    TEUCHOS_TEST_FOR_EXCEPTION(!fineLevel.IsAvailable("numDimensions", NoFactory::get()),
                               Exceptions::RuntimeError,
                               "StructuredRAPFactory: numDimensions was not provided on level 0.");
    TEUCHOS_TEST_FOR_EXCEPTION(!fineLevel.IsAvailable("lNodesPerDim", NoFactory::get()),
                               Exceptions::RuntimeError,
                               "StructuredRAPFactory: lNodesPerDim was not provided on level 0.");
    fineLevel.DeclareInput("numDimensions", NoFactory::get(), this);
    fineLevel.DeclareInput("lNodesPerDim", NoFactory::get(), this);
  } else {
    Input(fineLevel, "numDimensions");
    Input(fineLevel, "lNodesPerDim");
  }

  // call DeclareInput of all user-given transfer factories
  for (std::vector<RCP<const FactoryBase>>::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it)
    (*it)->CallDeclareInput(coarseLevel);

  hasDeclaredInput_ = true;
}

// Infer the fine-grid stencil of A from one representative interior node
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::FineStencilSpec
StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DetectFineStencil(
    const Matrix& A, const int numDimensions,
    const Teuchos::Array<LocalOrdinal>& lFineNodesPerDim) const {
  const std::string prefix = "StructuredRAPFactory::DetectFineStencil: ";

  TEUCHOS_TEST_FOR_EXCEPTION(numDimensions < 1 || numDimensions > 3, Exceptions::RuntimeError,
                             prefix << "the number of dimensions must be between one and three.");
  TEUCHOS_TEST_FOR_EXCEPTION(lFineNodesPerDim.size() < numDimensions, Exceptions::RuntimeError,
                             prefix << "insufficient fine-grid dimensions were supplied.");

  const auto rowMap = A.getRowMap();
  const auto colMap = A.getColMap();
  const auto comm   = rowMap->getComm();
  const auto graph  = A.getCrsGraph();
  TEUCHOS_TEST_FOR_EXCEPTION(graph.is_null(), Exceptions::RuntimeError,
                             prefix << "A does not expose a CRS graph.");

  const LO dofsPerNode = Teuchos::as<LO>(A.GetFixedBlockSize());
  TEUCHOS_TEST_FOR_EXCEPTION(dofsPerNode <= 0, Exceptions::RuntimeError,
                             prefix << "A must have a positive fixed block size.");

  Kokkos::Array<LO, 3> fineNodes{{LO(1), LO(1), LO(1)}};
  Kokkos::Array<LO, 3> interior{{LO(0), LO(0), LO(0)}};
  size_t numFineNodes   = 1;
  bool hasLocalInterior = true;
  // Select interior node in the local subdomain as a representative node for detecting the fine-grid stencil
  for (int dim = 0; dim < numDimensions; ++dim) {
    fineNodes[dim]   = lFineNodesPerDim[dim];
    hasLocalInterior = hasLocalInterior && fineNodes[dim] >= LO(3);
    if (fineNodes[dim] > LO(0)) {
      interior[dim] = fineNodes[dim] / LO(2);
      numFineNodes *= Teuchos::as<size_t>(fineNodes[dim]);
    } else {
      numFineNodes = 0;
    }
  }

  const bool localDimensionsMatch =
      rowMap->getLocalNumElements() ==
      numFineNodes * Teuchos::as<size_t>(dofsPerNode);
  int invalidLocalDimensions = localDimensionsMatch ? 0 : 1;
  int invalidDimensions      = 0;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1,
                     &invalidLocalDimensions, &invalidDimensions);
  TEUCHOS_TEST_FOR_EXCEPTION(
      invalidDimensions != 0, Exceptions::RuntimeError,
      prefix << "the structured local dimensions and fixed block size do not match A on every rank.");

  const LO interiorNode = interior[2] * fineNodes[0] * fineNodes[1] +
                          interior[1] * fineNodes[0] + interior[0];

  const auto localA            = A.getLocalMatrixDevice();
  const LO invalidLocalOrdinal = Teuchos::OrdinalTraits<LO>::invalid();
  std::vector<std::vector<LO>> ownedColumns(
      Teuchos::as<size_t>(dofsPerNode));
  bool hasCompleteLocalStencil = hasLocalInterior;
  if (hasLocalInterior) {
    const size_t firstRow     = Teuchos::as<size_t>(interiorNode * dofsPerNode);
    const auto detectedRowMap = Kokkos::subview(
        localA.graph.row_map,
        std::make_pair(firstRow,
                       firstRow + Teuchos::as<size_t>(dofsPerNode) + 1));
    const auto detectedRowMapHost = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), detectedRowMap);
    const size_t firstEntry = detectedRowMapHost(0);
    const size_t lastEntry =
        detectedRowMapHost(Teuchos::as<size_t>(dofsPerNode));
    const auto detectedEntries = Kokkos::subview(
        localA.graph.entries, std::make_pair(firstEntry, lastEntry));
    using device_type     = typename Node::device_type;
    using execution_space = typename device_type::execution_space;
    Kokkos::View<LO*, device_type> detectedOwnedColumns(
        Kokkos::ViewAllocateWithoutInitializing(
            "StructuredRAP: detected owned columns"),
        lastEntry - firstEntry);
    const auto rowLocalMap = rowMap->getLocalMap();
    const auto colLocalMap = colMap->getLocalMap();
    Kokkos::parallel_for(
        "StructuredRAP: map detected stencil columns",
        Kokkos::RangePolicy<execution_space>(0, lastEntry - firstEntry),
        KOKKOS_LAMBDA(const size_t entry) {
          const GO globalColumn =
              colLocalMap.getGlobalElement(detectedEntries(entry));
          detectedOwnedColumns(entry) =
              rowLocalMap.getLocalElement(globalColumn);
        });
    const auto detectedOwnedColumnsHost =
        Kokkos::create_mirror_view_and_copy(
            Kokkos::HostSpace(), detectedOwnedColumns);

    for (LO rowDof = 0; rowDof < dofsPerNode; ++rowDof) {
      const size_t rowBegin = detectedRowMapHost(Teuchos::as<size_t>(rowDof));
      const size_t rowEnd =
          detectedRowMapHost(Teuchos::as<size_t>(rowDof) + 1);
      if (rowBegin == rowEnd)
        hasCompleteLocalStencil = false;
      auto& rowColumns = ownedColumns[Teuchos::as<size_t>(rowDof)];
      rowColumns.reserve(rowEnd - rowBegin);
      for (size_t entry = rowBegin; entry < rowEnd; ++entry) {
        const LO ownedColumn =
            detectedOwnedColumnsHost(entry - firstEntry);
        rowColumns.push_back(ownedColumn);
        if (ownedColumn == invalidLocalOrdinal)
          hasCompleteLocalStencil = false;
      }
    }
  }

  const int myRank        = comm->getRank();
  const int numRanks      = comm->getSize();
  const int candidateRank = hasCompleteLocalStencil ? myRank : numRanks;
  int detectorRank        = numRanks;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &candidateRank, &detectorRank);
  TEUCHOS_TEST_FOR_EXCEPTION(
      detectorRank == numRanks, Exceptions::RuntimeError,
      prefix << "no MPI rank owns a local interior node whose complete stencil is locally owned; "
             << "increase the local structured subdomain size.");

  FineStencilSpec stencil;
  stencil.numDimensions = numDimensions;
  stencil.dofsPerNode   = dofsPerNode;

  // Only the elected rank performs detection
  if (myRank == detectorRank) {
    for (LO rowDof = 0; rowDof < dofsPerNode; ++rowDof) {
      const auto& rowColumns = ownedColumns[Teuchos::as<size_t>(rowDof)];
      const size_t rowLength = rowColumns.size();
      for (size_t entry = 0; entry < rowLength; ++entry) {
        const LO ownedColumn = rowColumns[entry];
        const LO columnNode  = ownedColumn / dofsPerNode;
        const LO columnDof   = ownedColumn % dofsPerNode;
        LO columnX, columnY, columnZ;
        StructuredRAPFactoryDetails::getLocalNodeIndices(
            columnNode, fineNodes, columnX, columnY, columnZ);
        const StencilOffset offset{
            Teuchos::as<int>(columnX - interior[0]),
            Teuchos::as<int>(columnY - interior[1]),
            Teuchos::as<int>(columnZ - interior[2])};

        bool knownOffset = false;
        for (const StencilOffset& known : stencil.stencilOffsets) {
          if (known.x == offset.x && known.y == offset.y && known.z == offset.z) {
            knownOffset = true;
            break;
          }
        }
        if (!knownOffset)
          stencil.stencilOffsets.push_back(offset);

        stencil.entries.push_back(FineStencilEntry{
            offset, rowDof, columnDof});
      }
    }
  }

  int counts[2] = {Teuchos::as<int>(stencil.stencilOffsets.size()),
                   Teuchos::as<int>(stencil.entries.size())};
  Teuchos::broadcast<int, int>(*comm, detectorRank, 2, counts);

  Teuchos::Array<int> offsetData(3 * counts[0]);
  Teuchos::Array<int> entryOffsetData(3 * counts[1]);
  Teuchos::Array<LO> entryData(2 * counts[1]);
  if (myRank == detectorRank) {
    for (int offset = 0; offset < counts[0]; ++offset) {
      offsetData[3 * offset]     = stencil.stencilOffsets[offset].x;
      offsetData[3 * offset + 1] = stencil.stencilOffsets[offset].y;
      offsetData[3 * offset + 2] = stencil.stencilOffsets[offset].z;
    }
    for (int entry = 0; entry < counts[1]; ++entry) {
      entryOffsetData[3 * entry]     = stencil.entries[entry].offset.x;
      entryOffsetData[3 * entry + 1] = stencil.entries[entry].offset.y;
      entryOffsetData[3 * entry + 2] = stencil.entries[entry].offset.z;
      entryData[2 * entry]           = stencil.entries[entry].rowDof;
      entryData[2 * entry + 1]       = stencil.entries[entry].columnDof;
    }
  }

  if (counts[0] > 0)
    Teuchos::broadcast<int, int>(*comm, detectorRank, 3 * counts[0],
                                 offsetData.getRawPtr());
  if (counts[1] > 0) {
    Teuchos::broadcast<int, int>(*comm, detectorRank, 3 * counts[1],
                                 entryOffsetData.getRawPtr());
    Teuchos::broadcast<int, LO>(*comm, detectorRank, 2 * counts[1],
                                entryData.getRawPtr());
  }

  if (myRank != detectorRank) {
    stencil.stencilOffsets.reserve(counts[0]);
    for (int offset = 0; offset < counts[0]; ++offset)
      stencil.stencilOffsets.push_back(StencilOffset{
          offsetData[3 * offset], offsetData[3 * offset + 1],
          offsetData[3 * offset + 2]});

    stencil.entries.reserve(counts[1]);
    for (int entry = 0; entry < counts[1]; ++entry)
      stencil.entries.push_back(FineStencilEntry{
          StencilOffset{entryOffsetData[3 * entry],
                        entryOffsetData[3 * entry + 1],
                        entryOffsetData[3 * entry + 2]},
          entryData[2 * entry], entryData[2 * entry + 1]});
  }

  TEUCHOS_TEST_FOR_EXCEPTION(stencil.stencilOffsets.empty(), Exceptions::RuntimeError,
                             prefix << "the representative interior rows are empty.");
  return stencil;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::StructuredGraphSpec
StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeriveCoarseRAPStencil(
    const FineStencilSpec& fineStencil, const int interpolationOrder) const {
  const std::string prefix = "StructuredRAPFactory::DeriveCoarseRAPStencil: ";

  TEUCHOS_TEST_FOR_EXCEPTION(fineStencil.numDimensions < 1 || fineStencil.numDimensions > 3,
                             Exceptions::RuntimeError,
                             prefix << "the number of dimensions must be between one and three.");
  TEUCHOS_TEST_FOR_EXCEPTION(fineStencil.dofsPerNode <= 0, Exceptions::RuntimeError,
                             prefix << "dofsPerNode must be positive.");
  TEUCHOS_TEST_FOR_EXCEPTION(fineStencil.stencilOffsets.empty(), Exceptions::RuntimeError,
                             prefix << "the fine-grid stencil is empty.");
  TEUCHOS_TEST_FOR_EXCEPTION(interpolationOrder < 0 || interpolationOrder > 1,
                             Exceptions::RuntimeError,
                             prefix << "interpolation order " << interpolationOrder
                                    << " is not supported; expected zero or one.");
  for (const StencilOffset& offset : fineStencil.stencilOffsets) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        offset.x < -1 || offset.x > 1 ||
            offset.y < -1 || offset.y > 1 ||
            offset.z < -1 || offset.z > 1,
        Exceptions::RuntimeError,
        prefix << "fine-grid stencil offset (" << offset.x << ", "
               << offset.y << ", " << offset.z
               << ") exceeds the supported radius one; disable "
                  "\"rap: prebuild coarse graph\" for wider stencils.");
  }

  StructuredGraphSpec coarseStencil;
  coarseStencil.numDimensions = fineStencil.numDimensions;
  coarseStencil.dofsPerNode   = fineStencil.dofsPerNode;

  if (interpolationOrder == 0) {
    coarseStencil.stencilOffsets = fineStencil.stencilOffsets;
    bool hasCenter               = false;
    for (const StencilOffset& offset : coarseStencil.stencilOffsets)
      hasCenter = hasCenter || (offset.x == 0 && offset.y == 0 && offset.z == 0);
    if (!hasCenter)
      coarseStencil.stencilOffsets.push_back(StencilOffset{0, 0, 0});
  } else {
    // The supports of neighboring piecewise-linear basis functions overlap.
    // For a radius-one fine stencil this produces the full radius-one coarse
    // stencil in every active dimension.
    const int minY = fineStencil.numDimensions > 1 ? -1 : 0;
    const int maxY = fineStencil.numDimensions > 1 ? 1 : 0;
    const int minZ = fineStencil.numDimensions > 2 ? -1 : 0;
    const int maxZ = fineStencil.numDimensions > 2 ? 1 : 0;
    for (int z = minZ; z <= maxZ; ++z)
      for (int y = minY; y <= maxY; ++y)
        for (int x = -1; x <= 1; ++x)
          coarseStencil.stencilOffsets.push_back(StencilOffset{x, y, z});
  }

  int minimum[3] = {coarseStencil.stencilOffsets[0].x,
                    coarseStencil.stencilOffsets[0].y,
                    coarseStencil.stencilOffsets[0].z};
  int maximum[3] = {minimum[0], minimum[1], minimum[2]};
  for (const StencilOffset& offset : coarseStencil.stencilOffsets) {
    const int components[3] = {offset.x, offset.y, offset.z};
    for (int dim = 0; dim < fineStencil.numDimensions; ++dim) {
      if (components[dim] < minimum[dim])
        minimum[dim] = components[dim];
      if (components[dim] > maximum[dim])
        maximum[dim] = components[dim];
    }
  }

  std::ostringstream description;
  for (int dim = 0; dim < fineStencil.numDimensions; ++dim) {
    if (dim != 0)
      description << "x";
    description << maximum[dim] - minimum[dim] + 1;
  }
  coarseStencil.description = description.str();
  return coarseStencil;
}

// Prebuild sparsity structure of coarse matrix
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetStructuredGraph(
    RCP<Matrix>& Ac, const RCP<Matrix>& P,
    const Teuchos::Array<LocalOrdinal>& lCoarseNodesPerDim,
    const StructuredGraphSpec& graphSpec) const {
  TEUCHOS_TEST_FOR_EXCEPTION(
      P->getDomainMap()->lib() != Xpetra::UseTpetra,
      Exceptions::RuntimeError,
      "StructuredRAPFactory::GetStructuredGraph requires the Tpetra backend.");

  Teuchos::Array<Tpetra::Experimental::StencilOffset> stencilOffsets;
  stencilOffsets.reserve(graphSpec.stencilOffsets.size());
  for (const StencilOffset& offset : graphSpec.stencilOffsets)
    stencilOffsets.push_back(
        Tpetra::Experimental::StencilOffset{offset.x, offset.y, offset.z});

  const auto tpetraRowMap = Xpetra::toTpetra(P->getDomainMap());
  const auto tpetraGraph =
      Tpetra::Experimental::createCrsGraphFromStencil(
          tpetraRowMap,
          lCoarseNodesPerDim(0, graphSpec.numDimensions),
          stencilOffsets(), graphSpec.dofsPerNode, Behavior::debug());
  RCP<CrsGraph> graph = Teuchos::rcp_const_cast<CrsGraph>(
      Xpetra::toXpetra(tpetraGraph));

  RCP<ParameterList> paramList = rcp(new ParameterList);
  paramList->set("No Nonlocal Changes", true);
  paramList->set("Optimize Storage", true);
  paramList->set("compute global constants", true);
  using values_type = typename Matrix::local_matrix_type::values_type;
  values_type values(
      Kokkos::ViewAllocateWithoutInitializing(
          "StructuredRAP: matrix values"),
      tpetraGraph->getLocalNumEntries());
  Ac = MatrixFactory::Build(graph, values, paramList);
  return;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  const bool doTranspose           = true;
  const bool doFillComplete        = true;
  const bool doOptimizeStorage     = true;
  const Teuchos::ParameterList& pL = GetParameterList();
  const bool useRAPDelegate =
      !pL.get<bool>("rap: prebuild coarse graph");

  TEUCHOS_TEST_FOR_EXCEPTION(
      !pL.get<bool>("transpose: use implicit"), Exceptions::RuntimeError,
      "StructuredRAPFactory requires \"transpose: use implicit\" = true because "
      "the prebuilt coarse graph assumes the Galerkin product P^T A P.");

  RCP<Matrix> Ac;

  TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_ == false, Exceptions::RuntimeError,
                             "MueLu::RAPFactory::Build(): CallDeclareInput has not been called before Build!");

  if (useRAPDelegate) {
    if (pL.get<std::string>("rap: triple product implementation") ==
        "structured") {
      GetOStream(Runtime0)
          << "StructuredRAP: Coarse graph prebuilding is disabled, "
             "redirecting the structured triple product to Xpetra."
          << std::endl;
    }
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
    const std::string labelstr = FormattingHelper::getColonLabel(coarseLevel.getObjectLabel());
    FactoryMonitor m(*this, "Computing Ac", coarseLevel);

    std::ostringstream levelstr;
    levelstr << coarseLevel.GetLevelID();

    TEUCHOS_TEST_FOR_EXCEPTION(
        pL.get<bool>("rap: triple product") == false, Exceptions::RuntimeError,
        "StructuredRAPFactory requires \"rap: triple product\" = true.");

    RCP<Matrix> A = Get<RCP<Matrix>>(fineLevel, "A");
    RCP<Matrix> P = Get<RCP<Matrix>>(coarseLevel, "P");
    const std::string tripleProductImplementation =
        pL.get<std::string>("rap: triple product implementation");
    TEUCHOS_TEST_FOR_EXCEPTION(
        tripleProductImplementation != "unstructured" && tripleProductImplementation != "structured",
        Exceptions::RuntimeError,
        "StructuredRAPFactory: unknown value \"" << tripleProductImplementation
                                                 << "\" for \"rap: triple product implementation\". "
                                                    "Valid values are \"unstructured\" and \"structured\".");

    const RCP<const FactoryBase> coarseDimensionsFactory =
        GetStructuredMetadataFactory("lCoarseNodesPerDim");
    const RCP<const FactoryBase> interpolationOrderFactory =
        GetStructuredMetadataFactory("structuredInterpolationOrder");
    const Teuchos::Array<LocalOrdinal> lCoarseNodesPerDim =
        fineLevel.Get<Teuchos::Array<LocalOrdinal>>(
            "lCoarseNodesPerDim", coarseDimensionsFactory.get());
    const int interpolationOrder =
        fineLevel.Get<int>("structuredInterpolationOrder",
                           interpolationOrderFactory.get());
    const int numDimensions =
        fineLevel.GetLevelID() == 0
            ? fineLevel.Get<int>("numDimensions", NoFactory::get())
            : Get<int>(fineLevel, "numDimensions");
    const Teuchos::Array<LocalOrdinal> lFineNodesPerDim =
        fineLevel.GetLevelID() == 0
            ? fineLevel.Get<Teuchos::Array<LocalOrdinal>>(
                  "lNodesPerDim", NoFactory::get())
            : Get<Teuchos::Array<LocalOrdinal>>(fineLevel, "lNodesPerDim");
    Teuchos::Array<int> structuredCoarseningRate;
    if (tripleProductImplementation == "structured") {
      // lCoarseNodesPerDim and structuredInterpolationOrder are produced by
      // StructuredAggregationFactory. Read the rate from that same producer
      // so the structured RAP path uses the existing
      // "aggregation: coarsening rate" setting instead of requiring a second,
      // potentially inconsistent RAP-specific parameter.
      const Factory* structuredFactory =
          dynamic_cast<const Factory*>(coarseDimensionsFactory.get());
      TEUCHOS_TEST_FOR_EXCEPTION(
          structuredFactory == nullptr, Exceptions::RuntimeError,
          "StructuredRAPFactory: the factory producing lCoarseNodesPerDim "
          "does not expose its aggregation parameters.");
      const ParameterList& structuredParameters =
          structuredFactory->GetParameterList();
      TEUCHOS_TEST_FOR_EXCEPTION(
          !structuredParameters.isType<std::string>(
              "aggregation: coarsening rate"),
          Exceptions::RuntimeError,
          "StructuredRAPFactory: the factory producing lCoarseNodesPerDim "
          "does not provide the string parameter \"aggregation: coarsening "
          "rate\".");
      try {
        structuredCoarseningRate = Teuchos::fromStringToArray<int>(
            structuredParameters.get<std::string>(
                "aggregation: coarsening rate"));
      } catch (const Teuchos::InvalidArrayStringRepresentation& e) {
        GetOStream(Errors, -1)
            << " *** \"aggregation: coarsening rate\" must be a string "
               "convertible into an array! *** "
            << std::endl;
        throw e;
      }
    }
    if (P.is_null()) {
      Ac = Teuchos::null;
      Set(coarseLevel, "A", Ac);
      return;
    }

    const FineStencilSpec fineStencil =
        DetectFineStencil(*A, numDimensions, lFineNodesPerDim);
    const StructuredGraphSpec graphSpec =
        DeriveCoarseRAPStencil(fineStencil, interpolationOrder);

    {
      RCP<ParameterList> RAPparams;
      {
        SubFactoryMonitor mGraph(*this, "Prebuilding coarse Ac graph", coarseLevel);

        RAPparams = rcp(new ParameterList);
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
          Ac->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());
        } else {
          // If reuse data is not available, prebuild the sparse graph from the
          // detected fine stencil and its derived coarse RAP stencil.
          GetOStream(Statistics1) << "StructuredRAP: Using " << graphSpec.description
                                  << " stencil with " << graphSpec.stencilOffsets.size()
                                  << " nodal entries." << std::endl;
          GetStructuredGraph(Ac, P, lCoarseNodesPerDim, graphSpec);
        }
      }

      // We *always* need global constants for the RAP, but not for the temps
      RAPparams->set("compute global constants: temporaries", RAPparams->get("compute global constants: temporaries", false));
      RAPparams->set("compute global constants", true);

      if (tripleProductImplementation == "unstructured") {
        SubFactoryMonitor m2(*this, "MxMxM: Xpetra P^T x A x P (implicit)", coarseLevel);
        Xpetra::TripleMatrixMultiply<SC, LO, GO, NO>::
            MultiplyRAP(*P, doTranspose, *A, !doTranspose, *P, !doTranspose, *Ac, doFillComplete,
                        doOptimizeStorage, labelstr + std::string("MueLu::Xpetra-P^T*A*P-implicit-") + levelstr.str(),
                        RAPparams);
      } else {
        SubFactoryMonitor m2(*this, "MxMxM: Structured P^T x A x P (implicit)", coarseLevel);
        Details::StructuredRAPKernel<SC, LO, GO, NO>::Compute(
            *A, *P, *Ac, fineStencil, interpolationOrder,
            lFineNodesPerDim, lCoarseNodesPerDim,
            structuredCoarseningRate, RAPparams);

        Ac->CreateView("stridedMaps", P, true, P, false);
      }

      {
        GetOStream(Statistics1) << "StructuredRAP: Ac nnz: local = "
                                << Ac->getLocalNumEntries()
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

        std::ostringstream oss;
        oss << "A_" << coarseLevel.GetLevelID();
        Ac->setObjectLabel(oss.str());
        Set(coarseLevel, "A", Ac);

        RAPparams->set("graph", Ac);
        Set(coarseLevel, "RAP reuse data", RAPparams);
      }
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
