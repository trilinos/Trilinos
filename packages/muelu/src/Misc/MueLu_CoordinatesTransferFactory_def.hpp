// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_COORDINATESTRANSFER_FACTORY_DEF_HPP
#define MUELU_COORDINATESTRANSFER_FACTORY_DEF_HPP

#include "Xpetra_ImportFactory.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_IO.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_CoordinatesTransferFactory_decl.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
CoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CoordinatesTransferFactory() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
CoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~CoordinatesTransferFactory() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> CoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Factory for coordinates generation");
  validParamList->set<RCP<const FactoryBase> >("Aggregates", Teuchos::null, "Factory for coordinates generation");
  validParamList->set<RCP<const FactoryBase> >("CoarseMap", Teuchos::null, "Generating factory of the coarse map");
  validParamList->set<bool>("structured aggregation", false, "Flag specifying that the geometric data is transferred for StructuredAggregationFactory");
  validParamList->set<bool>("aggregation coupled", false, "Flag specifying if the aggregation algorithm was used in coupled mode.");
  validParamList->set<bool>("Geometric", false, "Flag specifying that the coordinates are transferred for GeneralGeometricPFactory");
  validParamList->set<RCP<const FactoryBase> >("coarseCoordinates", Teuchos::null, "Factory for coarse coordinates generation");
  validParamList->set<RCP<const FactoryBase> >("gCoarseNodesPerDim", Teuchos::null, "Factory providing the global number of nodes per spatial dimensions of the mesh");
  validParamList->set<RCP<const FactoryBase> >("lCoarseNodesPerDim", Teuchos::null, "Factory providing the local number of nodes per spatial dimensions of the mesh");
  validParamList->set<RCP<const FactoryBase> >("numDimensions", Teuchos::null, "Factory providing the number of spatial dimensions of the mesh");
  validParamList->set<int>("write start", -1, "first level at which coordinates should be written to file");
  validParamList->set<int>("write end", -1, "last level at which coordinates should be written to file");
  validParamList->set<bool>("hybrid aggregation", false, "Flag specifying that hybrid aggregation data is transfered for HybridAggregationFactory");
  validParamList->set<RCP<const FactoryBase> >("aggregationRegionTypeCoarse", Teuchos::null, "Factory indicating what aggregation type is to be used on the coarse level of the region");
  validParamList->set<bool>("interface aggregation", false, "Flag specifying that interface aggregation data is transfered for HybridAggregationFactory");
  validParamList->set<RCP<const FactoryBase> >("coarseInterfacesDimensions", Teuchos::null, "Factory providing coarseInterfacesDimensions");
  validParamList->set<RCP<const FactoryBase> >("nodeOnCoarseInterface", Teuchos::null, "Factory providing nodeOnCoarseInterface");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  static bool isAvailableCoords = false;

  const ParameterList& pL = GetParameterList();
  if (pL.get<bool>("structured aggregation") == true) {
    if (pL.get<bool>("aggregation coupled") == true) {
      Input(fineLevel, "gCoarseNodesPerDim");
    }
    Input(fineLevel, "lCoarseNodesPerDim");
    Input(fineLevel, "numDimensions");
  } else if (pL.get<bool>("Geometric") == true) {
    Input(coarseLevel, "coarseCoordinates");
    Input(coarseLevel, "gCoarseNodesPerDim");
    Input(coarseLevel, "lCoarseNodesPerDim");
  } else if (pL.get<bool>("hybrid aggregation") == true) {
    Input(fineLevel, "aggregationRegionTypeCoarse");
    Input(fineLevel, "lCoarseNodesPerDim");
    Input(fineLevel, "numDimensions");
    if (pL.get<bool>("interface aggregation") == true) {
      Input(fineLevel, "coarseInterfacesDimensions");
      Input(fineLevel, "nodeOnCoarseInterface");
    }
  } else {
    if (coarseLevel.GetRequestMode() == Level::REQUEST)
      isAvailableCoords = coarseLevel.IsAvailable("Coordinates", this);

    if (isAvailableCoords == false) {
      Input(fineLevel, "Coordinates");
      Input(fineLevel, "Aggregates");
      Input(fineLevel, "CoarseMap");
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  using xdMV = Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>;

  GetOStream(Runtime0) << "Transferring coordinates" << std::endl;

  int numDimensions;
  RCP<xdMV> coarseCoords;
  RCP<xdMV> fineCoords;
  Array<GO> gCoarseNodesPerDir;
  Array<LO> lCoarseNodesPerDir;

  const ParameterList& pL = GetParameterList();

  if (pL.get<bool>("hybrid aggregation") == true) {
    std::string regionType = Get<std::string>(fineLevel, "aggregationRegionTypeCoarse");
    numDimensions          = Get<int>(fineLevel, "numDimensions");
    lCoarseNodesPerDir     = Get<Array<LO> >(fineLevel, "lCoarseNodesPerDim");
    Set<std::string>(coarseLevel, "aggregationRegionType", regionType);
    Set<int>(coarseLevel, "numDimensions", numDimensions);
    Set<Array<LO> >(coarseLevel, "lNodesPerDim", lCoarseNodesPerDir);

    if ((pL.get<bool>("interface aggregation") == true) && (regionType == "uncoupled")) {
      Array<LO> coarseInterfacesDimensions = Get<Array<LO> >(fineLevel, "coarseInterfacesDimensions");
      Array<LO> nodeOnCoarseInterface      = Get<Array<LO> >(fineLevel, "nodeOnCoarseInterface");
      Set<Array<LO> >(coarseLevel, "interfacesDimensions", coarseInterfacesDimensions);
      Set<Array<LO> >(coarseLevel, "nodeOnInterface", nodeOnCoarseInterface);
    }

  } else if (pL.get<bool>("structured aggregation") == true) {
    if (pL.get<bool>("aggregation coupled") == true) {
      gCoarseNodesPerDir = Get<Array<GO> >(fineLevel, "gCoarseNodesPerDim");
      Set<Array<GO> >(coarseLevel, "gNodesPerDim", gCoarseNodesPerDir);
    }
    lCoarseNodesPerDir = Get<Array<LO> >(fineLevel, "lCoarseNodesPerDim");
    Set<Array<LO> >(coarseLevel, "lNodesPerDim", lCoarseNodesPerDir);
    numDimensions = Get<int>(fineLevel, "numDimensions");
    Set<int>(coarseLevel, "numDimensions", numDimensions);

  } else if (pL.get<bool>("Geometric") == true) {
    coarseCoords       = Get<RCP<xdMV> >(coarseLevel, "coarseCoordinates");
    gCoarseNodesPerDir = Get<Array<GO> >(coarseLevel, "gCoarseNodesPerDim");
    lCoarseNodesPerDir = Get<Array<LO> >(coarseLevel, "lCoarseNodesPerDim");
    Set<Array<GO> >(coarseLevel, "gNodesPerDim", gCoarseNodesPerDir);
    Set<Array<LO> >(coarseLevel, "lNodesPerDim", lCoarseNodesPerDir);

    Set<RCP<xdMV> >(coarseLevel, "Coordinates", coarseCoords);

  } else {
    if (coarseLevel.IsAvailable("Coordinates", this)) {
      GetOStream(Runtime0) << "Reusing coordinates" << std::endl;
      return;
    }

    fineCoords               = Get<RCP<xdMV> >(fineLevel, "Coordinates");
    RCP<const Map> coarseMap = Get<RCP<const Map> >(fineLevel, "CoarseMap");

    // coarseMap is being used to set up the domain map of tentative P, and therefore, the row map of Ac
    // Therefore, if we amalgamate coarseMap, logical nodes in the coordinates vector would correspond to
    // logical blocks in the matrix
    LO blkSize = 1;
    if (rcp_dynamic_cast<const StridedMap>(coarseMap) != Teuchos::null)
      blkSize = rcp_dynamic_cast<const StridedMap>(coarseMap)->getFixedBlockSize();

    RCP<const Map> coarseCoordMap;
    RCP<const Map> uniqueMap = fineCoords->getMap();
    if (blkSize > 1) {
      // If the block size is greater than one, we need to create a coarse coordinate map
      // FIXME: The amalgamation should really be done on device.
      GO indexBase                     = coarseMap->getIndexBase();
      ArrayView<const GO> elementAList = coarseMap->getLocalElementList();
      size_t numElements               = elementAList.size() / blkSize;
      Array<GO> elementList(numElements);

      // Amalgamate the map
      for (LO i = 0; i < Teuchos::as<LO>(numElements); i++)
        elementList[i] = (elementAList[i * blkSize] - indexBase) / blkSize + indexBase;

      {
        SubFactoryMonitor sfm(*this, "MapFactory: coarseCoordMap", fineLevel);
        coarseCoordMap = MapFactory ::Build(coarseMap->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), elementList, indexBase, coarseMap->getComm());
      }
    } else {
      // If the block size is one, we can just use the coarse map for coordinates
      coarseCoordMap = coarseMap;
    }

    // Build the coarseCoords MultiVector
    coarseCoords = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(coarseCoordMap, fineCoords->getNumVectors());

    RCP<Aggregates> aggregates;
    bool aggregatesCrossProcessors;
    aggregates                = Get<RCP<Aggregates> >(fineLevel, "Aggregates");
    aggregatesCrossProcessors = aggregates->AggregatesCrossProcessors();

    // Create overlapped fine coordinates to reduce global communication
    RCP<xdMV> ghostedCoords = fineCoords;
    if (aggregatesCrossProcessors) {
      RCP<const Map> nonUniqueMap = aggregates->GetMap();
      RCP<const Import> importer  = ImportFactory::Build(uniqueMap, nonUniqueMap);

      ghostedCoords = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(nonUniqueMap, fineCoords->getNumVectors());
      ghostedCoords->doImport(*fineCoords, *importer, Xpetra::INSERT);
    }

    // The good news is that this graph has already been constructed for the
    // TentativePFactory and was cached in Aggregates. So this is a no-op.
    auto aggGraph = aggregates->GetGraph();
    auto numAggs  = aggGraph.numRows();

    auto fineCoordsView   = ghostedCoords->getDeviceLocalView(Xpetra::Access::ReadOnly);
    auto coarseCoordsView = coarseCoords->getDeviceLocalView(Xpetra::Access::OverwriteAll);

    // Fill in coarse coordinates
    {
      SubFactoryMonitor m2(*this, "AverageCoords", coarseLevel);

      const auto dim = ghostedCoords->getNumVectors();

      typename AppendTrait<decltype(fineCoordsView), Kokkos::RandomAccess>::type fineCoordsRandomView = fineCoordsView;
      for (size_t j = 0; j < dim; j++) {
        Kokkos::parallel_for(
            "MueLu:CoordinatesTransferF:Build:coord", Kokkos::RangePolicy<local_ordinal_type, execution_space>(0, numAggs),
            KOKKOS_LAMBDA(const LO i) {
              // A row in this graph represents all node ids in the aggregate
              // Therefore, averaging is very easy

              auto aggregate = aggGraph.rowConst(i);

              typename Teuchos::ScalarTraits<Scalar>::magnitudeType sum = 0.0;  // do not use Scalar here (Stokhos)
              for (size_t colID = 0; colID < static_cast<size_t>(aggregate.length); colID++)
                sum += fineCoordsRandomView(aggregate(colID), j);

              coarseCoordsView(i, j) = sum / aggregate.length;
            });
      }
    }

    Set<RCP<xdMV> >(coarseLevel, "Coordinates", coarseCoords);
  }

  int writeStart = pL.get<int>("write start"), writeEnd = pL.get<int>("write end");
  if (writeStart == 0 && fineLevel.GetLevelID() == 0 && writeStart <= writeEnd) {
    std::ostringstream buf;
    buf << fineLevel.GetLevelID();
    std::string fileName = "coordinates_before_rebalance_level_" + buf.str() + ".m";
    Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Write(fileName, *fineCoords);
  }
  if (writeStart <= coarseLevel.GetLevelID() && coarseLevel.GetLevelID() <= writeEnd) {
    std::ostringstream buf;
    buf << coarseLevel.GetLevelID();
    std::string fileName = "coordinates_before_rebalance_level_" + buf.str() + ".m";
    Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Write(fileName, *coarseCoords);
  }
}

}  // namespace MueLu

#endif  // MUELU_COORDINATESTRANSFER_FACTORY_DEF_HPP
