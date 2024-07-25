// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_BLOCKEDCOORDINATESTRANSFER_FACTORY_DEF_HPP
#define MUELU_BLOCKEDCOORDINATESTRANSFER_FACTORY_DEF_HPP

#include "Xpetra_ImportFactory.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_IO.hpp"

#include "MueLu_BlockedCoordinatesTransferFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> BlockedCoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Factory for coordinates generation");
  validParamList->set<RCP<const FactoryBase> >("CoarseMap", Teuchos::null, "Generating factory of the coarse map");
  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& /* fineLevel */, Level& coarseLevel) const {
  Input(coarseLevel, "CoarseMap");

  // Make sure the Level knows I need these sub-Factories
  const size_t numSubFactories = NumFactories();
  for (size_t i = 0; i < numSubFactories; i++) {
    const RCP<const FactoryBase>& myFactory = subFactories_[i];
    coarseLevel.DeclareInput("Coordinates", myFactory.getRawPtr(), this);
  }

  // call DeclareInput of all user-given transfer factories
  for (std::vector<RCP<const FactoryBase> >::const_iterator it = subFactories_.begin(); it != subFactories_.end(); ++it)
    (*it)->CallDeclareInput(coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& /* fineLevel */, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  typedef Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LO, GO, NO> dMV;
  typedef Xpetra::BlockedMultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LO, GO, NO> dBV;

  GetOStream(Runtime0) << "Transferring (blocked) coordinates" << std::endl;

  const size_t numSubFactories = NumFactories();
  std::vector<RCP<const Map> > subBlockMaps(numSubFactories);
  std::vector<RCP<dMV> > subBlockCoords(numSubFactories);

  if (coarseLevel.IsAvailable("Coordinates", this)) {
    GetOStream(Runtime0) << "Reusing coordinates" << std::endl;
    return;
  }

  // Get components
  for (size_t i = 0; i < numSubFactories; i++) {
    GetOStream(Runtime1) << "Generating Coordinates for block " << i << "/" << numSubFactories << std::endl;
    const RCP<const FactoryBase>& myFactory = subFactories_[i];
    myFactory->CallBuild(coarseLevel);
    subBlockCoords[i] = coarseLevel.Get<RCP<dMV> >("Coordinates", myFactory.get());
    subBlockMaps[i]   = subBlockCoords[i]->getMap();
  }

  // Blocked Map
  RCP<const BlockedMap> coarseCoordMapBlocked;

  {
    // coarseMap is being used to set up the domain map of tentative P, and therefore, the row map of Ac
    // Therefore, if we amalgamate coarseMap, logical nodes in the coordinates vector would correspond to
    // logical blocks in the matrix
    RCP<const BlockedMap> coarseMap = Get<RCP<const BlockedMap> >(coarseLevel, "CoarseMap");
    bool thyraMode                  = coarseMap->getThyraMode();

    ArrayView<const GO> elementAList = coarseMap->getFullMap()->getLocalElementList();

    LO blkSize = 1;
    if (rcp_dynamic_cast<const StridedMap>(coarseMap->getMap(0, thyraMode)) != Teuchos::null)
      blkSize = rcp_dynamic_cast<const StridedMap>(coarseMap->getMap(0, thyraMode))->getFixedBlockSize();

    for (size_t i = 1; i < numSubFactories; i++) {
      LO otherBlkSize = 1;
      if (rcp_dynamic_cast<const StridedMap>(coarseMap->getMap(i, thyraMode)) != Teuchos::null)
        otherBlkSize = rcp_dynamic_cast<const StridedMap>(coarseMap->getMap(i, thyraMode))->getFixedBlockSize();
      TEUCHOS_TEST_FOR_EXCEPTION(otherBlkSize != blkSize, Exceptions::RuntimeError, "BlockedCoordinatesTransferFactory: Subblocks have different Block sizes. This is not yet supported.");
    }

    GO indexBase       = coarseMap->getFullMap()->getIndexBase();
    size_t numElements = elementAList.size() / blkSize;
    Array<GO> elementList(numElements);

    // Amalgamate the map
    for (LO i = 0; i < Teuchos::as<LO>(numElements); i++)
      elementList[i] = (elementAList[i * blkSize] - indexBase) / blkSize + indexBase;

    RCP<const Map> coarseCoordMap = MapFactory::Build(coarseMap->getFullMap()->lib(),
                                                      Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), elementList, indexBase, coarseMap->getFullMap()->getComm());

    coarseCoordMapBlocked = rcp(new BlockedMap(coarseCoordMap, subBlockMaps, thyraMode));
  }

  // Build blocked coordinates vector
  RCP<dBV> bcoarseCoords = rcp(new dBV(coarseCoordMapBlocked, subBlockCoords));

  // Turn the blocked coordinates vector into an unblocked one
  RCP<dMV> coarseCoords = bcoarseCoords->Merge();
  Set<RCP<dMV> >(coarseLevel, "Coordinates", coarseCoords);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddFactory(const RCP<const FactoryBase>& factory) {
  subFactories_.push_back(factory);
}

}  // namespace MueLu

#endif  // MUELU_BLOCKEDCOORDINATESTRANSFER_FACTORY_DEF_HPP
