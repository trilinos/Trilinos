// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REBALANCEMAPFACTORY_DEF_HPP_
#define MUELU_REBALANCEMAPFACTORY_DEF_HPP_

#include "MueLu_RebalanceMapFactory_decl.hpp"

#include <Teuchos_Utils.hpp>

#include "MueLu_Exceptions.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> RebalanceMapFactory<LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("repartition: use subcommunicators");
#undef SET_VALID_ENTRY

  // Information about map that is to be rebalanced
  validParamList->set<std::string>("Map name", "", "Name of map to rebalanced.");
  validParamList->set<RCP<const FactoryBase> >("Map factory", MueLu::NoFactory::getRCP(), "Generating factory of map to be rebalanced.");

  // Importer object with rebalancing information
  validParamList->set<RCP<const FactoryBase> >("Importer", Teuchos::null, "Factory of the importer object used for the rebalancing");

  return validParamList;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceMapFactory<LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  const Teuchos::ParameterList &pL           = GetParameterList();
  std::string mapName                        = pL.get<std::string>("Map name");
  Teuchos::RCP<const FactoryBase> mapFactory = GetFactory("Map factory");
  currentLevel.DeclareInput(mapName, mapFactory.get(), this);

  Input(currentLevel, "Importer");
}  // DeclareInput()

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Level &level) const {
  FactoryMonitor m(*this, "Build", level);

  // Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

  // extract data from Level object
  const Teuchos::ParameterList &pL           = GetParameterList();
  std::string mapName                        = pL.get<std::string>("Map name");
  Teuchos::RCP<const FactoryBase> mapFactory = GetFactory("Map factory");

  RCP<const Import> rebalanceImporter = Get<RCP<const Import> >(level, "Importer");

  if (rebalanceImporter != Teuchos::null) {
    // input map (not rebalanced)
    RCP<const Map> map = level.Get<RCP<const Map> >(mapName, mapFactory.get());

    // create vector based on input map
    // Note, that the map can be a part only of the full map stored in rebalanceImporter.getSourceMap()
    RCP<Vector> v = VectorFactory::Build(map);
    v->putScalar(1.0);

    // create a new vector based on the full rebalanceImporter.getSourceMap()
    // import the partial map information to the full source map
    RCP<const Import> blowUpImporter = ImportFactory::Build(map, rebalanceImporter->getSourceMap());
    RCP<Vector> pv                   = VectorFactory::Build(rebalanceImporter->getSourceMap());
    pv->doImport(*v, *blowUpImporter, Xpetra::INSERT);

    // do rebalancing using rebalanceImporter
    RCP<Vector> ptv = VectorFactory::Build(rebalanceImporter->getTargetMap());
    ptv->doImport(*pv, *rebalanceImporter, Xpetra::INSERT);

    if (pL.get<bool>("repartition: use subcommunicators") == true)
      ptv->replaceMap(ptv->getMap()->removeEmptyProcesses());

    // reconstruct rebalanced partial map
    Teuchos::ArrayRCP<const Scalar> ptvData = ptv->getData(0);
    std::vector<GlobalOrdinal> localGIDs;  // vector with GIDs that are stored on current proc

    for (size_t k = 0; k < ptv->getLocalLength(); k++) {
      if (ptvData[k] == 1.0) {
        localGIDs.push_back(ptv->getMap()->getGlobalElement(k));
      }
    }

    const Teuchos::ArrayView<const GlobalOrdinal> localGIDs_view(&localGIDs[0], localGIDs.size());

    Teuchos::RCP<const Map> localGIDsMap = MapFactory::Build(
        map->lib(),
        Teuchos::OrdinalTraits<int>::invalid(),
        localGIDs_view,
        0, ptv->getMap()->getComm());  // use correct communicator here!

    // store rebalanced partial map using the same name and generating factory as the original map
    // in the level class
    level.Set(mapName, localGIDsMap, mapFactory.get());
  }
}  // Build()

}  // end namespace MueLu

#endif /* MUELU_REBALANCEMAPFACTORY_DEF_HPP_ */
