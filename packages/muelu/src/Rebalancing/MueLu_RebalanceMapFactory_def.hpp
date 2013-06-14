// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2013 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_REBALANCEMAPFACTORY_DEF_HPP_
#define MUELU_REBALANCEMAPFACTORY_DEF_HPP_

#include "MueLu_RebalanceMapFactory_decl.hpp"

#include <Teuchos_Utils.hpp>
#include <Teuchos_DefaultMpiComm.hpp> //TODO: fwd decl.
#include <Teuchos_OpaqueWrapper.hpp>  //TODO: fwd decl.

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

 /*template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
 RebalanceMapFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RebalanceMapFactory(std::string mapName, Teuchos::RCP<const FactoryBase> mapFact)
 : mapName_(mapName), mapFact_(mapFact)
 {}*/

 template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
 RCP<const ParameterList> RebalanceMapFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    //validParamList->set< RCP<const FactoryBase> >("A",                    Teuchos::null, "Factory of the matrix A");
    validParamList->set< std::string >           ("Map name"   , "", "Name of map to rebalanced.");
    validParamList->set< RCP<const FactoryBase> >("Map factory", MueLu::NoFactory::getRCP(), "Generating factory of map to be rebalanced.");

    validParamList->set< RCP<const FactoryBase> >("Importer",             Teuchos::null, "Factory of the importer object used for the rebalancing");

    return validParamList;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceMapFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & currentLevel) const {
    //Input(currentLevel, "A");
    const Teuchos::ParameterList & pL = GetParameterList();
    std::string mapName                        = pL.get<std::string> ("Map name");
    Teuchos::RCP<const FactoryBase> mapFactory = GetFactory          ("Map factory");
    currentLevel.DeclareInput(mapName,mapFactory.get(),this);

    Input(currentLevel, "Importer");
  } //DeclareInput()

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceMapFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &level) const {
    FactoryMonitor m(*this, "Build", level);

    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

    const Teuchos::ParameterList & pL = GetParameterList();
    std::string mapName = pL.get<std::string> ("Map name");
    Teuchos::RCP<const FactoryBase> mapFactory = GetFactory ("Map factory");

    RCP<const Import> rebalanceImporter = Get<RCP<const Import> >(level, "Importer");

    // input map (not rebalanced)
    RCP<const Map> map = level.Get< RCP<const Map> >(mapName,mapFactory.get());

    // create vector based on input map
    RCP<Vector> v = VectorFactory::Build(map);
    v->putScalar(1.0);

    //v->describe(*fos, Teuchos::VERB_EXTREME);

    // source map unique
    RCP<const Import> blowUpImporter = ImportFactory::Build(map, rebalanceImporter->getSourceMap());
    RCP<Vector> pv = VectorFactory::Build(rebalanceImporter->getSourceMap());
    pv->doImport(*v,*blowUpImporter,Xpetra::INSERT);

    //pv->describe(*fos, Teuchos::VERB_EXTREME);

    // do rebalancing
    RCP<Vector> ptv = VectorFactory::Build(rebalanceImporter->getTargetMap());
    ptv->doImport(*pv,*rebalanceImporter,Xpetra::INSERT);

    // reconstruct map




    Teuchos::ArrayRCP< const Scalar > ptvData = ptv->getData(0);
    std::vector<GlobalOrdinal> localGIDs;  // vector with GIDs that are stored on current proc

    for (size_t k = 0; k < ptv->getLocalLength(); k++) {
      if(ptvData[k] == 1.0) {
        localGIDs.push_back(ptv->getMap()->getGlobalElement(k));
      }
    }

    const Teuchos::ArrayView<const LocalOrdinal> localGIDs_view(&localGIDs[0],localGIDs.size());

    Teuchos::RCP<const Map> localGIDsMap = MapFactory::Build(
        map->lib(),
        Teuchos::OrdinalTraits<int>::invalid(),
        localGIDs_view,
        0, map->getComm());

    //localGIDsMap->describe(*fos, Teuchos::VERB_EXTREME);


    level.Set(mapName, localGIDsMap, mapFactory.get());
  } //Build()

} // end namespace MueLu

#endif /* MUELU_REBALANCEMAPFACTORY_DEF_HPP_ */
