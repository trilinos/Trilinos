// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REBALANCEACFACTORY_DEF_HPP
#define MUELU_REBALANCEACFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include "MueLu_RebalanceAcFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_RAPFactory.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RebalanceAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RebalanceAcFactory() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RebalanceAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~RebalanceAcFactory() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> RebalanceAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("repartition: use subcommunicators");
  SET_VALID_ENTRY("repartition: use subcommunicators in place");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A for rebalancing");
  validParamList->set<RCP<const FactoryBase> >("Importer", Teuchos::null, "Generating factory of the importer");
  validParamList->set<RCP<const FactoryBase> >("InPlaceMap", Teuchos::null, "Generating factory of the InPlaceMap");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level & /* fineLevel */, Level &coarseLevel) const {
  Input(coarseLevel, "A");
  const Teuchos::ParameterList &pL = GetParameterList();
  if (pL.isParameter("repartition: use subcommunicators in place") && pL.get<bool>("repartition: use subcommunicators in place") == true) {
    Input(coarseLevel, "InPlaceMap");
  } else
    Input(coarseLevel, "Importer");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level & /* fineLevel */, Level &coarseLevel) const {
  FactoryMonitor m(*this, "Computing Ac", coarseLevel);

  const Teuchos::ParameterList &pL = GetParameterList();
  RCP<Matrix> originalAc           = Get<RCP<Matrix> >(coarseLevel, "A");

  // This is a short-circuit for if we want to leave A where it is, but restrict its communicator
  // to something corresponding to a given map.  Maxwell1 is the prime customer for this.
  bool inPlace = pL.get<bool>("repartition: use subcommunicators in place");
  if (inPlace) {
    SubFactoryMonitor subM(*this, "Rebalancing existing Ac in-place", coarseLevel);
    RCP<const Map> newMap = Get<RCP<const Map> >(coarseLevel, "InPlaceMap");

    originalAc->removeEmptyProcessesInPlace(newMap);

    // The "in place" still leaves a dummy matrix here.  That needs to go
    if (newMap.is_null()) originalAc = Teuchos::null;

    Set(coarseLevel, "A", originalAc);
    return;
  }

  RCP<const Import> rebalanceImporter = Get<RCP<const Import> >(coarseLevel, "Importer");

  if (rebalanceImporter != Teuchos::null) {
    RCP<Matrix> rebalancedAc;
    {
      SubFactoryMonitor subM(*this, "Rebalancing existing Ac", coarseLevel);
      RCP<const Map> targetMap = rebalanceImporter->getTargetMap();

      ParameterList XpetraList;
      if (pL.get<bool>("repartition: use subcommunicators") == true) {
        GetOStream(Runtime0) << "Replacing maps with a subcommunicator" << std::endl;
        XpetraList.set("Restrict Communicator", true);
      }
      // NOTE: If the communicator is restricted away, Build returns Teuchos::null.
      XpetraList.set("Timer Label", "MueLu::RebalanceAc-" + Teuchos::toString(coarseLevel.GetLevelID()));
      {
        SubFactoryMonitor subM2(*this, "Rebalancing existing Ac: MatrixFactory::Build", coarseLevel);
        rebalancedAc = MatrixFactory::Build(originalAc, *rebalanceImporter, *rebalanceImporter, targetMap, targetMap, rcp(&XpetraList, false));
      }

      if (!rebalancedAc.is_null()) {
        rebalancedAc->SetFixedBlockSize(originalAc->GetFixedBlockSize());
        std::ostringstream oss;
        oss << "A_" << coarseLevel.GetLevelID();
        rebalancedAc->setObjectLabel(oss.str());
      }
      Set(coarseLevel, "A", rebalancedAc);
    }
    if (!rebalancedAc.is_null() && IsPrint(Statistics2)) {
      int oldRank = SetProcRankVerbose(rebalancedAc->getRowMap()->getComm()->getRank());

      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      params->set("printCommInfo", true);
      GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*rebalancedAc, "Ac (rebalanced)", params);

      SetProcRankVerbose(oldRank);
    }

  } else {
    // Ac already built by the load balancing process and no load balancing needed
    GetOStream(Runtime1) << "No rebalancing" << std::endl;
    Set(coarseLevel, "A", originalAc);
  }

  if (rebalanceFacts_.begin() != rebalanceFacts_.end()) {
    SubFactoryMonitor m2(*this, "Rebalance additional data", coarseLevel);

    // call Build of all user-given transfer factories
    for (std::vector<RCP<const FactoryBase> >::const_iterator it = rebalanceFacts_.begin(); it != rebalanceFacts_.end(); ++it) {
      GetOStream(Runtime0) << "RebalanceAc: call rebalance factory " << (*it).get() << ": " << (*it)->description() << std::endl;
      (*it)->CallBuild(coarseLevel);
    }
  }
}  // Build()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddRebalanceFactory(const RCP<const FactoryBase> &factory) {
  /*TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast,
                             "MueLu::RAPFactory::AddTransferFactory: Transfer factory is not derived from TwoLevelFactoryBase. "
                             "This is very strange. (Note: you can remove this exception if there's a good reason for)");
  TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_, Exceptions::RuntimeError, "MueLu::RAPFactory::AddTransferFactory: Factory is being added after we have already declared input");*/
  rebalanceFacts_.push_back(factory);
}  // AddRebalanceFactory()

}  // namespace MueLu

#endif  // MUELU_REBALANCEACFACTORY_DEF_HPP
