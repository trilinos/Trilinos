// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_BLOCKEDCOARSEMAPFACTORY_DEF_HPP_
#define MUELU_BLOCKEDCOARSEMAPFACTORY_DEF_HPP_

#include <Xpetra_MultiVector.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_BlockedCoarseMapFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_Factory.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> BlockedCoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase>>("Aggregates", Teuchos::null, "Generating factory for aggregates.");
  validParamList->set<RCP<const FactoryBase>>("Nullspace", Teuchos::null, "Generating factory for null space.");
  validParamList->set<RCP<const FactoryBase>>("CoarseMap", Teuchos::null, "Generating factory of previous coarse map. (must be set by user!).");

  // do we need this?
  validParamList->set<std::string>("Striding info", "{}", "Striding information");
  validParamList->set<LocalOrdinal>("Strided block id", -1, "Strided block id");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  this->Input(currentLevel, "Aggregates");
  this->Input(currentLevel, "Nullspace");

  // Get CoarseMap from previously defined block
  RCP<const FactoryBase> prevCoarseMapFact = this->GetFactory("CoarseMap");
  TEUCHOS_TEST_FOR_EXCEPTION(prevCoarseMapFact == Teuchos::null, Exceptions::RuntimeError, "MueLu::BlockedCoarseMapFactory::getDomainMapOffset: user did not specify CoarseMap of previous block. Do not forget to set the CoarseMap factory.");
  currentLevel.DeclareInput("CoarseMap", prevCoarseMapFact.get(), this);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  GlobalOrdinal domainGIDOffset = GetDomainGIDOffset(currentLevel);
  CoarseMapFactory::BuildCoarseMap(currentLevel, domainGIDOffset);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal BlockedCoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetDomainGIDOffset(
    Level &currentLevel) const {
  RCP<const FactoryBase> prevCoarseMapFact = this->GetFactory("CoarseMap");
  RCP<const Map> subPDomainMap             = currentLevel.Get<RCP<const Map>>("CoarseMap", prevCoarseMapFact.get());
  GlobalOrdinal maxGlobalIndex             = subPDomainMap->getMaxAllGlobalIndex();

  return maxGlobalIndex + Teuchos::ScalarTraits<GlobalOrdinal>::one();
}

}  // namespace MueLu

#endif /* MUELU_BLOCKEDCOARSEMAPFACTORY_DEF_HPP_ */
