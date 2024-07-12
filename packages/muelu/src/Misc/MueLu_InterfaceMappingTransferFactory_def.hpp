// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_INTERFACEMAPPINGTRANSFERFACTORY_DEF_HPP_
#define MUELU_INTERFACEMAPPINGTRANSFERFACTORY_DEF_HPP_

#include "MueLu_InterfaceMappingTransferFactory_decl.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> InterfaceMappingTransferFactory<LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  validParamList->set<RCP<const FactoryBase>>("CoarseDualNodeID2PrimalNodeID", Teuchos::null, "Generating factory of the CoarseDualNodeID2PrimalNodeID map");
  return validParamList;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceMappingTransferFactory<LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
  Input(fineLevel, "CoarseDualNodeID2PrimalNodeID");
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceMappingTransferFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
  Monitor m(*this, "Interface Mapping transfer factory");

  RCP<std::map<LocalOrdinal, LocalOrdinal>> coarseLagr2Dof = Get<RCP<std::map<LocalOrdinal, LocalOrdinal>>>(fineLevel, "CoarseDualNodeID2PrimalNodeID");
  Set(coarseLevel, "DualNodeID2PrimalNodeID", coarseLagr2Dof);
}

}  // namespace MueLu

#endif /* MUELU_INTERFACEMAPPINGTRANSFERFACTORY_DEF_HPP_ */
