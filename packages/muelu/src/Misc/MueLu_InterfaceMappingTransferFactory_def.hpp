#ifndef MUELU_INTERFACEMAPPINGTRANSFERFACTORY_DEF_HPP_
#define MUELU_INTERFACEMAPPINGTRANSFERFACTORY_DEF_HPP_

#include "MueLu_InterfaceMappingTransferFactory_decl.hpp"

namespace MueLu
{

template <class LocalOrdinal, class GlobalOrdinal, class Node>
InterfaceMappingTransferFactory<LocalOrdinal, GlobalOrdinal, Node>::InterfaceMappingTransferFactory()
{
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> InterfaceMappingTransferFactory<LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const
{
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    return validParamList;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceMappingTransferFactory<LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const
{
    fineLevel.DeclareInput("CoarseDualNodeID2PrimalNodeID", NoFactory::get());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceMappingTransferFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const
{
    Monitor m(*this, "Interface Mapping transfer factory");
    RCP<std::map<LocalOrdinal, LocalOrdinal>> coarseLagr2Dof = fineLevel.Get<RCP<std::map<LocalOrdinal, LocalOrdinal>>>("CoarseDualNodeID2PrimalNodeID", NoFactory::get());
    coarseLevel.Set<RCP<std::map<LocalOrdinal, LocalOrdinal>>>("DualNodeID2PrimalNodeID", coarseLagr2Dof, NoFactory::get());
}

} // namespace MueLu

#endif /* MUELU_INTERFACEMAPPINGTRANSFERFACTORY_DEF_HPP_ */
