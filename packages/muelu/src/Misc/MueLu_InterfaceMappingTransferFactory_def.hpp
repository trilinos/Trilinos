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
    fineLevel.DeclareInput("CoarseLagr2Dof", NoFactory::get(), this);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceMappingTransferFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const
{
    Monitor m(*this, "Interface Mapping transfer factory");
    RCP<std::map<GlobalOrdinal, GlobalOrdinal>> coarseLagr2Dof = fineLevel.Get<RCP<std::map<GlobalOrdinal, GlobalOrdinal>>>("CoarseLagr2Dof", NoFactory::get());
    coarseLevel.Set<RCP<std::map<GlobalOrdinal, GlobalOrdinal>>>("Lagr2Dof", coarseLagr2Dof, NoFactory::get());
}

} // namespace MueLu

#endif /* MUELU_INTERFACEMAPPINGTRANSFERFACTORY_DEF_HPP_ */