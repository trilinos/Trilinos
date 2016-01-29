/*
 * MueLu_TopRAPFactory_decl.hpp
 *
 *  Created on: Jan 25, 2016
 *      Author: tawiesn
 */

#ifndef PACKAGES_MUELU_SRC_MUECENTRAL_MUELU_TOPRAPFACTORY_DECL_HPP_
#define PACKAGES_MUELU_SRC_MUECENTRAL_MUELU_TOPRAPFACTORY_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"

//#include "MueLu_FactoryManager_fwd.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
//#include "MueLu_Hierarchy_fwd.hpp"
//#include "MueLu_HierarchyManager_fwd.hpp"

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class TopRAPFactory : public TwoLevelFactoryBase {
#undef MUELU_TOPRAPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    TopRAPFactory(RCP<const FactoryManagerBase> parentFactoryManager);
    TopRAPFactory(RCP<const FactoryManagerBase> parentFactoryManagerFine, RCP<const FactoryManagerBase> parentFactoryManagerCoarse);

    virtual ~TopRAPFactory();

    void DeclareInput(Level & fineLevel, Level & coarseLevel) const;

    void Build(Level & fineLevel, Level & coarseLevel) const;

  private:
    RCP<const FactoryBase> PFact_;
    RCP<const FactoryBase> RFact_;
    RCP<const FactoryBase> AcFact_;
  };

}

#define MUELU_TOPRAPFACTORY_SHORT
#endif /* PACKAGES_MUELU_SRC_MUECENTRAL_MUELU_TOPRAPFACTORY_DECL_HPP_ */
