/*
 * MueLu_TopSmootherFactory_decl.hpp
 *
 *  Created on: Jan 25, 2016
 *      Author: tawiesn
 */

#ifndef PACKAGES_MUELU_SRC_MUECENTRAL_MUELU_TOPSMOOTHERFACTORY_DECL_HPP_
#define PACKAGES_MUELU_SRC_MUECENTRAL_MUELU_TOPSMOOTHERFACTORY_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_FactoryManagerBase_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_SmootherBase_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_SmootherPrototype_fwd.hpp"
//#include "MueLu_TwoLevelFactoryBase.hpp"
//#include "MueLu_Hierarchy_fwd.hpp"
//#include "MueLu_HierarchyManager_fwd.hpp"

namespace MueLu {

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class TopSmootherFactory : public SingleLevelFactoryBase {  // TODO: inherit from SmootherFactoryBase ?
#undef MUELU_TOPSMOOTHERFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  TopSmootherFactory(RCP<const FactoryManagerBase> parentFactoryManager, const std::string& varName);

  virtual ~TopSmootherFactory();

  void DeclareInput(Level& level) const;

  void Build(Level& level) const;

 private:
  RCP<const FactoryBase> preSmootherFact_;
  RCP<const FactoryBase> postSmootherFact_;
};

}  // namespace MueLu

#define MUELU_TOPSMOOTHERFACTORY_SHORT
#endif /* PACKAGES_MUELU_SRC_MUECENTRAL_MUELU_TOPSMOOTHERFACTORY_DECL_HPP_ */
