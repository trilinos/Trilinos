#ifndef MUELU_FACTORYFACTORY_DECL_HPP
#define MUELU_FACTORYFACTORY_DECL_HPP

#include <string>

#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_Array.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FactoryFactory_fwd.hpp"

#include "MueLu_HierarchyFactory.hpp"

#include "MueLu_FactoryBase.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_Hierarchy_fwd.hpp"
#include "MueLu_FactoryManagerBase_fwd.hpp"
#include "MueLu_FactoryManager_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp" //TMP
#include "MueLu_TrilinosSmoother_fwd.hpp" //TMP
#include "MueLu_TrilinosSmoother.hpp" //TMP
#include "MueLu_SmootherFactory_fwd.hpp" //TMP
#include "MueLu_SmootherFactory.hpp" //TMP
#include "MueLu_TentativePFactory.hpp" //TMP
#include "MueLu_Exceptions.hpp" //TMP

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void, LocalOrdinal, Node>::SparseOps>
  class FactoryFactory {
#undef MUELU_FACTORYFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

    typedef std::map<std::string, RCP<const FactoryBase> > FactoryMap; // TODO: remove

  public:

    // Parameter List Parsing:
    // ---------
    //     <Parameter name="smootherFact0" type="string" value="TrilinosSmoother"/>
    //
    // or:
    //
    //     <ParameterList name="smootherFact1">
    //       <Parameter name="type" type="string" value="TrilinosSmoother"/>
    //       ...
    //     </ParameterList>
    //
    RCP<FactoryBase> BuildFactory(const Teuchos::ParameterEntry & param, const FactoryMap & factoryMapIn) const {
      // Find factory
      std::string type;
      Teuchos::ParameterList paramList;    
      if (!param.isList()) {
        type = Teuchos::getValue<std::string>(param);
      } else {
        paramList = Teuchos::getValue<Teuchos::ParameterList>(param);
        type = paramList.get<std::string>("type");
      } 
    
      // TODO: see how Teko handles this.
      if (type == "TrilinosSmoother") {
        return BuildTrilinosSmoother(paramList);
      }
      if (type == "TentativePFactory") {
        return BuildTentativePFactory(paramList);
      }

      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::HierarchyFactory:: unkown factory name.");

      return Teuchos::null;
    }
  
    //
    //
    //

    // FOLLOWING FUNCTIONS SHOULD LIVE WITH THE CORRESPONDING CLASS

    //
    //
    //

    //! TentativePFactory
    RCP<FactoryBase> BuildTentativePFactory(const Teuchos::ParameterList & paramList) const {
      //TODO: handle parameters   
      return rcp(new TentativePFactory());
    }

    //! TrilinosSmoother
    // Parameter List Parsing:
    //     <ParameterList name="smootherFact1">
    //       <Parameter name="type" type="string" value="TrilinosSmoother"/>
    //       <Parameter name="verbose" type="string" value="Warnings"/>
    //       <Parameter name="precond" type="string" value="RELAXATION"/>
    //       <ParameterList name="ParameterList">
    //       ...
    //       </ParameterList>
    //     </ParameterList>
    RCP<FactoryBase> BuildTrilinosSmoother(const Teuchos::ParameterList & paramList) const {
      if (paramList.begin() == paramList.end())
        return rcp(new SmootherFactory(rcp(new TrilinosSmoother())));
    
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("type") != "TrilinosSmoother", Exceptions::RuntimeError, "");

      std::string precond;           if(paramList.isParameter("precond"))       precond = paramList.get<std::string>("precond");
      std::string verbose;           if(paramList.isParameter("verbose"))       verbose = paramList.get<std::string>("verbose");
      Teuchos::ParameterList params; if(paramList.isParameter("ParameterList")) params  = paramList.get<Teuchos::ParameterList>("ParameterList");

      return rcp(new SmootherFactory(rcp(new TrilinosSmoother(precond, params))));
    }
    
  }; // class

} // namespace MueLu

#endif // MUELU_FACTORYFACTORY_DECL_HPP


