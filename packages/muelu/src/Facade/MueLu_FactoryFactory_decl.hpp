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

#include "MueLu_RAPFactory.hpp" //TMP
#include "MueLu_SaPFactory.hpp" //TMP
#include "MueLu_TrilinosSmoother.hpp" //TMP
#include "MueLu_SmootherFactory.hpp" //TMP
#include "MueLu_TentativePFactory.hpp" //TMP
#include "MueLu_UCAggregationFactory.hpp" //TMP
#include "MueLu_DirectSolver.hpp" //TMP
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
    //       <Parameter name="factory" type="string" value="TrilinosSmoother"/>
    //       ...
    //     </ParameterList>
    //
    RCP<const FactoryBase> BuildFactory(const Teuchos::ParameterEntry & param, const FactoryMap & factoryMapIn) const {

      // Find factory
      std::string factoryName;
      Teuchos::ParameterList paramList;
      if (!param.isList()) {
        factoryName = Teuchos::getValue<std::string>(param);
      } else {
        paramList = Teuchos::getValue<Teuchos::ParameterList>(param);
        factoryName = paramList.get<std::string>("factory");
      } 

      // TODO: see how Teko handles this.
      if (factoryName == "TentativePFactory") {
        return BuildTentativePFactory(paramList);
      }
      if (factoryName == "SaPFactory") {
        return BuildSaPFactory(paramList);
      }
      if (factoryName == "RAPFactory") {
        return BuildRAPFactory(paramList);
      }
      if (factoryName == "UCAggregationFactory") {
        return BuildUCAggregationFactory(paramList);
      }
      if (factoryName == "TrilinosSmoother") {
        return BuildTrilinosSmoother(paramList);
      }
      if (factoryName == "DirectSolver") {
        return BuildDirectSolver(paramList);
      }

      // Use a user defined factories (in <Factories> node)
      if (factoryMapIn.find(factoryName) != factoryMapIn.end()) {
        TEUCHOS_TEST_FOR_EXCEPTION(paramList.begin() != paramList.end(), Exceptions::RuntimeError, "MueLu::FactoryFactory:: parameters of factories defined in <Factories> node cannot be redefined.");
        return factoryMapIn.find(factoryName)->second;
      }

      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory:: unkown factory name : " << factoryName);

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
      return rcp(new TentativePFactory());
    }
    
    //! SaPFactory
    RCP<FactoryBase> BuildSaPFactory(const Teuchos::ParameterList & paramList) const {
      if (paramList.begin() == paramList.end()) // short-circuit. Use default parameters of constructor
        return rcp(new SaPFactory());

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "SaPFactory", Exceptions::RuntimeError, "");      

      RCP<SaPFactory> f = rcp(new SaPFactory());
      if (paramList.isParameter("DampingFactory")) f->SetDampingFactor(paramList.get<Scalar>("DampingFactory"));

      return f;
    }

    //! RaPFactory
    RCP<FactoryBase> BuildRAPFactory(const Teuchos::ParameterList & paramList) const {
      if (paramList.begin() == paramList.end()) // short-circuit. Use default parameters of constructor
        return rcp(new RAPFactory());

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "RAPFactory", Exceptions::RuntimeError, "");      

      return rcp(new RAPFactory());
    }

    //! UCAggregationFactory
    RCP<FactoryBase> BuildUCAggregationFactory(const Teuchos::ParameterList & paramList) const {
      if (paramList.begin() == paramList.end()) // short-circuit. Use default parameters of constructor
        return rcp(new UCAggregationFactory());

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "UCAggregationFactory", Exceptions::RuntimeError, "");      

      RCP<UCAggregationFactory> f = rcp(new UCAggregationFactory());

      if(paramList.isParameter("Ordering")) {
        std::string orderingStr = paramList.get<std::string>("Ordering");
        Ordering ordering;
        if (orderingStr == "Natural")
          ordering = NATURAL;
        else if (orderingStr == "Random")
          ordering = RANDOM;
        else if (orderingStr == "Graph")
          ordering = GRAPH;
        else TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory::BuildUCAggregationFactory()::Unknown Ordering type");

        f->SetOrdering(ordering);
      }

      if(paramList.isParameter("MaxNeighAlreadySelected")) {
        f->SetMaxNeighAlreadySelected(paramList.get<int>("MaxNeighAlreadySelected"));
      }
      
      if(paramList.isParameter("Phase3AggCreation")) {
        f->SetPhase3AggCreation(paramList.get<double>("Phase3AggCreation"));
      }

      if(paramList.isParameter("MinNodesPerAggregate")) {
        f->SetMinNodesPerAggregate(paramList.get<int>("MinNodesPerAggregate"));
      }

      return f;
    }

    //! TrilinosSmoother
    // Parameter List Parsing:
    //     <ParameterList name="smootherFact1">
    //       <Parameter name="factory" type="string" value="TrilinosSmoother"/>
    //       <Parameter name="verbose" type="string" value="Warnings"/>
    //       <Parameter name="type" type="string" value="RELAXATION"/>
    //       <ParameterList name="ParameterList">
    //       ...
    //       </ParameterList>
    //     </ParameterList>
    RCP<FactoryBase> BuildTrilinosSmoother(const Teuchos::ParameterList & paramList) const {
      if (paramList.begin() == paramList.end())
        return rcp(new SmootherFactory(rcp(new TrilinosSmoother())));

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "TrilinosSmoother", Exceptions::RuntimeError, "");

      std::string type;               if(paramList.isParameter("type"))         type = paramList.get<std::string>("type");
      // std::string verbose;         if(paramList.isParameter("verbose"))      verbose = paramList.get<std::string>("verbose");
      Teuchos::ParameterList params; if(paramList.isParameter("ParameterList")) params  = paramList.get<Teuchos::ParameterList>("ParameterList");

      return rcp(new SmootherFactory(rcp(new TrilinosSmoother(type, params))));
    }
    
    RCP<FactoryBase> BuildDirectSolver(const Teuchos::ParameterList & paramList) const {
      if (paramList.begin() == paramList.end())
        return rcp(new SmootherFactory(rcp(new DirectSolver())));
      
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "DirectSolver", Exceptions::RuntimeError, "");

      std::string type;              if(paramList.isParameter("type"))          type = paramList.get<std::string>("type");
      // std::string verbose;        if(paramList.isParameter("verbose"))       verbose = paramList.get<std::string>("verbose");
      Teuchos::ParameterList params; if(paramList.isParameter("ParameterList")) params  = paramList.get<Teuchos::ParameterList>("ParameterList");
      
      return rcp(new SmootherFactory(rcp(new DirectSolver(type, params))));
    }

  }; // class

} // namespace MueLu

#define MUELU_FACTORYFACTORY_SHORT
#endif // MUELU_FACTORYFACTORY_DECL_HPP

// TODO: handle factory parameters
// TODO: parameter validator
