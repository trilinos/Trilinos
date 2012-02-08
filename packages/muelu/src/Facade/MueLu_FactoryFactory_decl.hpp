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
    //       <Parameter name="class" type="string" value="TrilinosSmoother"/>
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
        type = paramList.get<std::string>("class");
      } 
    
      // TODO: see how Teko handles this.
      if (type == "TentativePFactory") {
        return BuildTentativePFactory(paramList);
      }
      if (type == "SaPFactory") {
        return BuildSaPFactory(paramList);
      }
      if (type == "RAPFactory") {
        return BuildRAPFactory(paramList);
      }
      if (type == "UCAggregationFactory") {
        return BuildUCAggregationFactory(paramList);
      }
      if (type == "TrilinosSmoother") {
        return BuildTrilinosSmoother(paramList);
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
      return rcp(new TentativePFactory());
    }
    
    //! SaPFactory
    RCP<FactoryBase> BuildSaPFactory(const Teuchos::ParameterList & paramList) const {
      if (paramList.begin() == paramList.end()) // short-circuit. Use default parameters of constructor
        return rcp(new SaPFactory());

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("class") != "SaPFactory", Exceptions::RuntimeError, "");      

      RCP<SaPFactory> f = rcp(new SaPFactory());
      if (paramList.isParameter("DampingFactory")) f->SetDampingFactor(paramList.get<Scalar>("DampingFactory"));

      return f;
    }

    //! RaPFactory
    RCP<FactoryBase> BuildRAPFactory(const Teuchos::ParameterList & paramList) const {
      if (paramList.begin() == paramList.end()) // short-circuit. Use default parameters of constructor
        return rcp(new RAPFactory());

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("class") != "RAPFactory", Exceptions::RuntimeError, "");      

      return rcp(new RAPFactory());
    }

    //! UCAggregationFactory
    RCP<FactoryBase> BuildUCAggregationFactory(const Teuchos::ParameterList & paramList) const {
      if (paramList.begin() == paramList.end()) // short-circuit. Use default parameters of constructor
        return rcp(new UCAggregationFactory());

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("class") != "UCAggregationFactory", Exceptions::RuntimeError, "");      

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
    //       <Parameter name="class" type="string" value="TrilinosSmoother"/>
    //       <Parameter name="verbose" type="string" value="Warnings"/>
    //       <Parameter name="precond" type="string" value="RELAXATION"/>
    //       <ParameterList name="ParameterList">
    //       ...
    //       </ParameterList>
    //     </ParameterList>
    RCP<FactoryBase> BuildTrilinosSmoother(const Teuchos::ParameterList & paramList) const {
      if (paramList.begin() == paramList.end())
        return rcp(new SmootherFactory(rcp(new TrilinosSmoother())));

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("class") != "TrilinosSmoother", Exceptions::RuntimeError, "");

      std::string precond;           if(paramList.isParameter("precond"))       precond = paramList.get<std::string>("precond");
      std::string verbose;           if(paramList.isParameter("verbose"))       verbose = paramList.get<std::string>("verbose");
      Teuchos::ParameterList params; if(paramList.isParameter("ParameterList")) params  = paramList.get<Teuchos::ParameterList>("ParameterList");

      return rcp(new SmootherFactory(rcp(new TrilinosSmoother(precond, params))));
    }
    
  }; // class

} // namespace MueLu

#define MUELU_FACTORYFACTORY_SHORT
#endif // MUELU_FACTORYFACTORY_DECL_HPP

// TODO: handle factory parameters
// TODO: parameter validator
