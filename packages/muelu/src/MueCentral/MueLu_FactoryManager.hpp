#ifndef MUELU_FACTORYMANAGER_HPP
#define MUELU_FACTORYMANAGER_HPP

#include <map>

#include <Xpetra_Operator.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_NoFactory.hpp" //TODO: remove
#include "MueLu_SmootherFactoryBase.hpp"
#include "MueLu_SmootherBase.hpp"

#include "MueLu_FactoryManagerBase.hpp"

// Headers for factories used by default:
#include "MueLu_SaPFactory.hpp"
//#include "MueLu_TentativePFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_ReUseFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class FactoryManager : public FactoryManagerBase {
#include "MueLu_UseShortNames.hpp"
      
  public:

    //@{

    //!
    FactoryManager(const RCP<const FactoryBase> PFact = Teuchos::null, const RCP<const FactoryBase> RFact = Teuchos::null, const RCP<const FactoryBase> AcFact = Teuchos::null)
      : PFact_(PFact), RFact_(RFact), AcFact_(AcFact)
    { }
    
    //! Destructor.
    virtual ~FactoryManager() { }

    //@}

    //@{ Get/Set functions.

    //! Get
    // Return ref because user also give ref to the Hierarchy.
    // Factory freed at the end of FillHierarchy() //->TODO
    virtual const FactoryBase & GetFactory2(const std::string & varName) const {
      // TODO: try/catch + better exception msg if not found
      return *factoryTable_.get(varName);
    }

    void SetFactory(const std::string & varName, const RCP<const FactoryBase> & factory) const { //TODO: remove const, remame SetFactory()
      // TODO: if (varName already exist) ...
      factoryTable_.put(varName, factory);
    }

    bool IsAvailable(const std::string & varName) const {
      return factoryTable_.containsKey(varName);
    }

    //@}

    //@{

    void SetSmootherFactory(const RCP<const FactoryBase> & smootherFact) {
      smootherFact_ = smootherFact;
    }
  
    void SetCoarsestSolverFactory(const RCP<const FactoryBase> & coarsestSolver) {
      coarsestSolverFact_ = coarsestSolver;
    }

    RCP<const FactoryBase> GetPFact()  const { return PFact_; }
    RCP<const FactoryBase> GetRFact()  const { return RFact_; }
    RCP<const FactoryBase> GetAcFact() const { return AcFact_; }

    RCP<const FactoryBase> GetSmootherFactory() const { return smootherFact_; }
    RCP<const FactoryBase> GetCoarsestSolverFactory() const { return coarsestSolverFact_; }

    //@}

    //@{ Get/Set functions.

    virtual const FactoryBase & GetFactory(const std::string & varName) const {
      if (! IsAvailable(varName)) {

        if (varName == "A")            return *NoFactory::get();
        if (varName == "P")            return *NoFactory::get();
        if (varName == "R")            return *NoFactory::get();

    	if (varName == "Nullspace")    return SetAndReturnDefaultFactory(varName, rcp(new NullspaceFactory()));
        if (varName == "Graph")        return SetAndReturnDefaultFactory(varName, rcp(new CoalesceDropFactory()));
        if (varName == "Aggregates")   return SetAndReturnDefaultFactory(varName, rcp(new UCAggregationFactory()));

        TEST_FOR_EXCEPTION(1, MueLu::Exceptions::RuntimeError, "FactoryManager::GetFactory(): No default factory available for building '"+varName+"'.");
      }

      return GetFactory2(varName);
    }

    //@}
    
  private:

    //! helper
    const FactoryBase & SetAndReturnDefaultFactory(const std::string & varName, const RCP<FactoryBase> factory) const {

      GetOStream(Warnings0, 0)  << "Warning: No factory have been specified for building '" << varName << "'." << std::endl;
      GetOStream(Warnings00, 0) << "         using default factory: ";
      { Teuchos::OSTab tab(getOStream(), 8); factory->describe(GetOStream(Warnings00), getVerbLevel()); }

      SetFactory(varName, factory);
      return GetFactory(varName); //TODO: replace by: return factory;
    }


  private:

    mutable Teuchos::Hashtable<std::string, RCP<const FactoryBase> > factoryTable_; //TODO: use std lib hashtable instead (Teuchos::Hashtable is deprecated)

    RCP<const FactoryBase> PFact_;
    RCP<const FactoryBase> RFact_;
    RCP<const FactoryBase> AcFact_;

    RCP<const FactoryBase> smootherFact_;
    RCP<const FactoryBase> coarsestSolverFact_;
    
  }; // class

} //namespace MueLu

#define MUELU_FACTORYMANAGER_SHORT
#endif // ifndef MUELU_FACTORYMANAGER_HPP

//TODO: clean up factory created at the end
//TODO: factoryTable_ must be cleaned at the end of hierarchy Populate() (because Hierarchy is not holding any factories after construction)
