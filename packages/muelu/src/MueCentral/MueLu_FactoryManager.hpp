#ifndef MUELU_FACTORYMANAGER_HPP
#define MUELU_FACTORYMANAGER_HPP

#include <map>

#include <Xpetra_Operator.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_SmootherFactoryBase.hpp" //TODO:remove
#include "MueLu_SmootherBase.hpp"

#include "MueLu_FactoryManagerBase.hpp"

// Headers for factories used by default:
#include "MueLu_NoFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
 
namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class FactoryManager : public FactoryManagerBase {
#include "MueLu_UseShortNames.hpp"
      
  public:

    //@{

    //!
    FactoryManager(const RCP<const FactoryBase> PFact = Teuchos::null, const RCP<const FactoryBase> RFact = Teuchos::null, const RCP<const FactoryBase> AcFact = Teuchos::null) { 
      if (PFact  != Teuchos::null) SetFactory("P", PFact);
      if (RFact  != Teuchos::null) SetFactory("R", RFact);
      if (AcFact != Teuchos::null) SetFactory("A", AcFact);
    }
    
    //! Destructor.
    virtual ~FactoryManager() { }

    //@}

    //@{ Get/Set functions.

    //! Set Factory
    void SetFactory(const std::string & varName, const RCP<const FactoryBase> & factory) {
      if (IsAvailable(varName, factoryTable_)) // TODO: too much warnings (for smoothers)
        GetOStream(Warnings1, 0) << "Warning: FactoryManager::SetFactory(): Changing an already defined factory for " << varName << std::endl;

      factoryTable_[varName] = factory;
    }

    //! Get Factory
    const RCP<const FactoryBase> & GetFactory(const std::string & varName) const {
      if (FactoryManager::IsAvailable(varName, factoryTable_))
        return factoryTable_.find(varName)->second; // == factoryTable_[varName] (operator std::map[] is not const)
      else 
        return GetDefaultFactory(varName);
    }

    //!
    const RCP<const FactoryBase> & GetDefaultFactory(const std::string & varName) const {
      if (IsAvailable(varName, defaultFactoryTable_)) {

        return defaultFactoryTable_[varName];

      }        else {
          
        //if (varName == "A")           return SetAndReturnDefaultFactory(varName, rcp(new RAPFactory())); will need some work
        if (varName == "A")             return SetAndReturnDefaultFactory(varName, NoFactory::getRCP());
        if (varName == "P")             return SetAndReturnDefaultFactory(varName, rcp(new TentativePFactory()));
        if (varName == "R")             return SetAndReturnDefaultFactory(varName, rcp(new TransPFactory()));

        if (varName == "Nullspace")     return SetAndReturnDefaultFactory(varName, rcp(new NullspaceFactory()));
        if (varName == "Graph")         return SetAndReturnDefaultFactory(varName, rcp(new CoalesceDropFactory()));
        if (varName == "Aggregates")    return SetAndReturnDefaultFactory(varName, rcp(new UCAggregationFactory()));

        if (varName == "PreSmoother")   return GetFactory("Smoother");
        if (varName == "PostSmoother")  return GetFactory("Smoother");

        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::FactoryManager::GetDefaultFactory(): No default factory available for building '"+varName+"'.");
      }
      
    }

    //@}

    void Clean() const { defaultFactoryTable_.clear(); }

  private:
    
    //
    // Helper functions
    //

    //! Add a factory to the default factory list and return it. This helper function is used by GetDefaultFactory()
    //TODO factory->setObjectLabel("Default " + varName + "Factory");

    const RCP<const FactoryBase> & SetAndReturnDefaultFactory(const std::string & varName, const RCP<const FactoryBase> & factory) const {
      TEUCHOS_TEST_FOR_EXCEPTION(factory == Teuchos::null, Exceptions::RuntimeError, "");

      GetOStream(Warnings0,  0) << "Warning: No factory have been specified for building '" << varName << "'." << std::endl;
      GetOStream(Warnings00, 0) << "        using default factory: ";
      { Teuchos::OSTab tab(getOStream(), 7); factory->describe(GetOStream(Warnings00), GetVerbLevel()); }

      defaultFactoryTable_[varName] = factory;

      return defaultFactoryTable_[varName];
    }

    //! Test if factoryTable_[varName] exists
    static bool IsAvailable(const std::string & varName, const std::map<std::string, RCP<const FactoryBase> > & factoryTable) {
      return factoryTable.find(varName) != factoryTable.end();
    }

    //
    // Data structures
    //

    // Note 1: we distinguish 'user defined factory' and 'default factory' to allow the desallocation of default factories separatly.
    // Note 2: defaultFactoryTable_ is mutable because default factories are only added to the list when they are requested to avoid allocation of unused factories.

    std::map<std::string, RCP<const FactoryBase> > factoryTable_;        // User defined factories

    mutable 
    std::map<std::string, RCP<const FactoryBase> > defaultFactoryTable_; // Default factories
    
  }; // class

} // namespace MueLu

#define MUELU_FACTORYMANAGER_SHORT
#endif // ifndef MUELU_FACTORYMANAGER_HPP

//TODO: add operator[]
//TODO: should we use a parameterList instead of a std::map? It might be useful to tag which factory have been used and report unused factory.
//TODO: add an option 'NoDefault' to check if we are using any default factory.
//TODO: use Teuchos::ConstNonConstObjectContainer to allow user to modify factories after a GetFactory()
