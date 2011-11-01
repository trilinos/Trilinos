#ifndef MUELU_FACTORYMANAGER_DECL_HPP
#define MUELU_FACTORYMANAGER_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

#include <map>

#include <Xpetra_Operator.hpp>

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
#include "MueLu_GaussSeidelSmoother.hpp"
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
    FactoryManager(const RCP<const FactoryBase> PFact = Teuchos::null, const RCP<const FactoryBase> RFact = Teuchos::null, const RCP<const FactoryBase> AcFact = Teuchos::null) ;
    
    //! Destructor.
    virtual ~FactoryManager() ;

    //@}

    //@{ Get/Set functions.

    //! Set Factory
    void SetFactory(const std::string & varName, const RCP<const FactoryBase> & factory) ;

    //! Get Factory
    const RCP<const FactoryBase> & GetFactory(const std::string & varName) const ;

    //!
    const RCP<const FactoryBase> & GetDefaultFactory(const std::string & varName) const ;

    //@}

    void Clean() const ;

  private:
    
    //
    // Helper functions
    //

    //! Add a factory to the default factory list and return it. This helper function is used by GetDefaultFactory()
    //TODO factory->setObjectLabel("Default " + varName + "Factory");

    const RCP<const FactoryBase> & SetAndReturnDefaultFactory(const std::string & varName, const RCP<const FactoryBase> & factory) const ;

    //! Test if factoryTable_[varName] exists
    static bool IsAvailable(const std::string & varName, const std::map<std::string, RCP<const FactoryBase> > & factoryTable) ;

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

//TODO: add operator[]
//TODO: should we use a parameterList instead of a std::map? It might be useful to tag which factory have been used and report unused factory.
//TODO: add an option 'NoDefault' to check if we are using any default factory.
//TODO: use Teuchos::ConstNonConstObjectContainer to allow user to modify factories after a GetFactory()

#define MUELU_FACTORYMANAGER_SHORT
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_FACTORYMANAGER_DECL_HPP
