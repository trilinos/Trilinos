#ifndef MUELU_FACTORYMANAGER_HPP
#define MUELU_FACTORYMANAGER_HPP

#include <map>

#include <Xpetra_Operator.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_NoFactory.hpp"
#include "MueLu_SmootherFactoryBase.hpp"
#include "MueLu_SmootherBase.hpp"

#include "MueLu_DefaultFactoryHandler.hpp"

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class FactoryManager  { // inherit from FactoryHandlerBase
#include "MueLu_UseShortNames.hpp"
      
  public:

    //@{

    //!
    FactoryManager(const RCP<const FactoryBase> PFact = Teuchos::null, const RCP<const FactoryBase> RFact = Teuchos::null, const RCP<const FactoryBase> AcFact = Teuchos::null)
      : PFact_(PFact), RFact_(RFact), AcFact_(AcFact), factoryManager_(rcp(new DefaultFactoryHandler()))
    { }
    
    //! Destructor.
    virtual ~FactoryManager() { }

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

    RCP<DefaultFactoryHandlerBase> GetFactoryManager() const { return factoryManager_; }

    //@}

  private:

    RCP<const FactoryBase> PFact_;
    RCP<const FactoryBase> RFact_;
    RCP<const FactoryBase> AcFact_;

    RCP<const FactoryBase> smootherFact_;
    RCP<const FactoryBase> coarsestSolverFact_;
    
    mutable RCP<DefaultFactoryHandlerBase> factoryManager_;

  }; // class

} //namespace MueLu

#define MUELU_FACTORYMANAGER_SHORT
#endif // ifndef MUELU_FACTORYMANAGER_HPP

//TODO: clean up factory created at the end
