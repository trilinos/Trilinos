#ifndef MUELU_HIERARCHY_HELPERS_HPP
#define MUELU_HIERARCHY_HELPERS_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_FactoryManagerBase.hpp"

namespace MueLu {

  //! An exception safe way to call the method 'Level::SetDefaultFactoryHandler()'
  class SetFactoryManager {
  public:

    //@{

    //!
    SetFactoryManager(Level & level, RCP<const FactoryManagerBase> & factoryManager)
      :  level_(level) 
    {
      level.SetDefaultFactoryHandler(factoryManager);
    }
    
    //! Destructor.
    virtual ~SetFactoryManager() { 
      level_.SetDefaultFactoryHandler(Teuchos::null);
    }

    //@}

  private:
    Level & level_;
  };

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class TopRAPFactory : public TwoLevelFactoryBase {
#include "MueLu_UseShortNames.hpp"

  public:

    TopRAPFactory(RCP<const FactoryManagerBase> parentFactoryManager, const RCP<const FactoryBase> PFact = Teuchos::null, const RCP<const FactoryBase> RFact = Teuchos::null, const RCP<const FactoryBase> AcFact = Teuchos::null)
      : parentFactoryManager_(parentFactoryManager), PFact_(PFact), RFact_(RFact), AcFact_(AcFact)
    { }
    
    virtual ~TopRAPFactory() { }

    void DeclareInput(Level & fineLevel, Level & coarseLevel) const {
      SetFactoryManager SFM1(coarseLevel, parentFactoryManager_);
      SetFactoryManager SFM2(fineLevel, parentFactoryManager_);
      
      if (PFact_  != Teuchos::null) coarseLevel.DeclareInput("P", PFact_.get());
      if (RFact_  != Teuchos::null) coarseLevel.DeclareInput("R", RFact_.get());
      if (AcFact_ != Teuchos::null) coarseLevel.DeclareInput("A", AcFact_.get());    
    }
    
    void Build(Level & fineLevel, Level & coarseLevel) const {
      SetFactoryManager SFM1(fineLevel,   parentFactoryManager_);
      SetFactoryManager SFM2(coarseLevel, parentFactoryManager_);

      if (PFact_ != Teuchos::null) {
        RCP<Operator> P = coarseLevel.Get<RCP<Operator> >("P", PFact_.get());
        coarseLevel.Set("P", P);
      }
      
      if (RFact_ != Teuchos::null) {
        RCP<Operator> R = coarseLevel.Get<RCP<Operator> >("R", RFact_.get());
        coarseLevel.Set("R", R);
      }
      
      if (AcFact_ != Teuchos::null) {
        RCP<Operator> Ac = coarseLevel.Get<RCP<Operator> >("A", AcFact_.get());
        coarseLevel.Set("A", Ac);
      }
    }
    
  private:
    mutable RCP<const FactoryManagerBase> parentFactoryManager_;
    RCP<const FactoryBase> PFact_;
    RCP<const FactoryBase> RFact_;
    RCP<const FactoryBase> AcFact_;
  };

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class TopSmootherFactory : public SingleLevelFactoryBase { //TODO: inherit from SmootherFactoryBase ?
#include "MueLu_UseShortNames.hpp"

  public:

    TopSmootherFactory(RCP<const FactoryManagerBase> parentFactoryManager, RCP<const FactoryBase> smootherFact)
      : parentFactoryManager_(parentFactoryManager), smootherFact_(smootherFact)
    { }

    virtual ~TopSmootherFactory() { }

    void DeclareInput(Level & level) const {
      SetFactoryManager SFM(level, parentFactoryManager_);
           
      if (smootherFact_ != Teuchos::null) {
        level.DeclareInput("PreSmoother",  smootherFact_.get());
        level.DeclareInput("PostSmoother", smootherFact_.get());
      }
    }

    void Build(Level & level) const {
      SetFactoryManager SFM(level, parentFactoryManager_);

      if (smootherFact_ != Teuchos::null) {
        smootherFact_->NewBuild(level);
          
        if (level.IsAvailable("PreSmoother", smootherFact_.get())) {
          RCP<SmootherBase> Pre  = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact_.get());
          level.Set("PreSmoother", Pre);
        }
          
        if (level.IsAvailable("PostSmoother", smootherFact_.get())) {
          RCP<SmootherBase> Post = level.Get<RCP<SmootherBase> >("PostSmoother", smootherFact_.get());
          level.Set("PostSmoother", Post);
        }
      }
    }

  private:
    mutable RCP<const FactoryManagerBase> parentFactoryManager_;
    RCP<const FactoryBase> smootherFact_;
  };

} // namespace MueLu
  
#define MUELU_HIERARCHY_HELPERS_SHORT
#endif //ifndef MUELU_HIERARCHY_HELPERS_HPP
