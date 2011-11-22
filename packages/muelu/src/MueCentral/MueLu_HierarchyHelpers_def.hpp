#ifndef MUELU_HIERARCHYHELPERS_DEF_HPP
#define MUELU_HIERARCHYHELPERS_DEF_HPP

#include <Xpetra_Operator.hpp>

#include "MueLu_HierarchyHelpers_decl.hpp"

#include "MueLu_SmootherBase.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TopRAPFactory(RCP<const FactoryManagerBase> parentFactoryManager)
    : factoryManager_(rcp( new InternalFactoryManager(parentFactoryManager))), PFact_(parentFactoryManager->GetFactory("P")), RFact_(parentFactoryManager->GetFactory("R")), AcFact_(parentFactoryManager->GetFactory("A"))
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~TopRAPFactory() { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & fineLevel, Level & coarseLevel) const {
    SetFactoryManager SFM2(fineLevel,   factoryManager_);
    SetFactoryManager SFM1(coarseLevel, factoryManager_);

    if (PFact_  != Teuchos::null) coarseLevel.DeclareInput("P", PFact_.get(), this);
    if (RFact_  != Teuchos::null) coarseLevel.DeclareInput("R", RFact_.get(), this);
    if (AcFact_ != Teuchos::null) coarseLevel.DeclareInput("A", AcFact_.get(), this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & fineLevel, Level & coarseLevel) const {
    SetFactoryManager SFM1(fineLevel,   factoryManager_);
    SetFactoryManager SFM2(coarseLevel, factoryManager_);

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

  //
  //
  //

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TopSmootherFactory(RCP<const FactoryManagerBase> parentFactoryManager, const std::string & varName)
  : factoryManager_(rcp( new InternalFactoryManager(parentFactoryManager))), smootherFact_(parentFactoryManager->GetFactory(varName))
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~TopSmootherFactory() { }
  
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & level) const {
    SetFactoryManager SFM(level, factoryManager_);

    if (smootherFact_ != Teuchos::null) {
      level.DeclareInput("PreSmoother",  smootherFact_.get(), this);
      level.DeclareInput("PostSmoother", smootherFact_.get(), this);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & level) const {
    typedef MueLu::SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> SmootherBase2; //TODO
 
    SetFactoryManager SFM(level, factoryManager_);

    if (smootherFact_ != Teuchos::null) {
      smootherFact_->CallBuild(level);

      if (level.IsAvailable("PreSmoother", smootherFact_.get())) {
        RCP<SmootherBase2> Pre  = level.Get<RCP<SmootherBase2> >("PreSmoother", smootherFact_.get());
        level.Set("PreSmoother", Pre);
      }

      if (level.IsAvailable("PostSmoother", smootherFact_.get())) {
        RCP<SmootherBase2> Post = level.Get<RCP<SmootherBase2> >("PostSmoother", smootherFact_.get());
        level.Set("PostSmoother", Post);
      }

    }
  }

} // namespace MueLu

// TODO: remove 'RCP' for TopRAPFactory::factoryManager_ and TopSmootherFactory::factoryManager_?

#define MUELU_HIERARCHY_HELPERS_SHORT
#endif // MUELU_HIERARCHYHELPERS_DEF_HPP
