#ifndef MUELU_HIERARCHYHELPERS_DEF_HPP
#define MUELU_HIERARCHYHELPERS_DEF_HPP

#include <Xpetra_Operator.hpp>

#include "MueLu_HierarchyHelpers_decl.hpp"

#include "MueLu_SmootherBase.hpp"

//TODO/FIXME: DeclareInput(, **this**) cannot be used here

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TopRAPFactory(RCP<const FactoryManagerBase> parentFactoryManager)
    : factoryManagerFine_(parentFactoryManager), factoryManagerCoarse_(parentFactoryManager), PFact_(parentFactoryManager->GetFactory("P")), RFact_(parentFactoryManager->GetFactory("R")), AcFact_(parentFactoryManager->GetFactory("A"))
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TopRAPFactory(RCP<const FactoryManagerBase> parentFactoryManagerFine, RCP<const FactoryManagerBase> parentFactoryManagerCoarse)
    :  factoryManagerFine_(parentFactoryManagerFine), factoryManagerCoarse_(parentFactoryManagerCoarse), PFact_(parentFactoryManagerCoarse->GetFactory("P")), RFact_(parentFactoryManagerCoarse->GetFactory("R")), AcFact_(parentFactoryManagerCoarse->GetFactory("A"))
  { }
  
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~TopRAPFactory() { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & fineLevel, Level & coarseLevel) const {
    SetFactoryManager SFM2(fineLevel,   factoryManagerFine_);
    SetFactoryManager SFM1(coarseLevel, factoryManagerCoarse_);

    if (PFact_  != Teuchos::null)                                       coarseLevel.DeclareInput("P", PFact_.get());
    if (RFact_  != Teuchos::null)                                       coarseLevel.DeclareInput("R", RFact_.get());
    if ((AcFact_ != Teuchos::null) && (AcFact_ != NoFactory::getRCP())) coarseLevel.DeclareInput("A", AcFact_.get());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & fineLevel, Level & coarseLevel) const {
    SetFactoryManager SFM1(fineLevel,   factoryManagerFine_);
    SetFactoryManager SFM2(coarseLevel, factoryManagerCoarse_);

    if (PFact_ != Teuchos::null) {
      RCP<Operator> P = coarseLevel.Get<RCP<Operator> >("P", PFact_.get());
      coarseLevel.Set("P", P, NoFactory::get());
      coarseLevel.AddKeepFlag("P", NoFactory::get(), MueLu::Final);       // FIXME2: Order of Remove/Add matter (data removed otherwise). Should do something about this
      coarseLevel.RemoveKeepFlag("P", NoFactory::get(), MueLu::UserData); // FIXME: This is a hack, I should change behavior of Level::Set() instead.
    }

    if (RFact_ != Teuchos::null) {
      RCP<Operator> R = coarseLevel.Get<RCP<Operator> >("R", RFact_.get());
      coarseLevel.Set("R", R, NoFactory::get());
      coarseLevel.AddKeepFlag("R", NoFactory::get(), MueLu::Final);
      coarseLevel.RemoveKeepFlag("R", NoFactory::get(), MueLu::UserData); // FIXME: This is a hack
    }

    if ((AcFact_ != Teuchos::null) && (AcFact_ != NoFactory::getRCP())) {
      RCP<Operator> Ac = coarseLevel.Get<RCP<Operator> >("A", AcFact_.get());
      coarseLevel.Set("A", Ac, NoFactory::get());
      coarseLevel.AddKeepFlag("A", NoFactory::get(), MueLu::Final);
      coarseLevel.RemoveKeepFlag("A", NoFactory::get(), MueLu::UserData); // FIXME: This is a hack
    }
  }

  //
  //
  //

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TopSmootherFactory(RCP<const FactoryManagerBase> parentFactoryManager, const std::string & varName)
  : factoryManager_(parentFactoryManager), smootherFact_(parentFactoryManager->GetFactory(varName))
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~TopSmootherFactory() { }
  
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & level) const {
    SetFactoryManager SFM(level, factoryManager_);

    if (smootherFact_ != Teuchos::null) {
      level.DeclareInput("PreSmoother",  smootherFact_.get());
      level.DeclareInput("PostSmoother", smootherFact_.get());
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & level) const {
    typedef MueLu::SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> SmootherBase2; //TODO
 
    SetFactoryManager SFM(level, factoryManager_);

    // Teuchos::null == skip
    if (smootherFact_ != Teuchos::null) {

      // Only call factory if at least one smoother is missing (mimic behavior of level.Get<> but level.Get<> cannot be used here as we don't know if the factory will produce both Pre and Post smoother)
      if (!level.IsAvailable("PreSmoother", smootherFact_.get()) || !level.IsAvailable("PostSmoother", smootherFact_.get())) {
        smootherFact_->CallBuild(level);
      }

      // Factory might or might not have created a pre smoother
      if (level.IsAvailable("PreSmoother", smootherFact_.get())) {
        RCP<SmootherBase2> Pre  = level.Get<RCP<SmootherBase2> >("PreSmoother", smootherFact_.get());
        level.Set("PreSmoother", Pre, NoFactory::get());
        level.AddKeepFlag("PreSmoother", NoFactory::get(), MueLu::Final);
        level.RemoveKeepFlag("PreSmoother", NoFactory::get(), MueLu::UserData); // FIXME: This is a hack
      }

      // Factory might or might not have created a post smoother
      if (level.IsAvailable("PostSmoother", smootherFact_.get())) {
        RCP<SmootherBase2> Post = level.Get<RCP<SmootherBase2> >("PostSmoother", smootherFact_.get());
        level.Set("PostSmoother", Post, NoFactory::get());
        level.AddKeepFlag("PostSmoother", NoFactory::get(), MueLu::Final);
        level.RemoveKeepFlag("PostSmoother", NoFactory::get(), MueLu::UserData); // FIXME: This is a hack
      }

    }
  }

} // namespace MueLu

// TODO: remove 'RCP' for TopRAPFactory::factoryManager_ and TopSmootherFactory::factoryManager_?

#define MUELU_HIERARCHY_HELPERS_SHORT
#endif // MUELU_HIERARCHYHELPERS_DEF_HPP
