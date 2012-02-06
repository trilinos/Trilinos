#ifndef MUELU_HIERARCHYHELPERS_DECL_HPP
#define MUELU_HIERARCHYHELPERS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_HierarchyHelpers_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_FactoryManagerBase.hpp"

#include "MueLu_Level_fwd.hpp"

// Warning: on TopRAPFactory and TopSmootherFactory constructors, Teuchos::null doesn't mean "default factory" but "no build"

namespace MueLu {

  //! An exception safe way to call the method 'Level::SetFactoryManager()'
  class SetFactoryManager {

  public:

    //@{

    //!
    SetFactoryManager(Level & level, const RCP<const FactoryManagerBase> & factoryManager)
      : level_(level), prevFactoryManager_(level.GetFactoryManager())
    {
      // set new factory manager
      level.SetFactoryManager(factoryManager);
    }

    //! Destructor.
    virtual ~SetFactoryManager() {
      // restore previous factory manager
      //FIXME      level_.SetFactoryManager(prevFactoryManager_);
    }

    //@}

  private:
    Level & level_;
    const RCP<const FactoryManagerBase> prevFactoryManager_; // save & restore previous factoryManager
  };

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class TopRAPFactory : public TwoLevelFactoryBase {
#include "MueLu_UseShortNames.hpp"

  public:

    TopRAPFactory(RCP<const FactoryManagerBase> parentFactoryManager);
    TopRAPFactory(RCP<const FactoryManagerBase> parentFactoryManagerFine, RCP<const FactoryManagerBase> parentFactoryManagerCoarse);

    virtual ~TopRAPFactory();

    void DeclareInput(Level & fineLevel, Level & coarseLevel) const;

    void Build(Level & fineLevel, Level & coarseLevel) const;

  private:
    RCP<const FactoryBase> PFact_;
    RCP<const FactoryBase> RFact_;
    RCP<const FactoryBase> AcFact_;
  };

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class TopSmootherFactory : public SingleLevelFactoryBase { //TODO: inherit from SmootherFactoryBase ?
#include "MueLu_UseShortNames.hpp"

  public:

    TopSmootherFactory(RCP<const FactoryManagerBase> parentFactoryManager, const std::string & varName);

    virtual ~TopSmootherFactory();

    void DeclareInput(Level & level) const;

    void Build(Level & level) const;

  private:
    RCP<const FactoryBase> smootherFact_;
  };

} // namespace MueLu

#define MUELU_HIERARCHYHELPERS_SHORT
#endif // MUELU_HIERARCHYHELPERS_DECL_HPP
