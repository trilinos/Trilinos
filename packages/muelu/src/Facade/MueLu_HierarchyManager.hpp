#ifndef MUELU_HIERARCHYMANAGER_DECL_HPP
#define MUELU_HIERARCHYMANAGER_DECL_HPP

#include <string>
#include <map>

#include <Teuchos_Array.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_HierarchyFactory.hpp"

#include "MueLu_Hierarchy.hpp"
#include "MueLu_Exceptions.hpp"

#include "MueLu_RAPFactory.hpp" //TMP

namespace MueLu {

  // This class stores the configuration of a Hierarchy.
  // The class also provides an algorithm to build a Hierarchy from the configuration.
  //
  // See also: FactoryManager
  //
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void, LocalOrdinal, Node>::SparseOps>
  class HierarchyManager : public HierarchyFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> { 
#undef MUELU_HIERARCHYMANAGER_SHORT
#include "MueLu_UseShortNames.hpp"
      
  public:

    //!
    HierarchyManager() { }

    //!
    virtual ~HierarchyManager() { }

    // Unused
    // void AddFactoryManager(int levelID, RCP<FactoryManagerBase>& manager) {
    //   if (levelManagers_.size() < levelID + 1) levelManagers_.resize(levelID + 1);
    //   levelManagers_[levelID] = manager;
    // }

    //!
    void AddFactoryManager(int startLevel, int numDesiredLevel, RCP<FactoryManagerBase>& manager) {
      const int lastLevel = startLevel + numDesiredLevel - 1;
      if (levelManagers_.size() < lastLevel + 1) levelManagers_.resize(lastLevel + 1);
      
      for(int iLevel = startLevel; iLevel <= lastLevel; iLevel++) {
        levelManagers_[iLevel] = manager;
      }
    }

    //!
    void SetFactoryManagerCoarsestLevel(RCP<FactoryManagerBase>& manager) {
      coarsestLevelManager_ = manager;
    }

    //!
    void CheckConfig() {
      for(int i=0; i<levelManagers_.size(); i++) {
        TEUCHOS_TEST_FOR_EXCEPTION(levelManagers_[i] == Teuchos::null, Exceptions::RuntimeError, "MueLu:HierarchyConfig::CheckConfig(): Undefined configuration for level:");
      }
    }

    //@{

    virtual RCP<Hierarchy> CreateHierarchy() const {
      return rcp(new Hierarchy());
    }
    
    //! Setup Hierarchy object
    virtual void SetupHierarchy(Hierarchy & H) const {
      H.Setup();
    }

    //@}

    typedef std::map<std::string, RCP<const FactoryBase> > FactoryMap;

  protected: //TODO: access function

    // Hierarchy parameters
    int numDesiredLevel_; 

  private:
    // Levels
    Array<RCP<FactoryManagerBase> > levelManagers_;        // one FactoryManager per level.
    RCP<FactoryManagerBase>         coarsestLevelManager_; // coarsest level manager

  }; // class HierarchyManager

} // namespace MueLu

#define MUELU_HIERARCHYMANAGER_SHORT
#endif // MUELU_HIERARCHYMANAGER_HPP

//TODO: split into _decl/_def
