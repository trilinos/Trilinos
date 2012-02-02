#ifndef MUELU_HIERARCHYFACTORY_DECL_HPP
#define MUELU_HIERARCHYFACTORY_DECL_HPP

#include <string>

#include <Teuchos_Array.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_HierarchyFactory_fwd.hpp"
#include "MueLu_HierarchyFactoryBase.hpp"

#include "MueLu_Hierarchy_fwd.hpp"
#include "MueLu_FactoryManagerBase_fwd.hpp"
#include "MueLu_FactoryManager_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp" //TMP

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class HierarchyFactory : public HierarchyFactoryBase {
#undef MUELU_HIERARCHYFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"
      
  public:

    //@{

    //!
    HierarchyFactory() { }

    //!
    HierarchyFactory(Teuchos::ParameterList & paramList);

    //!
    HierarchyFactory(const std::string & xmlFileName);
    
    //! Destructor.
    virtual ~HierarchyFactory() { }

    //@}

    void SetParameterList(Teuchos::ParameterList & paramList);

    //@{

    //!
    //    virtual RCP<const FactoryManagerBase> GetFactoryManager(int level = 0);

    //@}

    //@{

    //!
    RCP<Hierarchy> CreateHierarchy() const; // called created and not Build(), because just return an non-initialized object (MG setup not done)

    void SetupHierarchy(Hierarchy & H) const; // == h.Setup(this)

    //@}

  private:
    Array<FactoryManagerBase> config_; // store up to one FactoryManager per level.

  }; // class

} // namespace MueLu

#define MUELU_HIERARCHYFACTORY_SHORT
#endif // MUELU_HIERARCHYFACTORY_DECL_HPP

// TODO:
// - parameter list validator
// - SetParameterList
// - Set/Get directlty Level manager
// - build per level
// - comments
