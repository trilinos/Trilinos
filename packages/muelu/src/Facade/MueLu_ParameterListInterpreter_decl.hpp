#ifndef MUELU_PARAMETERLISTINTERPRETER_DECL_HPP
#define MUELU_PARAMETERLISTINTERPRETER_DECL_HPP

#include <string>

#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Xpetra_MultiVectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_HierarchyFactory_fwd.hpp"
#include "MueLu_HierarchyManager.hpp"

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class ParameterListInterpreter : public HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> {
#undef MUELU_PARAMETERLISTINTERPRETER_SHORT
#include "MueLu_UseShortNames.hpp"
      
  public:

    //@{

    //!
    ParameterListInterpreter() { }

    //!
    ParameterListInterpreter(Teuchos::ParameterList & paramList);

    //!
    ParameterListInterpreter(const std::string & xmlFileName);
    
    //! Destructor.
    virtual ~ParameterListInterpreter() { }

    //@}

    //@{

    void SetParameterList(const Teuchos::ParameterList & paramList);

    //@}

  private:
    typedef std::map<std::string, RCP<const FactoryBase> > FactoryMap; //TODO: remove this line

    void BuildFactoryMap(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn, FactoryMap & factoryMapOut) const;

    //@{ Operator configuration

    //! Setup Operator object
    virtual void SetupOperator(Operator & Op) const;

    //! Setup extra data
    virtual void SetupExtra(Hierarchy & H) const;

    // Operator configuration storage
    Teuchos::ParameterList operatorList_; //TODO: should it be stored in another format to avoid xml parsing in SetupOperator()?

    //@}

  }; // class

} // namespace MueLu

#define MUELU_PARAMETERLISTINTERPRETER_SHORT
#endif // MUELU_PARAMETERLISTINTERPRETER_DECL_HPP

// TODO:
// - parameter list validator
// - SetParameterList
// - Set/Get directly Level manager
// - build per level
// - comments/docs
// - use FactoryManager instead of FactoryMap
