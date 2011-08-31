#ifndef MUELU_NULLSPACEFACTORY_HPP
#define MUELU_NULLSPACEFACTORY_HPP

#include "Xpetra_Operator.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class NullspaceFactory : public SingleLevelFactoryBase {

#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    NullspaceFactory(RCP<FactoryBase> AFact = Teuchos::null) 
      : AFact_(AFact)
    { }

    //! Destructor
    virtual ~NullspaceFactory() {}

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) {
      currentLevel.Input("A", AFact_());
    }

    //@}

    bool Build(Level &currentLevel) const {
      RCP<Operator> A = currentLevel.NewGet< RCP<Operator> >("A", AFact_());

      //FIXME this doesn't check for the #dofs per node, or whether we have a blocked system

      RCP<MultiVector> nullspace = MultiVectorFactory::Build(A->getDomainMap(), 1);
      nullspace->putScalar(1.0);

      currentLevel.NewSet("Nullspace", nullspace, this);

      return true; //??

    } // Build

  private:
    //! A Factory
    RCP<FactoryBase> AFact_;

  }; //class NullspaceFactory

} //namespace MueLu

#define MUELU_NULLSPACEFACTORY_SHORT

#endif //ifndef MUELU_NULLSPACEFACTORY_HPP
