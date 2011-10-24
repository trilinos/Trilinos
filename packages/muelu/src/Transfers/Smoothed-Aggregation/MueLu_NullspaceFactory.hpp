#ifndef MUELU_NULLSPACEFACTORY_HPP
#define MUELU_NULLSPACEFACTORY_HPP

#include "Xpetra_Operator.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

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
    NullspaceFactory(RCP<const FactoryBase> AFact = Teuchos::null)
      : AFact_(AFact)
    { }

    //! Destructor
    virtual ~NullspaceFactory() {}

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const {
      currentLevel.DeclareInput("A", AFact_.get());
    }

    //@}

    void Build(Level &currentLevel) const {
      RCP<MultiVector> nullspace;
      
      Monitor m(*this, "Nullspace factory");

      if (currentLevel.IsAvailable("Nullspace")) {
        // When a fine nullspace have already been defined by user using Set("Nullspace", ...), we use it.
        nullspace = currentLevel.Get< RCP<MultiVector> >("Nullspace");
        GetOStream(Runtime1, 0) << "Use user-given nullspace: nullspace dimension=" << nullspace->getNumVectors() << std::endl;

      } else {
        
        RCP<Operator> A = currentLevel.Get< RCP<Operator> >("A", AFact_.get());

        //FIXME this doesn't check for the #dofs per node, or whether we have a blocked system

        nullspace = MultiVectorFactory::Build(A->getDomainMap(), 1);
        nullspace->putScalar(1.0);
        GetOStream(Runtime1, 0) << "Calculate nullspace: nullspace dimension=" << nullspace->getNumVectors() << std::endl;
      }

      if(currentLevel.GetLevelID() != 0)
        GetOStream(Warnings0, 0) << "NullspaceFactory::Build called for Level " << currentLevel.GetLevelID() << "! Check me!" << std::endl;

      currentLevel.Set("Nullspace", nullspace, this);

    } // Build

  private:

    //! A Factory
    RCP<const FactoryBase> AFact_;

  }; //class NullspaceFactory

} //namespace MueLu

#define MUELU_NULLSPACEFACTORY_SHORT

#endif //ifndef MUELU_NULLSPACEFACTORY_HPP
