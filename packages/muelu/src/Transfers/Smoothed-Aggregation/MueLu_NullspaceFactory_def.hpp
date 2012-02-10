#ifndef MUELU_NULLSPACEFACTORY_DEF_HPP
#define MUELU_NULLSPACEFACTORY_DEF_HPP

#include <Xpetra_Operator.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_NullspaceFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::NullspaceFactory(RCP<const FactoryBase> AFact)
    : AFact_(AFact)
  { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~NullspaceFactory() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    //GetOStream(Warnings1, 0) << "NullspaceFactory::DeclareInput: GetA by fac: " << AFact_.get() << std::endl;
    if (currentLevel.IsAvailable("Nullspace") == false && currentLevel.GetLevelID() == 0)
      currentLevel.DeclareInput("A", AFact_.get(),this);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const {
    RCP<MultiVector> nullspace;

    Monitor m(*this, "Nullspace factory");

    TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.GetLevelID() != 0, Exceptions::RuntimeError, "MueLu::NullspaceFactory::Build(): NullspaceFactory can be used for finest level (LevelID == 0) only.");

    if (currentLevel.IsAvailable("Nullspace")) {
      // When a fine nullspace have already been defined by user using Set("Nullspace", ...), we use it.
      nullspace = currentLevel.Get< RCP<MultiVector> >("Nullspace");
      GetOStream(Runtime1, 0) << "Use user-given nullspace: nullspace dimension=" << nullspace->getNumVectors() << std::endl;

    } else {
        
      RCP<Operator> A = currentLevel.Get< RCP<Operator> >("A", AFact_.get()); // no request since given by user

      //FIXME this doesn't check for the #dofs per node, or whether we have a blocked system

      nullspace = MultiVectorFactory::Build(A->getDomainMap(), 1);
      nullspace->putScalar(1.0);
      GetOStream(Runtime1, 0) << "Calculate nullspace: nullspace dimension=" << nullspace->getNumVectors() << std::endl;
    }

    currentLevel.Set("Nullspace", nullspace, this);

  } // Build

} //namespace MueLu

#endif // MUELU_NULLSPACEFACTORY_DEF_HPP
