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
  NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::NullspaceFactory(RCP<const FactoryBase> AFact, RCP<const FactoryBase> nullspaceFact)
    : AFact_(AFact), nullspaceFact_(nullspaceFact)
  { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~NullspaceFactory() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {

    // only request "A" in DeclareInput if
    // 1)there is not "Nullspace" is available in Level AND
    // 2) it is the finest level (i.e. LevelID == 0)
    if (currentLevel.IsAvailable("Nullspace") == false && currentLevel.GetLevelID() == 0)
      currentLevel.DeclareInput("A", AFact_.get(),this);

    if (currentLevel.GetLevelID() !=0) {
      currentLevel.DeclareInput("Nullspace", nullspaceFact_.get(),this);
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Nullspace factory", currentLevel);

    RCP<MultiVector> nullspace;

    //TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.GetLevelID() != 0, Exceptions::RuntimeError, "MueLu::NullspaceFactory::Build(): NullspaceFactory can be used for finest level (LevelID == 0) only.");

    if (currentLevel.GetLevelID() == 0) {

      if (currentLevel.IsAvailable("Nullspace")) {
        //FIXME: with the new version of Level::GetFactory(), this never happens.

        // When a fine nullspace have already been defined by user using Set("Nullspace", ...), we use it.
        nullspace = currentLevel.Get< RCP<MultiVector> >("Nullspace");
        GetOStream(Runtime1, 0) << "Use user-given nullspace: nullspace dimension=" << nullspace->getNumVectors() << std::endl;
        
      } else {
        
        RCP<Operator> A = currentLevel.Get< RCP<Operator> >("A", AFact_.get()); // no request since given by user
        
        LocalOrdinal numPDEs = A->GetFixedBlockSize();
        
        GetOStream(Runtime1, 0) << "Generating canonical nullspace: dimension = " << numPDEs << std::endl;
        nullspace = MultiVectorFactory::Build(A->getDomainMap(), numPDEs);
        
        for (int i=0; i<numPDEs; ++i) {
          ArrayRCP<Scalar> nsValues = nullspace->getDataNonConst(i);
          int numBlocks = nsValues.size() / numPDEs;
          for (int j=0; j< numBlocks; ++j) {
            nsValues[j*numPDEs + i] = 1.0;
          }
        }

      } 
    } else {

        nullspace = currentLevel.Get< RCP<MultiVector> >("Nullspace", nullspaceFact_.get());
       
      }

    currentLevel.Set("Nullspace", nullspace, this);
    
  } // Build

} //namespace MueLu

#endif // MUELU_NULLSPACEFACTORY_DEF_HPP
