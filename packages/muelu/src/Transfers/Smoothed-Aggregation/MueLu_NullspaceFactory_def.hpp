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
    : nspName_("Nullspace"), AFact_(AFact), nullspaceFact_(nullspaceFact)
  { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::NullspaceFactory(std::string nspName, RCP<const FactoryBase> nullspaceFact)
    : nspName_(nspName), AFact_(Teuchos::null), nullspaceFact_(nullspaceFact)
  {  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~NullspaceFactory() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {

    // only request "A" in DeclareInput if
    // 1) there is not nspName_ (e.g. "Nullspace") is available in Level AND
    // 2) it is the finest level (i.e. LevelID == 0)
    if (currentLevel.IsAvailable(nspName_) == false && currentLevel.GetLevelID() == 0)
      currentLevel.DeclareInput("A", AFact_.get(),this);

    if (currentLevel.GetLevelID() !=0) {

      // validate nullspaceFact_

      // 1) nullspaceFact_ must not be Teuchos::null, since the default factory for "Nullspace" is
      //    a NullspaceFactory
      // 2) nullspaceFact_ must be a TentativePFactory i.e. at least a TwoLevelFactoryBase derived object


      currentLevel.DeclareInput("Nullspace", nullspaceFact_.get(),this);
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Nullspace factory", currentLevel);

    RCP<MultiVector> nullspace;

    //TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.GetLevelID() != 0, Exceptions::RuntimeError, "MueLu::NullspaceFactory::Build(): NullspaceFactory can be used for finest level (LevelID == 0) only.");

    if (currentLevel.GetLevelID() == 0) {

      if (currentLevel.IsAvailable(nspName_)) {
        //FIXME: with the new version of Level::GetFactory(), this never happens.

        // When a fine nullspace have already been defined by user using Set("Nullspace", ...), we use it.
        nullspace = currentLevel.Get< RCP<MultiVector> >(nspName_);
        GetOStream(Runtime1, 0) << "Use user-given nullspace " << nspName_ << ": nullspace dimension=" << nullspace->getNumVectors() << std::endl;
      
      } else {
        // "Nullspace" (nspName_) is not available
        RCP<Operator> A = currentLevel.Get< RCP<Operator> >("A", AFact_.get()); // no request since given by user
        
        // determine numPDEs
        LocalOrdinal numPDEs = 1;
        if(A->IsView("stridedMaps")==true) {
          Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps"); // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
          TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap()) == Teuchos::null,Exceptions::BadCast,"MueLu::CoalesceFactory::Build: cast to strided row map failed.");
          numPDEs = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getFixedBlockSize();
          oldView = A->SwitchToView(oldView);
        }
        
        GetOStream(Runtime1, 0) << "Generating canonical nullspace: dimension = " << numPDEs << std::endl;
        nullspace = MultiVectorFactory::Build(A->getDomainMap(), numPDEs);
        
        for (int i=0; i<numPDEs; ++i) {
          ArrayRCP<Scalar> nsValues = nullspace->getDataNonConst(i);
          int numBlocks = nsValues.size() / numPDEs;
          for (int j=0; j< numBlocks; ++j) {
            nsValues[j*numPDEs + i] = 1.0;
          }
        }
      } // end if "Nullspace" not available
    } else {
        // on coarser levels always use "Nullspace" as variable name, since it is expected by
        // tentative P factory to be "Nullspace"

        nullspace = currentLevel.Get< RCP<MultiVector> >("Nullspace", nullspaceFact_.get());
       
      }

    // provide "Nullspace" variable on current level (used by TentativePFactory)
    currentLevel.Set("Nullspace", nullspace, this);
    
  } // Build

} //namespace MueLu

#endif // MUELU_NULLSPACEFACTORY_DEF_HPP
