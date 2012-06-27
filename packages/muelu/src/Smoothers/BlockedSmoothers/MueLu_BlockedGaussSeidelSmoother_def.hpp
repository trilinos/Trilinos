/*
 * MueLu_BlockedGaussSeidelSmoother_def.hpp
 *
 *  Created on: 30.01.2012
 *      Author: tobias
 */

#ifndef MUELU_BLOCKEDGAUSSSEIDELSMOOTHER_DEF_HPP_
#define MUELU_BLOCKEDGAUSSSEIDELSMOOTHER_DEF_HPP_

#include "Teuchos_ArrayViewDecl.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Operator.hpp>
#include <Xpetra_BlockedCrsOperator.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_BlockedGaussSeidelSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_HierarchyHelpers.hpp"
#include "MueLu_SmootherBase.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BlockedGaussSeidelSmoother(LocalOrdinal sweeps, Scalar omega, RCP<FactoryBase> AFact)
    : type_("blocked GaussSeidel"), nSweeps_(sweeps), omega_(omega), AFact_(AFact), A_(Teuchos::null)
  {
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~BlockedGaussSeidelSmoother() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager) {
    FactManager_.push_back(FactManager);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A", AFact_.get());

    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

    // loop over all factory managers for the subblocks of blocked operator A
    std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
    for(it = FactManager_.begin(); it!=FactManager_.end(); ++it) {
      SetFactoryManager currentSFM  (rcpFromRef(currentLevel),   *it);

      // request "Smoother" for current subblock row.
      currentLevel.DeclareInput("PreSmoother",(*it)->GetFactory("Smoother").get()); // TODO check me: what about postsmoother?
      currentLevel.DeclareInput("A",(*it)->GetFactory("A").get()); // request A for comparing maps (to setup a map of the real block rows and the ordering of given blocks in Inverse_)
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level &currentLevel) {
    //typedef Xpetra::BlockedCrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> BlockedCrsOOperator;

    FactoryMonitor m(*this, "Setup blocked Gauss-Seidel Smoother", currentLevel);
    if (SmootherPrototype::IsSetup() == true) this->GetOStream(Warnings0, 0) << "Warning: MueLu::BlockedGaussSeidelSmoother::Setup(): Setup() has already been called";

    // extract blocked operator A from current level
    A_ = currentLevel.Get< RCP<Operator> >("A", AFact_.get()); // A needed for extracting map extractors
    RCP<BlockedCrsOperator> bA = Teuchos::rcp_dynamic_cast<BlockedCrsOperator>(A_);
    TEUCHOS_TEST_FOR_EXCEPTION(bA==Teuchos::null, Exceptions::BadCast, "MueLu::BlockedPFactory::Build: input matrix A is not of type BlockedCrsOperator! error.");

    // plausibility check
    TEUCHOS_TEST_FOR_EXCEPTION(bA->Rows() != FactManager_.size(), Exceptions::RuntimeError, "MueLu::BlockedPFactory::Build: number of block rows of A does not match number of SubFactoryManagers. error.");
    TEUCHOS_TEST_FOR_EXCEPTION(bA->Cols() != FactManager_.size(), Exceptions::RuntimeError, "MueLu::BlockedPFactory::Build: number of block cols of A does not match number of SubFactoryManagers. error.");

    // store map extractors
    rangeMapExtractor_  = bA->getRangeMapExtractor();
    domainMapExtractor_ = bA->getDomainMapExtractor();

    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

    // TODO do setup -> called by SmootherFactory::BuildSmoother

    // loop over all factory managers for the subblocks of blocked operator A
    size_t bgsOrderingIndex = 0;
    //std::map<size_t,size_t> bgsOrderingIndex2blockRowIndex;
    std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
    for(it = FactManager_.begin(); it!=FactManager_.end(); ++it) {
      SetFactoryManager currentSFM  (rcpFromRef(currentLevel),   *it);

      // extract Smoother for current block row (BGS ordering)
      RCP<const SmootherBase> Smoo = currentLevel.Get< RCP<SmootherBase> >("PreSmoother",(*it)->GetFactory("Smoother").get());
      Inverse_.push_back(Smoo);

      // extract i-th diagonal block Aii -> determine block
      RCP<Operator> Aii = currentLevel.Get< RCP<Operator> >("A",(*it)->GetFactory("A").get());
      for(size_t i = 0; i<rangeMapExtractor_->NumMaps(); i++) {
        if(rangeMapExtractor_->getMap(i)->isSameAs(*(Aii->getRangeMap()))) {
          // map found: i is the true block row index
          // fill std::map with information bgsOrderingIndex -> i
          bgsOrderingIndex2blockRowIndex_[bgsOrderingIndex] = i;
          break;
        }
      }

      bgsOrderingIndex++;
    }

    SmootherPrototype::IsSetup(true);
  }

  // This function is equivalent to the std::map 'at' method introduced in C++11.
  // TODO: move to Utils
  template <class StdMapType>
  const typename StdMapType::mapped_type & at(const StdMapType& map, const typename StdMapType::key_type& x) {
    typename StdMapType::const_iterator it = map.find(x);
    TEUCHOS_TEST_FOR_EXCEPTION(it == map.end(), std::out_of_range, "MueLu::at(): element does not exist in the map");
    return it->second;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Apply(): Setup() has not been called");

    RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > residual = MultiVectorFactory::Build(B.getMap(), B.getNumVectors());
    RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tempres = MultiVectorFactory::Build(B.getMap(), B.getNumVectors());
    RCP<MultiVector> rcpX = Teuchos::rcpFromRef(X);


    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));


    // outer Richardson loop
    for (LocalOrdinal run = 0; run < nSweeps_; ++run) {

      // one BGS sweep

      //*fos << "BGS sweep: " << run << std::endl;

      // loop over all block rows
      for(size_t i = 0; i<Inverse_.size(); i++) {

        // calculate block residual r = B-A*X
        // note: A_ is the full blocked operator
        //*fos << "BGS sweep: " << run << ", r = B - A*X with i=" << i << std::endl;
        residual->update(1.0,B,0.0); // r = B
        A_->apply(X, *residual, Teuchos::NO_TRANS, -1.0, 1.0);

        // extract corresponding subvectors from X and residual
        size_t blockRowIndex = at(bgsOrderingIndex2blockRowIndex_, i); // == bgsOrderingIndex2blockRowIndex_.at(i) (only available since C++11)
        Teuchos::RCP<MultiVector> Xi = domainMapExtractor_->ExtractVector(rcpX, blockRowIndex);
        Teuchos::RCP<MultiVector> ri = rangeMapExtractor_->ExtractVector(residual, blockRowIndex);

        Teuchos::RCP<MultiVector> tXi = domainMapExtractor_->getVector(blockRowIndex, X.getNumVectors());

        // apply solver/smoother
        //*fos << "BGS sweep: " << run << ", x = A_ii^{-1} r with i = " << i << std::endl;
        Inverse_.at(i)->Apply(*tXi, *ri, false);

        // update vector
        //*fos << "BGS sweep: " << run << ", update x_" << i << std::endl;
        Xi->update(omega_,*tXi,1.0);  // X_{i+1} = X_i + omega \Delta X_i

        // update corresponding part of rhs and lhs
        //*fos << "BGS sweep: " << run << ", finish substep i=" << i << std::endl;
        domainMapExtractor_->InsertVector(Xi, blockRowIndex, rcpX); // TODO wrong! fix me
      }
    }

  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    return rcp( new BlockedGaussSeidelSmoother(*this) );
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << " Sweeps: " << nSweeps_ << " damping: " << omega_ << std::endl;
    }

    if (verbLevel & Debug) {
      std::map<size_t,size_t>::const_iterator itOrdering;
      for(itOrdering = bgsOrderingIndex2blockRowIndex_.begin(); itOrdering!=bgsOrderingIndex2blockRowIndex_.end(); itOrdering++) {
        std::cout << "block GaussSeidel ordering index: " << (*itOrdering).first << " -> block row in blocked A: " << (*itOrdering).second << std::endl;
      }
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
    }
  }

} // namespace MueLu



#endif /* MUELU_BLOCKEDGAUSSSEIDELSMOOTHER_DEF_HPP_ */
