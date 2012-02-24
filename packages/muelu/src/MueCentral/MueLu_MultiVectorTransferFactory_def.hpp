#ifndef MUELU_MULTIVECTORTRANSFER_FACTORY_DEF_HPP
#define MUELU_MULTIVECTORTRANSFER_FACTORY_DEF_HPP

#include "MueLu_MultiVectorTransferFactory_decl.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

// #include <Xpetra_Operator.hpp>

#include "MueLu_Level.hpp"
// #include "MueLu_Monitor.hpp"

namespace MueLu {

  // ----------------------------------------------------------------------------------------
  // Constructor
  // ----------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MultiVectorTransferFactory(
    std::string const & vectorName,
    std::string const & restrictionName,
    RCP<const FactoryBase> const &restrictionFact)
    : vectorName_(vectorName),
      restrictionName_(restrictionName),
      restrictionFact_(restrictionFact)
  { }

  // ----------------------------------------------------------------------------------------
  // Destructor
  // ----------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~MultiVectorTransferFactory() {}

  // ----------------------------------------------------------------------------------------
  // DeclareInput
  // ----------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    if (fineLevel.GetLevelID() == 0) fineLevel.DeclareInput(vectorName_, MueLu::NoFactory::get(), this);
    else                             fineLevel.DeclareInput(vectorName_, this                   , this);
    coarseLevel.DeclareInput(restrictionName_,restrictionFact_.get(),this);
    /*
    //FIXME ThreeLevels unit test dies
    fineLevel.DeclareInput(vectorName_, restrictionFact_.get(),this);
    coarseLevel.DeclareInput(restrictionName_,restrictionFact_.get(),this);
    */

  }

  // ----------------------------------------------------------------------------------------
  // Build
  // ----------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & fineLevel, Level &coarseLevel) const {

    //TEUCHOS_TEST_FOR_EXCEPTION(!fineLevel.IsAvailable(vectorName_,MueLu::NoFactory::get()), Exceptions::RuntimeError,
    //                    "MueLu::MultiVectorTransferFactory::Build(): vector '" + vectorName_ + "' is not available.");

    //RCP<MultiVector> vector  = fineLevel.Get<RCP<MultiVector> >(vectorName_,MueLu::NoFactory::get());
    RCP<MultiVector> vector;
    if (fineLevel.GetLevelID() == 0)
      vector  = fineLevel.Get<RCP<MultiVector> >(vectorName_,MueLu::NoFactory::get());
    else
      vector  = fineLevel.Get<RCP<MultiVector> >(vectorName_,this);
    RCP<Operator> transferOp = coarseLevel.Get<RCP<Operator> >(restrictionName_,restrictionFact_.get());

    /*
    //FIXME ThreeLevels unit test  dies
    RCP<MultiVector> vector  = fineLevel.Get<RCP<MultiVector> >(vectorName_,restrictionFact_.get());
    RCP<Operator> transferOp = coarseLevel.Get<RCP<Operator> >(restrictionName_,restrictionFact_.get());
    */

    RCP<MultiVector> result = MultiVectorFactory::Build(transferOp->getRangeMap(),1);
    GetOStream(Runtime0,0) << "Transferring multivector \"" << vectorName_ << "\"" << std::endl;
    transferOp->apply(*vector,*result);
    coarseLevel.Set<RCP<MultiVector> >(vectorName_,result,this);
    
  } //Build

  //TODO do we need methods to get name of MultiVector and or Transfer Operator?

} // namespace MueLu

#define MUELU_MULTIVECTORTRANSFER_FACTORY_SHORT
#endif // MUELU_MULTIVECTORTRANSFER_FACTORY_DEF_HPP
