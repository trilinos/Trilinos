#ifndef MUELU_MULTIVECTORTRANSFER_FACTORY_DEF_HPP
#define MUELU_MULTIVECTORTRANSFER_FACTORY_DEF_HPP

#include "MueLu_MultiVectorTransferFactory_decl.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

// #include <Xpetra_Operator.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

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

    FactoryMonitor m(*this, "MultiVectorTransferFactory", coarseLevel);

    //TEUCHOS_TEST_FOR_EXCEPTION(!fineLevel.IsAvailable(vectorName_,MueLu::NoFactory::get()), Exceptions::RuntimeError,
    //                    "MueLu::MultiVectorTransferFactory::Build(): vector '" + vectorName_ + "' is not available.");

    //RCP<MultiVector> vector  = fineLevel.Get<RCP<MultiVector> >(vectorName_,MueLu::NoFactory::get());

    RCP<MultiVector> vector;
    /*
    //FIXME JJH reenable this once we sort out dependencies of tranferring, permuting vectors
    if (fineLevel.GetLevelID() == 0) {
      vector  = fineLevel.Get<RCP<MultiVector> >(vectorName_,MueLu::NoFactory::get());
      std::cout << "MultiVectorTransferFactory::Build -- requesting " << vectorName_ << " from factory " << MueLu::NoFactory::get() << std::endl;
    } else {
      vector  = fineLevel.Get<RCP<MultiVector> >(vectorName_,this);
      std::cout << "MultiVectorTransferFactory::Build -- requesting " << vectorName_ << " from factory " << this << std::endl;
    }
    */

    if (vectorName_ == "Coordinates") { // very elegant!

      // Convert format to Xpetra::MultiVector
      if (fineLevel.IsAvailable("XCoordinates") && !fineLevel.IsAvailable("Coordinates")) {
        std::cout << "Converting coordinates from 3xArrayRCP to MultiVector" << std::endl;
        
        TEUCHOS_TEST_FOR_EXCEPTION(fineLevel.GetLevelID() != 0, Exceptions::RuntimeError, "??" << fineLevel.GetLevelID());

        Array< ArrayView<const double> > arrayOfPtrs;
        
        {
          ArrayRCP<SC> & coords = fineLevel.Get<ArrayRCP<SC> >("XCoordinates");
          arrayOfPtrs.push_back(coords());
        }
        
        if(fineLevel.IsAvailable("YCoordinates")) {
          ArrayRCP<SC> & coords = fineLevel.Get<ArrayRCP<SC> >("YCoordinates");
          arrayOfPtrs.push_back(coords());
        }
        
        if(fineLevel.IsAvailable("ZCoordinates")) {
          TEUCHOS_TEST_FOR_EXCEPTION(!fineLevel.IsAvailable("YCoordinates"), Exceptions::RuntimeError, "ZCoordinates specified but no YCoordinates");
          ArrayRCP<SC> & coords = fineLevel.Get<ArrayRCP<SC> >("ZCoordinates");
          arrayOfPtrs.push_back(coords());
        }
        
        RCP<const Map> map  = fineLevel.Get<RCP<Operator> >("A", NULL/*default A*/)->getRowMap(); //FIXME for multiple DOF per node!! I need the map of the coalesce graph
        RCP<MultiVector> coordinates = MultiVectorFactory::Build(map, arrayOfPtrs, arrayOfPtrs.size());
        
        fineLevel.Set("Coordinates", coordinates);
      }
    }

    //FIXME JJH and get rid of this
      vector  = fineLevel.Get<RCP<MultiVector> >(vectorName_,MueLu::NoFactory::get());
      //std::cout << "MultiVectorTransferFactory::Build -- requesting " << vectorName_ << " from factory " << MueLu::NoFactory::get() << std::endl;
    //FIXME JJH

    RCP<Operator> transferOp = coarseLevel.Get<RCP<Operator> >(restrictionName_,restrictionFact_.get());

    //FIXME debugging output

/*
    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    RCP<const Teuchos::Comm<int> > comm = transferOp->getRowMap()->getComm();
    for (int i=0; i<comm->getSize(); ++i) {
      if (comm->getRank() == i) {
        fos->setOutputToRootOnly(comm->getRank());
        *fos << "================== pid " << comm->getRank() << ": fine level in TentativePFactory =================" << std::endl;
        fineLevel.print(*fos,Teuchos::VERB_EXTREME);
        *fos << "================= end of level description =========================" << std::endl;
      }
      comm->barrier();
    }
    fos->setOutputToRootOnly(-1);
*/
    //FIXME end of debugging output

    /*
    //FIXME ThreeLevels unit test  dies
    RCP<MultiVector> vector  = fineLevel.Get<RCP<MultiVector> >(vectorName_,restrictionFact_.get());
    RCP<Operator> transferOp = coarseLevel.Get<RCP<Operator> >(restrictionName_,restrictionFact_.get());
    */

    RCP<MultiVector> result = MultiVectorFactory::Build(transferOp->getRangeMap(),vector->getNumVectors());
    GetOStream(Runtime0,0) << "Transferring multivector \"" << vectorName_ << "\"" << std::endl;
    /*
    //FIXME DEBUGGING INFO
    if (vectorName_ == "Coordinates") {
      RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      fos->setOutputToRootOnly(-1);
      transferOp->describe(*fos,Teuchos::VERB_EXTREME);
      *fos << "===================================================================" << std::endl;
      vector->describe(*fos,Teuchos::VERB_EXTREME);
    }
    */
    //FIXME
    transferOp->apply(*vector,*result);
    //coarseLevel.Set<RCP<MultiVector> >(vectorName_,result,this); //FIXME JJH
    coarseLevel.Set<RCP<MultiVector> >(vectorName_,result,NoFactory::get()); //FIXME JJH
    
  } //Build

  //TODO do we need methods to get name of MultiVector and or Transfer Operator?

} // namespace MueLu

#define MUELU_MULTIVECTORTRANSFER_FACTORY_SHORT
#endif // MUELU_MULTIVECTORTRANSFER_FACTORY_DEF_HPP
