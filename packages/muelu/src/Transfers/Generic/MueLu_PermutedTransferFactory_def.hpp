#ifndef MUELU_PERMUTEDTRANSFER_FACTORY_DEF_HPP
#define MUELU_PERMUTEDTRANSFER_FACTORY_DEF_HPP

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

#include "MueLu_PermutedTransferFactory_decl.hpp"
#include "MueLu_Utilities.hpp"

#include <Xpetra_Operator.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  PermutedTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PermutedTransferFactory(
    RCP<const FactoryBase> repartitionFact,
    RCP<const FactoryBase> initialAFact,
    RCP<const FactoryBase> initialTransferFact,
    TransferType     PorR,
    RCP<const FactoryBase> nullspaceFact,
    RCP<const FactoryBase> coordinateFact )
    : repartitionFact_(repartitionFact),
      initialAFact_(initialAFact),
      initialTransferFact_(initialTransferFact),
      PorR_(PorR),
      nullspaceFact_(nullspaceFact),
      coordinateFact_(coordinateFact)
  { }

  //-----------------------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  PermutedTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~PermutedTransferFactory() {}

  //-----------------------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PermutedTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {

    FactoryMonitor m(*this, "PermutedTransferFactory", coarseLevel);

    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    fos->setOutputToRootOnly(-1);
    GetOStream(Warnings0, 0) <<  "** In PermutedTransferFactory::DeclareInput **" << std::endl;
    //coarseLevel.print(*fos,Teuchos::VERB_EXTREME);


    coarseLevel.DeclareInput("A", initialAFact_.get(),this);
    if (PorR_ == MueLu::INTERPOLATION)
      coarseLevel.DeclareInput("P",initialTransferFact_.get(),this);
    else
      coarseLevel.DeclareInput("R",initialTransferFact_.get(),this);
    coarseLevel.DeclareInput("Permutation",repartitionFact_.get(),this);
  }

  //-----------------------------------------------------------------------------------------------------

  // TODO 1) Need to pass in to ctor generating factories for nullspace (TentativePFactory) and coordinates
  // TODO (MultiVectorTransferFactory).  DeclareInput should also call the DeclareInputs of these two guys.
  // TODO 2) Also add interface (ala RAPFactory) to register additional data/generating factories that must
  // TODO be permuted.  Call those DeclareInputs, as well?

  // partially done

  //-----------------------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PermutedTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {
    RCP<Operator> permMatrix;
    try {
      permMatrix = coarseLevel.Get< RCP<Operator> >("Permutation",repartitionFact_.get());
    }
    catch(Teuchos::ExceptionBase e) {
      std::string gridTransferType;
      GetOStream(Warnings0, 0) <<  "Skipping permuting of "
                               << ((PorR_ == MueLu::INTERPOLATION) ? "prolongator" : "restriction")
                               << ".  No permutation is available for the following reason:"
      << std::endl << e.what() << std::endl;
    }

    switch (PorR_) {

      case MueLu::INTERPOLATION:
        {
          RCP<Operator> originalP = coarseLevel.Get< RCP<Operator> >("P",initialTransferFact_.get());
          if (permMatrix != Teuchos::null) {
            GetOStream(Runtime0, 0) <<  "Permuting prolongator." << std::endl;
            RCP<Operator> permutedP = Utils::TwoMatrixMultiply(originalP,false,permMatrix,true); //P*transpose(perm)
            coarseLevel.Set< RCP<Operator> >("P",permutedP,this);
          } else {
            coarseLevel.Set< RCP<Operator> >("P",originalP,this);
          }
        }
        break;

      case MueLu::RESTRICTION:
        {
          //TODO how do we handle implicitly transposed restriction operators?
          RCP<Operator> originalR = coarseLevel.Get< RCP<Operator> >("R",initialTransferFact_.get());
          if (permMatrix != Teuchos::null) {
            GetOStream(Runtime0, 0) <<  "Permuting restriction." << std::endl;
            RCP<Operator> permutedR = Utils::TwoMatrixMultiply(permMatrix,false,originalR,false); //perm * R
            coarseLevel.Set< RCP<Operator> >("R",permutedR,this);
            //if (coarseLevel.IsAvailable("Coordinates",coordinateFact_.get()))  //FIXME JJH
            if (coarseLevel.IsAvailable("Coordinates")) //FIXME JJH
            {
              GetOStream(Runtime0, 0) <<  "Permuting coordinates." << std::endl;
              //RCP<MultiVector> coords  = coarseLevel.Get< RCP<MultiVector> >("Coordinates",coordinateFact_.get()); //FIXME JJH
              RCP<MultiVector> coords  = coarseLevel.Get< RCP<MultiVector> >("Coordinates"); //FIXME JJH
              RCP<MultiVector> permutedCoords  = MultiVectorFactory::Build(permMatrix->getRangeMap(),1);
              permMatrix->apply(*coords,*permutedCoords,Teuchos::NO_TRANS,1,0);
              coarseLevel.Set< RCP<MultiVector> >("Coordinates",permutedCoords); //FIXME JJH no generating factory specified
            }
            if (coarseLevel.IsAvailable("Nullspace")) {
              GetOStream(Runtime0, 0) <<  "Permuting nullspace." << std::endl;
              //RCP<MultiVector> nullspace  = coarseLevel.Get< RCP<MultiVector> >("Nullspace",nullspaceFact_.get());
              RCP<MultiVector> nullspace  = coarseLevel.Get< RCP<MultiVector> >("Nullspace");
              RCP<MultiVector> permutedNullspace  = MultiVectorFactory::Build(permMatrix->getRangeMap(),nullspace->getNumVectors());
              permMatrix->apply(*nullspace,*permutedNullspace,Teuchos::NO_TRANS,1,0);
              coarseLevel.Set< RCP<MultiVector> >("Nullspace",permutedNullspace); //FIXME no generating factory specified
            }
          } else {
            coarseLevel.Set< RCP<Operator> >("R",originalR,this);
          }
        }
        break;
    } //switch

  } //Build

} // namespace MueLu

#define MUELU_PERMUTEDTRANSFER_FACTORY_SHORT
#endif // MUELU_PERMUTEDTRANSFER_FACTORY_DEF_HPP
