#ifndef MUELU_PERMUTEDTRANSFER_FACTORY_DEF_HPP
#define MUELU_PERMUTEDTRANSFER_FACTORY_DEF_HPP

#include "MueLu_PermutedTransferFactory_decl.hpp"
#include "MueLu_Utilities.hpp"

// #include <Xpetra_Operator.hpp>

#include "MueLu_Level.hpp"
// #include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  PermutedTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PermutedTransferFactory(
    RCP<FactoryBase> repartitionFact,
    RCP<FactoryBase> initialAFact,
    RCP<FactoryBase> initialTransferFact,
    TransferType     PorR)
    : repartitionFact_(repartitionFact),
      initialAFact_(initialAFact),
      initialTransferFact_(initialTransferFact),
      PorR_(PorR)
  { }

  //-----------------------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  PermutedTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~PermutedTransferFactory() {}

  //-----------------------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PermutedTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    coarseLevel.DeclareInput("A", initialAFact_.get(),this);
    if (PorR_ == MueLu::INTERPOLATION)
      coarseLevel.DeclareInput("P",initialTransferFact_.get(),this);
    else
      coarseLevel.DeclareInput("R",initialTransferFact_.get(),this);
    coarseLevel.DeclareInput("Permutation",repartitionFact_.get(),this);
  }

  //-----------------------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PermutedTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {
    RCP<Operator> permMatrix;
    try {
      permMatrix = coarseLevel.Get< RCP<Operator> >("Permutation",repartitionFact_.get());
    }
    catch(Teuchos::ExceptionBase e) {
      GetOStream(Warnings0, 0) <<  "Skipping permuting of transfers. No permutation is available for the following reason: "
                               << e.what() << std::endl;
    }

    switch (PorR_) {

      case MueLu::INTERPOLATION:
        {
        RCP<Operator> originalP = coarseLevel.Get< RCP<Operator> >("P",initialTransferFact_.get());
        RCP<Operator> permutedP = Utils::TwoMatrixMultiply(originalP,false,permMatrix,false); //P*perm
        coarseLevel.Set< RCP<Operator> >("P",permutedP,this);
        }
        break;

      case MueLu::RESTRICTION:
        {
        //TODO how do we handle implicitly transposed restriction operators?
        RCP<Operator> originalR = coarseLevel.Get< RCP<Operator> >("R",initialTransferFact_.get());
        RCP<Operator> permutedR = Utils::TwoMatrixMultiply(permMatrix,true,originalR,false); //transpose(perm) * R
        coarseLevel.Set< RCP<Operator> >("R",permutedR,this);
        }
        break;
    } //switch

  } //Build

} // namespace MueLu

#define MUELU_PERMUTEDTRANSFER_FACTORY_SHORT
#endif // MUELU_PERMUTEDTRANSFER_FACTORY_DEF_HPP
