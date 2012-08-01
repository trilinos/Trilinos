/*
 * MueLu_VariableTransferFactory_def.hpp
 *
 *  Created on: Jul 30, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AGGSTATTRANSFERFACTORY_DEF_HPP_
#define MUELU_AGGSTATTRANSFERFACTORY_DEF_HPP_

#include "MueLu_AggStatTransferFactory_decl.hpp"

#include "MueLu_CheapAggregationAlgorithm.hpp"  // needed for NodeState enum

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  // ----------------------------------------------------------------------------------------
  // Constructor
  // ----------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  AggStatTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AggStatTransferFactory(
    std::string const & varName,
    RCP<const FactoryBase> const &genFact)
    : varName_(varName),
      genFact_(genFact)
  { }

  // ----------------------------------------------------------------------------------------
  // Destructor
  // ----------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  AggStatTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~AggStatTransferFactory() {}

  // ----------------------------------------------------------------------------------------
  // DeclareInput
  // ----------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AggStatTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    fineLevel.DeclareInput(varName_, genFact_.get(), this);
  }

  // ----------------------------------------------------------------------------------------
  // Build
  // ----------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AggStatTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & fineLevel, Level &coarseLevel) const {

    FactoryMonitor m(*this, "VariableTransferFactory", coarseLevel);

    // TODO find something smarter than distinction by variable name
    // for example one could use the underlaying type of the variable
    // Therefor we have to add the functionality to the Level class.
    // not sure we wanna do this. -> do decided later
    if(varName_ == "coarseAggStat") {
      Teuchos::ArrayRCP<NodeState> data = fineLevel.Get<Teuchos::ArrayRCP<NodeState> >(varName_,genFact_.get());
      coarseLevel.Set<Teuchos::ArrayRCP<NodeState> >(varName_, data, genFact_.get());
    }

  } //Build

} // namespace MueLu


#endif /* MUELU_AGGSTATTRANSFERFACTORY_DEF_HPP_ */
