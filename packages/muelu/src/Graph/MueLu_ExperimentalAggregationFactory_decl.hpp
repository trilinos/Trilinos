/*
 * MueLu_ExperimentalAggregationFactory_decl.hpp
 *
 *  Created on: Jul 26, 2012
 *      Author: wiesner
 */

#ifndef MUELU_EXPERIMENTALAGGREGATIONFACTORY_DECL_HPP_
#define MUELU_EXPERIMENTALAGGREGATIONFACTORY_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_ExperimentalAggregationFactory_fwd.hpp"

#include "MueLu_CheapAggregationAlgorithm_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"
#include "MueLu_Graph_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {
template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
class ExperimentalAggregationFactory : public SingleLevelFactoryBase {
#undef MUELU_EXPERIMENTALAGGREGATIONFACTORY_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  ExperimentalAggregationFactory(RCP<const FactoryBase> graphFact = Teuchos::null);

  //! Destructor.
  virtual ~ExperimentalAggregationFactory() { }

  //@}

  //! @name Set/get methods.
  //@{

  // Options algo1
  //void SetOrdering(Ordering ordering) { algo1_->SetOrdering(ordering); }
  void SetMaxNeighAlreadySelected(int maxNeighAlreadySelected) { /*algo1_->SetMaxNeighAlreadySelected(maxNeighAlreadySelected);*/ }
  //Ordering GetOrdering() const { return algo1_->GetOrdering(); }
  int GetMaxNeighAlreadySelected() const { return algo1_->GetMaxNeighAlreadySelected(); }

  // Options shared algo1 and algo2
  void SetMinNodesPerAggregate(int minNodesPerAggregate) { /*algo1_->SetMinNodesPerAggregate(minNodesPerAggregate);*/ }
  int GetMinNodesPerAggregate() const { return 0; /* algo1_->GetMinNodesPerAggregate();*/ }

  //@}

  //! Input
  //@{

  void DeclareInput(Level &currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*! @brief Build aggregates. */
  void Build(Level &currentLevel) const;

  //@}

private:

  //! Graph Factory
  RCP<const FactoryBase> graphFact_;

  //! Algorithms
  RCP<MueLu::CheapAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > algo1_;


}; // class ExperimentalAggregationFactory

}

#define MUELU_EXPERIMENTALAGGREGATIONFACTORY_SHORT
#endif /* MUELU_EXPERIMENTALAGGREGATIONFACTORY_DECL_HPP_ */
