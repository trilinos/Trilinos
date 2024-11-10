// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AGGREGATIONPHASE1ALGORITHM_DECL_HPP_
#define MUELU_AGGREGATIONPHASE1ALGORITHM_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_AggregationAlgorithmBase.hpp"
#include "MueLu_AggregationPhase1Algorithm_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_LWGraph.hpp"

namespace MueLu {
/*!
  @class AggregationPhase1Algorithm class.
  @brief Algorithm for coarsening a graph with uncoupled aggregation.

  @ingroup Aggregation

  ### Idea ###
  Phase 1 tries to build new aggregates which fulfill the user chosen aggregation
  criteria (i.e. minimum and maximum size of aggregates). Especially the chosen
  ordering for the input nodes may have some influence on the final aggregates.
  Phase 1 is the most important aggregation routine for building new aggregates.

  ### Parameters ###
  Parameter | Meaning
  ----------|--------
  aggregation: ordering | Ordering of graph nodes in which the nodes are processed for aggregation. The options are natural, random and graph.
  aggregation: max selected neighbors | Maximum number of neighbor nodes which have already been added to aggregates.
  aggregation: min agg size | minimum number of nodes which have to be in an aggregate.
  aggregation: max agg size | maximum allowed number of nodes in an aggregate

  ### Comments ###
  Only nodes with state READY are changed to AGGREGATED. Nodes with other states are not touched.
*/

template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class AggregationPhase1Algorithm : public MueLu::AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_AGGREGATIONPHASE1ALGORITHM_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  AggregationPhase1Algorithm(const RCP<const FactoryBase>& /* graphFact */ = Teuchos::null) {}

  //! Destructor.
  virtual ~AggregationPhase1Algorithm() {}

  //@}

  //! @name Aggregation methods.
  //@{

  /*! @brief Local aggregation. */

  void BuildAggregatesNonKokkos(const ParameterList& params, const LWGraph& graph, Aggregates& aggregates, typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType& aggStat, LO& numNonAggregatedNodes) const;

  void BuildAggregates(const Teuchos::ParameterList& params,
                       const LWGraph_kokkos& graph,
                       Aggregates& aggregates,
                       typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                       LO& numNonAggregatedNodes) const;

  void BuildAggregatesRandom(const LO maxAggSize,
                             const LWGraph_kokkos& graph,
                             Aggregates& aggregates,
                             typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                             LO& numNonAggregatedNodes) const;

  void BuildAggregatesDeterministic(const LO maxAggSize,
                                    const LWGraph_kokkos& graph,
                                    Aggregates& aggregates,
                                    typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                                    LO& numNonAggregatedNodes) const;
  //@}

  std::string description() const { return "Phase 1 (main)"; }

 private:
  /*! @brief Utility to take a list of integers and reorder them randomly (by using a local permutation).
    @param list On input, a bunch of integers. On output, the same integers in a different order
    that is determined randomly.
  */
  void RandomReorder(ArrayRCP<LO> list) const;

  /*! @brief Generate a random number in the range [min, max] */
  int RandomOrdinal(int min, int max) const;
};

}  // namespace MueLu

#define MUELU_AGGREGATIONPHASE1ALGORITHM_SHORT
#endif /* MUELU_AGGREGATIONPHASE1ALGORITHM_DECL_HPP_ */
