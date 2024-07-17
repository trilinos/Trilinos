// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AGGREGATIONPHASE2AALGORITHM_DECL_HPP_
#define MUELU_AGGREGATIONPHASE2AALGORITHM_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_AggregationAlgorithmBase.hpp"

#include "MueLu_AggregationPhase2aAlgorithm_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_LWGraph.hpp"

namespace MueLu {
/*!
  @class AggregationPhase2aAlgorithm class.
  @brief Among unaggregated points, see if we can make a reasonable size aggregate out of it.
  @ingroup Aggregation

  ### Idea ###
  Among unaggregated points, see if we can make a reasonable size
  aggregate out of it. We do this by looking at neighbors and seeing
  how many are unaggregated and on my processor. Loosely, base the
  number of new aggregates created on the percentage of unaggregated nodes.

  ### Parameters ###
  Parameter | Meaning
  ----------|--------
  aggregation: min agg size | minimum number of nodes which have to be in an aggregate.
  aggregation: max agg size | maximum allowed number of nodes in an aggregate

  ### Comments ###
  Only nodes with state READY are changed to AGGREGATED.

*/

template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class AggregationPhase2aAlgorithm : public MueLu::AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_AGGREGATIONPHASE2AALGORITHM_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  AggregationPhase2aAlgorithm(const RCP<const FactoryBase>& /* graphFact */ = Teuchos::null) {}

  //! Destructor.
  virtual ~AggregationPhase2aAlgorithm() {}

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

  void BuildAggregatesRandom(const Teuchos::ParameterList& params,
                             const LWGraph_kokkos& graph,
                             Aggregates& aggregates,
                             typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                             LO& numNonAggregatedNodes) const;

  void BuildAggregatesDeterministic(const Teuchos::ParameterList& params,
                                    const LWGraph_kokkos& graph,
                                    Aggregates& aggregates,
                                    typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                                    LO& numNonAggregatedNodes) const;
  //@}

  std::string description() const { return "Phase 2a (secondary)"; }
};

}  // namespace MueLu

#define MUELU_AGGREGATIONPHASE2AALGORITHM_SHORT

#endif /* MUELU_AGGREGATIONPHASE2AALGORITHM_DECL_HPP_ */
