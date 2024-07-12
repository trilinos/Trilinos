// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AGGREGATIONPHASE2BALGORITHM_DECL_HPP_
#define MUELU_AGGREGATIONPHASE2BALGORITHM_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_AggregationAlgorithmBase.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_LWGraph.hpp"

#include "MueLu_AggregationPhase2bAlgorithm_fwd.hpp"

namespace MueLu {
/*!
  @class AggregationPhase2bAlgorithm class.
  @brief Add leftovers to existing aggregates
  @ingroup Aggregation

  ### Idea ###
  In phase 2b non-aggregated nodes are added to existing aggregates.
  All neighbors of the unaggregated node are checked and the corresponding
  aggregate weight is increased. The unaggregated node is added to the aggregate
  with the best weight. A simple penalty strategy makes sure that the non-aggregated
  nodes are added to different aggregates.
  The routine runs twice to cover non-aggregate nodes which have a node distance
  of two to existing aggregates. Assuming that the node distance is not greater
  than 3 (the aggregate diameter size), running the algorithm only twice should
  be sufficient.

  ### Comments ###
  Only nodes with state READY are changed to AGGREGATED. There are no aggregation criteria considered. Especially the aggregation: max agg size criterion is ignored.
  This is not a problem, since after the previous aggregation phases one should not be able to build too large aggregates.
*/

template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class AggregationPhase2bAlgorithm : public MueLu::AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_AGGREGATIONPHASE2BALGORITHM_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  AggregationPhase2bAlgorithm(const RCP<const FactoryBase>& /* graphFact */ = Teuchos::null) {}

  //! Destructor.
  virtual ~AggregationPhase2bAlgorithm() {}

  //@}

  //! @name Aggregation methods.
  //@{

  /*! @brief Local aggregation. */

  void BuildAggregatesNonKokkos(const ParameterList& params, const LWGraph& graph, Aggregates& aggregates, typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType& aggStat, LO& numNonAggregatedNodes) const;

  void BuildAggregates(const ParameterList& params,
                       const LWGraph_kokkos& graph,
                       Aggregates& aggregates,
                       typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                       LO& numNonAggregatedNodes) const;

  void BuildAggregatesRandom(const ParameterList& params,
                             const LWGraph_kokkos& graph,
                             Aggregates& aggregates,
                             typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                             LO& numNonAggregatedNodes) const;

  void BuildAggregatesDeterministic(const ParameterList& params,
                                    const LWGraph_kokkos& graph,
                                    Aggregates& aggregates,
                                    typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                                    LO& numNonAggregatedNodes) const;
  //@}

  std::string description() const { return "Phase 2b (expansion)"; }
};

}  // namespace MueLu

#define MUELU_AGGREGATIONPHASE2BALGORITHM_SHORT

#endif /* MUELU_AGGREGATIONPHASE2BALGORITHM_DECL_HPP_ */
