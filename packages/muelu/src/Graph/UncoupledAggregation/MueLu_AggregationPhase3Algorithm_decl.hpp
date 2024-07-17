// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AGGREGATIONPHASE3ALGORITHM_DECL_HPP_
#define MUELU_AGGREGATIONPHASE3ALGORITHM_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_AggregationAlgorithmBase.hpp"
#include "MueLu_AggregationPhase3Algorithm_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_LWGraph.hpp"

namespace MueLu {
/*!
  @class AggregationPhase3Algorithm class.
  @brief Handle leftover nodes. Try to avoid singleton nodes
  @ingroup Aggregation

  ### Idea ###
  In phase 3 we try to stick unaggregated nodes into a neighboring aggregate.
  We try to avoid singletons: we first try to build a new aggregate containing
  all neighboring non-aggregated nodes. If we cannot build a new aggregate,
  we add the non-aggregated node to the first adjacent aggregate.
  Only if there is no adjacent aggregate, we create a singleton node aggregate.

  ### Comments ###
  Only nodes with state READY are changed to AGGREGATED.

*/

template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class AggregationPhase3Algorithm : public MueLu::AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_AGGREGATIONPHASE3ALGORITHM_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  AggregationPhase3Algorithm(const RCP<const FactoryBase>& /* graphFact */ = Teuchos::null) {}

  //! Destructor.
  virtual ~AggregationPhase3Algorithm() {}

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
  //@}

  std::string description() const { return "Phase 3 (cleanup)"; }
};

}  // namespace MueLu

#define MUELU_AGGREGATIONPHASE3ALGORITHM_SHORT

#endif /* MUELU_AGGREGATIONPHASE3ALGORITHM_DECL_HPP_ */
