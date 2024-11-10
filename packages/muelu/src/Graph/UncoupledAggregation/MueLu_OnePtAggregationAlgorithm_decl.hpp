// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * MueLu_OnePtAggregationAlgorithm_decl.hpp
 *
 *  Created on: Sep 18, 2012
 *      Author: Tobias Wiesner
 */

#ifndef MUELU_ONEPTAGGREGATIONALGORITHM_DECL_HPP_
#define MUELU_ONEPTAGGREGATIONALGORITHM_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_AggregationAlgorithmBase.hpp"
#include "MueLu_OnePtAggregationAlgorithm_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"

#include "MueLu_LWGraph.hpp"

namespace MueLu {
/*!
  @class OnePtAggregationAlgorithm class.
  @brief Algorithm for coarsening a graph with uncoupled aggregation.
  keep special marked nodes as singleton node aggregates over all multigrid levels

  @ingroup Aggregation

  ### Idea ###
  The user can mark some nodes as ONEPT to build some single node aggregates.
  This can be very useful for certain applications. We build single node aggregates
  for nodes with the state ONEPT. Then, the state is changed to ignored.
  The OnePtAggregationAlgorithm should run before the Phase1AggregationAlgorithm.

  ### Comments ###
  Only nodes with state ONEPT are changed to IGNORED.

*/

template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class OnePtAggregationAlgorithm : public MueLu::AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_ONEPTAGGREGATIONALGORITHM_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  OnePtAggregationAlgorithm(RCP<const FactoryBase> const& graphFact = Teuchos::null);

  //! Destructor.
  virtual ~OnePtAggregationAlgorithm() {}

  //@}

  //! @name Aggregation methods.
  //@{

  /*! @brief Local aggregation. */

  void BuildAggregatesNonKokkos(Teuchos::ParameterList const& params, LWGraph const& graph, Aggregates& aggregates, typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType& aggStat, LO& numNonAggregatedNodes) const;

  void BuildAggregates(Teuchos::ParameterList const& params,
                       LWGraph_kokkos const& graph,
                       Aggregates& aggregates,
                       typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                       LO& numNonAggregatedNodes) const;
  //@}

};  // class OnePtAggregationAlgorithm

}  // namespace MueLu

#define MUELU_ONEPTAGGREGATIONALGORITHM_SHORT
#endif /* MUELU_ONEPTAGGREGATIONALGORITHM_DECL_HPP_ */
