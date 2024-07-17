// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * MueLu_InterfaceAggregationAlgorithm_decl.hpp
 *
 *  Created on: Sep 18, 2012
 *      Author: Tobias Wiesner
 */

#ifndef MUELU_INTERFACEAGGREGATIONALGORITHM_DECL_HPP_
#define MUELU_INTERFACEAGGREGATIONALGORITHM_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_AggregationAlgorithmBase.hpp"
#include "MueLu_InterfaceAggregationAlgorithm_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"

#include "MueLu_LWGraph.hpp"

namespace MueLu {
/*!
  @class InterfaceAggregationAlgorithm class.
  @brief Algorithm for coarsening a graph with uncoupled aggregation.
  creates aggregates along an interface using specified root nodes.

  @ingroup Aggregation

  ### Idea ###
  The user can mark some nodes as INTERFACE to build aggregates across an interface.
  This can be very useful for certain applications. We build aggregates for nodes with
  the state INTERFACE. Then, the state is changed to AGGREGATED.
  The InterfaceAggregationAlgorithm should run before the Phase1AggregationAlgorithm.

*/

template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class InterfaceAggregationAlgorithm : public MueLu::AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_INTERFACEAGGREGATIONALGORITHM_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  InterfaceAggregationAlgorithm(RCP<const FactoryBase> const& graphFact = Teuchos::null);

  //! Destructor.
  virtual ~InterfaceAggregationAlgorithm() {}

  //@}

  //! @name Aggregation methods.
  //@{

  /*! @brief Local aggregation. */

  void BuildAggregatesNonKokkos(Teuchos::ParameterList const& params, LWGraph const& graph, Aggregates& aggregates, typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType& aggStat, LO& numNonAggregatedNodes) const;

  void BuildAggregates(const Teuchos::ParameterList& params,
                       const LWGraph_kokkos& graph,
                       Aggregates& aggregates,
                       typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                       LO& numNonAggregatedNodes) const;
  //@}

};  // class InterfaceAggregationAlgorithm

}  // namespace MueLu

#define MUELU_INTERFACEAGGREGATIONALGORITHM_SHORT
#endif /* MUELU_INTERFACEAGGREGATIONALGORITHM_DECL_HPP_ */
