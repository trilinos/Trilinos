// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_DECL_HPP_
#define MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_AggregationAlgorithmBase.hpp"
#include "MueLu_PreserveDirichletAggregationAlgorithm_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_LWGraph.hpp"

namespace MueLu {
/*!
  @class PreserveDirichletAggregationAlgorithm class.
  @brief Builds one-to-one aggregates for all Dirichlet boundary nodes. For some applications this might
         be necessary. (default = off)

  @ingroup Aggregation

  ### Idea ###
  Handles Dirichlet boundary nodes with the state Boundary.
  Depending on the boolean parameter "aggregation: preserve Dirichlet points" one-to-one aggregates
  with singleton nodes are built for all Dirichlet boundary nodes or the aggregates are just
  ignored (default behavior). The state of all boundary nodes (state = Boundary)
  is set to ignored. That means, that these nodes are not considered for further
  aggregation in the later aggregation phases.

  ### Parameters ###
  Parameter | Meaning
  ----------|--------
  aggregation: preserve Dirichlet points | Boolean parameter stating whether Dirichlet boundary nodes shall be aggregated in singleton aggregates (default: false).

  ### Comments ###
  Only nodes with state BOUNDARY are changed to IGNORED. No other nodes are touched.
*/

template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class PreserveDirichletAggregationAlgorithm : public MueLu::AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  PreserveDirichletAggregationAlgorithm(const RCP<const FactoryBase>& /* graphFact */ = Teuchos::null) {}

  //! Destructor.
  virtual ~PreserveDirichletAggregationAlgorithm() {}

  //@}

  //! @name Aggregation methods.
  //@{

  /*! @brief Local aggregation. */

  void BuildAggregatesNonKokkos(const Teuchos::ParameterList& params, const LWGraph& graph, Aggregates& aggregates, typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType& aggStat, LO& numNonAggregatedNodes) const;

  void BuildAggregates(const Teuchos::ParameterList& params,
                       const LWGraph_kokkos& graph,
                       Aggregates& aggregates,
                       typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                       LO& numNonAggregatedNodes) const;
  //@}

  std::string description() const { return "Phase - (Dirichlet)"; }

};  // class PreserveDirichletAggregationAlgorithm

}  // namespace MueLu

#define MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_SHORT

#endif /* MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_DECL_HPP_ */
