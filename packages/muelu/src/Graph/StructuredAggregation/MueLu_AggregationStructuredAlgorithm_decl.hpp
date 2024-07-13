// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AGGREGATIONSTRUCTUREDALGORITHM_DECL_HPP_
#define MUELU_AGGREGATIONSTRUCTUREDALGORITHM_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_AggregationAlgorithmBase.hpp"
#include "MueLu_AggregationStructuredAlgorithm_fwd.hpp"

#include "Xpetra_Vector.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_IndexManager_fwd.hpp"
#include "MueLu_IndexManager_kokkos_fwd.hpp"
#include "MueLu_LWGraph.hpp"
#include "MueLu_LWGraph_kokkos.hpp"

namespace MueLu {
/*!
  @class AggregationStructuredAlgorithm class.
  @brief Algorithm for coarsening a graph with structured aggregation.

  @ingroup Aggregation

  ### Idea ###
  Use the logical indexing of the mesh to obtain a very regular aggregation pattern and maintain
  lines and planes of the problem as they might be useful to the smoother.
  This algorithms is also very easy to parallelize on node due to its very regular and predictible
  memory access patern.

  ### Parameters ###
  Parameter | Meaning
  ----------|--------
  aggregation: coarsen | describe the coarsening rate to be used in each direction
*/

template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class AggregationStructuredAlgorithm : public MueLu::AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_AGGREGATIONSTRUCTUREDALGORITHM_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  using local_graph_type       = typename LWGraph_kokkos::local_graph_type;
  using non_const_row_map_type = typename local_graph_type::row_map_type::non_const_type;
  using size_type              = typename local_graph_type::size_type;
  using entries_type           = typename local_graph_type::entries_type;
  using device_type            = typename local_graph_type::device_type;
  using execution_space        = typename local_graph_type::device_type::execution_space;
  using memory_space           = typename local_graph_type::device_type::memory_space;

  using LOVectorView      = decltype(std::declval<LOVector>().getDeviceLocalView(Xpetra::Access::ReadWrite));
  using constIntTupleView = typename Kokkos::View<const int[3], device_type>;
  using constLOTupleView  = typename Kokkos::View<const LO[3], device_type>;
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  AggregationStructuredAlgorithm(const RCP<const FactoryBase>& /* graphFact */ = Teuchos::null) {}

  //! Destructor.
  virtual ~AggregationStructuredAlgorithm() {}

  //@}

  //! @name Aggregation methods.
  //@{

  /*! @brief Local aggregation. */

  void BuildAggregatesNonKokkos(const Teuchos::ParameterList& params, const LWGraph& graph,
                                Aggregates& aggregates,
                                typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType& aggStat,
                                LO& numNonAggregatedNodes) const;

  /*! @brief Local aggregation. */

  void BuildGraphOnHost(const LWGraph& graph, RCP<IndexManager>& geoData, const LO dofsPerNode,
                        RCP<CrsGraph>& myGraph, RCP<const Map>& coarseCoordinatesFineMap,
                        RCP<const Map>& coarseCoordinatesMap) const;

  /*! @brief Build aggregates object. */

  void BuildAggregates(const Teuchos::ParameterList& params,
                       const LWGraph_kokkos& graph,
                       Aggregates& aggregates,
                       typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                       LO& numNonAggregatedNodes) const;

  /*! @brief Build a CrsGraph instead of aggregates. */

  void BuildGraph(const LWGraph_kokkos& graph,
                  RCP<IndexManager_kokkos>& geoData,
                  const LO dofsPerNode,
                  RCP<CrsGraph>& myGraph) const;
  //@}

  std::string description() const { return "Aggretation: structured algorithm"; }

  struct fillAggregatesFunctor {
    IndexManager_kokkos geoData_;
    const int myRank_;
    Kokkos::View<unsigned*, device_type> aggStat_;
    LOVectorView vertex2AggID_;
    LOVectorView procWinner_;

    fillAggregatesFunctor(RCP<IndexManager_kokkos> geoData,
                          const int myRank,
                          Kokkos::View<unsigned*, device_type> aggStat,
                          LOVectorView vertex2AggID,
                          LOVectorView procWinner);

    KOKKOS_INLINE_FUNCTION
    void operator()(const LO nodeIdx, LO& lNumAggregatedNodes) const;

  };  // struct fillAggregatesFunctor

  struct computeGraphDataConstantFunctor {
    IndexManager_kokkos geoData_;
    const int numGhostedNodes_;
    const LO dofsPerNode_;
    constIntTupleView coarseRate_;
    constIntTupleView endRate_;
    constLOTupleView lFineNodesPerDir_;
    non_const_row_map_type rowPtr_;
    entries_type colIndex_;

    computeGraphDataConstantFunctor(RCP<IndexManager_kokkos> geoData,
                                    const LO numGhostedNodes, const LO dofsPerNode,
                                    constIntTupleView coarseRate, constIntTupleView endRate,
                                    constLOTupleView lFineNodesPerDir,
                                    non_const_row_map_type rowPtr, entries_type colIndex);

    KOKKOS_INLINE_FUNCTION
    void operator()(const LO nodeIdx) const;

  };  // struct computeGraphDataConstantFunctor

  struct computeGraphRowPtrFunctor {
    IndexManager_kokkos geoData_;
    const LO dofsPerNode_;
    const int numInterpolationPoints_;
    const LO numLocalRows_;
    constIntTupleView coarseRate_;
    constLOTupleView lFineNodesPerDir_;
    non_const_row_map_type rowPtr_;

    computeGraphRowPtrFunctor(RCP<IndexManager_kokkos> geoData,
                              const LO dofsPerNode,
                              const int numInterpolationPoints, const LO numLocalRows,
                              constIntTupleView coarseRate, constLOTupleView lFineNodesPerDir,
                              non_const_row_map_type rowPtr);

    KOKKOS_INLINE_FUNCTION
    void operator()(const LO rowIdx, GO& update, const bool final) const;
  };  // struct computeGraphRowPtrFunctor

  struct computeGraphDataLinearFunctor {
    IndexManager_kokkos geoData_;
    const int numDimensions_;
    const int numGhostedNodes_;
    const LO dofsPerNode_;
    const int numInterpolationPoints_;
    constIntTupleView coarseRate_;
    constIntTupleView endRate_;
    constLOTupleView lFineNodesPerDir_;
    constLOTupleView ghostedNodesPerDir_;
    non_const_row_map_type rowPtr_;
    entries_type colIndex_;

    computeGraphDataLinearFunctor(RCP<IndexManager_kokkos> geoData,
                                  const int numDimensions,
                                  const LO numGhostedNodes, const LO dofsPerNode,
                                  const int numInterpolationPoints,
                                  constIntTupleView coarseRate, constIntTupleView endRate,
                                  constLOTupleView lFineNodesPerDir,
                                  constLOTupleView ghostedNodesPerDir,
                                  non_const_row_map_type rowPtr, entries_type colIndex);

    KOKKOS_INLINE_FUNCTION
    void operator()(const LO nodeIdx) const;

  };  // struct computeGraphDataLinearFunctor

 private:
  void ComputeGraphDataConstant(const LWGraph& graph, RCP<IndexManager>& geoData,
                                const LO dofsPerNode, const int numInterpolationPoints,
                                ArrayRCP<size_t>& nnzOnRow, Array<size_t>& rowPtr,
                                Array<LO>& colIndex) const;

  void ComputeGraphDataLinear(const LWGraph& graph, RCP<IndexManager>& geoData,
                              const LO dofsPerNode, const int numInterpolationPoints,
                              ArrayRCP<size_t>& nnzOnRow, Array<size_t>& rowPtr,
                              Array<LO>& colIndex) const;
};

}  // namespace MueLu

#define MUELU_AGGREGATIONSTRUCTUREDALGORITHM_SHORT
#endif /* MUELU_AGGREGATIONSTRUCTUREDALGORITHM_DECL_HPP_ */
