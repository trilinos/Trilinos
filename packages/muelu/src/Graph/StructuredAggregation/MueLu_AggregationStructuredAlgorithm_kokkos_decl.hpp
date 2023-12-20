// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_AGGREGATIONSTRUCTUREDALGORITHM_KOKKOS_DECL_HPP
#define MUELU_AGGREGATIONSTRUCTUREDALGORITHM_KOKKOS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_AggregationAlgorithmBase_kokkos.hpp"
#include "MueLu_AggregationStructuredAlgorithm_kokkos_fwd.hpp"

#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_IndexManager_kokkos_fwd.hpp"
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
  All the parameters needed are passed to this class by the StructuredAggregationFactory class.
*/

template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class AggregationStructuredAlgorithm_kokkos : public MueLu::AggregationAlgorithmBase_kokkos<LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_AGGREGATIONSTRUCTUREDALGORITHM_KOKKOS_SHORT
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
  AggregationStructuredAlgorithm_kokkos() {}

  //! Destructor.
  virtual ~AggregationStructuredAlgorithm_kokkos() {}

  //@}

  //! @name Aggregation methods.
  //@{

  /*! @brief Build aggregates object. */

  void BuildAggregates(const Teuchos::ParameterList& params,
                       const LWGraph_kokkos& graph,
                       Aggregates& aggregates,
                       Kokkos::View<unsigned*, device_type>& aggStat,
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

};  // class AggregationStructuredAlgorithm_kokkos

}  // namespace MueLu

#define MUELU_AGGREGATIONSTRUCTUREDALGORITHM_KOKKOS_SHORT
#endif /* MUELU_AGGREGATIONSTRUCTUREDALGORITHM_DECL_HPP_ */
