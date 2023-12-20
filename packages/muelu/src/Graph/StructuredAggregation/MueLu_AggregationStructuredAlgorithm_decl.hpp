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
#ifndef MUELU_AGGREGATIONSTRUCTUREDALGORITHM_DECL_HPP_
#define MUELU_AGGREGATIONSTRUCTUREDALGORITHM_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_AggregationAlgorithmBase.hpp"
#include "MueLu_AggregationStructuredAlgorithm_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_IndexManager_fwd.hpp"
#include "MueLu_GraphBase.hpp"

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

  void BuildAggregates(const Teuchos::ParameterList& params, const GraphBase& graph,
                       Aggregates& aggregates, std::vector<unsigned>& aggStat,
                       LO& numNonAggregatedNodes) const;

  /*! @brief Local aggregation. */

  void BuildGraph(const GraphBase& graph, RCP<IndexManager>& geoData, const LO dofsPerNode,
                  RCP<CrsGraph>& myGraph, RCP<const Map>& coarseCoordinatesFineMap,
                  RCP<const Map>& coarseCoordinatesMap) const;
  //@}

  std::string description() const { return "Aggretation: structured algorithm"; }

 private:
  void ComputeGraphDataConstant(const GraphBase& graph, RCP<IndexManager>& geoData,
                                const LO dofsPerNode, const int numInterpolationPoints,
                                ArrayRCP<size_t>& nnzOnRow, Array<size_t>& rowPtr,
                                Array<LO>& colIndex) const;

  void ComputeGraphDataLinear(const GraphBase& graph, RCP<IndexManager>& geoData,
                              const LO dofsPerNode, const int numInterpolationPoints,
                              ArrayRCP<size_t>& nnzOnRow, Array<size_t>& rowPtr,
                              Array<LO>& colIndex) const;
};

}  // namespace MueLu

#define MUELU_AGGREGATIONSTRUCTUREDALGORITHM_SHORT
#endif /* MUELU_AGGREGATIONSTRUCTUREDALGORITHM_DECL_HPP_ */
