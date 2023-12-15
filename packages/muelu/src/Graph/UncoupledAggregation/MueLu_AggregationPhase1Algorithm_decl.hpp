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
#ifndef MUELU_AGGREGATIONPHASE1ALGORITHM_DECL_HPP_
#define MUELU_AGGREGATIONPHASE1ALGORITHM_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_AggregationAlgorithmBase.hpp"
#include "MueLu_AggregationPhase1Algorithm_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_GraphBase.hpp"

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

  void BuildAggregates(const ParameterList& params, const GraphBase& graph, Aggregates& aggregates, std::vector<unsigned>& aggStat, LO& numNonAggregatedNodes) const;
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
