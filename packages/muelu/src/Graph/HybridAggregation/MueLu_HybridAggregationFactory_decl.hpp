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
#ifndef MUELU_HYBRIDAGGREGATIONFACTORY_DECL_HPP_
#define MUELU_HYBRIDAGGREGATIONFACTORY_DECL_HPP_

#include <Xpetra_Map_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_HybridAggregationFactory_fwd.hpp"

#include "MueLu_AggregationAlgorithmBase.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_GraphBase_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Exceptions.hpp"

// Uncoupled Agg
#include "MueLu_InterfaceAggregationAlgorithm_fwd.hpp"
#include "MueLu_OnePtAggregationAlgorithm_fwd.hpp"
#include "MueLu_PreserveDirichletAggregationAlgorithm_fwd.hpp"

#include "MueLu_AggregationPhase1Algorithm_fwd.hpp"
#include "MueLu_AggregationPhase2aAlgorithm_fwd.hpp"
#include "MueLu_AggregationPhase2bAlgorithm_fwd.hpp"
#include "MueLu_AggregationPhase3Algorithm_fwd.hpp"

// Structured Agg
#include "MueLu_AggregationStructuredAlgorithm_fwd.hpp"

namespace MueLu {

/*!
    @class HybridAggregationFactory class.
    @brief Factory for building aggregates on meshes partly structured and partly unstructured.

    Factory for creating aggregates from partly structured grids. The hybrid aggregation method
    returns an aggregate structure used by prolongator factories.

    Internally, each node has a status which can be one of the following:

    Node status | Meaning
    ------------|---------
    READY       | Node is not aggregated and can be used for building a new aggregate or can be added to an existing aggregate.
    AGGREGATED  | Node is aggregated.
    IGNORED     | Node is not considered for aggregation (it may have been dropped or put into a singleton aggregate)
    BOUNDARY    | Node is a Dirichlet boundary node (with one or more Dirichlet boundary conditions).
    ONEPT       | The user forces the aggregation algorithm to treat the node as a singleton. Important: Do not forget to set aggregation: allow user-specified singletons to true! Otherwise Phase3 will just handle the ONEPT nodes and probably not build singletons

    @ingroup Aggregation

    ## Input/output of HybridAggregationFactory ##

    ### User parameters of HybridAggregationFactory ###
    Parameter | type | default | master.xml | validated | requested | description
    ----------|------|---------|:----------:|:---------:|:---------:|------------
     Graph    | Factory | null |            |   *       |    *      | Generating factory for variable 'Graph'
     DofsPerNode | Factory | null |   | * | * | Generating factory for variable 'DofsPerNode', usually the same as for 'Graph'
     OnePt aggregate map name  | string |  | | * | * | Name of input map for single node aggregates (default=''). Makes only sense if the parameter 'aggregation: allow user-specified singletons' is set to true.
     OnePt aggregate map factory | Factory | null |   | * | * | Generating factory of (DOF) map for single node aggregates.  Makes only sense if the parameter 'aggregation: allow user-specified singletons' is set to true.
     aggregation: max agg size | int | see master.xml | * | * |  | Maximum number of nodes per aggregate.
     aggregation: min agg size | int | see master.xml | * | * |  | Minimum number of nodes necessary to build a new aggregate.
     aggregation: max selected neighbors | int | see master.xml | * | * |  | Maximum number of neighbor nodes already in aggregate (needed in Phase1)
     aggregation: ordering | string | "natural" | * | * |  | Ordering of node aggregation (can be either "natural", "graph" or "random").
     aggregation: enable phase 1 | bool | true | * | * |   |Turn on/off phase 1 aggregation
     aggregation: enable phase 2a | bool | true | * | * |  |Turn on/off phase 2a aggregation
     aggregation: enable phase 2b | bool | true | * | * |  |Turn on/off phase 2b aggregation
     aggregation: enable phase 3 | bool | true | * | * |   |Turn on/off phase 3 aggregation
     aggregation: preserve Dirichlet points | bool | false | * | * |   | preserve Dirichlet points as singleton nodes (default=false, i.e., drop Dirichlet nodes during aggregation)
     aggregation: allow user-specified singletons | bool | false | * | * |  | Turn on/off OnePtAggregationAlgorithm (default=false)
     aggregationRegionType| Factory | null |  | * | * | Factory generating the type of aggregation to use on the region


    The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
    The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see HybridAggregationFactory::GetValidParameters).<br>
    The * in the @c requested column states that the data is requested as input with all dependencies (see HybridAggregationFactory::DeclareInput).

    ### Variables provided by HybridAggregationFactory ###

    After HybridAggregationFactory::Build the following data is available (if requested)

    Parameter | generated by | description
    ----------|--------------|------------
    | Aggregates   | HybridAggregationFactory   | Container class with aggregation information. See also Aggregates.
*/

template <class LocalOrdinal, class GlobalOrdinal, class Node>
class HybridAggregationFactory : public SingleLevelFactoryBase {
#undef MUELU_HYBRIDAGGREGATIONFACTORY_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  HybridAggregationFactory();

  //! Destructor.
  virtual ~HybridAggregationFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*! @brief Build aggregates. */
  void Build(Level& currentLevel) const;

  /*! @brief Specifically build aggregates along interfaces */
  void BuildInterfaceAggregates(Level& currentLevel, RCP<Aggregates> aggregates,
                                std::vector<unsigned>& aggStat, LO& numNonAggregatedNodes,
                                Array<LO> coarseRate) const;

  //@}

 private:
  //! aggregation algorithms
  // will be filled in Build routine
  mutable std::vector<RCP<MueLu::AggregationAlgorithmBase<LO, GO, Node> > > algos_;

  //! boolean flag: definition phase
  //! if true, the aggregation algorithms still can be set and changed.
  //! if false, no change in aggregation algorithms is possible any more
  mutable bool bDefinitionPhase_;

};  // class HybridAggregationFactory

}  // namespace MueLu

#define MUELU_HYBRIDAGGREGATIONFACTORY_SHORT
#endif /* MUELU_HYBRIDAGGREGATIONFACTORY_DECL_HPP_ */
