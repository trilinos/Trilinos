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
#ifndef MUELU_UNCOUPLEDAGGREGATIONFACTORY_KOKKOS_DECL_HPP
#define MUELU_UNCOUPLEDAGGREGATIONFACTORY_KOKKOS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_UncoupledAggregationFactory_kokkos_fwd.hpp"

#include "MueLu_Aggregates_kokkos_fwd.hpp"
#include "MueLu_AggregationAlgorithmBase_kokkos.hpp"
#include "MueLu_AggregationPhase1Algorithm_kokkos_fwd.hpp"
#include "MueLu_AggregationPhase2aAlgorithm_kokkos_fwd.hpp"
#include "MueLu_AggregationPhase2bAlgorithm_kokkos_fwd.hpp"
#include "MueLu_AggregationPhase3Algorithm_kokkos_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_IsolatedNodeAggregationAlgorithm_kokkos_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_OnePtAggregationAlgorithm_kokkos_fwd.hpp"
#include "MueLu_PreserveDirichletAggregationAlgorithm_kokkos_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

namespace MueLu {

/*!
    @class UncoupledAggregationFactory class.
    @brief Factory for building uncoupled aggregates.

    Factory for creating uncoupled aggregates from the amalgamated graph of A. The uncoupled aggregation method
    uses several aggregation phases which put together all nodes into aggregates.

    ## Aggregation phases ##
    AggregationAlgorithm | Short description
    ---------------------|------------------
    PreserveDirichletAggregationAlgorithm |  Handle Dirichlet nodes. Decide whether to drop/ignore them in the aggregation or keep them as singleton nodes.
    OnePtAggregationAlgorithm | Special handling for nodes with status ONEPT. A user can mark special nodes for singleton aggregates or a user-specified handling. This aggregation phase has to be switched on by the user if necessary (default = off).
    AggregationPhase1Algorithm | Build new aggregates
    AggregationPhase2aAlgorithm | Build aggregates of reasonable size from leftover nodes
    AggregationPhase2bAlgorithm | Add leftover nodes to existing aggregates
    AggregationPhase3Algorithm | Handle leftover nodes. Try to avoid singletons
    IsolatedNodeAggregationAlgorithm | Drop/ignore leftover nodes

    Internally, each node has a status which can be one of the following:

    Node status | Meaning
    ------------|---------
    READY       | Node is not aggregated and can be used for building a new aggregate or can be added to an existing aggregate.
    AGGREGATED  | Node is aggregated.
    IGNORED     | Node is not considered for aggregation (it may have been dropped or put into a singleton aggregate)
    BOUNDARY    | Node is a Dirichlet boundary node (with one or more Dirichlet boundary conditions).
    ONEPT       | The user forces the aggregation algorithm to treat the node as a singleton. Important: Do not forget to set aggregation: allow user-specified singletons to true! Otherwise Phase3 will just handle the ONEPT nodes and probably not build singletons

    @ingroup Aggregation

    ## Input/output of UncoupledAggregationFactory ##

    ### User parameters of UncoupledAggregationFactory ###
    Parameter | type | default | master.xml | validated | requested | description
    ----------|------|---------|:----------:|:---------:|:---------:|------------
     Graph              | Factory | null |   | * | * | Generating factory of the graph of A
     DofsPerNode        | Factory | null |   | * | * | Generating factory for variable 'DofsPerNode', usually the same as for 'Graph'
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


    The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
    The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see UncoupledAggregationFactory::GetValidParameters).<br>
    The * in the @c requested column states that the data is requested as input with all dependencies (see UncoupledAggregationFactory::DeclareInput).

    ### Variables provided by UncoupledAggregationFactory ###

    After UncoupledAggregationFactory::Build the following data is available (if requested)

    Parameter | generated by | description
    ----------|--------------|------------
    | Aggregates   | UncoupledAggregationFactory   | Container class with aggregation information. See also Aggregates.
*/

  template<class LocalOrdinal = DefaultLocalOrdinal,
           class GlobalOrdinal = DefaultGlobalOrdinal,
           class Node = DefaultNode>
  class UncoupledAggregationFactory_kokkos : public SingleLevelFactoryBase {
#undef MUELU_UNCOUPLEDAGGREGATIONFACTORY_KOKKOS_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    UncoupledAggregationFactory_kokkos();

    //! Destructor.
    virtual ~UncoupledAggregationFactory_kokkos() { }

    RCP<const ParameterList> GetValidParameterList() const;

    //@}

    //! @name Set/get methods.
    //@{

    // Options shared by all aggregation algorithms

    // deprecated
    void SetOrdering(const std::string& ordering) {
      SetParameter("aggregation: ordering", ParameterEntry(ordering));
    }
    // deprecated
    void SetMaxNeighAlreadySelected(int maxNeighAlreadySelected) {
      SetParameter("aggregation: max selected neighbors", ParameterEntry(Teuchos::as<LocalOrdinal>(maxNeighAlreadySelected))); // revalidate
    }
    // deprecated
    void SetMinNodesPerAggregate(int minNodesPerAggregate) {
      SetParameter("aggregation: min agg size", ParameterEntry(Teuchos::as<LocalOrdinal>(minNodesPerAggregate))); // revalidate
    }
    // set information about 1-node aggregates (map name and generating factory)
    void SetOnePtMapName(const std::string name, Teuchos::RCP<const FactoryBase> mapFact) {
      SetParameter("OnePt aggregate map name", ParameterEntry(std::string(name))); // revalidate
      SetFactory("OnePt aggregate map factory",mapFact);
    }

    // deprecated
    const std::string& GetOrdering() const {
      const ParameterList& pL = GetParameterList();
      return pL.get<std::string>("aggregation: ordering");
    }
    // deprecated
    int GetMaxNeighAlreadySelected() const {
      const ParameterList& pL = GetParameterList();
      return Teuchos::as<int>(pL.get<LocalOrdinal>("aggregation: max selected neighbors"));
    }
    // deprecated
    int GetMinNodesPerAggregate() const {
      const ParameterList& pL = GetParameterList();
      return Teuchos::as<int>(pL.get<LocalOrdinal>("aggregation: min agg size"));
    }

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    //! @name Build methods.
    //@{

    /*! @brief Build aggregates. */
    void Build(Level &currentLevel) const;

    //@}

    //! @name Definition methods
    //@{

    /*! @brief Append a new aggregation algorithm to list of aggregation algorithms */
    //void Append(const RCP<MueLu::AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node> > & alg);

    /*! @brief Remove all aggregation algorithms from list */
    //void ClearAggregationAlgorithms() { algos_.clear(); }
    //@}

  private:

    //! aggregation algorithms
    // will be filled in Build routine
    mutable std::vector<RCP<MueLu::AggregationAlgorithmBase_kokkos<LocalOrdinal, GlobalOrdinal, Node> > > algos_;

    //! boolean flag: definition phase
    //! if true, the aggregation algorithms still can be set and changed.
    //! if false, no change in aggregation algorithms is possible any more
    mutable bool bDefinitionPhase_;

  }; // class UncoupledAggregationFactory_kokkos

}

#define MUELU_UNCOUPLEDAGGREGATIONFACTORY_KOKKOS_SHORT
#endif // MUELU_UNCOUPLEDAGGREGATIONFACTORY_KOKKOS_DECL_HPP
