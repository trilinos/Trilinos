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
#ifndef MUELU_STRUCTUREDAGGREGATIONFACTORY_DECL_HPP
#define MUELU_STRUCTUREDAGGREGATIONFACTORY_DECL_HPP

// #include <Xpetra_Map_fwd.hpp>
// #include <Xpetra_CrsGraph_fwd.hpp>
// #include <Xpetra_CrsGraphFactory.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_StructuredAggregationFactory_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_AggregationStructuredAlgorithm_fwd.hpp"

namespace MueLu {

/*!
    @class StructuredAggregationFactory class.
    @brief Factory for building aggregates on structured grids.

    Factory for creating aggregates from grid structure of the problem. The structured aggregation
    method can return an aggregate structure or a geometric structure used by prolongator factories.

    Internally, each node has a status which can be one of the following:

    Node status | Meaning
    ------------|---------
    READY       | Node is not aggregated and can be used for building a new aggregate or can be added to an existing aggregate.
    AGGREGATED  | Node is aggregated.
    IGNORED     | Node is not considered for aggregation (it may have been dropped or put into a singleton aggregate)
    BOUNDARY    | Node is a Dirichlet boundary node (with one or more Dirichlet boundary conditions).
    ONEPT       | The user forces the aggregation algorithm to treat the node as a singleton. Important: Do not forget to set aggregation: allow user-specified singletons to true! Otherwise Phase3 will just handle the ONEPT nodes and probably not build singletons

    @ingroup Aggregation

    ## Input/output of StructuredAggregationFactory ##

    ### User parameters of StructuredAggregationFactory ###
    Parameter | type | default | master.xml | validated | requested | description
    ----------|------|---------|:----------:|:---------:|:---------:|------------
     DofsPerNode        | Factory | null |   | * | * | Generating factory for variable 'DofsPerNode', usually the same as for 'Graph'
     OnePt aggregate map name  | string |  | | * | * | Name of input map for single node aggregates (default=''). Makes only sense if the parameter 'aggregation: allow user-specified singletons' is set to true.
     OnePt aggregate map factory | Factory | null |   | * | * | Generating factory of (DOF) map for single node aggregates.  Makes only sense if the parameter 'aggregation: allow user-specified singletons' is set to true.


    The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
    The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see StructuredAggregationFactory::GetValidParameters).<br>
    The * in the @c requested column states that the data is requested as input with all dependencies (see StructuredAggregationFactory::DeclareInput).

    ### Variables provided by StructuredAggregationFactory ###

    After StructuredAggregationFactory::Build the following data is available (if requested)

    Parameter | generated by | description
    ----------|--------------|------------
    | Aggregates   | StructuredAggregationFactory   | Container class with aggregation information. See also Aggregates.
*/

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class StructuredAggregationFactory : public SingleLevelFactoryBase {
#undef MUELU_STRUCTUREDAGGREGATIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  StructuredAggregationFactory();

  //! Destructor.
  virtual ~StructuredAggregationFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! @name Set/get methods.
  //@{
  // set information about 1-node aggregates (map name and generating factory)
  void SetOnePtMapName(const std::string name, Teuchos::RCP<const FactoryBase> mapFact) {
    SetParameter("OnePt aggregate map name", ParameterEntry(std::string(name)));  // revalidate
    SetFactory("OnePt aggregate map factory", mapFact);
  }

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*! @brief Build aggregates. */
  void Build(Level& currentLevel) const;

  //@}

 private:
  //! boolean flag: definition phase
  //! if true, the aggregation algorithms still can be set and changed.
  //! if false, no change in aggregation algorithms is possible any more
  mutable bool bDefinitionPhase_;

};  // class StructuredAggregationFactory

}  // namespace MueLu

#define MUELU_STRUCTUREDAGGREGATIONFACTORY_SHORT
#endif /* MUELU_STRUCTUREDAGGREGATIONFACTORY_DECL_HPP */
