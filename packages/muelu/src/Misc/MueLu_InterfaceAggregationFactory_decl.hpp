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
#ifndef MUELU_INTERFACEAGGREGATIONFACTORY_DECL_HPP_
#define MUELU_INTERFACEAGGREGATIONFACTORY_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

#include "MueLu_Exceptions.hpp"
namespace MueLu
{

/*!
    @class InterfaceAggregationFactory class.
    @brief Factory for building aggregates for Lagrange multipliers.

    ## Context, assumptions, and use cases ##

    This factory is intended to be used for saddle-point systems of surface-coupled problems,
    where constraints are enforced with Lagrange multipliers.
    In addition to the primal unknowns, Lagrange multipliers are considered as dual unknowns.
    The presence of Lagrange multipliers make this a primal/dual problem.

    It is assumed that each primal slave-side interface node (carrying primal unknowns) is replicated
    with a dual node carrying the dual unknowns.
    Furthermore, the number of degrees of freedom per dual node is required to constant for all dual nodes.

    ## Idea ##

    This factory will generate aggregates for the dual nodes such that they "match" the aggregates of their primal counterpart.
    Instead of performing an actual aggregation procedure on the dual nodes,
    we grep the existing primal aggregates and use a user-given mapping of dual-to-primal node IDs
    to create the dual aggregates.

    @ingroup Aggregation

    ## Input/output of InterfaceAggregationFactory ##

    ### User parameters of InterfaceAggregationFactory ###
    Parameter | type | default | master.xml | validated | requested | description
    ----------|------|---------|:----------:|:---------:|:---------:|------------
    A                       | Factory | null |   | * | * | Generating factory of the matrix A
    Aggregates              | Factory | null |   | * | * | Generating factory of the aggregates (of type "Aggregates" produced, e.g., by the UncoupledAggregationFactory)
    DualNodeID2PrimalNodeID | Factory | null |   | * | * | Generating factory of the fine dual-to-primal node mapping

    The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
    The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see InterfaceAggregationFactory::GetValidParameters).<br>
    The * in the @c requested column states that the data is requested as input with all dependencies (see InterfaceAggregationFactory::DeclareInput).

    The mapping of dual-to-primal node IDs, \c DualNodeID2PrimalNodeID, is of data type \c std::map<LocalOrdinal,LocalOrdinal>.
    The 'key' refers to the local ID of the dual node, while the 'value' represents the local ID of its primal counterpart.

    ### Variables provided by InterfaceAggregationFactory ###

    After InterfaceAggregationFactory::Build the following data is available (if requested)

    Parameter | generated by | description
    ----------|--------------|------------
    | Aggregates                    | InterfaceAggregationFactory | Aggregates of "dual nodes" carrying Lagrange multipliers in surface-coupled problems with primal and dual variables.
    | CoarseDualNodeID2PrimalNodeID | InterfaceAggregationFactory | Coarsened mapping of dual node IDs two primal node IDs.
*/

template <class Scalar = DefaultScalar,
          class LocalOrdinal = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node = DefaultNode>
class InterfaceAggregationFactory : public SingleLevelFactoryBase
{
#undef MUELU_INTERFACEAGGREGATIONFACTORY_SHORT
#include "MueLu_UseShortNamesScalar.hpp"

public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  InterfaceAggregationFactory() = default;

  //! Destructor.
  ~InterfaceAggregationFactory() {}

  RCP<const ParameterList> GetValidParameterList() const override;

  //@}

  //! Input
  //@{

  void DeclareInput(Level &currentLevel) const override;

  //@}

  //! @name Build methods.
  //@{

  /*! @brief Build aggregates. */
  void Build(Level &currentLevel) const override;

  //@}
};

} // namespace MueLu

#define MUELU_INTERFACEAGGREGATIONFACTORY_SHORT
#endif /* MUELU_INTERFACEAGGREGATIONFACTORY_DECL_HPP_ */
