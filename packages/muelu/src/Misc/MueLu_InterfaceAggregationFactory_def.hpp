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
#ifndef MUELU_INTERFACEAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_INTERFACEAGGREGATIONFACTORY_DEF_HPP_

#include "MueLu_InterfaceAggregationFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu
{

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
InterfaceAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::InterfaceAggregationFactory()
{
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> InterfaceAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const
{
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase>>("A", null, "Generating factory of A");
  validParamList->set<RCP<const FactoryBase>>("Aggregates", null, "Generating factory of the Aggregates (for block 0,0)");

  validParamList->set<LocalOrdinal>("number of DOFs per node", Teuchos::ScalarTraits<LocalOrdinal>::one(), "Number of DOFs per node");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const
{
  Input(currentLevel, "A");
  Input(currentLevel, "Aggregates");

  currentLevel.DeclareInput("Lagr2Dof", NoFactory::get(), this);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const
{
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
  typedef Aggregates<LocalOrdinal, GlobalOrdinal, Node> Aggregates;
  FactoryMonitor m(*this, "Build", currentLevel);

  RCP<Matrix> A = Get<RCP<Matrix>>(currentLevel, "A");
  RCP<Aggregates> aggs00 = Get< RCP<Aggregates> >(currentLevel, "Aggregates");
  ArrayRCP<LocalOrdinal> vertex2AggIdInput = aggs00->GetVertex2AggId()->getDataNonConst(0);
  ArrayRCP<LocalOrdinal> procWinnerInput = aggs00->GetProcWinner()->getDataNonConst(0);

  RCP<std::map<GlobalOrdinal, GlobalOrdinal>> lagr2dof = currentLevel.Get<RCP<std::map<GlobalOrdinal, GlobalOrdinal>>>("Lagr2Dof", NoFactory::get());

  RCP<const Map> aggDomainMap = A->getDomainMap();
  RCP<const Map> aggRangeMap = A->getRangeMap();
  const size_t myRank = A->getRangeMap()->getComm()->getRank();

  RCP<const Map> aggVertexMap;

  if (aggDomainMap->getGlobalNumElements() == lagr2dof->size()) // True if one dof per node
    aggVertexMap = aggDomainMap;
  else
  {
    size_t n_dof_per_lag = (aggDomainMap->getGlobalNumElements()) / lagr2dof->size();
    LocalOrdinal numVertices = lagr2dof->size();
    GlobalOrdinal indexBase = aggDomainMap->getIndexBase();
    auto comm = aggDomainMap->getComm();
    std::vector<GlobalOrdinal> myVertices = {};

    for (size_t i = 0; i < aggDomainMap->getNodeNumElements(); i += n_dof_per_lag)
      myVertices.push_back((aggDomainMap->getGlobalElement(i) - indexBase) / n_dof_per_lag + indexBase);

    aggVertexMap = rcp(new const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(numVertices, myVertices, indexBase, comm));
  }

  RCP<Aggregates> aggregates = rcp(new Aggregates(aggVertexMap));
  aggregates->setObjectLabel("IA");
  aggregates->AggregatesCrossProcessors(false);

  ArrayRCP<LocalOrdinal> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
  ArrayRCP<LocalOrdinal> procWinner = aggregates->GetProcWinner()->getDataNonConst(0);

  RCP<std::map<GlobalOrdinal, GlobalOrdinal>> coarseLagr2Dof = rcp(new std::map<GlobalOrdinal, GlobalOrdinal>());
  RCP<std::map<GlobalOrdinal, GlobalOrdinal>> coarseDof2Lagr = rcp(new std::map<GlobalOrdinal, GlobalOrdinal>());

  LocalOrdinal numAggId = 0;

  LocalOrdinal n_local_lag_nodes = lagr2dof->size();

  // Loop over the local "nodes" of the block 11
  for (LocalOrdinal local_lag_node_id = 0; local_lag_node_id < n_local_lag_nodes; ++local_lag_node_id)
  {
    // Find the aggregate of block 00 associated with the current node of the block 11
    auto currentAggIdInput = vertex2AggIdInput[(*lagr2dof)[local_lag_node_id]];

    // Test if the current aggregate of block 00 has no associated aggregate of block 11
    if (coarseDof2Lagr->count(currentAggIdInput) == 0)
    {
      // Associate a new aggregate of block 11 to the current aggregate of block 00
      (*coarseDof2Lagr)[currentAggIdInput] = numAggId;
      (*coarseLagr2Dof)[numAggId] = currentAggIdInput;
      ++numAggId;
    }

    // Fill the block 11 aggregate information
    vertex2AggId[local_lag_node_id] = (*coarseDof2Lagr)[currentAggIdInput];
    procWinner[local_lag_node_id] = myRank;
  }

  aggregates->SetNumAggregates(numAggId);
  Set(currentLevel, "Aggregates", aggregates);

  currentLevel.Set<RCP<std::map<GlobalOrdinal, GlobalOrdinal>>>("CoarseLagr2Dof", coarseLagr2Dof, NoFactory::get());

  GetOStream(Statistics1) << aggregates->description() << std::endl;
}
} //namespace MueLu

#endif /* MUELU_INTERFACEAGGREGATIONFACTORY_DEF_HPP_ */
