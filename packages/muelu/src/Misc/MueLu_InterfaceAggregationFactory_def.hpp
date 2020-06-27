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
#include "Xpetra_MapFactory.hpp"

namespace MueLu
{

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> InterfaceAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const
{
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Generating factory of A");
  validParamList->set<RCP<const FactoryBase>>("Aggregates", Teuchos::null, "Generating factory of the Aggregates (for block 0,0)");
  validParamList->set<RCP<const FactoryBase>>("DualNodeID2PrimalNodeID", Teuchos::null, "Generating factory of the DualNodeID2PrimalNodeID map");

  validParamList->set<LocalOrdinal>("number of DOFs per dual node", Teuchos::ScalarTraits<LocalOrdinal>::one(), "Number of DOFs per dual node");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const
{
  Input(currentLevel, "A");
  Input(currentLevel, "Aggregates");

  if (currentLevel.GetLevelID() == 0)
  {
    if(currentLevel.IsAvailable("DualNodeID2PrimalNodeID", NoFactory::get())) {
      currentLevel.DeclareInput("DualNodeID2PrimalNodeID", NoFactory::get(), this);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.IsAvailable("DualNodeID2PrimalNodeID", NoFactory::get()),
                                 Exceptions::RuntimeError,
                                 "DualNodeID2PrimalNodeID was not provided by the user on level 0!");
    }
  }
  else
  {
    Input(currentLevel, "DualNodeID2PrimalNodeID");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const
{
  using Map = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using MapFactory = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>;
  using Aggregates = Aggregates<LocalOrdinal, GlobalOrdinal, Node>;
  using Dual2Primal_type = std::map<LocalOrdinal, LocalOrdinal>;

  const char prefix[] = "MueLu::InterfaceAggregationFactory::Build: ";
  const ParameterList &pL = GetParameterList();

  FactoryMonitor m(*this, "Build", currentLevel);

  RCP<Matrix> A = Get<RCP<Matrix>>(currentLevel, "A");
  RCP<Aggregates> aggs00 = Get<RCP<Aggregates>>(currentLevel, "Aggregates");
  ArrayRCP<LocalOrdinal> vertex2AggIdInput = aggs00->GetVertex2AggId()->getDataNonConst(0);

  ArrayRCP<LocalOrdinal> procWinnerInput = aggs00->GetProcWinner()->getDataNonConst(0);

  RCP<Dual2Primal_type> lagr2dof;

  if (currentLevel.GetLevelID() == 0)
    lagr2dof = currentLevel.Get<RCP<Dual2Primal_type>>("DualNodeID2PrimalNodeID", NoFactory::get());
  else
    lagr2dof = Get<RCP<Dual2Primal_type>>(currentLevel, "DualNodeID2PrimalNodeID");

  const LocalOrdinal nDOFPerDualNode = pL.get<LocalOrdinal>("number of DOFs per dual node");

  RCP<const Map> aggRangeMap = A->getRangeMap();
  const size_t myRank = aggRangeMap->getComm()->getRank();

  LocalOrdinal globalNumDualNodes = aggRangeMap->getGlobalNumElements() / nDOFPerDualNode;
  LocalOrdinal numDualNodes = aggRangeMap->getNodeNumElements() / nDOFPerDualNode;

  TEUCHOS_TEST_FOR_EXCEPTION(numDualNodes != LocalOrdinal(lagr2dof->size()), std::runtime_error, prefix << " MueLu requires the range map and the DualNodeID2PrimalNodeID map to be compatible.");

  RCP<const Map> aggVertexMap;

  if (nDOFPerDualNode == 1)
    aggVertexMap = aggRangeMap;
  else
  {
    GlobalOrdinal indexBase = aggRangeMap->getIndexBase();
    auto comm = aggRangeMap->getComm();
    std::vector<GlobalOrdinal> myDualNodes = {};

    for (size_t i = 0; i < aggRangeMap->getNodeNumElements(); i += nDOFPerDualNode)
      myDualNodes.push_back((aggRangeMap->getGlobalElement(i) - indexBase) / nDOFPerDualNode + indexBase);

    aggVertexMap = MapFactory::Build(aggRangeMap->lib(), globalNumDualNodes, myDualNodes, indexBase, comm);
  }

  RCP<Aggregates> aggregates = rcp(new Aggregates(aggVertexMap));
  aggregates->setObjectLabel("IA");
  aggregates->AggregatesCrossProcessors(false);

  ArrayRCP<LocalOrdinal> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
  ArrayRCP<LocalOrdinal> procWinner = aggregates->GetProcWinner()->getDataNonConst(0);

  RCP<Dual2Primal_type> coarseLagr2Dof = rcp(new Dual2Primal_type());
  RCP<Dual2Primal_type> coarseDof2Lagr = rcp(new Dual2Primal_type());

  LocalOrdinal numAggId = 0;

  // Loop over the local "nodes" of the block 11
  for (LocalOrdinal localDualNodeID = 0; localDualNodeID < numDualNodes; ++localDualNodeID)
  {
    // Get local ID of the primal node associated to the current dual node
    LocalOrdinal localPrimalNodeID = (*lagr2dof)[localDualNodeID];

    // Find the aggregate of block 00 associated with the current primal node
    LocalOrdinal currentAggIdInput = vertex2AggIdInput[localPrimalNodeID];

    // Test if the current aggregate of block 00 has no associated aggregate of block 11
    if (coarseDof2Lagr->count(currentAggIdInput) == 0)
    {
      // Associate a new aggregate of block 11 to the current aggregate of block 00
      (*coarseDof2Lagr)[currentAggIdInput] = numAggId;
      (*coarseLagr2Dof)[numAggId] = currentAggIdInput;
      ++numAggId;
    }

    // Fill the block 11 aggregate information
    vertex2AggId[localDualNodeID] = (*coarseDof2Lagr)[currentAggIdInput];
    procWinner[localDualNodeID] = myRank;
  }

  aggregates->SetNumAggregates(numAggId);
  Set(currentLevel, "Aggregates", aggregates);
  Set(currentLevel, "CoarseDualNodeID2PrimalNodeID", coarseLagr2Dof);
  GetOStream(Statistics1) << aggregates->description() << std::endl;
}
} //namespace MueLu

#endif /* MUELU_INTERFACEAGGREGATIONFACTORY_DEF_HPP_ */
