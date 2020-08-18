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

#include "Xpetra_Map.hpp"
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
} // GetValidParameterList()

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
} // DeclareInput

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const
{
  using Map = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using MapFactory = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>;
  using Aggregates = Aggregates<LocalOrdinal, GlobalOrdinal, Node>;
  using Dual2Primal_type = std::map<LocalOrdinal, LocalOrdinal>;

  const char prefix[] = "MueLu::InterfaceAggregationFactory::Build: ";

  FactoryMonitor m(*this, "Build", currentLevel);

  // Access the input data
  const ParameterList &pL = GetParameterList();
  RCP<const Matrix> A = Get<RCP<Matrix>>(currentLevel, "A");
  const LocalOrdinal numDofsPerDualNode = pL.get<LocalOrdinal>("number of DOFs per dual node");
  RCP<const Aggregates> primalAggregates = Get<RCP<Aggregates>>(currentLevel, "Aggregates");
  ArrayRCP<const LocalOrdinal> primalVertex2AggId = primalAggregates->GetVertex2AggId()->getData(0);

  // Get the user-prescribed mapping of dual to primal node IDs
  RCP<Dual2Primal_type> mapNodesDualToPrimal;
  if (currentLevel.GetLevelID() == 0)
    mapNodesDualToPrimal = currentLevel.Get<RCP<Dual2Primal_type>>("DualNodeID2PrimalNodeID", NoFactory::get());
  else
    mapNodesDualToPrimal = Get<RCP<Dual2Primal_type>>(currentLevel, "DualNodeID2PrimalNodeID");

  RCP<const Map> operatorRangeMap = A->getRangeMap();
  const size_t myRank = operatorRangeMap->getComm()->getRank();

  LocalOrdinal globalNumDualNodes = operatorRangeMap->getGlobalNumElements() / numDofsPerDualNode;
  LocalOrdinal localNumDualNodes = operatorRangeMap->getNodeNumElements() / numDofsPerDualNode;

  TEUCHOS_TEST_FOR_EXCEPTION(localNumDualNodes != Teuchos::as<LocalOrdinal>(mapNodesDualToPrimal->size()),
      std::runtime_error, prefix << " MueLu requires the range map and the DualNodeID2PrimalNodeID map to be compatible.");

  RCP<const Map> dualNodeMap = Teuchos::null;
  if (numDofsPerDualNode == 1)
    dualNodeMap = operatorRangeMap;
  else
  {
    GlobalOrdinal indexBase = operatorRangeMap->getIndexBase();
    auto comm = operatorRangeMap->getComm();
    std::vector<GlobalOrdinal> myDualNodes = {};

    for (size_t i = 0; i < operatorRangeMap->getNodeNumElements(); i += numDofsPerDualNode)
      myDualNodes.push_back((operatorRangeMap->getGlobalElement(i) - indexBase) / numDofsPerDualNode + indexBase);

    dualNodeMap = MapFactory::Build(operatorRangeMap->lib(), globalNumDualNodes, myDualNodes, indexBase, comm);
  }

  RCP<Aggregates> dualAggregates = rcp(new Aggregates(dualNodeMap));
  dualAggregates->setObjectLabel("InterfaceAggregation");

  // Copy setting from primal aggregates, as we copy the interface part of primal aggregates anyways
  dualAggregates->AggregatesCrossProcessors(primalAggregates->AggregatesCrossProcessors());

  ArrayRCP<LocalOrdinal> dualVertex2AggId = dualAggregates->GetVertex2AggId()->getDataNonConst(0);
  ArrayRCP<LocalOrdinal> dualProcWinner = dualAggregates->GetProcWinner()->getDataNonConst(0);

  RCP<Dual2Primal_type> coarseMapNodesDualToPrimal = rcp(new Dual2Primal_type());
  RCP<Dual2Primal_type> coarseMapNodesPrimalToDual = rcp(new Dual2Primal_type());

  LocalOrdinal numLocalDualAggregates = 0;

  /* Loop over the local dual nodes and
   *
   * - assign dual nodes to dual aggregates
   * - recursively coarsen the dual-to-primal node mapping
   */
  LocalOrdinal localPrimalNodeID = - Teuchos::ScalarTraits<LocalOrdinal>::one();
  LocalOrdinal currentPrimalAggId = - Teuchos::ScalarTraits<LocalOrdinal>::one();
  for (LocalOrdinal localDualNodeID = 0; localDualNodeID < localNumDualNodes; ++localDualNodeID)
  {
    // Get local ID of the primal node associated to the current dual node
    localPrimalNodeID = (*mapNodesDualToPrimal)[localDualNodeID];

    // Find the primal aggregate that owns the current primal node
    currentPrimalAggId = primalVertex2AggId[localPrimalNodeID];

    // Test if the current primal aggregate has no associated dual aggregate, yet.
    // Create new dual aggregate, if necessary.
    if (coarseMapNodesPrimalToDual->count(currentPrimalAggId) == 0)
    {
      // Associate a new dual aggregate w/ the current primal aggregate
      (*coarseMapNodesPrimalToDual)[currentPrimalAggId] = numLocalDualAggregates;
      (*coarseMapNodesDualToPrimal)[numLocalDualAggregates] = currentPrimalAggId;
      ++numLocalDualAggregates;
    }

    // Fill the dual aggregate
    dualVertex2AggId[localDualNodeID] = (*coarseMapNodesPrimalToDual)[currentPrimalAggId];
    dualProcWinner[localDualNodeID] = myRank;
  }

  // Store dual aggregeate data as well as coarsening information
  dualAggregates->SetNumAggregates(numLocalDualAggregates);
  Set(currentLevel, "Aggregates", dualAggregates);
  Set(currentLevel, "CoarseDualNodeID2PrimalNodeID", coarseMapNodesDualToPrimal);
  GetOStream(Statistics1) << dualAggregates->description() << std::endl;
} // Build

} // namespace MueLu

#endif /* MUELU_INTERFACEAGGREGATIONFACTORY_DEF_HPP_ */
