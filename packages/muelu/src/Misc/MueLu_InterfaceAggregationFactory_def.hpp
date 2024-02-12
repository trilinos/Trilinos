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

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_StridedMap.hpp>

#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_InterfaceAggregationFactory_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> InterfaceAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Generating factory of A (matrix block related to dual DOFs)");
  validParamList->set<RCP<const FactoryBase>>("Aggregates", Teuchos::null, "Generating factory of the Aggregates (for block 0,0)");

  validParamList->set<std::string>("Dual/primal mapping strategy", "vague",
                                   "Strategy to represent mapping between dual and primal quantities [node-based, dof-based]");

  validParamList->set<RCP<const FactoryBase>>("DualNodeID2PrimalNodeID", Teuchos::null,
                                              "Generating factory of the DualNodeID2PrimalNodeID map as input data in a Moertel-compatible std::map<LO,LO> to map local IDs of dual nodes to local IDs of primal nodes");
  validParamList->set<LocalOrdinal>("number of DOFs per dual node", -Teuchos::ScalarTraits<LocalOrdinal>::one(),
                                    "Number of DOFs per dual node");

  validParamList->set<RCP<const FactoryBase>>("Primal interface DOF map", Teuchos::null,
                                              "Generating factory of the primal DOF row map of slave side of the coupling surface");

  return validParamList;
}  // GetValidParameterList()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  Input(currentLevel, "A");  // matrix block of dual variables
  Input(currentLevel, "Aggregates");

  const ParameterList &pL = GetParameterList();
  TEUCHOS_TEST_FOR_EXCEPTION(pL.get<std::string>("Dual/primal mapping strategy") == "vague", Exceptions::InvalidArgument,
                             "Strategy for dual/primal mapping not selected. Please select one of the available strategies.")
  if (pL.get<std::string>("Dual/primal mapping strategy") == "node-based") {
    if (currentLevel.GetLevelID() == 0) {
      TEUCHOS_TEST_FOR_EXCEPTION(!currentLevel.IsAvailable("DualNodeID2PrimalNodeID", NoFactory::get()),
                                 Exceptions::RuntimeError, "DualNodeID2PrimalNodeID was not provided by the user on level 0!");

      currentLevel.DeclareInput("DualNodeID2PrimalNodeID", NoFactory::get(), this);
    } else {
      Input(currentLevel, "DualNodeID2PrimalNodeID");
    }
  } else if (pL.get<std::string>("Dual/primal mapping strategy") == "dof-based") {
    if (currentLevel.GetLevelID() == 0)
      currentLevel.DeclareInput("Primal interface DOF map", NoFactory::get(), this);
    else
      Input(currentLevel, "Primal interface DOF map");
  } else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::InvalidArgument, "Unknown strategy for dual/primal mapping.")

}  // DeclareInput

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
  const std::string prefix = "MueLu::InterfaceAggregationFactory::Build: ";

  FactoryMonitor m(*this, "Build", currentLevel);

  // Call a specialized build routine based on the format of user-given input
  const ParameterList &pL         = GetParameterList();
  const std::string parameterName = "Dual/primal mapping strategy";
  if (pL.get<std::string>(parameterName) == "node-based")
    BuildBasedOnNodeMapping(prefix, currentLevel);
  else if (pL.get<std::string>(parameterName) == "dof-based")
    BuildBasedOnPrimalInterfaceDofMap(prefix, currentLevel);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::InvalidArgument,
                               "MueLu::InterfaceAggregationFactory::Builld(): Unknown strategy for dual/primal mapping. Set a valid value for the parameter \"" << parameterName << "\".")
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildBasedOnNodeMapping(const std::string &prefix,
                                                                                                     Level &currentLevel) const {
  using Dual2Primal_type = std::map<LocalOrdinal, LocalOrdinal>;

  const ParameterList &pL = GetParameterList();

  RCP<const Matrix> A                   = Get<RCP<Matrix>>(currentLevel, "A");
  const LocalOrdinal numDofsPerDualNode = pL.get<LocalOrdinal>("number of DOFs per dual node");
  TEUCHOS_TEST_FOR_EXCEPTION(numDofsPerDualNode < Teuchos::ScalarTraits<LocalOrdinal>::one(), Exceptions::InvalidArgument,
                             "Number of dual DOFs per node < 0 (default value). Specify a valid \"number of DOFs per dual node\" in the parameter list for the InterfaceAggregationFactory.");

  RCP<const Aggregates> primalAggregates          = Get<RCP<Aggregates>>(currentLevel, "Aggregates");
  ArrayRCP<const LocalOrdinal> primalVertex2AggId = primalAggregates->GetVertex2AggId()->getData(0);

  // Get the user-prescribed mapping of dual to primal node IDs
  RCP<Dual2Primal_type> mapNodesDualToPrimal;
  if (currentLevel.GetLevelID() == 0)
    mapNodesDualToPrimal = currentLevel.Get<RCP<Dual2Primal_type>>("DualNodeID2PrimalNodeID", NoFactory::get());
  else
    mapNodesDualToPrimal = Get<RCP<Dual2Primal_type>>(currentLevel, "DualNodeID2PrimalNodeID");

  RCP<const Map> operatorRangeMap = A->getRangeMap();
  const size_t myRank             = operatorRangeMap->getComm()->getRank();

  LocalOrdinal globalNumDualNodes = operatorRangeMap->getGlobalNumElements() / numDofsPerDualNode;
  LocalOrdinal localNumDualNodes  = operatorRangeMap->getLocalNumElements() / numDofsPerDualNode;

  TEUCHOS_TEST_FOR_EXCEPTION(localNumDualNodes != Teuchos::as<LocalOrdinal>(mapNodesDualToPrimal->size()),
                             std::runtime_error, prefix << " MueLu requires the range map and the DualNodeID2PrimalNodeID map to be compatible.");

  RCP<const Map> dualNodeMap = Teuchos::null;
  if (numDofsPerDualNode == 1)
    dualNodeMap = operatorRangeMap;
  else {
    GlobalOrdinal indexBase                = operatorRangeMap->getIndexBase();
    auto comm                              = operatorRangeMap->getComm();
    std::vector<GlobalOrdinal> myDualNodes = {};

    for (size_t i = 0; i < operatorRangeMap->getLocalNumElements(); i += numDofsPerDualNode)
      myDualNodes.push_back((operatorRangeMap->getGlobalElement(i) - indexBase) / numDofsPerDualNode + indexBase);

    dualNodeMap = MapFactory::Build(operatorRangeMap->lib(), globalNumDualNodes, myDualNodes, indexBase, comm);
  }
  TEUCHOS_TEST_FOR_EXCEPTION(localNumDualNodes != Teuchos::as<LocalOrdinal>(dualNodeMap->getLocalNumElements()),
                             std::runtime_error, prefix << " Local number of dual nodes given by user is incompatible to the dual node map.");

  RCP<Aggregates> dualAggregates = rcp(new Aggregates(dualNodeMap));
  dualAggregates->setObjectLabel("InterfaceAggregation");

  // Copy setting from primal aggregates, as we copy the interface part of primal aggregates anyways
  dualAggregates->AggregatesCrossProcessors(primalAggregates->AggregatesCrossProcessors());

  ArrayRCP<LocalOrdinal> dualVertex2AggId = dualAggregates->GetVertex2AggId()->getDataNonConst(0);
  ArrayRCP<LocalOrdinal> dualProcWinner   = dualAggregates->GetProcWinner()->getDataNonConst(0);

  RCP<Dual2Primal_type> coarseMapNodesDualToPrimal = rcp(new Dual2Primal_type());
  RCP<Dual2Primal_type> coarseMapNodesPrimalToDual = rcp(new Dual2Primal_type());

  LocalOrdinal numLocalDualAggregates = 0;

  /* Loop over the local dual nodes and
   *
   * - assign dual nodes to dual aggregates
   * - recursively coarsen the dual-to-primal node mapping
   */
  LocalOrdinal localPrimalNodeID  = -Teuchos::ScalarTraits<LocalOrdinal>::one();
  LocalOrdinal currentPrimalAggId = -Teuchos::ScalarTraits<LocalOrdinal>::one();
  for (LocalOrdinal localDualNodeID = 0; localDualNodeID < localNumDualNodes; ++localDualNodeID) {
    // Get local ID of the primal node associated to the current dual node
    localPrimalNodeID = (*mapNodesDualToPrimal)[localDualNodeID];

    // Find the primal aggregate that owns the current primal node
    currentPrimalAggId = primalVertex2AggId[localPrimalNodeID];

    // Test if the current primal aggregate has no associated dual aggregate, yet.
    // Create new dual aggregate, if necessary.
    if (coarseMapNodesPrimalToDual->count(currentPrimalAggId) == 0) {
      // Associate a new dual aggregate w/ the current primal aggregate
      (*coarseMapNodesPrimalToDual)[currentPrimalAggId]     = numLocalDualAggregates;
      (*coarseMapNodesDualToPrimal)[numLocalDualAggregates] = currentPrimalAggId;
      ++numLocalDualAggregates;
    }

    // Fill the dual aggregate
    dualVertex2AggId[localDualNodeID] = (*coarseMapNodesPrimalToDual)[currentPrimalAggId];
    dualProcWinner[localDualNodeID]   = myRank;
  }

  // Store dual aggregeate data as well as coarsening information
  dualAggregates->SetNumAggregates(numLocalDualAggregates);
  Set(currentLevel, "Aggregates", dualAggregates);
  Set(currentLevel, "CoarseDualNodeID2PrimalNodeID", coarseMapNodesDualToPrimal);
  GetOStream(Statistics1) << dualAggregates->description() << std::endl;
}  // BuildBasedOnNodeMapping

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildBasedOnPrimalInterfaceDofMap(
    const std::string &prefix, Level &currentLevel) const {
  const GlobalOrdinal GO_ZERO = Teuchos::ScalarTraits<LocalOrdinal>::zero();
  const GlobalOrdinal GO_ONE  = Teuchos::ScalarTraits<GlobalOrdinal>::one();

  // filled with striding information from A01
  LocalOrdinal numDofsPerDualNode   = 0;
  LocalOrdinal numDofsPerPrimalNode = 0;

  // Grab the off-diagonal block (0,1) from the global blocked operator
  RCP<const Matrix> A01                           = Get<RCP<Matrix>>(currentLevel, "A");
  RCP<const Aggregates> primalAggregates          = Get<RCP<Aggregates>>(currentLevel, "Aggregates");
  ArrayRCP<const LocalOrdinal> primalVertex2AggId = primalAggregates->GetVertex2AggId()->getData(0);

  auto comm        = A01->getRowMap()->getComm();
  const int myRank = comm->getRank();

  RCP<const Map> primalInterfaceDofRowMap = Teuchos::null;
  if (currentLevel.GetLevelID() == 0) {
    // Use NoFactory, since the fine level asks for user data
    primalInterfaceDofRowMap = currentLevel.Get<RCP<const Map>>("Primal interface DOF map", NoFactory::get());
  } else {
    primalInterfaceDofRowMap = Get<RCP<const Map>>(currentLevel, "Primal interface DOF map");
  }
  TEUCHOS_ASSERT(!primalInterfaceDofRowMap.is_null());

  if (A01->IsView("stridedMaps") && rcp_dynamic_cast<const StridedMap>(A01->getRowMap("stridedMaps")) != Teuchos::null) {
    auto stridedRowMap   = rcp_dynamic_cast<const StridedMap>(A01->getRowMap("stridedMaps"));
    auto stridedColMap   = rcp_dynamic_cast<const StridedMap>(A01->getColMap("stridedMaps"));
    numDofsPerPrimalNode = Teuchos::as<LocalOrdinal>(stridedRowMap->getFixedBlockSize());
    numDofsPerDualNode   = Teuchos::as<LocalOrdinal>(stridedColMap->getFixedBlockSize());

    if (numDofsPerPrimalNode != numDofsPerDualNode) {
      GetOStream(Warnings) << "InterfaceAggregation attempts to work with "
                           << numDofsPerPrimalNode << " primal DOFs per node and " << numDofsPerDualNode << " dual DOFs per node."
                           << "Be careful! Algorithm is not well-tested, if number of primal and dual DOFs per node differ." << std::endl;
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION(numDofsPerPrimalNode == 0, Exceptions::RuntimeError,
                             "InterfaceAggregationFactory could not extract the number of primal DOFs per node from striding information. At least, make sure that StridedMap information has actually been provided.");
  TEUCHOS_TEST_FOR_EXCEPTION(numDofsPerDualNode == 0, Exceptions::RuntimeError,
                             "InterfaceAggregationFactory could not extract the number of dual DOFs per node from striding information. At least, make sure that StridedMap information has actually been provided.");

  /* Determine block information for primal block
   *
   * primalDofOffset: global offset of primal DOF GIDs (usually is zero (default))
   * primalBlockDim: block dim for fixed size blocks
   * - is 2 or 3 (for 2d or 3d problems) on the finest level (# displacement dofs per node)
   * - is 3 or 6 (for 2d or 3d problems) on coarser levels (# nullspace vectors)
   */
  GlobalOrdinal primalDofOffset = GO_ZERO;
  LocalOrdinal primalBlockDim   = numDofsPerPrimalNode;

  /* Determine block information for Lagrange multipliers
   *
   * dualDofOffset: usually > zero (set by domainOffset for Ptent11Fact)
   * dualBlockDim:
   * - is primalBlockDim (for 2d or 3d problems) on the finest level (1 Lagrange multiplier per
   *   displacement dof)
   * - is 2 or 3 (for 2d or 3d problems) on coarser levels (same as on finest level, whereas there
   *   are 3 or 6 displacement dofs per node)
   */
  GlobalOrdinal dualDofOffset = A01->getColMap()->getMinAllGlobalIndex();
  LocalOrdinal dualBlockDim   = numDofsPerDualNode;

  // Generate global replicated mapping "lagrNodeId -> dispNodeId"
  RCP<const Map> dualDofMap    = A01->getDomainMap();
  GlobalOrdinal gMaxDualNodeId = AmalgamationFactory::DOFGid2NodeId(
      dualDofMap->getMaxAllGlobalIndex(), dualBlockDim, dualDofOffset, dualDofMap->getIndexBase());
  GlobalOrdinal gMinDualNodeId = AmalgamationFactory::DOFGid2NodeId(
      dualDofMap->getMinAllGlobalIndex(), dualBlockDim, dualDofOffset, dualDofMap->getIndexBase());

  GetOStream(Runtime1) << "  Dual DOF map: index base = " << dualDofMap->getIndexBase()
                       << ", block dim = " << dualBlockDim
                       << ", gid offset = " << dualDofOffset
                       << std::endl;

  GetOStream(Runtime1) << "  [primal / dual] DOFs per node = [" << numDofsPerPrimalNode
                       << "/" << numDofsPerDualNode << "]" << std::endl;

  // Generate locally replicated vector for mapping dual node IDs to primal node IDs
  Array<GlobalOrdinal> dualNodeId2primalNodeId(gMaxDualNodeId - gMinDualNodeId + 1, -GO_ONE);
  Array<GlobalOrdinal> local_dualNodeId2primalNodeId(gMaxDualNodeId - gMinDualNodeId + 1, -GO_ONE);

  // Generate locally replicated vector for mapping dual node IDs to primal aggregate ID
  Array<GlobalOrdinal> dualNodeId2primalAggId(gMaxDualNodeId - gMinDualNodeId + 1, -GO_ONE);
  Array<GlobalOrdinal> local_dualNodeId2primalAggId(gMaxDualNodeId - gMinDualNodeId + 1, -GO_ONE);

  Array<GlobalOrdinal> dualDofId2primalDofId(primalInterfaceDofRowMap->getGlobalNumElements(), -GO_ONE);
  Array<GlobalOrdinal> local_dualDofId2primalDofId(primalInterfaceDofRowMap->getGlobalNumElements(), -GO_ONE);

  // Fill mapping of Lagrange Node IDs to displacement aggregate IDs
  const size_t numMyPrimalInterfaceDOFs = primalInterfaceDofRowMap->getLocalNumElements();
  for (size_t r = 0; r < numMyPrimalInterfaceDOFs; r += numDofsPerPrimalNode) {
    GlobalOrdinal gPrimalRowId = primalInterfaceDofRowMap->getGlobalElement(r);

    if (A01->getRowMap()->isNodeGlobalElement(gPrimalRowId))  // Remove this if?
    {
      const LocalOrdinal lPrimalRowId   = A01->getRowMap()->getLocalElement(gPrimalRowId);
      const GlobalOrdinal gPrimalNodeId = AmalgamationFactory::DOFGid2NodeId(gPrimalRowId, primalBlockDim, primalDofOffset, primalInterfaceDofRowMap->getIndexBase());
      const LocalOrdinal lPrimalNodeId  = lPrimalRowId / numDofsPerPrimalNode;
      const LocalOrdinal primalAggId    = primalVertex2AggId[lPrimalNodeId];

      const GlobalOrdinal gDualDofId = A01->getColMap()->getGlobalElement(r);

      const GlobalOrdinal gDualNodeId = AmalgamationFactory::DOFGid2NodeId(gDualDofId, dualBlockDim, dualDofOffset, 0);

      if (local_dualNodeId2primalNodeId[gDualNodeId - gMinDualNodeId] == -GO_ONE) {
        local_dualNodeId2primalNodeId[gDualNodeId - gMinDualNodeId] = gPrimalNodeId;
        local_dualNodeId2primalAggId[gDualNodeId - gMinDualNodeId]  = primalAggId;
      } else {
        GetOStream(Warnings) << "PROC: " << myRank << " gDualNodeId " << gDualNodeId << " is already connected to primal nodeId "
                             << local_dualNodeId2primalNodeId[gDualNodeId - gMinDualNodeId]
                             << ". Ignore new dispNodeId: " << gPrimalNodeId << std::endl;
      }
    }
  }

  const int dualNodeId2primalNodeIdSize = Teuchos::as<int>(local_dualNodeId2primalNodeId.size());
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, dualNodeId2primalNodeIdSize,
                     &local_dualNodeId2primalNodeId[0], &dualNodeId2primalNodeId[0]);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, dualNodeId2primalNodeIdSize,
                     &local_dualNodeId2primalAggId[0], &dualNodeId2primalAggId[0]);

  // build node map for dual variables
  // generate "artificial nodes" for lagrange multipliers
  // the node map is also used for defining the Aggregates for the lagrange multipliers
  std::vector<GlobalOrdinal> dualNodes;
  for (size_t r = 0; r < A01->getDomainMap()->getLocalNumElements(); r++) {
    // determine global Lagrange multiplier row Dof
    // generate a node id using the grid, lagr_blockdim and lagr_offset // todo make sure, that
    // nodeId is unique and does not interfer with the displacement nodes
    GlobalOrdinal gDualDofId  = A01->getDomainMap()->getGlobalElement(r);
    GlobalOrdinal gDualNodeId = AmalgamationFactory::DOFGid2NodeId(gDualDofId, dualBlockDim, dualDofOffset, 0);
    dualNodes.push_back(gDualNodeId);
  }

  // remove all duplicates
  dualNodes.erase(std::unique(dualNodes.begin(), dualNodes.end()), dualNodes.end());

  // define node map for Lagrange multipliers
  Teuchos::RCP<const Map> dualNodeMap = MapFactory::Build(A01->getRowMap()->lib(),
                                                          Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), dualNodes, A01->getRowMap()->getIndexBase(), comm);

  // Build aggregates using the lagrange multiplier node map
  Teuchos::RCP<Aggregates> dualAggregates = Teuchos::rcp(new Aggregates(dualNodeMap));
  dualAggregates->setObjectLabel("UC (dual variables)");

  // extract aggregate data structures to fill
  Teuchos::ArrayRCP<LocalOrdinal> dualVertex2AggId = dualAggregates->GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> dualProcWinner   = dualAggregates->GetProcWinner()->getDataNonConst(0);

  // loop over local lagrange multiplier node ids
  LocalOrdinal nLocalAggregates = 0;
  std::map<GlobalOrdinal, LocalOrdinal> primalAggId2localDualAggId;
  for (size_t lDualNodeID = 0; lDualNodeID < dualNodeMap->getLocalNumElements(); ++lDualNodeID) {
    const GlobalOrdinal gDualNodeId = dualNodeMap->getGlobalElement(lDualNodeID);
    const GlobalOrdinal primalAggId = dualNodeId2primalAggId[gDualNodeId - gMinDualNodeId];
    if (primalAggId2localDualAggId.count(primalAggId) == 0)
      primalAggId2localDualAggId[primalAggId] = nLocalAggregates++;
    dualVertex2AggId[lDualNodeID] = primalAggId2localDualAggId[primalAggId];
    dualProcWinner[lDualNodeID]   = myRank;
  }

  const LocalOrdinal fullblocksize    = numDofsPerDualNode;
  const GlobalOrdinal offset          = A01->getColMap()->getMinAllGlobalIndex();
  const LocalOrdinal blockid          = -1;
  const LocalOrdinal nStridedOffset   = 0;
  const LocalOrdinal stridedblocksize = fullblocksize;

  RCP<Array<LO>> rowTranslation = rcp(new Array<LO>());
  RCP<Array<LO>> colTranslation = rcp(new Array<LO>());
  const size_t numMyDualNodes   = dualNodeMap->getLocalNumElements();
  for (size_t lDualNodeID = 0; lDualNodeID < numMyDualNodes; ++lDualNodeID) {
    for (LocalOrdinal dof = 0; dof < numDofsPerDualNode; ++dof) {
      rowTranslation->push_back(lDualNodeID);
      colTranslation->push_back(lDualNodeID);
    }
  }

  TEUCHOS_ASSERT(A01->isFillComplete());

  RCP<AmalgamationInfo> dualAmalgamationInfo = rcp(new AmalgamationInfo(rowTranslation, colTranslation,
                                                                        A01->getDomainMap(), A01->getDomainMap(), A01->getDomainMap(),
                                                                        fullblocksize, offset, blockid, nStridedOffset, stridedblocksize));

  dualAggregates->SetNumAggregates(nLocalAggregates);
  dualAggregates->AggregatesCrossProcessors(primalAggregates->AggregatesCrossProcessors());

  if (dualAggregates->AggregatesCrossProcessors())
    GetOStream(Runtime1) << "Interface aggregates cross processor boundaries." << std::endl;
  else
    GetOStream(Runtime1) << "Interface aggregates do not cross processor boundaries." << std::endl;

  currentLevel.Set("Aggregates", dualAggregates, this);
  currentLevel.Set("UnAmalgamationInfo", dualAmalgamationInfo, this);

}  // BuildBasedOnPrimalInterfaceDofMap

}  // namespace MueLu

#endif /* MUELU_INTERFACEAGGREGATIONFACTORY_DEF_HPP_ */
