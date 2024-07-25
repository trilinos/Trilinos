// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_USERAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_USERAGGREGATIONFACTORY_DEF_HPP_

#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapFactory.hpp>

#include "MueLu_UserAggregationFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> UserAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  // input parameters
  validParamList->set<std::string>("filePrefix", "", "The data is read from files of this name: <filePrefix><level>_<PID>.<fileExt>");
  validParamList->set<std::string>("fileExt", "", "The data is read from files of this name: <filePrefix><level>_<PID>.<fileExt>");

  return validParamList;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UserAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& /* currentLevel */) const {}

/**
 * The function reads aggregate information from a file.
 * The file structure is the following:
 *  * line 1 : <number of aggregates>
 *  * line 2+: <aggregate size> <node 1 (root node) GID> <node 2 GID> ... <node last GID>
 */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UserAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  const ParameterList& pL = GetParameterList();

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  const int myRank                    = comm->getRank();

  std::string fileName = pL.get<std::string>("filePrefix") + toString(currentLevel.GetLevelID()) + "_" + toString(myRank) + "." + pL.get<std::string>("fileExt");
  std::ifstream ifs(fileName.c_str());
  TEUCHOS_TEST_FOR_EXCEPTION(!ifs.good(), Exceptions::RuntimeError, "Cannot read data from \"" << fileName << "\"");

  LO numVertices, numAggregates;
  ifs >> numVertices;
  TEUCHOS_TEST_FOR_EXCEPTION(!ifs.good(), Exceptions::RuntimeError, "Cannot read data from \"" << fileName << "\"");
  ifs >> numAggregates;
  TEUCHOS_TEST_FOR_EXCEPTION(numVertices <= 0, Exceptions::InvalidArgument, "Number of vertices   must be > 0");
  TEUCHOS_TEST_FOR_EXCEPTION(numAggregates <= 0, Exceptions::InvalidArgument, "Number of aggregates must be > 0");

  Xpetra::UnderlyingLib lib = currentLevel.lib();
  const int indexBase       = 0;
  RCP<Map> map              = MapFactory::Build(lib, numVertices, indexBase, comm);

  RCP<Aggregates> aggregates = rcp(new Aggregates(map));
  aggregates->setObjectLabel("User");

  aggregates->SetNumAggregates(numAggregates);

  Teuchos::ArrayRCP<LO> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LO> procWinner   = aggregates->GetProcWinner()->getDataNonConst(0);

  for (LO i = 0; i < numAggregates; i++) {
    int aggSize = 0;
    ifs >> aggSize;

    std::vector<LO> list(aggSize);
    for (int k = 0; k < aggSize; k++) {
      // FIXME: File contains GIDs, we need LIDs
      // for now, works on a single processor
      ifs >> list[k];
    }

    // Mark first node as root node for the aggregate
    aggregates->SetIsRoot(list[0]);

    // Fill vertex2AggId and procWinner structure with information
    for (int k = 0; k < aggSize; k++) {
      vertex2AggId[list[k]] = i;
      procWinner[list[k]]   = myRank;
    }
  }

  // FIXME: do the proper check whether aggregates cross interprocessor boundary
  aggregates->AggregatesCrossProcessors(false);

  Set(currentLevel, "Aggregates", aggregates);

  GetOStream(Statistics0) << aggregates->description() << std::endl;
}

}  // namespace MueLu

#endif /* MUELU_USERAGGREGATIONFACTORY_DEF_HPP_ */
