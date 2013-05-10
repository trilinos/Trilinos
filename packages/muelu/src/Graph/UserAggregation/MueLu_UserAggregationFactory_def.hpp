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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_USERAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_USERAGGREGATIONFACTORY_DEF_HPP_

#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapFactory.hpp>

#include "MueLu_UserAggregationFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<const ParameterList> UserAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  // input parameters
  validParamList->set<std::string>("filePrefix", "", "The data is read from files of this name: <filePrefix><level>_<PID>.<fileExt>");
  validParamList->set<std::string>("fileExt",    "", "The data is read from files of this name: <filePrefix><level>_<PID>.<fileExt>");

  return validParamList;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void UserAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
}

/**
 * The function reads aggregate information from a file.
 * The file structure is the following:
 *  * line 1 : <number of aggregates>
 *  * line 2+: <aggregate size> <node 1 (root node) GID> <node 2 GID> ... <node last GID>
 */
template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void UserAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  const ParameterList & pL = GetParameterList();

  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  const int myRank = comm->getRank();

  std::string fileName = pL.get<std::string>("filePrefix") + toString(currentLevel.GetLevelID()) + "_" + toString(myRank) + "." + pL.get<std::string>("fileExt");
  std::ifstream ifs(fileName.c_str());
  if (!ifs.good())
    throw Exceptions::RuntimeError("Cannot read data from \"" + fileName + "\"");

  LO numVertices, numAggregates;
  ifs >> numVertices >> numAggregates;

  std::cout << "numVertices = " << numVertices << std::endl;

  // FIXME: what is the map?
  Xpetra::UnderlyingLib  lib       = Xpetra::UseEpetra;
  const int              indexBase = 0;
  RCP<Map> map = MapFactory::Build(lib, numVertices, indexBase, comm);

  RCP<Aggregates> aggregates = rcp(new Aggregates(map));
  aggregates->setObjectLabel("User");

  aggregates->SetNumAggregates(numAggregates);

  Teuchos::ArrayRCP<LO> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LO> procWinner   = aggregates->GetProcWinner()  ->getDataNonConst(0);

  for (LO i = 0; i < numAggregates; i++) {
    int aggSize = 0;
    ifs >> aggSize;

    std::vector<LO> list(aggSize);
    for (unsigned int k = 0; k < aggSize; k++) {
      // FIXME: File contains GIDs, we need LIDs
      // for now, works on a single processor
      ifs >> list[k];
    }

    // Mark first node as root node for the aggregate
    aggregates->SetIsRoot(list[0]);

    // Fill vertex2AggId and procWinner structure with information
    for (unsigned int k = 0; k < aggSize; k++) {
      vertex2AggId[list[k]] = i;
      procWinner  [list[k]] = myRank;
    }
  }

  // FIXME: do the proper check whether aggregates cross interprocessor boundary
  aggregates->AggregatesCrossProcessors(false);

  Set(currentLevel, "Aggregates", aggregates);

  GetOStream(Statistics0, 0) << aggregates->description() << std::endl;
}

} //namespace MueLu

#endif /* MUELU_USERAGGREGATIONFACTORY_DEF_HPP_ */
