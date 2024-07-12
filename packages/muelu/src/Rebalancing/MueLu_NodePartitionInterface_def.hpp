// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_NODEPARTITIONINTERFACE_DEF_HPP
#define MUELU_NODEPARTITIONINTERFACE_DEF_HPP

#include <sstream>
#include <set>

#include "MueLu_NodePartitionInterface_decl.hpp"
#if defined(HAVE_MPI)
#include <Teuchos_Utils.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_MasterList.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
NodePartitionInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NodePartitionInterface() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> NodePartitionInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("repartition: node id");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Factory of the matrix A");
  validParamList->set<RCP<const FactoryBase> >("number of partitions", Teuchos::null, "Instance of RepartitionHeuristicFactory.");
  validParamList->set<RCP<const FactoryBase> >("Node Comm", Teuchos::null, "Generating factory of the node level communicator");

  // We don't really need this, but we might need it on coarser levels
  validParamList->set<RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Factory of the coordinates");
  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void NodePartitionInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
  Input(currentLevel, "number of partitions");
  Input(currentLevel, "Node Comm");
  Input(currentLevel, "Coordinates");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void NodePartitionInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& level) const {
  FactoryMonitor m(*this, "Build", level);
  RCP<Matrix> A         = Get<RCP<Matrix> >(level, "A");
  RCP<const Map> rowMap = A->getRowMap();

  int numParts = Get<int>(level, "number of partitions");
  if (numParts == 1 || numParts == -1) {
    // Single processor, decomposition is trivial: all zeros
    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, true);
    Set(level, "Partition", decomposition);
    return;
  }

  // Let us repartition nodally
  RCP<const Teuchos::Comm<int> > NodeComm = Get<RCP<const Teuchos::Comm<int> > >(level, "Node Comm");
  TEUCHOS_TEST_FOR_EXCEPTION(NodeComm.is_null(), Exceptions::RuntimeError, "MueLu::NodePartitionInterface::Build(): NodeComm is null.");

  // Get the rank (in current comm) of rank 0 in my NodeComm
  int nodeZeroRank = A->getMap()->getComm()->getRank();
  Teuchos::broadcast<int, int>(*NodeComm, 0, Teuchos::inOutArg(nodeZeroRank));

  // A "Partition" from a *Interface is supposed to be is a vector of length # of my rows with the partition number to which the unknown is assigned
  // BUT, since we're bypassing remap for NodePartition, we'll return a *rank* of the guy who gets each unknown (which is what remap outputs).
  RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, false);
  decomposition->putScalar(Teuchos::as<GO>(nodeZeroRank));

  Set(level, "Partition", decomposition);
}

}  // namespace MueLu

#endif  // if defined(HAVE_MPI)

#endif  // MUELU_NODEPARTITIONINTERFACE_DEF_HPP
