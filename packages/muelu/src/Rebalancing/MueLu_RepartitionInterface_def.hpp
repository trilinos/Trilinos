/*
 * MueLu_RepartitionInterface_def.hpp
 *
 *  Created on: 5 Sep 2013
 *      Author: wiesner
 */

#ifndef MUELU_REPARTITIONINTERFACE_DEF_HPP_
#define MUELU_REPARTITIONINTERFACE_DEF_HPP_

#include "MueLu_RepartitionInterface_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> RepartitionInterface<LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Factory of the matrix A");
  validParamList->set<RCP<const FactoryBase> >("number of partitions", Teuchos::null, "Instance of RepartitionHeuristicFactory.");
  validParamList->set<RCP<const FactoryBase> >("AmalgamatedPartition", Teuchos::null, "(advanced) Factory generating the AmalgamatedPartition (e.g. an IsorropiaInterface)");

  return validParamList;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void RepartitionInterface<LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  Input(currentLevel, "A");
  Input(currentLevel, "number of partitions");
  Input(currentLevel, "AmalgamatedPartition");
}  // DeclareInput()

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void RepartitionInterface<LocalOrdinal, GlobalOrdinal, Node>::Build(Level &level) const {
  FactoryMonitor m(*this, "Build", level);

  RCP<Matrix> A                                       = Get<RCP<Matrix> >(level, "A");
  RCP<Xpetra::Vector<GO, LO, GO, NO> > amalgPartition = Get<RCP<Xpetra::Vector<GO, LO, GO, NO> > >(level, "AmalgamatedPartition");
  int numParts                                        = Get<int>(level, "number of partitions");

  RCP<const Map> rowMap = A->getRowMap();

  // standard case: use matrix info and amalgamated rebalancing info to create "Partition" vector
  RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();

  // Short cut: if we only need one partition, then create a dummy partition vector
  if (numParts == 1 || numParts == -1) {
    // Single processor, decomposition is trivial: all zeros
    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, true);
    Set(level, "Partition", decomposition);
    return;
  } /* else if (numParts == -1) {
     // No repartitioning
     RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Teuchos::null; //Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, true);
     //decomposition->putScalar(Teuchos::as<Scalar>(comm->getRank()));
     Set(level, "Partition", decomposition);
     return;
   }*/

  ArrayRCP<GO> amalgPartitionData = amalgPartition->getDataNonConst(0);
  RCP<const Map> nodeMap          = amalgPartition->getMap();

  // extract amalgamation information from matrix A
  LO blockdim         = 1;         // block dim for fixed size blocks
  LO blockid          = -1;        // block id in strided map
  LO nStridedOffset   = 0;         // DOF offset for strided block id "blockid" (default = 0)
  LO stridedblocksize = blockdim;  // size of strided block id "blockid" (default = fullblocksize, only if blockid!=-1 stridedblocksize <= fullblocksize)

  // 1) check for blocking/striding information
  //    fill above variables
  if (A->IsView("stridedMaps") &&
      Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
    Xpetra::viewLabel_t oldView  = A->SwitchToView("stridedMaps");  // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
    RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap());
    TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null, Exceptions::BadCast, "MueLu::RepartitionInterface::Build: cast to strided row map failed.");
    blockdim = strMap->getFixedBlockSize();
    blockid  = strMap->getStridedBlockId();
    if (blockid > -1) {
      std::vector<size_t> stridingInfo = strMap->getStridingData();
      for (size_t j = 0; j < Teuchos::as<size_t>(blockid); j++)
        nStridedOffset += stridingInfo[j];
      stridedblocksize = Teuchos::as<LocalOrdinal>(stridingInfo[blockid]);

    } else {
      stridedblocksize = blockdim;
    }
    oldView = A->SwitchToView(oldView);
    // GetOStream(Statistics0, -1) << "RepartitionInterface::Build():" << " found blockdim=" << blockdim << " from strided maps (blockid=" << blockid << ", strided block size=" << stridedblocksize << "). offset=" << offset << std::endl;
  } else
    GetOStream(Statistics0, -1) << "RepartitionInterface::Build(): no striding information available. Use blockdim=1 with offset=0" << std::endl;

  // vector which stores final (unamalgamated) repartitioning
  RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, false);
  ArrayRCP<GO> decompEntries                         = decomposition->getDataNonConst(0);

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<int>(nodeMap->getLocalNumElements()) * stridedblocksize != Teuchos::as<int>(rowMap->getLocalNumElements()), Exceptions::RuntimeError, "Inconsistency between nodeMap and dofMap: we are supporting block maps only. No support for general strided maps, yet!");

  // RCP<std::map<GO,std::vector<GO> > > nodegid2dofgids = amalgInfo->GetGlobalAmalgamationParams();

  // fill vector with information about partitioning
  // TODO: we assume simple block maps here
  // TODO: adapt this to usage of nodegid2dofgids
  for (size_t i = 0; i < nodeMap->getLocalNumElements(); i++) {
    // not fully sure about this. We're filling local ids in the decomposition vector with
    // the results stored in array. The decomposition vector is created using the rowMap of A

    // transform local node id to global node id.
    // GO gNodeId = nodeMap->getGlobalElement(i);

    // extract global DOF ids that belong to gNodeId
    /*std::vector<GlobalOrdinal> DOFs = (*nodegid2dofgids)[gNodeId];
    for(size_t j=0; j<stridedblocksize; j++) {
      decompEntries[i*stridedblocksize + j] = myRank;
    }*/
    for (LO j = 0; j < stridedblocksize /*DOFs.size()*/; j++) {
      // transform global DOF ids to local DOF ids using rowMap
      // note: The vector decomposition is based on rowMap
      // LO lDofId = rowMap->getLocalElement(DOFs[j]);     // -> i doubt that we need this!

      // put the same domain id to all DOFs of the same node
      decompEntries[i * stridedblocksize + j] = amalgPartitionData[i];
      // decompEntries[lDofId] = amalgPartitionData[i];
    }
  }

  Set(level, "Partition", decomposition);

}  // Build()

}  // namespace MueLu

#endif /* MUELU_REPARTITIONINTERFACE_DEF_HPP_ */
