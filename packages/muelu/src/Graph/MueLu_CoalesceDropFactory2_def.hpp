/*
 * MueLu_CoalesceDropFactory2_def.hpp
 *
 *  Created on: 10.07.2012
 *      Author: tobias
 */

#ifndef MUELU_COALESCEDROPFACTORY2_DEF_HPP_
#define MUELU_COALESCEDROPFACTORY2_DEF_HPP_

#include <Xpetra_Operator.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_StridedMap.hpp>

#include "MueLu_CoalesceDropFactory2_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_Monitor.hpp"
#include <MueLu_AmalgamationInfo.hpp>

namespace MueLu {

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
CoalesceDropFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CoalesceDropFactory2(RCP<const FactoryBase> AFact)
: AFact_(AFact)
  {

  }

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void CoalesceDropFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
  currentLevel.DeclareInput("A", AFact_.get(), this);
}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void CoalesceDropFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const {
  FactoryMonitor m(*this, "CoalesceDropFactory2", currentLevel);

  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  RCP<Operator> A = currentLevel.Get< RCP<Operator> >("A", AFact_.get());

  LocalOrdinal blockdim = 1; // block dim for fixed size blocks
  GlobalOrdinal offset = 0;  // global offset of dof gids

  // 1) check for blocking/striding information
  if(A->IsView("stridedMaps") &&
      Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
    Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps"); // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap()) == Teuchos::null,Exceptions::BadCast,"MueLu::CoalesceFactory::Build: cast to strided row map failed.");
    blockdim = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getFixedBlockSize();
    offset   = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getOffset();
    oldView = A->SwitchToView(oldView);
    GetOStream(Debug, 0) << "CoalesceDropFactory::Build():" << " found blockdim=" << blockdim << " from strided maps. offset=" << offset << std::endl;
  } else std::cout << "TODO: CoalesceDropFactory2: fix me: no striding information available." << std::endl;
  // TODO: maybe no striding information on coarser levels -> misuse nullspace vector?


  // 2) prepare maps for amalgamated graph of A and
  //    setup unamalgamation information

  RCP<std::vector<GlobalOrdinal> > gNodeIds; // contains global node ids on current proc
  gNodeIds = Teuchos::rcp(new std::vector<GlobalOrdinal>);
  gNodeIds->empty();

  // in nodegid2dofgids_ for each node on the current proc a vector of
  // the corresponding DOFs gids is stored.
  // The map contains all nodes the current proc has connections to (including
  // nodes that are stored on other procs when there are off-diagonal entries in A)
  nodegid2dofgids_ = Teuchos::rcp(new std::map<GlobalOrdinal,std::vector<GlobalOrdinal> >);

  // extract information from overlapping column map of A
  GlobalOrdinal cnt_amalRows = 0;
  for(LocalOrdinal i=0; i<Teuchos::as<LocalOrdinal>(A->getColMap()->getNodeNumElements());i++) {
    // get global DOF id
    GlobalOrdinal gDofId = A->getColMap()->getGlobalElement(i);

    // translate DOFGid to node id
    GlobalOrdinal gNodeId = DOFGid2NodeId(gDofId, A, blockdim, offset);

    // gblockid -> gDofId/lDofId
    if(nodegid2dofgids_->count(gNodeId) == 0) {

      // current column DOF gDofId belongs to a node that has not been added
      // to nodeid2dofgids_ yet. Do it now and add ALL DOFs of node gNodeId to
      // unamalgamation information.
      // Note: we use offset and blockdim, ie. information from strided maps indirectly
      std::vector<GlobalOrdinal> DOFs; DOFs.reserve(blockdim);
      for(LocalOrdinal k=0; k<blockdim; k++) {
        DOFs.push_back(offset + gNodeId*blockdim + k);
      }

      (*nodegid2dofgids_)[gNodeId] = DOFs; // store this vector? // TODO check me

      if(A->getRowMap()->isNodeGlobalElement(gDofId)) {
        gNodeIds->push_back(gNodeId);
        cnt_amalRows++; // new local block row in amalgamated matrix graph
      }
    }
  }

  // inter processor communication: sum up number of block ids
  GlobalOrdinal num_blockids = 0;
  Teuchos::reduceAll<int,GlobalOrdinal>(*(A->getRowMap()->getComm()),Teuchos::REDUCE_SUM, cnt_amalRows, Teuchos::ptr(&num_blockids) );
  GetOStream(Debug, 0) << "CoalesceDropFactory::SetupAmalgamationData()" << " # of amalgamated blocks=" << num_blockids << std::endl;

  // 3) generate row map for amalgamated matrix (graph of A)
  //    with same distribution over all procs as row map of A
  Teuchos::ArrayRCP<GlobalOrdinal> arr_gNodeIds = Teuchos::arcp( gNodeIds );
  Teuchos::RCP<Map> nodeMap = MapFactory::Build(A->getRowMap()->lib(), num_blockids, arr_gNodeIds(), A->getRowMap()->getIndexBase(), A->getRowMap()->getComm());
  GetOStream(Debug, 0) << "CoalesceDropFactory: nodeMap " << nodeMap->getNodeNumElements() << "/" << nodeMap->getGlobalNumElements() << " elements" << std::endl;

  // 4) create graph of amalgamated matrix
  RCP<CrsGraph> crsGraph = CrsGraphFactory::Build(nodeMap, 10, Xpetra::DynamicProfile);

  // 5) do amalgamation. generate graph of amalgamated matrix
  for(LocalOrdinal row=0; row<Teuchos::as<LocalOrdinal>(A->getRowMap()->getNodeNumElements()); row++) {
    // get global DOF id
    GlobalOrdinal grid = A->getRowMap()->getGlobalElement(row);

    // translate grid to nodeid
    GlobalOrdinal nodeId = DOFGid2NodeId(grid, A, blockdim, offset);

    size_t nnz = A->getNumEntriesInLocalRow(row);
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> vals;
    A->getLocalRowView(row, indices, vals);
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::CoalesceFactory::Amalgamate: number of nonzeros not equal to number of indices? Error.");

    RCP<std::vector<GlobalOrdinal> > cnodeIds = Teuchos::rcp(new std::vector<GlobalOrdinal>);  // global column block ids
    LocalOrdinal realnnz = 0;
    for(LocalOrdinal col=0; col<Teuchos::as<LocalOrdinal>(nnz); col++) {
      TEUCHOS_TEST_FOR_EXCEPTION(A->getColMap()->isNodeLocalElement(indices[col])==false,Exceptions::RuntimeError, "MueLu::CoalesceFactory::Amalgamate: Problem with columns. Error.");
      GlobalOrdinal gcid = A->getColMap()->getGlobalElement(indices[col]); // global column id

      if(vals[col]!=0.0) {
        GlobalOrdinal cnodeId = DOFGid2NodeId(gcid, A, blockdim, offset);
        cnodeIds->push_back(cnodeId);
        realnnz++; // increment number of nnz in matrix row
      }
    }

    Teuchos::ArrayRCP<GlobalOrdinal> arr_cnodeIds = Teuchos::arcp( cnodeIds );


    TEUCHOS_TEST_FOR_EXCEPTION(crsGraph->getRowMap()->isNodeGlobalElement(nodeId)==false,Exceptions::RuntimeError, "MueLu::CoalesceFactory::Amalgamate: global row id does not belong to current proc. Error.");
    crsGraph->insertGlobalIndices(nodeId, arr_cnodeIds());
  }
  // fill matrix graph
  crsGraph->fillComplete(nodeMap,nodeMap);

  // 6) create MueLu Graph object
  RCP<Graph> graph = rcp(new Graph(crsGraph, "amalgamated graph of A"));

  // store information in Graph object for unamalgamation of vectors
  // TODO remove this
  graph->SetAmalgamationParams(nodegid2dofgids_);

  // store (un)amalgamation information on current level
  RCP<AmalgamationInfo> amalgamationData = rcp(new AmalgamationInfo());
  amalgamationData->SetAmalgamationParams(nodegid2dofgids_);
  currentLevel.Set("UnAmalgamationInfo", amalgamationData, this);

  // 7) store results in Level
  currentLevel.Set("DofsPerNode", blockdim, this);
  currentLevel.Set("Graph", graph, this);
}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
GlobalOrdinal CoalesceDropFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DOFGid2NodeId(GlobalOrdinal gid, const RCP<Operator>& A, LocalOrdinal blockSize, const GlobalOrdinal offset) const {
  //GetOStream(Runtime0, 0) << "fixed block size..." << std::endl;
  GlobalOrdinal globalblockid = ((GlobalOrdinal) gid - offset) / blockSize;
  return globalblockid;
}

}

#endif /* MUELU_COALESCEDROPFACTORY2_DEF_HPP_ */
