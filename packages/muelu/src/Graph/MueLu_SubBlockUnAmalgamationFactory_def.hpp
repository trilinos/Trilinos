#ifndef MUELU_SUBBLOCKUNAMALGAMATIONFACTORY_DEF_HPP
#define MUELU_SUBBLOCKUNAMALGAMATIONFACTORY_DEF_HPP

#include <Xpetra_Operator.hpp>
#include <Xpetra_BlockedCrsOperator.hpp>

#include "MueLu_SubBlockUnAmalgamationFactory.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {
  
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SubBlockUnAmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SubBlockUnAmalgamationFactory(RCP<const FactoryBase> AFact)
  : AFact_(AFact)
  {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SubBlockUnAmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
   
    currentLevel.DeclareInput("A", AFact_.get(),    this); // sub-block from blocked A

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SubBlockUnAmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const
  {  
    FactoryMonitor m(*this, "SubBlockUnAmalgamationFactory", currentLevel);
    
    RCP<Operator> A = currentLevel.Get< RCP<Operator> >("A", AFact_.get());

    LocalOrdinal  fullblocksize = 1;         // block dim for fixed size blocks
    GlobalOrdinal offset = 0;          // global offset of dof gids
    LocalOrdinal blockid = -1;         // block id in strided map
    LocalOrdinal nStridedOffset = 0;   // DOF offset for strided block id "blockid" (default = 0)
    LocalOrdinal stridedblocksize = fullblocksize; // size of strided block id "blockid" (default = fullblocksize, only if blockid!=-1 stridedblocksize <= fullblocksize)

    // 1) check for blocking/striding information
    if(A->IsView("stridedMaps") &&
	Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
      Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps"); // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
      TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap()) == Teuchos::null,Exceptions::BadCast,"MueLu::CoalesceFactory::Build: cast to strided row map failed.");
      fullblocksize = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getFixedBlockSize(); // TODO shorten code
      offset   = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getOffset();
      blockid  = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getStridedBlockId();
      if (blockid > -1) {
	std::vector<size_t> stridingInfo = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getStridingData();
	for(size_t j=0; j<Teuchos::as<size_t>(blockid); j++) {
	  nStridedOffset += stridingInfo[j];
	}
	stridedblocksize = Teuchos::as<LocalOrdinal>(stridingInfo[blockid]);
      } else {
	stridedblocksize = fullblocksize;
      }
      oldView = A->SwitchToView(oldView);
      GetOStream(Debug, 0) << "SubBlockUnAmalagamationFactory::Build():" << " found fullblocksize=" << fullblocksize << " from strided maps. offset=" << offset << std::endl;
      /*std::cout << "fullblocksize: " << fullblocksize << std::endl;
      std::cout << "offset: " << offset << std::endl;
      std::cout << "blockid: " << blockid << std::endl;
      std::cout << "nStridedOffset: " << nStridedOffset << std::endl;
      std::cout << "stridedblocksize: " << stridedblocksize << std::endl;*/
    } else GetOStream(Debug, 0) << "SubBlockUnAmalagamationFactory::Build(): no striding information available. Use blockdim=1 with offset=0" << std::endl;
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
    GlobalOrdinal cnt_amalRows = 0; // counts number of nodes (rows in amalgamated matrix) on current proc
    for(LocalOrdinal i=0; i<Teuchos::as<LocalOrdinal>(A->getColMap()->getNodeNumElements());i++) {
      // get global DOF id
      GlobalOrdinal gDofId = A->getColMap()->getGlobalElement(i);

      // translate DOFGid to node id
      GlobalOrdinal gNodeId = DOFGid2NodeId(gDofId, A, fullblocksize, offset);

      // gblockid -> gDofId/lDofId
      if(nodegid2dofgids_->count(gNodeId) == 0) {

	// current column DOF gDofId belongs to a node that has not been added
	// to nodeid2dofgids_ yet. Do it now and add ALL DOFs of node gNodeId to
	// unamalgamation information.
	// Note: we use offset and fullblocksize, ie. information from strided maps indirectly
	std::vector<GlobalOrdinal> DOFs; 
	
	DOFs.reserve(stridedblocksize);
	for(LocalOrdinal k=0; k<stridedblocksize; k++) {
	  DOFs.push_back(offset + gNodeId*fullblocksize + nStridedOffset + k);
	}


	(*nodegid2dofgids_)[gNodeId] = DOFs;

	if(A->getRowMap()->isNodeGlobalElement(gDofId)) {
	  gNodeIds->push_back(gNodeId);
	  cnt_amalRows++; // new local block row in amalgamated matrix graph
	}
      }
    }

    // store (un)amalgamation information on current level
    RCP<AmalgamationInfo> amalgamationData = rcp(new AmalgamationInfo());
    amalgamationData->SetAmalgamationParams(nodegid2dofgids_);
    amalgamationData->SetNodeGIDVector(gNodeIds);
    amalgamationData->SetNumberOfNodes(cnt_amalRows);
    currentLevel.Set("UnAmalgamationInfo", amalgamationData, this);
  }
  
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const GlobalOrdinal SubBlockUnAmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DOFGid2NodeId(GlobalOrdinal gid, const RCP<Operator>& A, LocalOrdinal blockSize, const GlobalOrdinal offset) {
    GlobalOrdinal globalblockid = ((GlobalOrdinal) gid - offset) / blockSize;
    return globalblockid;
  }
  
} //namespace MueLu

#endif /* MUELU_SUBBLOCKUNAMALGAMATIONFACTORY_DEF_HPP */

