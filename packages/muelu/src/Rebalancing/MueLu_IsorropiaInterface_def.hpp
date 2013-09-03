/*
 * MueLu_IsorropiaInterface_def.hpp
 *
 *  Created on: Jun 10, 2013
 *      Author: tobias
 */

#ifndef MUELU_ISORROPIAINTERFACE_DEF_HPP_
#define MUELU_ISORROPIAINTERFACE_DEF_HPP_

#include "MueLu_IsorropiaInterface_decl.hpp"
#if defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MPI)

#include <Teuchos_Utils.hpp>
#include <Teuchos_DefaultMpiComm.hpp> //TODO: fwd decl.
#include <Teuchos_OpaqueWrapper.hpp>  //TODO: fwd decl.

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsGraphFactory.hpp>

#ifdef HAVE_MUELU_ISORROPIA
#include <Isorropia_Exception.hpp>


#ifdef HAVE_MUELU_EPETRA
#include <Xpetra_EpetraCrsGraph.hpp>
#include <Epetra_CrsGraph.h>
#include <Isorropia_EpetraPartitioner.hpp>
#endif

#ifdef HAVE_MUELU_TPETRA
#include <Xpetra_TpetraCrsGraph.hpp>
#include <Tpetra_CrsGraph.hpp>
#ifdef HAVE_ISORROPIA_TPETRA
#include <Isorropia_TpetraPartitioner.hpp>
#endif // HAVE_ISORROPIA_TPETRA
#endif
#endif // ENDIF HAVE_MUELU_ISORROPIA

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

 template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
 RCP<const ParameterList> IsorropiaInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A",                    Teuchos::null, "Factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("number of partitions", Teuchos::null, "(advanced) Factory computing the number of partition.");

    return validParamList;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void IsorropiaInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "number of partitions");
  } //DeclareInput()

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void IsorropiaInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &level) const {
    FactoryMonitor m(*this, "Build", level);

    RCP<Matrix>      A             = Get< RCP<Matrix> >     (level, "A");
    GO               numPartitions = Get<GO>                (level, "number of partitions");

    RCP<const Map> rowMap        = A->getRowMap();
    RCP<const Map> colMap = A->getColMap();

    if (numPartitions == 1) {
      // Running on one processor, so decomposition is the trivial one, all zeros.
      RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, true);
      Set(level, "Partition", decomposition);
      return;
    }

    // ok, reconstruct graph information of matrix A
    // Note, this is the non-rebalanced matrix A and we need the graph
    // of the non-rebalanced matrix A. We cannot make use of the "Graph"
    // that is/will be built for the aggregates later for several reasons
    // 1) the "Graph" for the aggregates is meant to be based on the rebalanced matrix A
    // 2) we cannot call a CoalesceDropFactory::Build here since this is very complicated and
    //    completely messes up the whole factory chain
    // 3) CoalesceDropFactory is meant to provide some minimal Graph information for the aggregation
    //    (LWGraph), but here we need the full CrsGraph for Isorropia as input

    // That is, why this code is somewhat redundant to CoalesceDropFactory


    LO blockdim = 1;                          // block dim for fixed size blocks
    GO indexBase = rowMap->getIndexBase();    // index base of maps
    GO offset    = 0;
    LO blockid          = -1;  // block id in strided map
    LO nStridedOffset   = 0;   // DOF offset for strided block id "blockid" (default = 0)
    LO stridedblocksize = blockdim; // size of strided block id "blockid" (default = fullblocksize, only if blockid!=-1 stridedblocksize <= fullblocksize)

    // 1) check for blocking/striding information
    //    fill above variables
    if(A->IsView("stridedMaps") &&
       Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
      Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps"); // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
      RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap());
      TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null,Exceptions::BadCast,"MueLu::IsorropiaInterface::Build: cast to strided row map failed.");
      blockdim = strMap->getFixedBlockSize();
      offset   = strMap->getOffset();
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
      GetOStream(Statistics0, -1) << "IsorropiaInterface::Build():" << " found blockdim=" << blockdim << " from strided maps (blockid=" << blockid << ", strided block size=" << stridedblocksize << "). offset=" << offset << std::endl;
    } else GetOStream(Statistics0, -1) << "IsorropiaInterface::Build(): no striding information available. Use blockdim=1 with offset=0" << std::endl;

    // 2) build (un)amalgamation information

    // gNodeIds contains a list of all global node ids on current processor
    RCP<std::vector<GlobalOrdinal> > gNodeIds;
    gNodeIds = Teuchos::rcp(new std::vector<GlobalOrdinal>);
    gNodeIds->empty();

    // in nodegid2dofgids_ for each node on the current proc a vector of
    // the corresponding DOFs gids is stored.
    // The map contains all nodes the current proc has connections to (including
    // nodes that are stored on other procs when there are off-diagonal entries in A)
    std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > nodegid2dofgids;

    GlobalOrdinal cnt_amalRows = 0; // counts number of nodes (rows in amalgamated matrix) on current proc
    LocalOrdinal nColEle = Teuchos::as<LocalOrdinal>(colMap->getNodeNumElements());
    for (LocalOrdinal i = 0; i < nColEle; i++) {
      // get global DOF id
      GlobalOrdinal gDofId = colMap->getGlobalElement(i);

      // translate DOFGid to node id
      GlobalOrdinal gNodeId = AmalgamationFactory::DOFGid2NodeId(gDofId, blockdim, offset, indexBase);

      // gblockid -> gDofId/lDofId
      if (nodegid2dofgids.count(gNodeId) == 0) {

        // current column DOF gDofId belongs to a node that has not been added
        // to nodeid2dofgids_ yet. Do it now and add ALL DOFs of node gNodeId to
        // unamalgamation information.
        // Note: we use offset and fullblocksize, ie. information from strided maps indirectly
        std::vector<GlobalOrdinal> DOFs;

        DOFs.reserve(stridedblocksize);
        for (LocalOrdinal k = 0; k < stridedblocksize; k++) {
          // here, the assumption is, that the node map has the same indexBase as the dof map
          //                            this is the node map index base                    this is the dof map index base
          GO gDofIndex = offset + (gNodeId-indexBase)*blockdim + nStridedOffset + k + indexBase;
          if (colMap->isNodeGlobalElement(gDofIndex))
            DOFs.push_back(gDofIndex);
        }

        nodegid2dofgids[gNodeId] = DOFs;

        if (rowMap->isNodeGlobalElement(gDofId)) {
          gNodeIds->push_back(gNodeId);
          cnt_amalRows++; // new local block row in amalgamated matrix graph
        }
      }
    }

    // inter processor communication: sum up number of block ids
    GO num_blockids = 0;
    Teuchos::reduceAll<int,GO>(*(A->getRowMap()->getComm()),Teuchos::REDUCE_SUM, cnt_amalRows, Teuchos::ptr(&num_blockids) );
    GetOStream(Statistics0, -1) << "IsorropiaInterface::SetupAmalgamationData()" << " # of amalgamated blocks=" << num_blockids << std::endl;

    // 3) generate row map for amalgamated matrix (graph of A)
    //    with same distribution over all procs as row map of A

    Teuchos::ArrayRCP<GO> arr_gNodeIds = Teuchos::arcp( gNodeIds );
    Teuchos::RCP<Map> nodeMap = MapFactory::Build(A->getRowMap()->lib(), num_blockids, arr_gNodeIds(), indexBase, A->getRowMap()->getComm()); // note: nodeMap has same indexBase as row map of A (=dof map)
    GetOStream(Statistics0, -1) << "IsorropiaInterface: nodeMap " << nodeMap->getNodeNumElements() << "/" << nodeMap->getGlobalNumElements() << " local/global elements" << std::endl;

    // 4) create graph of amalgamated matrix
    RCP<CrsGraph> crsGraph = CrsGraphFactory::Build(nodeMap, 10, Xpetra::DynamicProfile);

    // 5) do amalgamation. generate graph of amalgamated matrix
    for(LO row=0; row<Teuchos::as<LO>(A->getRowMap()->getNodeNumElements()); row++) {
      // get global DOF id
      GO grid = rowMap->getGlobalElement(row);

      // translate grid to nodeid
      GO nodeId = AmalgamationFactory::DOFGid2NodeId(grid, blockdim, offset, indexBase);

      size_t nnz = A->getNumEntriesInLocalRow(row);
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const SC> vals;
      A->getLocalRowView(row, indices, vals);
      //TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::CoalesceFactory::Amalgamate: number of nonzeros not equal to number of indices? Error.");

      RCP<std::vector<GO> > cnodeIds = Teuchos::rcp(new std::vector<GO>);  // global column block ids
      LO realnnz = 0;
      for(LO col=0; col<Teuchos::as<LO>(nnz); col++) {
        //TEUCHOS_TEST_FOR_EXCEPTION(A->getColMap()->isNodeLocalElement(indices[col])==false,Exceptions::RuntimeError, "MueLu::CoalesceFactory::Amalgamate: Problem with columns. Error.");
        GO gcid = colMap->getGlobalElement(indices[col]); // global column id

        if(vals[col]!=0.0) {
          GO cnodeId = AmalgamationFactory::DOFGid2NodeId(gcid, blockdim, offset, indexBase);
          cnodeIds->push_back(cnodeId);
          realnnz++; // increment number of nnz in matrix row
        }
      }

      Teuchos::ArrayRCP<GO> arr_cnodeIds = Teuchos::arcp( cnodeIds );

      //TEUCHOS_TEST_FOR_EXCEPTION(crsGraph->getRowMap()->isNodeGlobalElement(nodeId)==false,Exceptions::RuntimeError, "MueLu::CoalesceFactory::Amalgamate: global row id does not belong to current proc. Error.");
      if(arr_cnodeIds.size() > 0 )
        crsGraph->insertGlobalIndices(nodeId, arr_cnodeIds());
    }
    // fill matrix graph
    crsGraph->fillComplete(nodeMap,nodeMap);

#ifdef HAVE_MUELU_ISORROPIA

    // prepare parameter list for Isorropia
    Teuchos::ParameterList paramlist;
    std::stringstream ss; ss << numPartitions;
    paramlist.set("NUM PARTS",ss.str());


    /*STRUCTURALLY SYMMETRIC [NO/yes] (is symmetrization required?)
    PARTITIONING METHOD [block/cyclic/random/rcb/rib/hsfc/graph/HYPERGRAPH]
    NUM PARTS [int k] (global number of parts)
    IMBALANCE TOL [float tol] (1.0 is perfect balance)
    BALANCE OBJECTIVE [ROWS/nonzeros]
    */
    Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
    sublist.set("LB_APPROACH", "PARTITION");



#ifdef HAVE_MUELU_EPETRA
    RCP< Xpetra::EpetraCrsGraph> epCrsGraph = Teuchos::rcp_dynamic_cast<Xpetra::EpetraCrsGraph>(crsGraph);
    if(epCrsGraph != Teuchos::null) {
      RCP< const Epetra_CrsGraph> epetraCrsGraph = epCrsGraph->getEpetra_CrsGraph();

      RCP<Isorropia::Epetra::Partitioner> isoPart = Teuchos::rcp(new Isorropia::Epetra::Partitioner(epetraCrsGraph,paramlist));

      int size = 0;
      const int* array = NULL;
      isoPart->extractPartsView(size,array);

      // consistency check
      TEUCHOS_TEST_FOR_EXCEPTION(size*stridedblocksize/*blockdim*/ != Teuchos::as<int>(rowMap->getNodeNumElements()), Exceptions::RuntimeError, "length of array returned from extractPartsView does not match local length of rowMap");

      RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, false);
      ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);

      // fill vector with information about partitioning
      // TODO: we assume simple block maps here
      // TODO: adapt this to usage of nodegid2dofgids
      for(int i = 0; i<size; i++) {
        // not fully sure about this. We're filling local ids in the decomposition vector with
        // the results stored in array. The decomposition vector is created using the rowMap of A
#if 0
        for (LO j=0; j<stridedblocksize/*blockdim*/; j++) {
          decompEntries[i*stridedblocksize/*blockdim*/ + j] = Teuchos::as<GO>(array[i]);
        }
#else
        // transform local node id to global node id.
        GO gNodeId = nodeMap->getGlobalElement(i);

        // extract global DOF ids that belong to gNodeId
        std::vector<GlobalOrdinal> DOFs = nodegid2dofgids[gNodeId];

        for(size_t j=0; j<DOFs.size(); j++) {
          // transform global DOF ids to local DOF ids using rowMap
          // note: The vector decomposition is based on rowMap
          //LO lDofId = rowMap->getLocalElement(DOFs[j]);

          // put the same domain id to all DOFs of the same node
          decompEntries[i*stridedblocksize + j] = Teuchos::as<GO>(array[i]);
        }
#endif

      }

      Set(level, "Partition", decomposition);
    }
#endif // ENDIF HAVE_MUELU_EPETRA

#ifdef HAVE_MUELU_TPETRA
#ifdef HAVE_MUELU_INST_DOUBLE_INT_INT

    RCP< Xpetra::TpetraCrsGraph<LO, GO, Node, LocalMatOps> > tpCrsGraph = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsGraph<LO, GO, Node, LocalMatOps> >(crsGraph);
    if(tpCrsGraph != Teuchos::null) {
#ifdef HAVE_ISORROPIA_TPETRA
      RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpetraCrsGraph = tpCrsGraph->getTpetra_CrsGraph();
      RCP<Isorropia::Tpetra::Partitioner<Node> > isoPart = rcp(new Isorropia::Tpetra::Partitioner<Node>(tpetraCrsGraph, paramlist));

      int size = 0;
      const int* array = NULL;
      isoPart->extractPartsView(size,array);

      TEUCHOS_TEST_FOR_EXCEPTION(size*stridedblocksize/*blockdim*/ != Teuchos::as<int>(rowMap->getNodeNumElements()), Exceptions::RuntimeError, "length of array returned from extractPartsView does not match local length of rowMap");

      RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, false);
      ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);

      // fill vector with information about partitioning
      // TODO: we assume simple block maps here
      // TODO: adapt this to usage of nodegid2dofgids
      for(int i = 0; i<size; i++) {
       for (LO j=0; j<stridedblocksize/*blockdim*/; j++) {
          decompEntries[i*stridedblocksize/*blockdim*/ + j] = Teuchos::as<GO>(array[i]);
        }
      }

      Set(level, "Partition", decomposition);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "Tpetra is not enabled for Isorropia. Recompile Isorropia with Tpetra support.");
#endif // ENDIF HAVE_ISORROPIA_TPETRA
    }
#else
    TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "Isorropia is an interface to Zoltan which only has support for LO=GO=int and SC=double.");
#endif // ENDIF HAVE_MUELU_INST_DOUBLE_INT_INT
#endif // ENDIF HAVE_MUELU_TPETRA
#endif // HAVE_MUELU_ISORROPIA

    // throw a more helpful error message if something failed
    //TEUCHOS_TEST_FOR_EXCEPTION(!level.IsAvailable("Partition"), Exceptions::RuntimeError, "IsorropiaInterface::Build : no \'Partition\' vector available on level. Isorropia failed to build a partition of the non-repartitioned graph of A. Please make sure, that Isorropia is correctly compiled (Epetra/Tpetra).");

  } //Build()



} //namespace MueLu

#endif //if defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MPI)


#endif /* MUELU_ISORROPIAINTERFACE_DEF_HPP_ */
