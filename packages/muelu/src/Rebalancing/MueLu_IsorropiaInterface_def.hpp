/*
 * MueLu_IsorropiaInterface_def.hpp
 *
 *  Created on: Jun 10, 2013
 *      Author: tobias
 */

#ifndef MUELU_ISORROPIAINTERFACE_DEF_HPP_
#define MUELU_ISORROPIAINTERFACE_DEF_HPP_

#include "MueLu_IsorropiaInterface_decl.hpp"

#include <Teuchos_Utils.hpp>
//#include <Teuchos_DefaultMpiComm.hpp> //TODO: fwd decl.
//#include <Teuchos_OpaqueWrapper.hpp>  //TODO: fwd decl.

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_BlockedMap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

#ifdef HAVE_MUELU_ISORROPIA
#include <Isorropia_Exception.hpp>

#ifdef HAVE_MUELU_EPETRA
#include <Xpetra_EpetraCrsGraph.hpp>
#include <Epetra_CrsGraph.h>
#include <Isorropia_EpetraPartitioner.hpp>
#endif

#include <Xpetra_TpetraCrsGraph.hpp>
#endif  // ENDIF HAVE_MUELU_ISORROPIA

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_AmalgamationInfo.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> IsorropiaInterface<LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Factory of the matrix A");
  validParamList->set<RCP<const FactoryBase> >("number of partitions", Teuchos::null, "Instance of RepartitionHeuristicFactory.");
  validParamList->set<RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory of UnAmalgamationInfo");

  return validParamList;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void IsorropiaInterface<LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
  Input(currentLevel, "number of partitions");
  Input(currentLevel, "UnAmalgamationInfo");
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void IsorropiaInterface<LocalOrdinal, GlobalOrdinal, Node>::Build(Level& level) const {
  FactoryMonitor m(*this, "Build", level);
  typedef Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node> BlockMap;

  RCP<Matrix> A                  = Get<RCP<Matrix> >(level, "A");
  RCP<AmalgamationInfo> amalInfo = Get<RCP<AmalgamationInfo> >(level, "UnAmalgamationInfo");
  GO numParts                    = Get<int>(level, "number of partitions");

  RCP<const Map> rowMap = A->getRowMap();
  RCP<const Map> colMap = A->getColMap();

  if (numParts == 1 || numParts == -1) {
    // Running on one processor, so decomposition is the trivial one, all zeros.
    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, true);
    Set(level, "AmalgamatedPartition", decomposition);
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

  LO blockdim       = 1;                       // block dim for fixed size blocks
  GO indexBase      = rowMap->getIndexBase();  // index base of maps
  GO offset         = 0;
  LO blockid        = -1;  // block id in strided map
  LO nStridedOffset = 0;   // DOF offset for strided block id "blockid" (default = 0)
  // LO stridedblocksize = blockdim; // size of strided block id "blockid" (default = fullblocksize, only if blockid!=-1 stridedblocksize <= fullblocksize)

  // 1) check for blocking/striding information
  //    fill above variables
  if (A->IsView("stridedMaps") &&
      Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
    Xpetra::viewLabel_t oldView  = A->SwitchToView("stridedMaps");  // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
    RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap());
    TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null, Exceptions::BadCast, "MueLu::IsorropiaInterface::Build: cast to strided row map failed.");
    blockdim = strMap->getFixedBlockSize();
    offset   = strMap->getOffset();
    blockid  = strMap->getStridedBlockId();
    if (blockid > -1) {
      std::vector<size_t> stridingInfo = strMap->getStridingData();
      for (size_t j = 0; j < Teuchos::as<size_t>(blockid); j++)
        nStridedOffset += stridingInfo[j];
      // stridedblocksize = Teuchos::as<LocalOrdinal>(stridingInfo[blockid]);

    }  // else {
    //  stridedblocksize = blockdim;
    //}
    oldView = A->SwitchToView(oldView);
    // GetOStream(Statistics0) << "IsorropiaInterface::Build():" << " found blockdim=" << blockdim << " from strided maps (blockid=" << blockid << ", strided block size=" << stridedblocksize << "). offset=" << offset << std::endl;
  } else
    GetOStream(Statistics0) << "IsorropiaInterface::Build(): no striding information available. Use blockdim=1 with offset=0" << std::endl;

  // 2) get row map for amalgamated matrix (graph of A)
  //    with same distribution over all procs as row map of A
  RCP<const Map> nodeMap         = amalInfo->getNodeRowMap();
  RCP<const BlockedMap> bnodeMap = Teuchos::rcp_dynamic_cast<const BlockedMap>(nodeMap);
  if (!bnodeMap.is_null()) nodeMap = bnodeMap->getMap();

  GetOStream(Statistics0) << "IsorropiaInterface:Build(): nodeMap " << nodeMap->getLocalNumElements() << "/" << nodeMap->getGlobalNumElements() << " elements" << std::endl;

  // 3) create graph of amalgamated matrix
  RCP<CrsGraph> crsGraph = CrsGraphFactory::Build(nodeMap, A->getLocalMaxNumRowEntries() * blockdim);

  // 4) do amalgamation. generate graph of amalgamated matrix
  for (LO row = 0; row < Teuchos::as<LO>(A->getRowMap()->getLocalNumElements()); row++) {
    // get global DOF id
    GO grid = rowMap->getGlobalElement(row);

    // translate grid to nodeid
    // JHU 2019-20-May this is identical to AmalgamationFactory::DOFGid2NodeId(), and is done
    // to break a circular dependence when libraries are built statically
    GO nodeId = (grid - offset - indexBase) / blockdim + indexBase;

    size_t nnz = A->getNumEntriesInLocalRow(row);
    Teuchos::ArrayView<const LO> indices;
    Teuchos::ArrayView<const SC> vals;
    A->getLocalRowView(row, indices, vals);

    RCP<std::vector<GO> > cnodeIds = Teuchos::rcp(new std::vector<GO>);  // global column block ids
    LO realnnz                     = 0;
    for (LO col = 0; col < Teuchos::as<LO>(nnz); col++) {
      GO gcid = colMap->getGlobalElement(indices[col]);  // global column id

      if (vals[col] != 0.0) {
        // JHU 2019-20-May this is identical to AmalgamationFactory::DOFGid2NodeId(), and is done
        // to break a circular dependence when libraries are built statically
        GO cnodeId = (gcid - offset - indexBase) / blockdim + indexBase;
        cnodeIds->push_back(cnodeId);
        realnnz++;  // increment number of nnz in matrix row
      }
    }

    Teuchos::ArrayRCP<GO> arr_cnodeIds = Teuchos::arcp(cnodeIds);

    if (arr_cnodeIds.size() > 0)
      crsGraph->insertGlobalIndices(nodeId, arr_cnodeIds());
  }
  // fill matrix graph
  crsGraph->fillComplete(nodeMap, nodeMap);

#ifdef HAVE_MPI
#ifdef HAVE_MUELU_ISORROPIA

  // prepare parameter list for Isorropia
  Teuchos::ParameterList paramlist;
  paramlist.set("NUM PARTS", toString(numParts));

  /*STRUCTURALLY SYMMETRIC [NO/yes] (is symmetrization required?)
  PARTITIONING METHOD [block/cyclic/random/rcb/rib/hsfc/graph/HYPERGRAPH]
  NUM PARTS [int k] (global number of parts)
  IMBALANCE TOL [float tol] (1.0 is perfect balance)
  BALANCE OBJECTIVE [ROWS/nonzeros]
  */
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  sublist.set("LB_APPROACH", "PARTITION");

#ifdef HAVE_MUELU_EPETRA
  RCP<Xpetra::EpetraCrsGraphT<GO, Node> > epCrsGraph = Teuchos::rcp_dynamic_cast<Xpetra::EpetraCrsGraphT<GO, Node> >(crsGraph);
  if (epCrsGraph != Teuchos::null) {
    RCP<const Epetra_CrsGraph> epetraCrsGraph = epCrsGraph->getEpetra_CrsGraph();

    RCP<Isorropia::Epetra::Partitioner> isoPart = Teuchos::rcp(new Isorropia::Epetra::Partitioner(epetraCrsGraph, paramlist));

    int size         = 0;
    const int* array = NULL;
    isoPart->extractPartsView(size, array);

    TEUCHOS_TEST_FOR_EXCEPTION(size != Teuchos::as<int>(nodeMap->getLocalNumElements()), Exceptions::RuntimeError, "length of array returned from extractPartsView does not match local length of rowMap");

    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(nodeMap, false);
    ArrayRCP<GO> decompEntries                         = decomposition->getDataNonConst(0);

    // fill vector with amalgamated information about partitioning
    for (int i = 0; i < size; i++) {
      decompEntries[i] = Teuchos::as<GO>(array[i]);
    }

    Set(level, "AmalgamatedPartition", decomposition);
  }
#endif  // ENDIF HAVE_MUELU_EPETRA

#ifdef HAVE_MUELU_INST_DOUBLE_INT_INT
  RCP<Xpetra::TpetraCrsGraph<LO, GO, Node> > tpCrsGraph = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsGraph<LO, GO, Node> >(crsGraph);
  TEUCHOS_TEST_FOR_EXCEPTION(tpCrsGraph != Teuchos::null, Exceptions::RuntimeError, "Tpetra is not supported with Isorropia.");
#else
  TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "Isorropia is an interface to Zoltan which only has support for LO=GO=int and SC=double.");
#endif  // ENDIF HAVE_MUELU_INST_DOUBLE_INT_INT
#endif  // HAVE_MUELU_ISORROPIA
#else   // if we don't have MPI

  // Running on one processor, so decomposition is the trivial one, all zeros.
  RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, true);
  Set(level, "AmalgamatedPartition", decomposition);

#endif  // HAVE_MPI
  // throw a more helpful error message if something failed
  // TEUCHOS_TEST_FOR_EXCEPTION(!level.IsAvailable("AmalgamatedPartition"), Exceptions::RuntimeError, "IsorropiaInterface::Build : no \'Partition\' vector available on level. Isorropia failed to build a partition of the non-repartitioned graph of A. Please make sure, that Isorropia is correctly compiled (Epetra/Tpetra).");

}  // Build()

}  // namespace MueLu

//#endif //if defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MPI)

#endif /* MUELU_ISORROPIAINTERFACE_DEF_HPP_ */
