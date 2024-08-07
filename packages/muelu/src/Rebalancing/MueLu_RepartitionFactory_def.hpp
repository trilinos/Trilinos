// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REPARTITIONFACTORY_DEF_HPP
#define MUELU_REPARTITIONFACTORY_DEF_HPP

#include <algorithm>
#include <iostream>
#include <sstream>

#include "MueLu_RepartitionFactory_decl.hpp"  // TMP JG NOTE: before other includes, otherwise I cannot test the fwd declaration in _def

#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Details_MpiTypeTraits.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include "MueLu_Utilities.hpp"

#include "MueLu_CloneRepartitionInterface.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RepartitionFactory() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~RepartitionFactory() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("repartition: print partition distribution");
  SET_VALID_ENTRY("repartition: remap parts");
  SET_VALID_ENTRY("repartition: remap num values");
  SET_VALID_ENTRY("repartition: remap accept partition");
  SET_VALID_ENTRY("repartition: node repartition level");
  SET_VALID_ENTRY("repartition: save importer");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Factory of the matrix A");
  validParamList->set<RCP<const FactoryBase> >("number of partitions", Teuchos::null, "Instance of RepartitionHeuristicFactory.");
  validParamList->set<RCP<const FactoryBase> >("Partition", Teuchos::null, "Factory of the partition");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
  Input(currentLevel, "number of partitions");
  Input(currentLevel, "Partition");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  const Teuchos::ParameterList& pL = GetParameterList();
  // Access parameters here to make sure that we set the parameter entry flag to "used" even in case of short-circuit evaluation.
  // TODO (JG): I don't really know if we want to do this.
  bool remapPartitions = pL.get<bool>("repartition: remap parts");

  // TODO: We only need a CrsGraph. This class does not have to be templated on Scalar types.
  RCP<Matrix> A = Get<RCP<Matrix> >(currentLevel, "A");
  if (A == Teuchos::null) {
    Set<RCP<const Import> >(currentLevel, "Importer", Teuchos::null);
    return;
  }
  RCP<const Map> rowMap     = A->getRowMap();
  GO indexBase              = rowMap->getIndexBase();
  Xpetra::UnderlyingLib lib = rowMap->lib();

  RCP<const Teuchos::Comm<int> > origComm = rowMap->getComm();
  RCP<const Teuchos::Comm<int> > comm     = origComm;

  int myRank   = comm->getRank();
  int numProcs = comm->getSize();

  RCP<const Teuchos::MpiComm<int> > tmpic = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
  TEUCHOS_TEST_FOR_EXCEPTION(tmpic == Teuchos::null, Exceptions::RuntimeError, "Cannot cast base Teuchos::Comm to Teuchos::MpiComm object.");
  RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm = tmpic->getRawMpiComm();

  /////
  int numPartitions = Get<int>(currentLevel, "number of partitions");

  // ======================================================================================================
  // Construct decomposition vector
  // ======================================================================================================
  RCP<GOVector> decomposition = Get<RCP<GOVector> >(currentLevel, "Partition");

  // check which factory provides "Partition"
  if (remapPartitions == true && Teuchos::rcp_dynamic_cast<const CloneRepartitionInterface>(GetFactory("Partition")) != Teuchos::null) {
    // if "Partition" is provided by a CloneRepartitionInterface class we have to switch of remapPartitions
    // as we can assume the processor Ids in Partition to be the expected ones. If we would do remapping we
    // would get different processors for the different blocks which screws up matrix-matrix multiplication.
    remapPartitions = false;
  }

  // check special cases
  if (numPartitions == 1) {
    // Trivial case: decomposition is the trivial one, all zeros. We skip the call to Zoltan_Interface
    // (this is mostly done to avoid extra output messages, as even if we didn't skip there is a shortcut
    // in Zoltan[12]Interface).
    // TODO: We can probably skip more work in this case (like building all extra data structures)
    GetOStream(Runtime0) << "Only one partition: Skip call to the repartitioner." << std::endl;
  } else if (numPartitions == -1) {
    // No repartitioning necessary: decomposition should be Teuchos::null
    GetOStream(Runtime0) << "No repartitioning necessary: partitions were left unchanged by the repartitioner" << std::endl;
    Set<RCP<const Import> >(currentLevel, "Importer", Teuchos::null);
    return;
  }

  // If we're doing node away, we need to be sure to get the mapping to the NodeComm's rank 0.
  const int nodeRepartLevel = pL.get<int>("repartition: node repartition level");
  if (currentLevel.GetLevelID() == nodeRepartLevel) {
    // NodePartitionInterface returns the *ranks* of the guy who gets the info, not the *partition number*
    // In a sense, we've already done remap here.

    // FIXME: We need a low-comm import construction
    remapPartitions = false;
  }

  // ======================================================================================================
  // Remap if necessary
  // ======================================================================================================
  // From a user perspective, we want user to not care about remapping, thinking of it as only a performance feature.
  // There are two problems, however.
  // (1) Next level aggregation depends on the order of GIDs in the vector, if one uses "natural" or "random" orderings.
  //     This also means that remapping affects next level aggregation, despite the fact that the _set_ of GIDs for
  //     each partition is the same.
  // (2) Even with the fixed order of GIDs, the remapping may influence the aggregation for the next-next level.
  //     Let us consider the following example. Lets assume that when we don't do remapping, processor 0 would have
  //     GIDs {0,1,2}, and processor 1 GIDs {3,4,5}, and if we do remapping processor 0 would contain {3,4,5} and
  //     processor 1 {0,1,2}. Now, when we run repartitioning algorithm on the next level (say Zoltan1 RCB), it may
  //     be dependent on whether whether it is [{0,1,2}, {3,4,5}] or [{3,4,5}, {0,1,2}]. Specifically, the tie-breaking
  //     algorithm can resolve these differently. For instance, running
  //         mpirun -np 5 ./MueLu_ScalingTestParamList.exe --xml=easy_sa.xml --nx=12 --ny=12 --nz=12
  //     with
  //         <ParameterList name="MueLu">
  //           <Parameter name="coarse: max size"                type="int"      value="1"/>
  //           <Parameter name="repartition: enable"             type="bool"     value="true"/>
  //           <Parameter name="repartition: min rows per proc"  type="int"      value="2"/>
  //           <ParameterList name="level 1">
  //             <Parameter name="repartition: remap parts"      type="bool"     value="false/true"/>
  //           </ParameterList>
  //         </ParameterList>
  //     produces different repartitioning for level 2.
  //     This different repartitioning may then escalate into different aggregation for the next level.
  //
  // We fix (1) by fixing the order of GIDs in a vector by sorting the resulting vector.
  // Fixing (2) is more complicated.
  // FIXME: Fixing (2) in Zoltan may not be enough, as we may use some arbitration in MueLu,
  // for instance with CoupledAggregation. What we really need to do is to use the same order of processors containing
  // the same order of GIDs. To achieve that, the newly created subcommunicator must be conforming with the order. For
  // instance, if we have [{0,1,2}, {3,4,5}], we create a subcommunicator where processor 0 gets rank 0, and processor 1
  // gets rank 1. If, on the other hand, we have [{3,4,5}, {0,1,2}], we assign rank 1 to processor 0, and rank 0 to processor 1.
  // This rank permutation requires help from Epetra/Tpetra, both of which have no such API in place.
  // One should also be concerned that if we had such API in place, rank 0 in subcommunicator may no longer be rank 0 in
  // MPI_COMM_WORLD, which may lead to issues for logging.
  if (remapPartitions) {
    SubFactoryMonitor m1(*this, "DeterminePartitionPlacement", currentLevel);

    bool acceptPartition = pL.get<bool>("repartition: remap accept partition");
    bool allSubdomainsAcceptPartitions;
    int localNumAcceptPartition = acceptPartition;
    int globalNumAcceptPartition;
    MueLu_sumAll(comm, localNumAcceptPartition, globalNumAcceptPartition);
    GetOStream(Statistics2) << "Number of ranks that accept partitions: " << globalNumAcceptPartition << std::endl;
    if (globalNumAcceptPartition < numPartitions) {
      GetOStream(Warnings0) << "Not enough ranks are willing to accept a partition, allowing partitions on all ranks." << std::endl;
      acceptPartition               = true;
      allSubdomainsAcceptPartitions = true;
    } else if (numPartitions > numProcs) {
      // We are trying to repartition to a larger communicator.
      allSubdomainsAcceptPartitions = true;
    } else {
      allSubdomainsAcceptPartitions = false;
    }

    DeterminePartitionPlacement(*A, *decomposition, numPartitions, acceptPartition, allSubdomainsAcceptPartitions);
  }

  // ======================================================================================================
  // Construct importer
  // ======================================================================================================
  // At this point, the following is true:
  //  * Each processors owns 0 or 1 partitions
  //  * If a processor owns a partition, that partition number is equal to the processor rank
  //  * The decomposition vector contains the partitions ids that the corresponding GID belongs to

  ArrayRCP<const GO> decompEntries;
  if (decomposition->getLocalLength() > 0)
    decompEntries = decomposition->getData(0);

#ifdef HAVE_MUELU_DEBUG
  // Test range of partition ids
  int incorrectRank = -1;
  for (int i = 0; i < decompEntries.size(); i++)
    if (decompEntries[i] >= numProcs || decompEntries[i] < 0) {
      incorrectRank = myRank;
      break;
    }

  int incorrectGlobalRank = -1;
  MueLu_maxAll(comm, incorrectRank, incorrectGlobalRank);
  TEUCHOS_TEST_FOR_EXCEPTION(incorrectGlobalRank > -1, Exceptions::RuntimeError, "pid " + Teuchos::toString(incorrectGlobalRank) + " encountered a partition number is that out-of-range");
#endif

  Array<GO> myGIDs;
  myGIDs.reserve(decomposition->getLocalLength());

  // Step 0: Construct mapping
  //    part number -> GIDs I own which belong to this part
  // NOTE: my own part GIDs are not part of the map
  typedef std::map<GO, Array<GO> > map_type;
  map_type sendMap;
  for (LO i = 0; i < decompEntries.size(); i++) {
    GO id  = decompEntries[i];
    GO GID = rowMap->getGlobalElement(i);

    if (id == myRank)
      myGIDs.push_back(GID);
    else
      sendMap[id].push_back(GID);
  }
  decompEntries = Teuchos::null;

  if (IsPrint(Statistics2)) {
    GO numLocalKept = myGIDs.size(), numGlobalKept, numGlobalRows = A->getGlobalNumRows();
    MueLu_sumAll(comm, numLocalKept, numGlobalKept);
    GetOStream(Statistics2) << "Unmoved rows: " << numGlobalKept << " / " << numGlobalRows << " (" << 100 * Teuchos::as<double>(numGlobalKept) / numGlobalRows << "%)" << std::endl;
  }

  int numSend = sendMap.size(), numRecv;

  // Arrayify map keys
  Array<GO> myParts(numSend), myPart(1);
  int cnt   = 0;
  myPart[0] = myRank;
  for (typename map_type::const_iterator it = sendMap.begin(); it != sendMap.end(); it++)
    myParts[cnt++] = it->first;

  // Step 1: Find out how many processors send me data
  // partsIndexBase starts from zero, as the processors ids start from zero
  {
    SubFactoryMonitor m1(*this, "Mapping Step 1", currentLevel);
    GO partsIndexBase       = 0;
    RCP<Map> partsIHave     = MapFactory ::Build(lib, Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), myParts(), partsIndexBase, comm);
    RCP<Map> partsIOwn      = MapFactory ::Build(lib, numProcs, myPart(), partsIndexBase, comm);
    RCP<Export> partsExport = ExportFactory::Build(partsIHave, partsIOwn);

    RCP<GOVector> partsISend    = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(partsIHave);
    RCP<GOVector> numPartsIRecv = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(partsIOwn);
    if (numSend) {
      ArrayRCP<GO> partsISendData = partsISend->getDataNonConst(0);
      for (int i = 0; i < numSend; i++)
        partsISendData[i] = 1;
    }
    (numPartsIRecv->getDataNonConst(0))[0] = 0;

    numPartsIRecv->doExport(*partsISend, *partsExport, Xpetra::ADD);
    numRecv = (numPartsIRecv->getData(0))[0];
  }

  // Step 2: Get my GIDs from everybody else
  MPI_Datatype MpiType = Teuchos::Details::MpiTypeTraits<GO>::getType();
  int msgTag           = 12345;  // TODO: use Comm::dup for all internal messaging

  // Post sends
  Array<MPI_Request> sendReqs(numSend);
  cnt = 0;
  for (typename map_type::iterator it = sendMap.begin(); it != sendMap.end(); it++)
    MPI_Isend(static_cast<void*>(it->second.getRawPtr()), it->second.size(), MpiType, Teuchos::as<GO>(it->first), msgTag, *rawMpiComm, &sendReqs[cnt++]);

  map_type recvMap;
  size_t totalGIDs = myGIDs.size();
  for (int i = 0; i < numRecv; i++) {
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, msgTag, *rawMpiComm, &status);

    // Get rank and number of elements from status
    int fromRank = status.MPI_SOURCE, count;
    MPI_Get_count(&status, MpiType, &count);

    recvMap[fromRank].resize(count);
    MPI_Recv(static_cast<void*>(recvMap[fromRank].getRawPtr()), count, MpiType, fromRank, msgTag, *rawMpiComm, &status);

    totalGIDs += count;
  }

  // Do waits on send requests
  if (numSend) {
    Array<MPI_Status> sendStatuses(numSend);
    MPI_Waitall(numSend, sendReqs.getRawPtr(), sendStatuses.getRawPtr());
  }

  // Merge GIDs
  myGIDs.reserve(totalGIDs);
  for (typename map_type::const_iterator it = recvMap.begin(); it != recvMap.end(); it++) {
    int offset = myGIDs.size(), len = it->second.size();
    if (len) {
      myGIDs.resize(offset + len);
      memcpy(myGIDs.getRawPtr() + offset, it->second.getRawPtr(), len * sizeof(GO));
    }
  }
  // NOTE 2: The general sorting algorithm could be sped up by using the knowledge that original myGIDs and all received chunks
  // (i.e. it->second) are sorted. Therefore, a merge sort would work well in this situation.
  std::sort(myGIDs.begin(), myGIDs.end());

  // Step 3: Construct importer
  RCP<Map> newRowMap;
  {
    SubFactoryMonitor m1(*this, "Map construction", currentLevel);
    newRowMap = MapFactory ::Build(lib, rowMap->getGlobalNumElements(), myGIDs(), indexBase, origComm);
  }

  RCP<const Import> rowMapImporter;

  RCP<const BlockedMap> blockedRowMap = Teuchos::rcp_dynamic_cast<const BlockedMap>(rowMap);

  {
    SubFactoryMonitor m1(*this, "Import construction", currentLevel);
    // Generate a nonblocked rowmap if we need one
    if (blockedRowMap.is_null())
      rowMapImporter = ImportFactory::Build(rowMap, newRowMap);
    else {
      rowMapImporter = ImportFactory::Build(blockedRowMap->getMap(), newRowMap);
    }
  }

  // If we're running BlockedCrs we should chop up the newRowMap into a newBlockedRowMap here (and do likewise for importers)
  if (!blockedRowMap.is_null()) {
    SubFactoryMonitor m1(*this, "Blocking newRowMap and Importer", currentLevel);
    RCP<const BlockedMap> blockedTargetMap = MueLu::UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GeneratedBlockedTargetMap(*blockedRowMap, *rowMapImporter);

    // NOTE: This code qualifies as "correct but not particularly performant"  If this needs to be sped up, we can probably read data from the existing importer to
    // build sub-importers rather than generating new ones ex nihilo
    size_t numBlocks = blockedRowMap->getNumMaps();
    std::vector<RCP<const Import> > subImports(numBlocks);

    for (size_t i = 0; i < numBlocks; i++) {
      RCP<const Map> source = blockedRowMap->getMap(i);
      RCP<const Map> target = blockedTargetMap->getMap(i);
      subImports[i]         = ImportFactory::Build(source, target);
    }
    Set(currentLevel, "SubImporters", subImports);
  }

  Set(currentLevel, "Importer", rowMapImporter);

  // Importer saving
  bool save_importer = pL.get<bool>("repartition: save importer");
  if (save_importer) {
    currentLevel.Set("Importer", rowMapImporter, NoFactory::get());
    currentLevel.AddKeepFlag("Importer", NoFactory::get(), MueLu::Final);
    currentLevel.RemoveKeepFlag("Importer", NoFactory::get(), MueLu::UserData);  // FIXME: This is a hack
  }
  // ======================================================================================================
  // Print some data
  // ======================================================================================================
  if (!rowMapImporter.is_null() && IsPrint(Statistics2)) {
    // int oldRank = SetProcRankVerbose(rebalancedAc->getRowMap()->getComm()->getRank());
    GetOStream(Statistics2) << PerfUtils::PrintImporterInfo(rowMapImporter, "Importer for rebalancing");
    // SetProcRankVerbose(oldRank);
  }
  if (pL.get<bool>("repartition: print partition distribution") && IsPrint(Statistics2)) {
    // Print the grid of processors
    GetOStream(Statistics2) << "Partition distribution over cores (ownership is indicated by '+')" << std::endl;

    char amActive = (myGIDs.size() ? 1 : 0);
    std::vector<char> areActive(numProcs, 0);
    MPI_Gather(&amActive, 1, MPI_CHAR, &areActive[0], 1, MPI_CHAR, 0, *rawMpiComm);

    int rowWidth = std::min(Teuchos::as<int>(ceil(sqrt(numProcs))), 100);
    for (int proc = 0; proc < numProcs; proc += rowWidth) {
      for (int j = 0; j < rowWidth; j++)
        if (proc + j < numProcs)
          GetOStream(Statistics2) << (areActive[proc + j] ? "+" : ".");
        else
          GetOStream(Statistics2) << " ";

      GetOStream(Statistics2) << "      " << proc << ":" << std::min(proc + rowWidth, numProcs) - 1 << std::endl;
    }
  }

}  // Build

//----------------------------------------------------------------------
template <typename T, typename W>
struct Triplet {
  T i, j;
  W v;
};
template <typename T, typename W>
static bool compareTriplets(const Triplet<T, W>& a, const Triplet<T, W>& b) {
  return (a.v > b.v);  // descending order
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DeterminePartitionPlacement(const Matrix& A, GOVector& decomposition, GO numPartitions, bool willAcceptPartition, bool allSubdomainsAcceptPartitions) const {
  RCP<const Map> rowMap = A.getRowMap();

  RCP<const Teuchos::Comm<int> > comm = rowMap->getComm()->duplicate();
  int numProcs                        = comm->getSize();

  RCP<const Teuchos::MpiComm<int> > tmpic = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
  TEUCHOS_TEST_FOR_EXCEPTION(tmpic == Teuchos::null, Exceptions::RuntimeError, "Cannot cast base Teuchos::Comm to Teuchos::MpiComm object.");
  RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm = tmpic->getRawMpiComm();

  const Teuchos::ParameterList& pL = GetParameterList();

  // maxLocal is a constant which determins the number of largest edges which are being exchanged
  // The idea is that we do not want to construct the full bipartite graph, but simply a subset of
  // it, which requires less communication. By selecting largest local edges we hope to achieve
  // similar results but at a lower cost.
  const int maxLocal = pL.get<int>("repartition: remap num values");
  const int dataSize = 2 * maxLocal;

  ArrayRCP<GO> decompEntries;
  if (decomposition.getLocalLength() > 0)
    decompEntries = decomposition.getDataNonConst(0);

  // Step 1: Sort local edges by weight
  // Each edge of a bipartite graph corresponds to a triplet (i, j, v) where
  //   i: processor id that has some piece of part with part_id = j
  //   j: part id
  //   v: weight of the edge
  // We set edge weights to be the total number of nonzeros in rows on this processor which
  // correspond to this part_id. The idea is that when we redistribute matrix, this weight
  // is a good approximation of the amount of data to move.
  // We use two maps, original which maps a partition id of an edge to the corresponding weight,
  // and a reverse one, which is necessary to sort by edges.
  std::map<GO, GO> lEdges;
  if (willAcceptPartition)
    for (LO i = 0; i < decompEntries.size(); i++)
      lEdges[decompEntries[i]] += A.getNumEntriesInLocalRow(i);

  // Reverse map, so that edges are sorted by weight.
  // This results in multimap, as we may have edges with the same weight
  std::multimap<GO, GO> revlEdges;
  for (typename std::map<GO, GO>::const_iterator it = lEdges.begin(); it != lEdges.end(); it++)
    revlEdges.insert(std::make_pair(it->second, it->first));

  // Both lData and gData are arrays of data which we communicate. The data is stored
  // in pairs, so that data[2*i+0] is the part index, and data[2*i+1] is the corresponding edge weight.
  // We do not store processor id in data, as we can compute that by looking on the offset in the gData.
  Array<GO> lData(dataSize, -1), gData(numProcs * dataSize);
  int numEdges = 0;
  for (typename std::multimap<GO, GO>::reverse_iterator rit = revlEdges.rbegin(); rit != revlEdges.rend() && numEdges < maxLocal; rit++) {
    lData[2 * numEdges + 0] = rit->second;  // part id
    lData[2 * numEdges + 1] = rit->first;   // edge weight
    numEdges++;
  }

  // Step 2: Gather most edges
  // Each processors contributes maxLocal edges by providing maxLocal pairs <part id, weight>, which is of size dataSize
  MPI_Datatype MpiType = Teuchos::Details::MpiTypeTraits<GO>::getType();
  MPI_Allgather(static_cast<void*>(lData.getRawPtr()), dataSize, MpiType, static_cast<void*>(gData.getRawPtr()), dataSize, MpiType, *rawMpiComm);

  // Step 3: Construct mapping

  // Construct the set of triplets
  Teuchos::Array<Triplet<int, int> > gEdges(numProcs * maxLocal);
  Teuchos::Array<bool> procWillAcceptPartition(numProcs, allSubdomainsAcceptPartitions);
  size_t k = 0;
  for (LO i = 0; i < gData.size(); i += 2) {
    int procNo = i / dataSize;  // determine the processor by its offset (since every processor sends the same amount)
    GO part    = gData[i + 0];
    GO weight  = gData[i + 1];
    if (part != -1) {  // skip nonexistent edges
      gEdges[k].i                     = procNo;
      gEdges[k].j                     = part;
      gEdges[k].v                     = weight;
      procWillAcceptPartition[procNo] = true;
      k++;
    }
  }
  gEdges.resize(k);

  // Sort edges by weight
  // NOTE: compareTriplets is actually a reverse sort, so the edges weight is in decreasing order
  std::sort(gEdges.begin(), gEdges.end(), compareTriplets<int, int>);

  // Do matching
  std::map<int, int> match;
  Teuchos::Array<char> matchedRanks(numProcs, 0), matchedParts(numPartitions, 0);
  int numMatched = 0;
  for (typename Teuchos::Array<Triplet<int, int> >::const_iterator it = gEdges.begin(); it != gEdges.end(); it++) {
    GO rank = it->i;
    GO part = it->j;
    if (matchedRanks[rank] == 0 && matchedParts[part] == 0) {
      matchedRanks[rank] = 1;
      matchedParts[part] = 1;
      match[part]        = rank;
      numMatched++;
    }
  }
  GetOStream(Statistics1) << "Number of unassigned partitions before cleanup stage: " << (numPartitions - numMatched) << " / " << numPartitions << std::endl;

  // Step 4: Assign unassigned partitions if necessary.
  // We do that through desperate matching for remaining partitions:
  // We select the lowest rank that can still take a partition.
  // The reason it is done this way is that we don't need any extra communication, as we don't
  // need to know which parts are valid.
  if (numPartitions - numMatched > 0) {
    Teuchos::Array<char> partitionCounts(numPartitions, 0);
    for (typename std::map<int, int>::const_iterator it = match.begin(); it != match.end(); it++)
      partitionCounts[it->first] += 1;
    for (int part = 0, matcher = 0; part < numPartitions; part++) {
      if (partitionCounts[part] == 0) {
        // Find first non-matched rank that accepts partitions
        while (matchedRanks[matcher] || !procWillAcceptPartition[matcher])
          matcher++;

        match[part] = matcher++;
        numMatched++;
      }
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION(numMatched != numPartitions, Exceptions::RuntimeError, "MueLu::RepartitionFactory::DeterminePartitionPlacement: Only " << numMatched << " partitions out of " << numPartitions << " got assigned to ranks.");

  // Step 5: Permute entries in the decomposition vector
  for (LO i = 0; i < decompEntries.size(); i++)
    decompEntries[i] = match[decompEntries[i]];
}

}  // namespace MueLu

#endif  // ifdef HAVE_MPI

#endif  // MUELU_REPARTITIONFACTORY_DEF_HPP
