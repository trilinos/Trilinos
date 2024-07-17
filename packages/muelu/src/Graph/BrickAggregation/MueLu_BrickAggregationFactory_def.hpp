// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_BRICKAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_BRICKAGGREGATIONFACTORY_DEF_HPP_

#include "MueLu_BrickAggregationFactory_decl.hpp"
#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#endif
#include <Teuchos_OrdinalTraits.hpp>

#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_Aggregates.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_LWGraph.hpp"

#include "MueLu_LWGraph.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("aggregation: brick x size");
  SET_VALID_ENTRY("aggregation: brick y size");
  SET_VALID_ENTRY("aggregation: brick z size");
  SET_VALID_ENTRY("aggregation: brick x Dirichlet");
  SET_VALID_ENTRY("aggregation: brick y Dirichlet");
  SET_VALID_ENTRY("aggregation: brick z Dirichlet");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory for matrix");
  validParamList->set<RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Generating factory for coordinates");
  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
  Input(currentLevel, "Coordinates");
}

// The current implementation cannot deal with bricks larger than 3x3(x3) in
// parallel. The reason is that aggregation infrastructure in place has
// major drawbacks.
//
// Aggregates class is constructed with a help of a provided map, either
// taken from a graph, or provided directly. This map is usually taken to be
// a column map of a matrix. The reason for that is that if we have an
// overlapped aggregation, we want the processor owning aggregates to store
// agg id for all nodes in this aggregate. If we used row map, there would
// be no way for the processor to know whether there are some other nodes on
// a different processor which belong to its aggregate. On the other hand,
// using column map allows both vertex2AggId and procWinner arrays in
// Aggregates class to store some extra data, such as whether nodes belonging
// to a different processor belong to this processor aggregate.
//
// The drawback of this is that it stores only overlap=1 data. For aggressive
// coarsening, such a brick aggregation with a large single dimension of
// brick, it could happen that we need to know depth two or more extra nodes
// in the other processor subdomain.
//
// Another issue is that we may have some implicit connection between
// aggregate map and maps of A used in the construction of a tentative
// prolongator.
//
// Another issue is that it seems that some info is unused or not required.
// Specifically, it seems that if a node belongs to an aggregate on a
// different processor, we don't actually need to set vertex2AggId and
// procWinner, despite the following comment in
// Aggregates decl:
//      vertex2AggId[k] gives a local id
//      corresponding to the aggregate to which
//      local id k has been assigned.  While k
//      is the local id on my processor (MyPID)
//      vertex2AggId[k] is the local id on the
//      processor which actually owns the
//      aggregate. This owning processor has id
//      given by procWinner[k].
// It is possible that that info is only used during arbitration in
// CoupledAggregationFactory.
//
// The steps that we need to do to resolve this issue:
//   - Break the link between maps in TentativePFactory, allowing any maps in Aggregates
//   - Allow Aggregates to construct their own maps, if necessary, OR
//   - construct aggregates based on row map
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  typedef Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> MultiVector_d;

  const ParameterList& pL   = GetParameterList();
  RCP<MultiVector_d> coords = Get<RCP<MultiVector_d> >(currentLevel, "Coordinates");
  RCP<Matrix> A             = Get<RCP<Matrix> >(currentLevel, "A");
  RCP<const Map> rowMap     = A->getRowMap();
  RCP<const Map> colMap     = A->getColMap();
  GO GO_INVALID             = Teuchos::OrdinalTraits<GO>::invalid();

  RCP<const Teuchos::Comm<int> > comm = rowMap->getComm();
  int numProcs                        = comm->getSize();
  int myRank                          = comm->getRank();

  int numPoints = colMap->getLocalNumElements();

  bx_ = pL.get<int>("aggregation: brick x size");
  by_ = pL.get<int>("aggregation: brick y size");
  bz_ = pL.get<int>("aggregation: brick z size");

  dirichletX_ = pL.get<bool>("aggregation: brick x Dirichlet");
  dirichletY_ = pL.get<bool>("aggregation: brick y Dirichlet");
  dirichletZ_ = pL.get<bool>("aggregation: brick z Dirichlet");
  if (dirichletX_) GetOStream(Runtime0) << "Dirichlet boundaries in the x direction" << std::endl;
  if (dirichletY_) GetOStream(Runtime0) << "Dirichlet boundaries in the y direction" << std::endl;
  if (dirichletZ_) GetOStream(Runtime0) << "Dirichlet boundaries in the z direction" << std::endl;

  if (numProcs > 1) {
    // TODO: deal with block size > 1  (see comments above)
    // TEUCHOS_TEST_FOR_EXCEPTION(bx_ > 3 || by_ > 3 || bz_ > 3, Exceptions::RuntimeError, "Currently cannot deal with brick size > 3");
  }

  RCP<MultiVector_d> overlappedCoords = coords;
  RCP<const Import> importer          = ImportFactory::Build(coords->getMap(), colMap);
  if (!importer.is_null()) {
    overlappedCoords = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(colMap, coords->getNumVectors());
    overlappedCoords->doImport(*coords, *importer, Xpetra::INSERT);
  }

  // Setup misc structures
  // Logically, we construct enough data to query topological information of a rectangular grid
  Setup(comm, overlappedCoords, colMap);

  GetOStream(Runtime0) << "Using brick size: " << bx_
                       << (nDim_ > 1 ? "x " + toString(by_) : "")
                       << (nDim_ > 2 ? "x " + toString(bz_) : "") << std::endl;

  // Build the graph
  BuildGraph(currentLevel, A);

  // Construct aggregates
  RCP<Aggregates> aggregates = rcp(new Aggregates(colMap));
  aggregates->setObjectLabel("Brick");

  ArrayRCP<LO> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
  ArrayRCP<LO> procWinner   = aggregates->GetProcWinner()->getDataNonConst(0);

  // In the first pass, we set a mapping from a vertex to aggregate global id. We deal with a structured
  // rectangular mesh, therefore we know the structure of aggregates. For each vertex we can tell exactly
  // which aggregate it belongs to.
  // If we determine that the aggregate does not belong to us (i.e. the root vertex does not belong to this
  // processor, or is outside and we lost "" arbitration), we record the global aggregate id in order to
  // fetch the local info from the processor owning the aggregate. This is required for aggregates, as it
  // uses the local aggregate ids of the owning processor.
  std::set<GO> myAggGIDs, remoteAggGIDs;
  for (LO LID = 0; LID < numPoints; LID++) {
    GO aggGID = getAggGID(LID);
    //      printf("[%d] (%d,%d,%d) => agg %d\n",LID,(int)(*xMap_)[x_[LID]],nDim_ > 1 ? (int)(*yMap_)[y_[LID]] : -1,nDim_ > 2 ? (int)(*zMap_)[z_[LID]] : -1,(int)aggGID);
    if (aggGID == GO_INVALID) continue;
    //      printf("[%d] getRoot = %d\n",(int)LID,(int)getRoot(LID));

    if ((revMap_.find(getRoot(LID)) != revMap_.end()) && rowMap->isNodeGlobalElement(colMap->getGlobalElement(revMap_[getRoot(LID)]))) {
      // Root of the brick aggregate containing GID (<- LID) belongs to us
      vertex2AggId[LID] = aggGID;
      myAggGIDs.insert(aggGID);

      if (isRoot(LID))
        aggregates->SetIsRoot(LID);
      //	printf("[%d] initial vertex2AggId = %d\n",(int)LID,(int)vertex2AggId[LID]);
    } else {
      remoteAggGIDs.insert(aggGID);
    }
  }
  size_t numAggregates = myAggGIDs.size();
  size_t numRemote     = remoteAggGIDs.size();
  aggregates->SetNumAggregates(numAggregates);

  std::map<GO, LO> AggG2L;   // Map: Agg GID -> Agg LID (possibly on a different processor)
  std::map<GO, int> AggG2R;  // Map: Agg GID -> processor rank owning aggregate

  Array<GO> myAggGIDsArray(numAggregates), remoteAggGIDsArray(numRemote);

  // Fill in the maps for aggregates that we own
  size_t ind = 0;
  for (typename std::set<GO>::const_iterator it = myAggGIDs.begin(); it != myAggGIDs.end(); it++) {
    AggG2L[*it] = ind;
    AggG2R[*it] = myRank;

    myAggGIDsArray[ind++] = *it;
  }

  // The map is a convenient way to fetch remote local indices from global indices.
  RCP<Map> aggMap = MapFactory::Build(rowMap->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                      myAggGIDsArray, 0, comm);

  ind = 0;
  for (typename std::set<GO>::const_iterator it = remoteAggGIDs.begin(); it != remoteAggGIDs.end(); it++)
    remoteAggGIDsArray[ind++] = *it;

  // Fetch the required aggregate local ids and ranks
  Array<int> remoteProcIDs(numRemote);
  Array<LO> remoteLIDs(numRemote);
  aggMap->getRemoteIndexList(remoteAggGIDsArray, remoteProcIDs, remoteLIDs);

  // Fill in the maps for aggregates that we don't own but which have some of our vertices
  for (size_t i = 0; i < numRemote; i++) {
    AggG2L[remoteAggGIDsArray[i]] = remoteLIDs[i];
    AggG2R[remoteAggGIDsArray[i]] = remoteProcIDs[i];
  }

  // Remap aggregate GIDs to LIDs and set up owning processors
  for (LO LID = 0; LID < numPoints; LID++) {
    if (revMap_.find(getRoot(LID)) != revMap_.end() && rowMap->isNodeGlobalElement(colMap->getGlobalElement(revMap_[getRoot(LID)]))) {
      GO aggGID = vertex2AggId[LID];
      if (aggGID != MUELU_UNAGGREGATED) {
        vertex2AggId[LID] = AggG2L[aggGID];
        procWinner[LID]   = AggG2R[aggGID];
      }
    }
  }

  GO numGlobalRemote;
  MueLu_sumAll(comm, as<GO>(numRemote), numGlobalRemote);
  aggregates->AggregatesCrossProcessors(numGlobalRemote);

  Set(currentLevel, "Aggregates", aggregates);

  GetOStream(Statistics1) << aggregates->description() << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Setup(const RCP<const Teuchos::Comm<int> >& comm, const RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> >& coords, const RCP<const Map>& /* map */) const {
  nDim_ = coords->getNumVectors();

  x_    = coords->getData(0);
  xMap_ = Construct1DMap(comm, x_);
  nx_   = xMap_->size();

  ny_ = 1;
  if (nDim_ > 1) {
    y_    = coords->getData(1);
    yMap_ = Construct1DMap(comm, y_);
    ny_   = yMap_->size();
  }

  nz_ = 1;
  if (nDim_ > 2) {
    z_    = coords->getData(2);
    zMap_ = Construct1DMap(comm, z_);
    nz_   = zMap_->size();
  }

  for (size_t ind = 0; ind < coords->getLocalLength(); ind++) {
    GO i = (*xMap_)[(coords->getData(0))[ind]], j = 0, k = 0;
    if (nDim_ > 1)
      j = (*yMap_)[(coords->getData(1))[ind]];
    if (nDim_ > 2)
      k = (*zMap_)[(coords->getData(2))[ind]];

    revMap_[k * ny_ * nx_ + j * nx_ + i] = ind;
  }

  // Get the number of aggregates in each direction, correcting for Dirichlet
  int xboost = dirichletX_ ? 1 : 0;
  int yboost = dirichletY_ ? 1 : 0;
  int zboost = dirichletZ_ ? 1 : 0;
  naggx_     = (nx_ - 2 * xboost) / bx_ + ((nx_ - 2 * xboost) % bx_ ? 1 : 0);

  if (nDim_ > 1)
    naggy_ = (ny_ - 2 * yboost) / by_ + ((ny_ - 2 * yboost) % by_ ? 1 : 0);
  else
    naggy_ = 1;

  if (nDim_ > 2)
    naggz_ = (nz_ - 2 * zboost) / bz_ + ((nz_ - 2 * zboost) % bz_ ? 1 : 0);
  else
    naggz_ = 1;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<typename BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::container>
BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Construct1DMap(const RCP<const Teuchos::Comm<int> >& comm,
                   const ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& x) const {
  int n = x.size();

  // Step 1: Create a local vector with unique coordinate points
  RCP<container> gMap = rcp(new container);
  for (int i = 0; i < n; i++)
    (*gMap)[x[i]] = 0;

#ifdef HAVE_MPI
  // Step 2: exchange coordinates
  // NOTE: we assume the coordinates are double, or double compatible
  // That means that for complex case, we assume that all imaginary parts are zeros
  int numProcs = comm->getSize();
  if (numProcs > 1) {
    RCP<const Teuchos::MpiComm<int> > dupMpiComm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm->duplicate());

    MPI_Comm rawComm = (*dupMpiComm->getRawMpiComm())();

    int sendCnt = gMap->size(), cnt = 0, recvSize;
    Array<int> recvCnt(numProcs), Displs(numProcs);
    Array<double> sendBuf, recvBuf;

    sendBuf.resize(sendCnt);
    for (typename container::const_iterator cit = gMap->begin(); cit != gMap->end(); cit++)
      sendBuf[cnt++] = Teuchos::as<double>(STS::real(cit->first));

    MPI_Allgather(&sendCnt, 1, MPI_INT, recvCnt.getRawPtr(), 1, MPI_INT, rawComm);
    Displs[0] = 0;
    for (int i = 0; i < numProcs - 1; i++)
      Displs[i + 1] = Displs[i] + recvCnt[i];
    recvSize = Displs[numProcs - 1] + recvCnt[numProcs - 1];
    recvBuf.resize(recvSize);
    MPI_Allgatherv(sendBuf.getRawPtr(), sendCnt, MPI_DOUBLE, recvBuf.getRawPtr(), recvCnt.getRawPtr(), Displs.getRawPtr(), MPI_DOUBLE, rawComm);

    for (int i = 0; i < recvSize; i++)
      (*gMap)[as<SC>(recvBuf[i])] = 0;
  }
#endif

  GO cnt = 0;
  for (typename container::iterator it = gMap->begin(); it != gMap->end(); it++)
    it->second = cnt++;

  return gMap;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isRoot(LocalOrdinal LID) const {
  int i, j, k;
  getIJK(LID, i, j, k);

  return (k * ny_ * nx_ + j * nx_ + i) == getRoot(LID);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isDirichlet(LocalOrdinal LID) const {
  bool boundary = false;
  int i, j, k;
  getIJK(LID, i, j, k);
  if (dirichletX_ && (i == 0 || i == nx_ - 1))
    boundary = true;
  if (nDim_ > 1 && dirichletY_ && (j == 0 || j == ny_ - 1))
    boundary = true;
  if (nDim_ > 2 && dirichletZ_ && (k == 0 || k == nz_ - 1))
    boundary = true;

  return boundary;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRoot(LocalOrdinal LID) const {
  if (isDirichlet(LID))
    return Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();

  int aggI, aggJ, aggK;
  getAggIJK(LID, aggI, aggJ, aggK);
  int xboost = dirichletX_ ? 1 : 0;
  int yboost = dirichletY_ ? 1 : 0;
  int zboost = dirichletZ_ ? 1 : 0;

  int i = xboost + aggI * bx_ + (bx_ - 1) / 2;
  int j = (nDim_ > 1) ? yboost + aggJ * by_ + (by_ - 1) / 2 : 0;
  int k = (nDim_ > 2) ? zboost + aggK * bz_ + (bz_ - 1) / 2 : 0;

  return k * ny_ * nx_ + j * nx_ + i;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getIJK(LocalOrdinal LID, int& i, int& j, int& k) const {
  i = (*xMap_)[x_[LID]];
  j = (nDim_ > 1) ? (*yMap_)[y_[LID]] : 0;
  k = (nDim_ > 2) ? (*zMap_)[z_[LID]] : 0;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getAggIJK(LocalOrdinal LID, int& i, int& j, int& k) const {
  int xboost = dirichletX_ ? 1 : 0;
  int yboost = dirichletY_ ? 1 : 0;
  int zboost = dirichletZ_ ? 1 : 0;
  int pointI, pointJ, pointK;
  getIJK(LID, pointI, pointJ, pointK);
  i = (pointI - xboost) / bx_;

  if (nDim_ > 1)
    j = (pointJ - yboost) / by_;
  else
    j = 0;

  if (nDim_ > 2)
    k = (pointK - zboost) / bz_;
  else
    k = 0;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getAggGID(LocalOrdinal LID) const {
  bool boundary = false;

  int i, j, k;
  getIJK(LID, i, j, k);
  int ii, jj, kk;
  getAggIJK(LID, ii, jj, kk);

  if (dirichletX_ && (i == 0 || i == nx_ - 1)) boundary = true;
  if (nDim_ > 1 && dirichletY_ && (j == 0 || j == ny_ - 1)) boundary = true;
  if (nDim_ > 2 && dirichletZ_ && (k == 0 || k == nz_ - 1)) boundary = true;

  /*
  if(boundary)
    printf("[%d] coord = (%d,%d,%d) {%d,%d,%d} agg = (%d,%d,%d) {%d,%d,%d} => agg %s\n",LID,i,j,k,nx_,ny_,nz_,ii,jj,kk,naggx_,naggy_,naggz_,"BOUNDARY");
  else
    printf("[%d] coord = (%d,%d,%d) {%d,%d,%d} agg = (%d,%d,%d) {%d,%d,%d} => agg %d\n",LID,i,j,k,nx_,ny_,nz_,ii,jj,kk,naggx_,naggy_,naggz_,kk*naggy_*naggx_ + jj*naggx_ + ii);
  */

  if (boundary)
    return Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();
  else
    return Teuchos::as<GlobalOrdinal>(kk * naggy_ * naggx_) + Teuchos::as<GlobalOrdinal>(jj * naggx_) + ii;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildGraph(Level& currentLevel, const RCP<Matrix>& A) const {
  // TODO: Currently only works w/ 1 DOF per node
  double dirichletThreshold = 0.0;

  if (bx_ > 1 && (nDim_ <= 1 || by_ > 1) && (nDim_ <= 2 || bz_ > 1)) {
    FactoryMonitor m(*this, "Generating Graph (trivial)", currentLevel);
    /*** Case 1: Use the matrix is the graph ***/
    // Bricks are of non-trivial size in all active dimensions
    RCP<LWGraph> graph = rcp(new LWGraph(A->getCrsGraph(), "graph of A"));
    auto boundaryNodes = MueLu::Utilities<SC, LO, GO, NO>::DetectDirichletRows_kokkos_host(*A, dirichletThreshold);
    graph->SetBoundaryNodeMap(boundaryNodes);

    if (GetVerbLevel() & Statistics1) {
      GO numLocalBoundaryNodes  = 0;
      GO numGlobalBoundaryNodes = 0;
      for (size_t i = 0; i < boundaryNodes.size(); ++i)
        if (boundaryNodes(i))
          numLocalBoundaryNodes++;
      RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
      MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
      GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
    }
    Set(currentLevel, "DofsPerNode", 1);
    Set(currentLevel, "Graph", graph);
    Set(currentLevel, "Filtering", false);
  } else {
    FactoryMonitor m(*this, "Generating Graph", currentLevel);
    /*** Case 2: Dropping required ***/
    // There is at least one active dimension in which we are not coarsening.
    // Those connections need to be dropped
    bool drop_x = (bx_ == 1);
    bool drop_y = (nDim_ > 1 && by_ == 1);
    bool drop_z = (nDim_ > 2 && bz_ == 1);

    typename LWGraph::row_type::non_const_type rows("rows", A->getLocalNumRows() + 1);
    typename LWGraph::entries_type::non_const_type columns("columns", A->getLocalNumEntries());

    size_t N = A->getRowMap()->getLocalNumElements();

    // FIXME: Do this on the host because indexing functions are host functions
    auto G      = A->getLocalMatrixHost().graph;
    auto rowptr = G.row_map;
    auto colind = G.entries;

    int ct  = 0;
    rows(0) = 0;
    for (size_t row = 0; row < N; row++) {
      // NOTE: Assumes that the first part of the colmap is the rowmap
      int ir, jr, kr;
      LO row2 = A->getColMap()->getLocalElement(A->getRowMap()->getGlobalElement(row));
      getIJK(row2, ir, jr, kr);

      for (size_t cidx = rowptr[row]; cidx < rowptr[row + 1]; cidx++) {
        int ic, jc, kc;
        LO col = colind[cidx];
        getIJK(col, ic, jc, kc);

        if ((row2 != col) && ((drop_x && ir != ic) || (drop_y && jr != jc) || (drop_z && kr != kc))) {
          // Drop it
          //            printf("[%4d] DROP row = (%d,%d,%d) col = (%d,%d,%d)\n",(int)row,ir,jr,kr,ic,jc,kc);
        } else {
          // Keep it
          //            printf("[%4d] KEEP row = (%d,%d,%d) col = (%d,%d,%d)\n",(int)row,ir,jr,kr,ic,jc,kc);
          columns(ct) = col;
          ct++;
        }
      }
      rows(row + 1) = ct;
    }  // end for

    RCP<LWGraph> graph = rcp(new LWGraph(rows, columns, A->getRowMap(), A->getColMap(), "thresholded graph of A"));

    auto boundaryNodes = MueLu::Utilities<SC, LO, GO, NO>::DetectDirichletRows_kokkos_host(*A, dirichletThreshold);
    graph->SetBoundaryNodeMap(boundaryNodes);

    if (GetVerbLevel() & Statistics1) {
      GO numLocalBoundaryNodes  = 0;
      GO numGlobalBoundaryNodes = 0;
      for (size_t i = 0; i < boundaryNodes.size(); ++i)
        if (boundaryNodes(i))
          numLocalBoundaryNodes++;
      RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
      MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
      GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
    }
    Set(currentLevel, "DofsPerNode", 1);
    Set(currentLevel, "Graph", graph);
    Set(currentLevel, "Filtering", true);
  }  // end else

}  // end BuildGraph

}  // namespace MueLu

#endif /* MUELU_BRICKAGGREGATIONFACTORY_DEF_HPP_ */
