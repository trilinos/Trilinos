// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Create a Tpetra::MultiVector, and time the following:
//   1.  Build a multivector with contiguous global ids.
//   2.  Migrate.
//   3.  Divide procs into two groups and build new a new multivector
//       in each group with the non-contiguous global ids of the migrated data.
//   4.  Migrate the new multivector in the subgroup.
//
// Repeat this using Epetra, and also Zoltan1.
//
// This mimics the recursive bisection algorithm in Zoltan2.
//
// Because global ordinals are "int" in the this test, Zoltan1 must be
// configured with cmake option ZOLTAN_ENABLE_UINT_IDS.

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Import.hpp>

#include <Epetra_Vector.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Import.h>

#include <zoltan.h>

#include <Zoltan2_Util.hpp>

#include <string>
#include <sstream>
#include <iostream>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::TimeMonitor;
using Teuchos::Time;

#ifdef HAVE_MPI

typedef int TPETRA_LO;
typedef int TPETRA_GO;
typedef unsigned int ZOLTAN_LO;
typedef unsigned int ZOLTAN_GO;
#define MPI_ZOLTAN_GO_TYPE MPI_UNSIGNED
typedef double TPETRA_SCALAR;

RCP<Time> tmvBuild;
RCP<Time> tmvMigrate;
RCP<Time> tmvBuildN;
RCP<Time> tmvMigrateN;

RCP<Time> emvBuild;
RCP<Time> emvMigrate;
RCP<Time> emvBuildN;
RCP<Time> emvMigrateN;

RCP<Time> ztnBuild;
RCP<Time> ztnMigrate;
RCP<Time> ztnBuildN;
RCP<Time> ztnMigrateN;

using Teuchos::MpiComm;
using Teuchos::Comm;

enum testConstants {COORDDIM=3};

static void usage(char *arg[]){
  std::cout << "Usage:" << std::endl;
  std::cout << arg[0] << " {num coords}" << std::endl;
}

///////////////////////////////////////////////
// Generate global Ids for different mappings
///////////////////////////////////////////////
TPETRA_LO numSequentialGlobalIds(TPETRA_GO numGlobalCoords, int nprocs, int rank)
{
  TPETRA_LO share = numGlobalCoords / nprocs;
  TPETRA_LO extra = numGlobalCoords % nprocs;
  TPETRA_LO numLocalCoords = share;
  if (rank < extra)
    numLocalCoords++;

  return numLocalCoords;
}
template <typename T>
void roundRobinGlobalIds(T numGlobalCoords, int nprocs, int rank,
  T *&gids)
{
  // Local number of coordinates does not change.
  T share = numGlobalCoords / nprocs;
  T extra = numGlobalCoords % nprocs;
  T numLocalCoords = share;
  if (static_cast<T> (rank) < extra)
    numLocalCoords++;

  gids = new T [numLocalCoords];
  if (!gids)
    throw std::bad_alloc();

  T next = 0;
  for (T i=rank; i < numGlobalCoords; i+=nprocs)
    gids[next++] = i;

  return;
}
template <typename T>
void subGroupGloballyIncreasingIds(T numGlobalCoords,
  int nprocs, int rank, T *&gids)
{
  int numProcsLeftHalf = nprocs / 2;
  T share = numGlobalCoords / nprocs;
  T extra = numGlobalCoords % nprocs;
  T numCoordsLeftHalf = 0;
  T firstIdx = 0, endIdx = 0;

  T endP = ((numProcsLeftHalf > rank) ?  numProcsLeftHalf : rank);

  for (T p=0; p < endP ; p++){
    T numLocalCoords = share;
    if (p < extra)
      numLocalCoords++;

    if (p < static_cast<T> (rank))
      firstIdx += numLocalCoords;

    if (p < static_cast<T> (numProcsLeftHalf))
      numCoordsLeftHalf += numLocalCoords;
  }

  endIdx = firstIdx + share;
  if (rank < int(extra))
    endIdx++;

  if (rank >= numProcsLeftHalf){
    firstIdx -= numCoordsLeftHalf;
    endIdx -= numCoordsLeftHalf;
  }

  int firstProc=0, endProc=0;

  if (rank < numProcsLeftHalf){
    firstProc = 0;
    endProc = numProcsLeftHalf;
  }
  else{
    firstProc = numProcsLeftHalf;
    endProc = nprocs;
  }

  int numProcsInMyHalf = endProc - firstProc;

  // Picture the round robin global ids as a matrix where the
  // columns are the processes and row values are global ids.
  // Row values start at 0 in the upper left corner and increase
  // to nprocs-1 in the upper right corner.  The values in
  // the second row are nprocs greater than the first row,
  // and so on.
  //
  // The processes were divided into two halves, represented
  // by a vertical line through the matrix dividing the
  // processes in the left half from the processes in the
  // right half.
  //
  // Now we want to enumerate the global ids in my half
  // in increasing order.

  T numLocalCoords = endIdx - firstIdx;
  gids = new T [numLocalCoords];
  if (!gids)
    throw std::bad_alloc();

  T firstRow = firstIdx / numProcsInMyHalf;
  T firstCol = (firstIdx % numProcsInMyHalf) + firstProc;
  int next = 0;

  for (T row = firstRow; static_cast<T> (next) < numLocalCoords; row++){
    T firstGid = row * nprocs + firstCol;
    for (T col = firstCol; static_cast<T> (next) < numLocalCoords && col < static_cast<T> (endProc); col++){
      gids[next++] = firstGid++;
    }
    firstCol = firstProc;
  }
}

void timeEpetra(TPETRA_GO numGlobalCoords, const RCP<const MpiComm<int> > &comm, bool);
void
timeTpetra (const TPETRA_GO numGlobalCoords,
            const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
            const bool doMemory);
void timeZoltan(ZOLTAN_GO numGlobalCoords, bool);

///////////////////////////////////////////////
// Main
///////////////////////////////////////////////

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > genComm = Tpetra::getDefaultComm();
  Teuchos::RCP<const MpiComm<int> > comm =
    rcp_dynamic_cast<const MpiComm<int> >(genComm);

  int rank = genComm->getRank();

  if (narg < 2){
    if (rank == 0)
      usage(arg);
    return 1;
  }

  TPETRA_GO numGlobalCoords = 0;
  std::string theArg(arg[1]);
  std::istringstream iss(theArg);
  iss >> numGlobalCoords;
  if (numGlobalCoords < genComm->getSize()){
    if (rank == 0)
      usage(arg);
    return 1;
  }

  if (rank == 0)
    std::cout << numGlobalCoords << " coordinates" << std::endl;

  tmvBuild = TimeMonitor::getNewTimer("CONSEC build Tpetra");
  tmvMigrate = TimeMonitor::getNewTimer("CONSEC migrate Tpetra");
  tmvBuildN = TimeMonitor::getNewTimer("!CONSEC build Tpetra");
  tmvMigrateN = TimeMonitor::getNewTimer("!CONSEC migrate Tpetra");

  ztnBuild = TimeMonitor::getNewTimer("CONSEC build Zoltan1");
  ztnMigrate = TimeMonitor::getNewTimer("CONSEC migrate Zoltan1");
  ztnBuildN = TimeMonitor::getNewTimer("!CONSEC build Zoltan1");
  ztnMigrateN = TimeMonitor::getNewTimer("!CONSEC migrate Zoltan1");

  emvBuild = TimeMonitor::getNewTimer("CONSEC build Epetra");
  emvMigrate = TimeMonitor::getNewTimer("CONSEC migrate Epetra");
  emvBuildN = TimeMonitor::getNewTimer("!CONSEC build Epetra");
  emvMigrateN = TimeMonitor::getNewTimer("!CONSEC migrate Epetra");

  TimeMonitor::zeroOutTimers();

  int ntests = 3;

  // Test with Zoltan_Comm and Zoltan_DataDirectory

  for (int i=0; i < ntests; i++){
    if (rank == 0)
      std::cout << "Zoltan test " << i+1 << std::endl;
    timeZoltan(numGlobalCoords, i==ntests-1);
  }

  // Test with Epetra_MultiVector

  for (int i=0; i < ntests; i++){
    if (rank == 0)
      std::cout << "Epetra test " << i+1 << std::endl;
    timeEpetra(numGlobalCoords, comm, i==ntests-1);
  }

  // Test with Tpetra::MultiVector

  for (int i=0; i < ntests; i++){
    if (rank == 0)
      std::cout << "Tpetra test " << i+1 << std::endl;
    timeTpetra(numGlobalCoords, comm, i==ntests-1);
  }

  // Report

  TimeMonitor::summarize();

  if (comm->getRank() == 0)
    std::cout << "PASS" << std::endl;
}

void
timeTpetra (const TPETRA_GO numGlobalCoords,
            const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
            const bool doMemory)
{
  using Teuchos::arcp;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Tpetra::Map<TPETRA_LO, TPETRA_GO> map_type;
  typedef Tpetra::MultiVector<TPETRA_SCALAR, TPETRA_LO, TPETRA_GO> MV;
  typedef ArrayView<const TPETRA_SCALAR> coordList_t;

  const int nprocs = comm->getSize ();
  const int rank = comm->getRank ();

  ///////////// Step 1 //////////////////////////////////
  // Create a MV with contiguous global IDs

  const TPETRA_LO numLocalCoords = numSequentialGlobalIds(numGlobalCoords, nprocs, rank);

  RCP<const map_type> tmap;
  RCP<MV> mvector;
  TPETRA_SCALAR* coords = NULL;
  {
    Teuchos::TimeMonitor timeMon (*tmvBuild);

    tmap = rcp (new map_type (numGlobalCoords, numLocalCoords, 0, comm));

    coords = new TPETRA_SCALAR [COORDDIM * numLocalCoords];
    memset (coords, 0, sizeof(TPETRA_SCALAR) * numLocalCoords * COORDDIM);

    coordList_t *avList = new coordList_t [COORDDIM];
    TPETRA_LO offset = 0;

    for (int dim = 0; dim < COORDDIM; ++dim) {
      avList[dim] = coordList_t(coords + offset, numLocalCoords);
      offset += numLocalCoords;
    }

    ArrayRCP<const coordList_t> vectors = arcp (avList, 0, COORDDIM);
    mvector = rcp (new MV (tmap, vectors.view (0, COORDDIM), COORDDIM));
  }

  if (rank == 0 && doMemory) {
    const long nkb = Zoltan2::getProcessKilobytes ();
    std::cout << "Create mvector 1: " << nkb << std::endl;
  }

  ///////////// Step 2 //////////////////////////////////
  // Migrate the MV.

  ArrayRCP<const TPETRA_GO> newGidArray;
  {
    TPETRA_GO *newGids = NULL;
    roundRobinGlobalIds<TPETRA_GO> (numGlobalCoords, nprocs, rank, newGids);
    newGidArray = arcp<const TPETRA_GO> (newGids, 0, numLocalCoords, true);
  }

  RCP<const map_type> newTmap;
  RCP<Tpetra::Import<TPETRA_LO, TPETRA_GO> > importer;
  RCP<MV> newMvector;
  {
    Teuchos::TimeMonitor timeMon (*tmvMigrate);

    newTmap = rcp (new map_type (numGlobalCoords, newGidArray.view(0, numLocalCoords), 0, comm));
    importer = rcp (new Tpetra::Import<TPETRA_LO, TPETRA_GO> (tmap, newTmap));
    newMvector = rcp (new MV (newTmap, COORDDIM, true));

    newMvector->doImport (*mvector, *importer, Tpetra::INSERT);
    mvector = newMvector;
  }

  delete [] coords;

  if (rank == 0 && doMemory) {
    const long nkb = Zoltan2::getProcessKilobytes ();
    std::cout << "Create mvector 2: " << nkb << std::endl;
  }

  ///////////// Step 3 //////////////////////////////////
  // Divide processes into two halves.

  RCP<Comm<int> > subComm;
  {
    int groupSize = 0;
    int leftHalfNumProcs = nprocs / 2;
    int *myHalfProcs = NULL;

    if (rank < leftHalfNumProcs){
      groupSize = leftHalfNumProcs;
      myHalfProcs = new int [groupSize];
      for (int i=0; i < groupSize; i++)
        myHalfProcs[i] = i;
    }
    else {
      groupSize = nprocs - leftHalfNumProcs;
      myHalfProcs = new int [groupSize];
      int firstNum = leftHalfNumProcs;
      for (int i=0; i < groupSize; i++)
        myHalfProcs[i] = firstNum++;
    }

    ArrayView<const int> idView(myHalfProcs, groupSize);
    subComm = comm->createSubcommunicator (idView);
    delete [] myHalfProcs;
  }

  // Divide the multivector into two.  Each process group is creating
  // a multivector with non-contiguous global ids.  For one group,
  // base gid is not 0.

  size_t globalSize = Teuchos::OrdinalTraits<size_t>::invalid ();
  RCP<map_type> subMap;
  RCP<MV> subMvector;
  {
    Teuchos::TimeMonitor timeMon (*tmvBuildN);

    ArrayView<const TPETRA_GO> gidList = mvector->getMap ()->getLocalElementList ();
    subMap = rcp (new map_type (globalSize, gidList, 0, subComm));
    globalSize = subMap->getGlobalNumElements ();

    // Get a view of the block of rows to copy.
    RCP<MV> tmp = mvector->offsetViewNonConst (subMap, 0);
    // Create a new multivector to hold the group's rows.
    subMvector = rcp (new MV (subMap, mvector->getNumVectors ()));
    // Copy the view into the new multivector.
    Tpetra::deep_copy (*subMvector, *tmp);
  }

  ///////////// Step 4 //////////////////////////////////
  // Each subgroup migrates the sub-multivector so the
  // global Ids are increasing with process rank.

  TPETRA_GO *increasingGids = NULL;
  subGroupGloballyIncreasingIds<TPETRA_GO> (numGlobalCoords, nprocs,
                                     rank, increasingGids);

  ArrayRCP<const TPETRA_GO> incrGidArray (increasingGids, 0, numLocalCoords, true);

  RCP<const map_type> newSubMap;
  RCP<Tpetra::Import<TPETRA_LO, TPETRA_GO> > subImporter;
  RCP<MV> newSubMvector;
  {
    Teuchos::TimeMonitor timeMon (*tmvMigrateN);

    newSubMap =
      rcp (new map_type (globalSize, incrGidArray.view (0, numLocalCoords),
                         0, subComm));
    subImporter = rcp (new Tpetra::Import<TPETRA_LO, TPETRA_GO> (subMap, newSubMap));
    newSubMvector = rcp (new MV (newSubMap, COORDDIM, true));
    newSubMvector->doImport (*subMvector, *subImporter, Tpetra::INSERT);
    mvector = newSubMvector;
  }
}

void timeEpetra(TPETRA_GO numGlobalCoords, const RCP<const MpiComm<int> > &comm,
  bool doMemory)
{
  RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > commPtr =
    comm->getRawMpiComm();

  RCP<Epetra_MpiComm> ecomm = rcp(new Epetra_MpiComm((*commPtr)()));

  int nprocs = comm->getSize();
  int rank = comm->getRank();

  ///////////// Step 1 //////////////////////////////////
  // Create a MV with contiguous global IDs

  TPETRA_LO numLocalCoords = numSequentialGlobalIds(numGlobalCoords, nprocs, rank);

  emvBuild->start();
  emvBuild->incrementNumCalls();

  RCP<Epetra_BlockMap> emap = rcp(new Epetra_BlockMap(numGlobalCoords,
    numLocalCoords, 1, 0, *ecomm));

  TPETRA_SCALAR *coords = new TPETRA_SCALAR [COORDDIM * numLocalCoords];
  memset(coords, 0, sizeof(TPETRA_SCALAR) * numLocalCoords * COORDDIM);

  RCP<Epetra_MultiVector> mvector = rcp(new Epetra_MultiVector(
    View, *emap, coords, 1, COORDDIM));

  emvBuild->stop();

  if (rank==0 && doMemory){
    long nkb = Zoltan2::getProcessKilobytes();
    std::cout << "Create mvector 1: " << nkb << std::endl;;
  }

  ///////////// Step 2 //////////////////////////////////
  // Migrate the MV.

  TPETRA_GO *newGids = NULL;
  roundRobinGlobalIds<TPETRA_GO>(numGlobalCoords, nprocs, rank, newGids);

  emvMigrate->start();
  emvMigrate->incrementNumCalls();

  RCP<Epetra_BlockMap> newMap = rcp(new Epetra_BlockMap(numGlobalCoords,
    numLocalCoords, newGids, 1, 0, *ecomm));

  RCP<Epetra_Import> importer = rcp(new Epetra_Import(*newMap, *emap));

  RCP<Epetra_MultiVector> newMvector = rcp(new Epetra_MultiVector(
    *newMap, COORDDIM));

  newMvector->Import(*mvector, *importer, Insert);

  mvector = newMvector;

  emvMigrate->stop();

  delete [] coords;

  if (rank==0 && doMemory){
    long nkb = Zoltan2::getProcessKilobytes();
    std::cout << "Create mvector 2: " << nkb << std::endl;;
  }

  ///////////// Step 3 //////////////////////////////////
  // Divide processes into two halves.

  int groupSize = 0;
  int leftHalfNumProcs = nprocs / 2;
  int *myHalfProcs = NULL;

  if (rank < leftHalfNumProcs){
    groupSize = leftHalfNumProcs;
    myHalfProcs = new int [groupSize];
    for (int i=0; i < groupSize; i++)
      myHalfProcs[i] = i;
  }
  else {
    groupSize = nprocs - leftHalfNumProcs;
    myHalfProcs = new int [groupSize];
    int firstNum = leftHalfNumProcs;
    for (int i=0; i < groupSize; i++)
      myHalfProcs[i] = firstNum++;
  }

  ArrayView<const int> idView(myHalfProcs, groupSize);
  // TODO - memory leak in createSubcommunicator.
  RCP<Comm<int> > newComm = comm->createSubcommunicator(idView);
  RCP<MpiComm<int> > genSubComm = rcp_dynamic_cast<MpiComm<int> >(newComm);

  commPtr = genSubComm->getRawMpiComm();

  RCP<Epetra_MpiComm> subComm = rcp(new Epetra_MpiComm((*commPtr)()));

  delete [] myHalfProcs;

  // Divide the multivector into two.  Each process group is creating
  // a multivector with non-contiguous global ids.  For one group,
  // base gid is not 0.

  emvBuildN->start();
  emvBuildN->incrementNumCalls();

  RCP<Epetra_BlockMap> subMap = rcp(new Epetra_BlockMap(-1,
   numLocalCoords, newGids, 1, 0, *subComm));

  TPETRA_SCALAR **avSubList = new TPETRA_SCALAR * [COORDDIM];

  for (int dim=0; dim < COORDDIM; dim++)
    (*mvector)(dim)->ExtractView(avSubList + dim);

  RCP<Epetra_MultiVector> subMvector = rcp(new Epetra_MultiVector(
    View, *subMap, avSubList, COORDDIM));

  mvector = subMvector;

  delete [] avSubList;
  delete [] newGids;

  emvBuildN->stop();

  ///////////// Step 4 //////////////////////////////////
  // Each subgroup migrates the sub-multivector so the
  // global Ids are increasing with process rank.

  TPETRA_GO *increasingGids = NULL;
  subGroupGloballyIncreasingIds<TPETRA_GO>(numGlobalCoords, nprocs, rank,
    increasingGids);

  emvMigrateN->start();
  emvMigrateN->incrementNumCalls();

  RCP<Epetra_BlockMap> newSubMap = rcp(new Epetra_BlockMap(-1,
    numLocalCoords, increasingGids, 1, 0, *subComm));

  RCP<Epetra_Import> subImporter = rcp(new Epetra_Import(
    *subMap, *newSubMap));

  RCP<Epetra_MultiVector> newSubMvector = rcp(new Epetra_MultiVector(
    *newSubMap, COORDDIM));

  newSubMvector->Import(*subMvector, *subImporter, Insert);

  mvector = newSubMvector;

  emvMigrateN->stop();

  delete [] increasingGids;
}

void timeZoltan(ZOLTAN_GO numGlobalCoords,
  bool doMemory)
{
  int nprocs, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  ///////////// Step 1 //////////////////////////////////
  // Create a global data directory with contiguous global IDs.
  // (We don't need this, but it is analygous to a Tpetra::Map.)

  ZOLTAN_LO numLocalCoords = numSequentialGlobalIds(numGlobalCoords, nprocs, rank);

  ZOLTAN_GO offset=0;
  MPI_Exscan(&numLocalCoords, &offset, 1,
    MPI_ZOLTAN_GO_TYPE, MPI_SUM, MPI_COMM_WORLD);

  ZOLTAN_GO *gids = new ZOLTAN_GO [numLocalCoords];
  for (ZOLTAN_LO i=0; i < numLocalCoords; i++){
    gids[i] = offset++;
  }

  ztnBuild->start();
  ztnBuild->incrementNumCalls();

  struct Zoltan_DD_Struct *dd = NULL;
  int rc = Zoltan_DD_Create(&dd, MPI_COMM_WORLD, 1, 1, 0, numLocalCoords, 0);
  if (rc != ZOLTAN_OK)
    exit(1);

  rc = Zoltan_DD_Update(dd, gids, NULL, NULL, NULL, numLocalCoords);

  // Create an array of coordinates associated with the global Ids.

  TPETRA_SCALAR *coords = new TPETRA_SCALAR [COORDDIM * numLocalCoords];
  memset(coords, 0, sizeof(TPETRA_SCALAR) * numLocalCoords * COORDDIM);

  ztnBuild->stop();

  if (rank==0 && doMemory){
    long nkb = Zoltan2::getProcessKilobytes();
    std::cout << "Create mvector 1: " << nkb << std::endl;;
  }

  Zoltan_DD_Destroy(&dd);

  ///////////// Step 2 //////////////////////////////////
  // Migrate the array of coordinates.

  ZOLTAN_GO *newGids = NULL;
  roundRobinGlobalIds<ZOLTAN_GO>(numGlobalCoords, nprocs, rank, newGids);

  ztnMigrate->start();
  ztnMigrate->incrementNumCalls();

  struct Zoltan_DD_Struct *ddNew = NULL;  // new "map"
  rc = Zoltan_DD_Create(&ddNew, MPI_COMM_WORLD, 1, 1, 0, numLocalCoords, 0);
  if (rc != ZOLTAN_OK)
    exit(1);

  rc = Zoltan_DD_Update(ddNew, newGids, NULL, NULL, NULL, numLocalCoords);
  if (rc != ZOLTAN_OK)
    exit(1);

  int *procOwners = new int [numLocalCoords];  // procs to get my data
  rc = Zoltan_DD_Find(ddNew, gids, NULL, NULL, NULL,
    numLocalCoords, procOwners);
  if (rc != ZOLTAN_OK)
    exit(1);

  Zoltan_DD_Destroy(&ddNew);

  struct Zoltan_Comm_Obj *commPlan = NULL;  // global communication plan
  int tag = 10000;
  int numReceive = 0;

  rc = Zoltan_Comm_Create(&commPlan, numLocalCoords, procOwners, MPI_COMM_WORLD,
    tag, &numReceive);
  if (rc != ZOLTAN_OK)
    exit(1);

  TPETRA_SCALAR *newCoords = new TPETRA_SCALAR [COORDDIM * numReceive];

  tag = 11000;

  // To prevent compile warnings or errors
  void *x = static_cast<void *>(coords);
  char *charCoords = static_cast<char *>(x);
  x = static_cast<void *>(newCoords);
  char *charNewCoords = static_cast<char *>(x);

  rc = Zoltan_Comm_Do(commPlan, tag, charCoords,
    sizeof(TPETRA_SCALAR)*COORDDIM, charNewCoords);

  if (rc != ZOLTAN_OK)
    exit(1);

  ztnMigrate->stop();

  Zoltan_Comm_Destroy(&commPlan);
  delete [] coords;
  delete [] gids;

  if (rank==0 && doMemory){
    long nkb = Zoltan2::getProcessKilobytes();
    std::cout << "Create mvector 2: " << nkb << std::endl;;
  }

  ///////////// Step 3 //////////////////////////////////
  // Divide processes into two halves.

  int groupSize = 0;
  int leftHalfNumProcs = nprocs / 2;
  int *myHalfProcs = NULL;

  if (rank < leftHalfNumProcs){
    groupSize = leftHalfNumProcs;
    myHalfProcs = new int [groupSize];
    for (int i=0; i < groupSize; i++)
      myHalfProcs[i] = i;
  }
  else {
    groupSize = nprocs - leftHalfNumProcs;
    myHalfProcs = new int [groupSize];
    for (int i=0; i < groupSize; i++)
      myHalfProcs[i] = i + leftHalfNumProcs;
  }

  MPI_Group group, subGroup;
  MPI_Comm subComm;

  MPI_Comm_group(MPI_COMM_WORLD, &group);
  MPI_Group_incl(group, groupSize, myHalfProcs, &subGroup);
  MPI_Comm_create(MPI_COMM_WORLD, subGroup, &subComm);
  MPI_Group_free(&subGroup);

  delete [] myHalfProcs;

  // Create global data directories for our sub groups.
  // (Analogous to creating the new MultiVectors in Tpetra.)

  ztnBuildN->start();
  ztnBuildN->incrementNumCalls();

  struct Zoltan_DD_Struct *ddSub = NULL;  // subgroup "map"
  rc = Zoltan_DD_Create(&ddSub, subComm, 1, 1, 0, numLocalCoords, 0);
  if (rc != ZOLTAN_OK)
    exit(1);

  rc = Zoltan_DD_Update(ddSub, newGids, NULL, NULL, NULL, numLocalCoords);
  if (rc != ZOLTAN_OK)
    exit(1);

  ztnBuildN->stop();

  Zoltan_DD_Destroy(&ddSub);

  ///////////// Step 4 //////////////////////////////////
  // Each subgroup migrates the sub-arrays so the
  // global Ids are again increasing with process rank.

  ZOLTAN_GO *increasingGids = NULL;
  subGroupGloballyIncreasingIds<ZOLTAN_GO>(
    numGlobalCoords, nprocs, rank, increasingGids);

  // Global "map" corresponding to new contiguous ids.

  ztnMigrateN->start();
  ztnMigrateN->incrementNumCalls();

  struct Zoltan_DD_Struct *ddNewSub = NULL;
  rc = Zoltan_DD_Create(&ddNewSub, subComm, 1, 1, 0, numLocalCoords, 0);
  if (rc != ZOLTAN_OK)
    exit(1);

  rc = Zoltan_DD_Update(ddNewSub, increasingGids, NULL, NULL, NULL,
    numLocalCoords);

  // Which processes gets my current coordinates in new map?

  rc = Zoltan_DD_Find(ddNewSub, newGids, NULL, NULL, NULL, numLocalCoords, procOwners);
  if (rc != ZOLTAN_OK)
    exit(1);

  delete [] newGids;

  Zoltan_DD_Destroy(&ddNewSub);

  struct Zoltan_Comm_Obj *subCommPlan = NULL;  // global communication plan
  tag = 12000;

  rc = Zoltan_Comm_Create(&subCommPlan, numLocalCoords, procOwners, subComm,
    tag, &numReceive);
  if (rc != ZOLTAN_OK)
    exit(1);

  delete [] procOwners;

  TPETRA_SCALAR *newContigCoords = new TPETRA_SCALAR [COORDDIM * numReceive];

  tag = 13000;
  // To prevent compile warnings or errors
  x = static_cast<void *>(newContigCoords);
  char *charNewContigCoords = static_cast<char *>(x);

  rc = Zoltan_Comm_Do(subCommPlan, tag, charNewCoords,
    sizeof(TPETRA_SCALAR)*COORDDIM, charNewContigCoords);
  if (rc != ZOLTAN_OK)
    exit(1);

  ztnMigrateN->stop();

  delete [] newCoords;
  delete [] newContigCoords;
  delete [] increasingGids;
  Zoltan_Comm_Destroy(&subCommPlan);
}
#else
int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > genComm = Tpetra::getDefaultComm();

  if (genComm->getRank() == 0){
    std::cout << "Test not run because MPI is not available." << std::endl;
    std::cout << "PASS" << std::endl;
  }
  return 0;
}
#endif  // HAVE_MPI
