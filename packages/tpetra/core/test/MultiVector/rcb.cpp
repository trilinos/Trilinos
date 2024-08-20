// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Ensure that if CUDA and KokkosCompat are enabled, then only the .cu
// version of this file will actually be compiled.
#include <Tpetra_ConfigDefs.hpp>

// This test comes from zoltan2/test/temp/multivectorTest.cpp.
// It mimics the recursive bisection algorithm in Zoltan2.
//
// Create a Tpetra::MultiVector, and time the following:
//   1.  Build a multivector with contiguous global ids.
//   2.  Migrate.
//   3.  Divide procs into two groups and build new a new multivector
//       in each group with the non-contiguous global ids of the migrated data.
//   4.  Migrate the new multivector in the subgroup.
//
// It only makes sense to run this test in an MPI build.  Otherwise,
// no data migration is necessary.
//
// mfh 29 Jul 2014: I adapted this test from Zoltan2.  They do some
// things that are not idiomatic Tpetra.  I wanted to get this test in
// the repository because it exercises some interesting compatibility
// features of the Kokkos refactor version of Tpetra.  In particular,
// I had to change the test (see comments close to the call to
// offsetViewNonConst) to fix incorrect use of MultiVector::getData.

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Import.hpp>

#include <string>
#include <sstream>
#include <iostream>

typedef Tpetra::Map<>::local_ordinal_type LO;
typedef Tpetra::Map<>::global_ordinal_type GO;
typedef Tpetra::MultiVector<>::scalar_type Scalar;

Teuchos::RCP<Teuchos::Time> tmvBuild;
Teuchos::RCP<Teuchos::Time> tmvMigrate;
Teuchos::RCP<Teuchos::Time> tmvBuildN;
Teuchos::RCP<Teuchos::Time> tmvMigrateN;

enum testConstants {COORDDIM=3};

static void usage (char *argv[]) {
  std::cout << "Usage:" << std::endl;
  std::cout << argv[0] << " {num coords}" << std::endl;
}

///////////////////////////////////////////////
// Generate global Ids for different mappings
///////////////////////////////////////////////
LO numSequentialGlobalIds (const GO numGlobalCoords,
                           const int nprocs, const int rank)
{
  LO share = numGlobalCoords / nprocs;
  LO extra = numGlobalCoords % nprocs;
  LO numLocalCoords = share;
  if (rank < extra) {
    numLocalCoords++;
  }
  return numLocalCoords;
}

template <typename T>
Teuchos::ArrayRCP<T>
roundRobinGlobalIds (const T numGlobalCoords,
                     const int nprocs,
                     const int rank)
{
  // Local number of coordinates does not change.
  T share = numGlobalCoords / nprocs;
  T extra = numGlobalCoords % nprocs;
  T numLocalCoords = share;
  if (static_cast<T> (rank) < extra) {
    numLocalCoords++;
  }

  Teuchos::ArrayRCP<T> gids (numLocalCoords);

  T next = 0;
  for (T i=rank; i < numGlobalCoords; i+=nprocs) {
    gids[next++] = i;
  }

  return gids;
}

template <typename T>
Teuchos::ArrayRCP<T>
subGroupGloballyIncreasingIds (const T numGlobalCoords,
                               const int nprocs,
                               const int rank)
{
  const int numProcsLeftHalf = nprocs / 2;
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
  if (rank < static_cast<T> (extra))
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
  Teuchos::ArrayRCP<T> gids (numLocalCoords);

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

  return gids;
}

void
timeTpetra (const GO numGlobalCoords,
            const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
            const bool doMemory);

///////////////////////////////////////////////
// Main
///////////////////////////////////////////////

int
main (int argc, char *argv[])
{
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using std::cerr;
  using std::cout;
  using std::endl;

  Teuchos::GlobalMPISession session (&argc, &argv, NULL);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
  const int rank = comm->getRank ();

  GO numGlobalCoords = 0;
  if (argc < 2) {
    numGlobalCoords = comm->getSize () * 100;
  }
  else {
    std::string theArg(argv[1]);
    std::istringstream iss(theArg);
    iss >> numGlobalCoords;
    if (numGlobalCoords < comm->getSize ()) {
      if (rank == 0) {
        usage (argv);
      }
      return EXIT_FAILURE;
    }
  }

  if (rank == 0) {
    cerr << numGlobalCoords << " coordinates" << endl;
  }

  tmvBuild = TimeMonitor::getNewTimer("CONSEC build Tpetra");
  tmvMigrate = TimeMonitor::getNewTimer("CONSEC migrate Tpetra");
  tmvBuildN = TimeMonitor::getNewTimer("!CONSEC build Tpetra");
  tmvMigrateN = TimeMonitor::getNewTimer("!CONSEC migrate Tpetra");

  TimeMonitor::zeroOutTimers ();

  int ntests = 3;

  // Test with Tpetra::MultiVector
  for (int i=0; i < ntests; i++) {
    if (rank == 0) {
      cerr << "Tpetra test " << i+1 << endl;
    }
    timeTpetra (numGlobalCoords, comm, i==ntests-1);
  }

  TimeMonitor::summarize (); // report timing results
  if (comm->getRank() == 0) {
    cout << "End Result: TEST PASSED" << endl;
  }
  return EXIT_SUCCESS;
}

void
timeTpetra (const GO numGlobalCoords,
            const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
            const bool doMemory)
{
  using Teuchos::arcp;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cerr;
  using std::endl;
  typedef Tpetra::Map<LO, GO> map_type;
  typedef Tpetra::MultiVector<Scalar, LO, GO> MV;
  typedef ArrayView<const Scalar> coordList_t;

  const int nprocs = comm->getSize ();
  const int rank = comm->getRank ();

  ///////////// Step 1 //////////////////////////////////
  // Create a MV with contiguous global IDs

  if (rank == 0) {
    cerr << "Step 1: Create a MultiVector with contiguous GIDs" << endl;
  }
  const LO numLocalCoords = numSequentialGlobalIds (numGlobalCoords, nprocs, rank);

  RCP<const map_type> tmap;
  RCP<MV> mvector;
  Scalar* coords = NULL;
  {
    Teuchos::TimeMonitor timeMon (*tmvBuild);

    tmap = rcp (new map_type (numGlobalCoords, numLocalCoords, 0, comm));

    coords = new Scalar [COORDDIM * numLocalCoords];
    memset (coords, 0, sizeof(Scalar) * numLocalCoords * COORDDIM);

    coordList_t *avList = new coordList_t [COORDDIM];
    LO offset = 0;

    for (int dim = 0; dim < COORDDIM; ++dim) {
      avList[dim] = coordList_t(coords + offset, numLocalCoords);
      offset += numLocalCoords;
    }

    ArrayRCP<const coordList_t> vectors = arcp (avList, 0, COORDDIM);
    mvector = rcp (new MV (tmap, vectors.view (0, COORDDIM), COORDDIM));
  }

  ///////////// Step 2 //////////////////////////////////
  // Migrate the MV.

  if (rank == 0) {
    cerr << "Step 2: Migrate the MultiVector" << endl;
  }

  ArrayRCP<const GO> newGidArray =
    roundRobinGlobalIds<GO> (numGlobalCoords, nprocs, rank);

  RCP<const map_type> newTmap;
  RCP<Tpetra::Import<LO, GO> > importer;
  RCP<MV> newMvector;
  {
    Teuchos::TimeMonitor timeMon (*tmvMigrate);

    newTmap = rcp (new map_type (numGlobalCoords, newGidArray.view(0, numLocalCoords), 0, comm));
    importer = rcp (new Tpetra::Import<LO, GO> (tmap, newTmap));
    newMvector = rcp (new MV (newTmap, COORDDIM, true));

    newMvector->doImport (*mvector, *importer, Tpetra::INSERT);
    mvector = newMvector;
  }

  delete [] coords;

  ///////////// Step 3 //////////////////////////////////
  // Divide processes into two halves.

  if (rank == 0) {
    cerr << "Step 3: Divide processes into two halves" << endl;
  }

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

    //const size_t localSize = mvector->getLocalLength ();
    ArrayView<const GO> gidList = mvector->getMap ()->getLocalElementList ();
    subMap = rcp (new map_type (globalSize, gidList, 0, subComm));
    globalSize = subMap->getGlobalNumElements ();

    // FIXME (mfh 29 Jul 2014) The commented-out code below is unsafe,
    // because it gets non-reference-counted views of the
    // MultiVector's data.  That amounts to a dangling reference,
    // given ArrayRCP semantics.  Furthermore, the code is
    // unnecessarily complicated, since it can be replaced with
    // offsetViewNonConst and a deep copy.  We do this below.

    // Teuchos::Array<ArrayView<const Scalar> > avSublist (COORDDIM);
    // //coordList_t *avSubList = new coordList_t [COORDDIM];
    // for (int dim = 0; dim < COORDDIM; ++dim) {
    //   avSublist[dim] = mvector->getData (dim).view (0, localSize);
    // }
    // // ArrayRCP<const ArrayView<const Scalar> > subVectors =
    // //   arcp (avSubList, 0, COORDDIM);
    // //subMvector = rcp (new MV (subMap, subVectors.view (0, COORDDIM), COORDDIM));
    // subMvector = rcp (new MV (subMap, avSublist (), COORDDIM));

    RCP<MV> tmp = mvector->offsetViewNonConst (subMap, 0); // to copy
    subMvector = rcp (new MV (subMap, mvector->getNumVectors ()));
    Tpetra::deep_copy (*subMvector, *tmp);
  }

  if (rank == 0) {
    cerr << "Step 4: migrate the sub-multivector" << endl;
  }

  ///////////// Step 4 //////////////////////////////////
  // Each subgroup migrates the sub-multivector so the
  // global Ids are increasing with process rank.

  ArrayRCP<const GO> incrGidArray =
    subGroupGloballyIncreasingIds<GO> (numGlobalCoords, nprocs, rank);

  RCP<const map_type> newSubMap;
  RCP<Tpetra::Import<LO, GO> > subImporter;
  RCP<MV> newSubMvector;
  {
    Teuchos::TimeMonitor timeMon (*tmvMigrateN);

    newSubMap = rcp (new map_type (globalSize, incrGidArray.view (0, numLocalCoords), 0, subComm));
    subImporter = rcp (new Tpetra::Import<LO, GO> (subMap, newSubMap));
    newSubMvector = rcp (new MV (newSubMap, COORDDIM, true));
    newSubMvector->doImport (*subMvector, *subImporter, Tpetra::INSERT);
    mvector = newSubMvector;
  }

  if (rank == 0) {
    cerr << "Done" << endl;
  }
}

