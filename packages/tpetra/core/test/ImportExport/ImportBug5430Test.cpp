// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This is a test which manifests Bug 5430, a bug in the Tpetra Import
// constructor.  It builds a test problem concurrently using Epetra
// and Tpetra objects.  The Epetra test works, but the Tpetra test
// results in incorrect communication.  That manifests as incorrectly
// imported values (visible as -1 values) in the target vector of the
// Import.
//
// The original test was written by Bill Spotz, Sandia National
// Laboratories, wfspotz@sandia.gov.  It was expanded by Mark Hoemmen,
// Sandia National Laboratories, mhoemme@sandia.gov.

//
// The DIMENSIONS macro can equal "2" or "3".  The bug manifests only
// in the 3-D case, not in the 2-D case.
//
#define DIMENSIONS 3

#include <Tpetra_ConfigDefs.hpp>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

// Epetra includes
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Vector.hpp"

// mfh 16 Dec 2011: Teuchos::GlobalMPISession doesn't seem to work
// with this example; the first call to a Teuchos::Comm method makes
// MPI kill the program.  I'm not sure why this is.
#ifdef HAVE_MPI
#  include "mpi.h" // Needed for Epetra comparison
#else
#  error "This example requires building with MPI."
#endif // HAVE_MPI

namespace {

bool
countFailures (const Teuchos::RCP<const Teuchos::Comm<int> >& teuchosComm,
               // Epetra widgets
               const Epetra_Map& epetraOwnedMap,
               const Epetra_Vector& epetraOwnedVector,
               const Epetra_Map& epetraOverlapMap,
               Epetra_Vector& epetraOverlapVector,
               // Tpetra widgets
               const Teuchos::RCP<Tpetra::Map<int> >& tpetraOwnedMap,
               const Tpetra::Vector<double,int>& tpetraOwnedVector,
               const Teuchos::RCP<Tpetra::Map<int> >& tpetraOverlapMap,
               Tpetra::Vector<double,int>& tpetraOverlapVector,
               // options
               const bool verbose)
{
  using Teuchos::ArrayRCP;
  using Teuchos::ptr;
  using Teuchos::RCP;
  using Teuchos::reduceAll;
  using std::cout;
  using std::endl;

  const int pid = teuchosComm->getRank();

  // (Const) view of the Tpetra Vector's data, for comparison with the
  // Epetra Vector's data.
  ArrayRCP<const double> tpetraOverlapArray = tpetraOverlapVector.getData(0);

  // Loop over all elements of the overlap data and verify that they
  // equal the corresponding global IDs.  If things are working
  // correctly, nothing should get printed out.  What happens is that
  // the Epetra objects behave as expected, but the Tpetra data is
  // incorrect in the overlap regions.
  int localEpetraFailureCount = 0;
  int localTpetraFailureCount = 0;
  std::ostringstream localOut;
  for (int lid = 0; lid < tpetraOverlapMap->getLocalElementList().size(); ++lid) {
    int gid = tpetraOverlapMap->getGlobalElement(lid);
    if (epetraOverlapVector[lid] != gid) {
      localOut << "Process " << pid << ": epetraOverlapVector[" << lid << "] = "
               << epetraOverlapVector[lid] << " (should equal " << gid << ")"
               << endl;
      ++localEpetraFailureCount;
    }
    if (tpetraOverlapArray[lid] != gid) {
      localOut << "Process " << pid << ": tpetraOverlapArray[" << lid << "] = "
           << tpetraOverlapArray[lid] << " (should equal " << gid << ")"
           << endl;
      ++localTpetraFailureCount;
    }
  }

  // In verbose mode, print out the collected error messages from the
  // different processes.  We assume here that any MPI process can
  // print to stdout.
  if (verbose) {
    const int numProcs = teuchosComm->getSize();
    for (int p = 0; p < numProcs; ++p) {
      if (p == pid) {
        cout << "Process " << p << ":" << endl;
        if (localTpetraFailureCount > 0) {
          cout << localOut.str() << endl;
        } else {
          cout << "No errors" << endl;
        }
        cout << std::flush; // Write an endl and flush the output stream.
      }
      // Do a few barriers to ensure that output to stdout completed.
      teuchosComm->barrier();
      teuchosComm->barrier();
      teuchosComm->barrier();
    }
  }

  // Get global failure counts.
  int globalEpetraFailureCount = 0;
  int globalTpetraFailureCount = 0;
  reduceAll (*teuchosComm, Teuchos::REDUCE_SUM, localEpetraFailureCount,
             ptr (&globalEpetraFailureCount));
  reduceAll (*teuchosComm, Teuchos::REDUCE_SUM, localTpetraFailureCount,
             ptr (&globalTpetraFailureCount));

  // Report the number of failures over all processes.
  Teuchos::oblackholestream blackHole;
  std::ostream& globalOut = (pid == 0) ? std::cout : blackHole;
  bool success = true;
  if (globalEpetraFailureCount > 0) {
    globalOut << "Number of Epetra failures: " << globalEpetraFailureCount << endl;
    success = false;
  }
  if (globalTpetraFailureCount > 0) {
    globalOut << "Number of Tpetra failures: " << globalTpetraFailureCount << endl;
    success = false;
  }
  return success;
}

} // namespace (anonymous)

bool
testMain (int argc, char *argv[])
{
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::CommandLineProcessor;
  using Teuchos::FancyOStream;
  using Teuchos::getFancyOStream;
  using Teuchos::OSTab;
  using Teuchos::ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using std::cout;
  using std::endl;

  bool success = true; // May be changed by tests

  Teuchos::oblackholestream blackHole;

  //
  // Construct communicators, and verify that we are on 4 processors.
  //

  // Construct a Teuchos Comm object.
  RCP<const Comm<int> > teuchosComm = Teuchos::DefaultComm<int>::getComm();
  const int numProcs = teuchosComm->getSize();
  const int pid = teuchosComm->getRank();
  RCP<FancyOStream> pOut =
    getFancyOStream (rcpFromRef ((pid == 0) ? std::cout : blackHole));
  FancyOStream& out = *pOut;
  // Verify that we are on four processors (which manifests the bug).
  if (teuchosComm->getSize() != 4) {
    out << "This test must be run on four processors.  Exiting ..." << endl;
    return false;
  }

  // We also need an Epetra Comm, so that we can compare Tpetra and
  // Epetra results.
  Epetra_MpiComm epetraComm (MPI_COMM_WORLD);

  //
  // Default values of command-line options.
  //
  bool verbose = false;
  bool printEpetra = false;
  bool printTpetra = false;
  CommandLineProcessor cmdp (false,true);
  //
  // Set command-line options.
  //
  cmdp.setOption ("verbose", "quiet", &verbose, "Print verbose output.");
  // Epetra and Tpetra output will ask the Maps and Import objects to
  // print themselves in distributed, maximally verbose fashion.  It's
  // best to turn on either Epetra or Tpetra, but not both.  Then you
  // can compare their output side by side.
  cmdp.setOption ("printEpetra", "dontPrintEpetra", &printEpetra,
                  "Print Epetra output (in verbose mode only).");
  cmdp.setOption ("printTpetra", "dontPrintTpetra", &printTpetra,
                  "Print Tpetra output (in verbose mode only).");
  // Parse command-line options.
  if (cmdp.parse (argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return false;
  }

  if (verbose) {
    out << "Running test on " << numProcs << " process"
        << (numProcs != 1 ? "es" : "") << "." << endl;
  }

  // The maps for this problem are derived from a 3D structured mesh.
  // In this example, the dimensions are 4x4x2 and there are 2
  // processors assigned to the first dimension and 2 processors
  // assigned to the second dimension, with no parallel decomposition
  // along the third dimension.  The "owned" arrays represent the
  // one-to-one map, with each array representing a 2x2x2 slice.  If
  // DIMENSIONS == 2, then only the first 4 values will be used,
  // representing a 2x2(x1) slice.
  int owned0[8] = { 0, 1, 4, 5,16,17,20,21};
  int owned1[8] = { 2, 3, 6, 7,18,19,22,23};
  int owned2[8] = { 8, 9,12,13,24,25,28,29};
  int owned3[8] = {10,11,14,15,26,27,30,31};

  // The "overlap" arrays represent the map with communication
  // elements, with each array representing a 3x3x2 slice.  If
  // DIMENSIONS == 2, then only the first 9 values will be used,
  // representing a 3x3(x1) slice.
  int overlap0[18] = {0,1,2,4, 5, 6, 8, 9,10,16,17,18,20,21,22,24,25,26};
  int overlap1[18] = {1,2,3,5, 6, 7, 9,10,11,17,18,19,21,22,23,25,26,27};
  int overlap2[18] = {4,5,6,8, 9,10,12,13,14,20,21,22,24,25,26,28,29,30};
  int overlap3[18] = {5,6,7,9,10,11,13,14,15,21,22,23,25,26,27,29,30,31};

  // Construct the owned and overlap maps for both Epetra and Tpetra.
  int* owned;
  int* overlap;
  if (pid == 0) {
    owned   = owned0;
    overlap = overlap0;
  }
  else if (pid == 1) {
    owned   = owned1;
    overlap = overlap1;
  }
  else if (pid == 2) {
    owned   = owned2;
    overlap = overlap2;
  }
  else {
    owned   = owned3;
    overlap = overlap3;
  }

#if DIMENSIONS == 2
  int ownedSize   = 4;
  int overlapSize = 9;
#elif DIMENSIONS == 3
  int ownedSize   =  8;
  int overlapSize = 18;
#endif

  // Create the two Epetra Maps.  Source for the Import is the owned
  // map; target for the Import is the overlap map.
  Epetra_Map epetraOwnedMap (  -1, ownedSize,   owned,   0, epetraComm);
  Epetra_Map epetraOverlapMap (-1, overlapSize, overlap, 0, epetraComm);

  if (verbose && printEpetra) {
    // Have the Epetra_Map objects describe themselves.
    //
    // Epetra_BlockMap::Print() takes an std::ostream&, and expects
    // all MPI processes to be able to write to it.  (The method
    // handles its own synchronization.)
    out << "Epetra owned map:" << endl;
    epetraOwnedMap.Print (std::cout);
    out << "Epetra overlap map:" << endl;
    epetraOverlapMap.Print (std::cout);
  }

  // Create the two Tpetra Maps.  The "invalid" global element count
  // input tells Tpetra::Map to compute the global number of elements
  // itself.
  const int invalid = Teuchos::OrdinalTraits<int>::invalid();
  RCP<Tpetra::Map<int> > tpetraOwnedMap =
    rcp (new Tpetra::Map<int> (invalid, ArrayView<int> (owned, ownedSize),
                               0, teuchosComm));
  tpetraOwnedMap->setObjectLabel ("Owned Map");
  RCP<Tpetra::Map<int> > tpetraOverlapMap =
    rcp (new Tpetra::Map<int> (invalid, ArrayView<int> (overlap, overlapSize),
                               0, teuchosComm));
  tpetraOverlapMap->setObjectLabel ("Overlap Map");

  // In verbose mode, have the Tpetra::Map objects describe themselves.
  if (verbose && printTpetra) {
    Teuchos::EVerbosityLevel verb = Teuchos::VERB_EXTREME;

    // Tpetra::Map::describe() takes a FancyOStream, but expects all
    // MPI processes to be able to write to it.  (The method handles
    // its own synchronization.)
    RCP<FancyOStream> globalOut = getFancyOStream (rcpFromRef (std::cout));
    out << "Tpetra owned map:" << endl;
    {
      OSTab tab (globalOut);
      tpetraOwnedMap->describe (*globalOut, verb);
    }
    out << "Tpetra overlap map:" << endl;
    {
      OSTab tab (globalOut);
      tpetraOverlapMap->describe (*globalOut, verb);
    }
  }

  // Use the owned and overlap maps to construct an importer for both
  // Epetra and Tpetra.
  Epetra_Import       epetraImporter (epetraOverlapMap, epetraOwnedMap  );
  Tpetra::Import<int> tpetraImporter (tpetraOwnedMap  , tpetraOverlapMap);

  // In verbose mode, have the Epetra_Import object describe itself.
  if (verbose && printEpetra) {
    out << "Epetra importer:" << endl;
    // The importer's Print() method takes an std::ostream& and plans
    // to write to it on all MPI processes (handling synchronization
    // itself).
    epetraImporter.Print (std::cout);
    out << endl;
  }

  // In verbose mode, have the Tpetra::Import object describe itself.
  if (verbose && printTpetra) {
    out << "Tpetra importer:" << endl;
    // The importer doesn't implement Teuchos::Describable.  It wants
    // std::cout and plans to write to it on all MPI processes (with
    // its own synchronization).
    tpetraImporter.print (std::cout);
    out << endl;
  }

  // Construct owned and overlap vectors for both Epetra and Tpetra.
  Epetra_Vector epetraOwnedVector   (epetraOwnedMap  );
  Epetra_Vector epetraOverlapVector (epetraOverlapMap);
  Tpetra::Vector<double,int> tpetraOwnedVector   (tpetraOwnedMap  );
  Tpetra::Vector<double,int> tpetraOverlapVector (tpetraOverlapMap);

  // The test is as follows: initialize the owned and overlap vectors
  // with global IDs in the owned regions.  Initialize the overlap
  // vectors to equal -1 in the overlap regions.  Then perform a
  // communication from the owned vectors to the overlap vectors.  The
  // resulting overlap vectors should have global IDs everywhere and
  // all of the -1 values should be overwritten.

  // Initialize.  We cannot assign directly to the Tpetra Vectors;
  // instead, we extract nonconst views and assign to those.  The
  // results aren't guaranteed to be committed to the vector unless
  // the views are released (by assigning Teuchos::null to them).
  epetraOverlapVector.PutScalar(-1);
  tpetraOverlapVector.putScalar(-1);
  ArrayRCP<double> tpetraOwnedArray   = tpetraOwnedVector.getDataNonConst(0);
  ArrayRCP<double> tpetraOverlapArray = tpetraOverlapVector.getDataNonConst(0);
  for (int owned_lid = 0;
       owned_lid < tpetraOwnedMap->getLocalElementList().size();
       ++owned_lid) {
    int gid         = tpetraOwnedMap->getGlobalElement(owned_lid);
    int overlap_lid = tpetraOverlapMap->getLocalElement(gid);
    epetraOwnedVector[owned_lid]     = gid;
    epetraOverlapVector[overlap_lid] = gid;
    tpetraOwnedArray[owned_lid]      = gid;
    tpetraOverlapArray[overlap_lid]  = gid;
  }
  // Make sure that the changes to the Tpetra Vector were committed,
  // by releasing the nonconst views.
  tpetraOwnedArray = Teuchos::null;
  tpetraOverlapArray = Teuchos::null;

  // Test the Epetra and Tpetra Import.
  if (verbose) {
    out << "Testing Import from owned Map to overlap Map:" << endl << endl;
  }
  epetraOverlapVector.Import(  epetraOwnedVector, epetraImporter, Insert);
  tpetraOverlapVector.doImport(tpetraOwnedVector, tpetraImporter,
                               Tpetra::INSERT);
  // Check the Import results.
  success = countFailures (teuchosComm, epetraOwnedMap, epetraOwnedVector,
                           epetraOverlapMap, epetraOverlapVector,
                           tpetraOwnedMap, tpetraOwnedVector,
                           tpetraOverlapMap, tpetraOverlapVector, verbose);

  const bool testOtherDirections = false;
  if (testOtherDirections) {
    //
    // Reinitialize the Tpetra vectors and test whether Export works.
    //
    tpetraOverlapVector.putScalar(-1);
    tpetraOwnedArray   = tpetraOwnedVector.getDataNonConst(0);
    tpetraOverlapArray = tpetraOverlapVector.getDataNonConst(0);
    for (int owned_lid = 0;
         owned_lid < tpetraOwnedMap->getLocalElementList().size();
         ++owned_lid)
      {
        int gid         = tpetraOwnedMap->getGlobalElement(owned_lid);
        int overlap_lid = tpetraOverlapMap->getLocalElement(gid);
        tpetraOwnedArray[owned_lid]      = gid;
        tpetraOverlapArray[overlap_lid]  = gid;
      }
    // Make sure that the changes to the Tpetra Vector were committed,
    // by releasing the nonconst views.
    tpetraOwnedArray = Teuchos::null;
    tpetraOverlapArray = Teuchos::null;

    // Make a Tpetra Export object, and test the export.
    Tpetra::Export<int> tpetraExporter1 (tpetraOwnedMap, tpetraOverlapMap);
    if (verbose) {
      out << "Testing Export from owned Map to overlap Map:" << endl << endl;
    }
    tpetraOverlapVector.doExport (tpetraOwnedVector, tpetraExporter1,
                                  Tpetra::INSERT);

    // Check the Export results.
    success = countFailures (teuchosComm, epetraOwnedMap, epetraOwnedVector,
                             epetraOverlapMap, epetraOverlapVector,
                             tpetraOwnedMap, tpetraOwnedVector,
                             tpetraOverlapMap, tpetraOverlapVector, verbose);
    //
    // Reinitialize the Tpetra vectors and see what Import in the
    // other direction does.
    //
    tpetraOverlapVector.putScalar(-1);
    tpetraOwnedArray   = tpetraOwnedVector.getDataNonConst(0);
    tpetraOverlapArray = tpetraOverlapVector.getDataNonConst(0);
    for (int owned_lid = 0;
         owned_lid < tpetraOwnedMap->getLocalElementList().size();
         ++owned_lid)
      {
        int gid         = tpetraOwnedMap->getGlobalElement(owned_lid);
        int overlap_lid = tpetraOverlapMap->getLocalElement(gid);
        tpetraOwnedArray[owned_lid]      = gid;
        tpetraOverlapArray[overlap_lid]  = gid;
      }
    // Make sure that the changes to the Tpetra Vector were committed,
    // by releasing the nonconst views.
    tpetraOwnedArray = Teuchos::null;
    tpetraOverlapArray = Teuchos::null;

    if (verbose) {
      out << "Testing Import from overlap Map to owned Map:" << endl << endl;
    }
    Tpetra::Import<int> tpetraImporter2 (tpetraOverlapMap, tpetraOwnedMap);
    tpetraOwnedVector.doImport (tpetraOverlapVector, tpetraImporter2,
                                Tpetra::INSERT);
    // Check the Import results.
    success = countFailures (teuchosComm, epetraOwnedMap, epetraOwnedVector,
                             epetraOverlapMap, epetraOverlapVector,
                             tpetraOwnedMap, tpetraOwnedVector,
                             tpetraOverlapMap, tpetraOverlapVector, verbose);
  } // if testOtherDirections

  return success;
}

int
main (int argc, char *argv[])
{
  MPI_Init (&argc, &argv);
  int myRank = 0;
  MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
  const bool success = testMain (argc, argv);
  MPI_Finalize ();

  if (myRank == 0) {
    std::cout << "End Result: TEST " << (success ? "PASSED" : "FAILED") << std::endl;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}






