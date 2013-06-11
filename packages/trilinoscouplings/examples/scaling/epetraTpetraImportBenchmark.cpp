// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

// TrilinosCouplings includes
#include <TrilinosCouplings_config.h>

#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#ifdef EPETRA_MPI
#  include <Epetra_MpiComm.h>
#  include <mpi.h>
#else
#  include <Epetra_SerialComm.h>
#endif

#include <Tpetra_Import.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Tpetra::global_size_t;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::CommandLineProcessor;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Time;
using Teuchos::TimeMonitor;

typedef double ST;
typedef int LO;
typedef int GO; // So that Epetra and Tpetra can use the same GID lists
typedef Kokkos::SerialNode NT;

// Create a new timer with the given name if it hasn't already been
// created, else get the previously created timer with that name.
RCP<Time> getTimer (const std::string& timerName) {
  RCP<Time> timer =
    TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }
  return timer;
}

void
benchmarkTpetraImport (ArrayView<const GO> srcGlobalElts,
                       ArrayView<const GO> destGlobalElts,
                       const GO indexBase,
                       RCP<const Comm<int> > comm,
                       RCP<NT> node,
                       const int numMapCreateTrials,
                       const int numImportCreateTrials,
                       const int numVectorCreateTrials,
                       const int numImportExecTrials)
{
  typedef Tpetra::Import<LO, GO, NT> import_type;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::Vector<ST, LO, GO, NT> vector_type;
  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid ();

  TEUCHOS_TEST_FOR_EXCEPTION(
    numMapCreateTrials < 1 && numImportCreateTrials > 0, std::invalid_argument,
    "numMapCreateTrials must be > 0 if numImportCreateTrials > 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numImportCreateTrials < 1 && numImportExecTrials > 0, std::invalid_argument,
    "numImportCreateTrials must be > 0 if numImportExecTrials > 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numVectorCreateTrials < 1 && numMapCreateTrials > 0, std::invalid_argument,
    "numVectorCreateTrials must be > 0 if numMapCreateTrials > 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numVectorCreateTrials < 1 && numImportExecTrials > 0, std::invalid_argument,
    "numVectorCreateTrials must be > 0 if numImportExecTrials > 0.");

  RCP<Time> mapCreateTimer = getTimer ("Tpetra: Map: Create");
  RCP<Time> importCreateTimer = getTimer ("Tpetra: Import: Create");
  RCP<Time> vectorCreateTimer = getTimer ("Tpetra: Vector: Create");
  RCP<Time> importExecTimer = getTimer ("Tpetra: Import: Execute");

  RCP<map_type> srcMap, destMap;
  {
    TimeMonitor timeMon (*mapCreateTimer);
    for (int k = 0; k < numMapCreateTrials; ++k) {
      srcMap = rcp (new map_type (INVALID, srcGlobalElts, indexBase, comm, node));
      destMap = rcp (new map_type (INVALID, destGlobalElts, indexBase, comm, node));
    }
  }
  RCP<import_type> import;
  {
    TimeMonitor timeMon (*importCreateTimer);
    for (int k = 0; k < numImportCreateTrials; ++k) {
      import = rcp (new import_type (srcMap, destMap));
    }
  }
  RCP<vector_type> srcVec, destVec;
  {
    TimeMonitor timeMon (*vectorCreateTimer);
    // This benchmarks both vector creation and vector destruction.
    for (int k = 0; k < numVectorCreateTrials; ++k) {
      srcVec = rcp (new vector_type (srcMap));
      destVec = rcp (new vector_type (destMap));
    }
  }
  {
    TimeMonitor timeMon (*importExecTimer);
    for (int k = 0; k < numImportExecTrials; ++k) {
      destVec->doImport (*srcVec, *import, Tpetra::ADD);
    }
  }
}

void
benchmarkEpetraImport (ArrayView<const int> srcGlobalElts,
                       ArrayView<const int> destGlobalElts,
                       const int indexBase,
                       const Epetra_Comm& comm,
                       const int numMapCreateTrials,
                       const int numImportCreateTrials,
                       const int numVectorCreateTrials,
                       const int numImportExecTrials)
{
  const int INVALID = -1;

  TEUCHOS_TEST_FOR_EXCEPTION(
    numMapCreateTrials < 1 && numImportCreateTrials > 0, std::invalid_argument,
    "numMapCreateTrials must be > 0 if numImportCreateTrials > 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numImportCreateTrials < 1 && numImportExecTrials > 0, std::invalid_argument,
    "numImportCreateTrials must be > 0 if numImportExecTrials > 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numVectorCreateTrials < 1 && numMapCreateTrials > 0, std::invalid_argument,
    "numVectorCreateTrials must be > 0 if numMapCreateTrials > 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numVectorCreateTrials < 1 && numImportExecTrials > 0, std::invalid_argument,
    "numVectorCreateTrials must be > 0 if numImportExecTrials > 0.");

  RCP<Time> mapCreateTimer = getTimer ("Epetra: Map: Create");
  RCP<Time> importCreateTimer = getTimer ("Epetra: Import: Create");
  RCP<Time> vectorCreateTimer = getTimer ("Epetra: Vector: Create");
  RCP<Time> importExecTimer = getTimer ("Epetra: Import: Execute");

  RCP<Epetra_Map> srcMap, destMap;
  {
    TimeMonitor timeMon (*mapCreateTimer);
    for (int k = 0; k < numMapCreateTrials; ++k) {
      const int srcNumElts = as<int> (srcGlobalElts.size ());
      const int* const srcElts = srcGlobalElts.getRawPtr ();
      srcMap = rcp (new Epetra_Map (INVALID, srcNumElts, srcElts, indexBase, comm));

      const int destNumElts = as<int> (destGlobalElts.size ());
      const int* const destElts = destGlobalElts.getRawPtr ();
      destMap = rcp (new Epetra_Map (INVALID, destNumElts, destElts, indexBase, comm));
    }
  }
  RCP<Epetra_Import> import;
  {
    TimeMonitor timeMon (*importCreateTimer);
    for (int k = 0; k < numImportCreateTrials; ++k) {
      import = rcp (new Epetra_Import (*srcMap, *destMap));
    }
  }
  RCP<Epetra_Vector> srcVec, destVec;
  {
    TimeMonitor timeMon (*vectorCreateTimer);
    // This benchmarks both vector creation and vector destruction.
    for (int k = 0; k < numVectorCreateTrials; ++k) {
      srcVec = rcp (new Epetra_Vector (*srcMap));
      destVec = rcp (new Epetra_Vector (*destMap));
    }
  }
  {
    TimeMonitor timeMon (*importExecTimer);
    for (int k = 0; k < numImportExecTrials; ++k) {
      (void) destVec->Import (*srcVec, *import, Add);
    }
  }
}

void
createGidLists (Array<GO>& srcGlobalElts,
                Array<GO>& destGlobalElts,
                const int numProcs,
                const int myRank,
                const int numEltsPerProc,
                const GO indexBase)
{
  const int overlap = (myRank == 0 || myRank == numProcs-1) ? 1 : 2;

  srcGlobalElts.resize (numEltsPerProc);
  destGlobalElts.resize (numEltsPerProc + overlap);

  const GO myStartSrcGid = indexBase + myRank * numEltsPerProc;
  const GO myEndSrcGid = myStartSrcGid + numEltsPerProc - 1; // inclusive

  // Put the GIDs in reverse order, to simulate a noncontiguously
  // ordered Map.  Neither Epetra nor Tpetra optimize for this case.
  for (int k = 0; k < numEltsPerProc; ++k) {
    const GO curGid = myStartSrcGid + as<GO> (k);
    TEUCHOS_TEST_FOR_EXCEPTION(
      curGid < indexBase, std::logic_error,
      "curGid = " << curGid << " < indexBase = " << indexBase << ".");
    srcGlobalElts[numEltsPerProc - k - 1] = curGid;
    destGlobalElts[numEltsPerProc - k - 1] = curGid;
  }

  // 1-D Poisson overlap.
  if (myRank == 0) {
    const GO rightGid = myEndSrcGid + 1;
    TEUCHOS_TEST_FOR_EXCEPTION(
      rightGid < indexBase, std::logic_error,
      "At left proc, right boundary GID = " << rightGid
      << " < indexBase = " << indexBase << ".");
    destGlobalElts[numEltsPerProc] = rightGid;
  } else if (myRank == numProcs-1) {
    const GO leftGid = myStartSrcGid - 1;
    TEUCHOS_TEST_FOR_EXCEPTION(
      leftGid < indexBase, std::logic_error,
      "At right proc, left boundary GID = " << leftGid
      << " < indexBase = " << indexBase << ".");
    destGlobalElts[numEltsPerProc] = leftGid;
  } else {
    const GO leftGid = myStartSrcGid - 1;
    TEUCHOS_TEST_FOR_EXCEPTION(
      leftGid < indexBase, std::logic_error,
      "At middle proc, left boundary GID = " << leftGid
      << " < indexBase = " << indexBase << ".");
    destGlobalElts[numEltsPerProc] = leftGid;

    const GO rightGid = myEndSrcGid + 1;
    TEUCHOS_TEST_FOR_EXCEPTION(
      rightGid < indexBase, std::logic_error,
      "At middle proc, right boundary GID = " << rightGid
      << " < indexBase = " << indexBase << ".");
    destGlobalElts[numEltsPerProc+1] = rightGid;
  }
}


int main (int argc, char* argv[]) {
  using std::cout;
  using std::endl;
  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);

  RCP<const Teuchos::Comm<int> > tpetraComm;
#ifdef EPETRA_MPI
  Epetra_MpiComm epetraComm (MPI_COMM_WORLD);
  tpetraComm = rcp (new Teuchos::MpiComm<int> (MPI_COMM_WORLD));
#else
  Epetra_SerialComm epetraComm;
  tpetraComm = rcp (new Teuchos::SerialComm<int>);
#endif // EPETRA_MPI
  RCP<NT> node = Kokkos::Details::getNode<NT> ();

  const int numProcs = tpetraComm->getSize ();
  const int myRank = tpetraComm->getRank ();
  const int indexBase = 1; // of interest to Sierra

  // Benchmark parameters
  int numEltsPerProc = 100000;
  int numTrials = 100;
  bool runEpetra = true;
  bool runTpetra = true;

  CommandLineProcessor cmdp;
  cmdp.setOption ("numEltsPerProc", &numEltsPerProc, "Number of global indices "
                  "owned by each process");
  cmdp.setOption ("numTrials", &numTrials, "Number of times to repeat each "
                  "operation in a timing loop");
  cmdp.setOption ("runEpetra", "noEpetra", &runEpetra,
                  "Whether to run the Epetra benchmark");
  cmdp.setOption ("runTpetra", "noTpetra", &runTpetra,
                  "Whether to run the Tpetra benchmark");
  const CommandLineProcessor::EParseCommandLineReturn parseResult =
    cmdp.parse (argc, argv);
  if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED) {
    // The user specified --help at the command line to print help
    // with command-line arguments.  We printed help already, so quit
    // with a happy return code.
    return EXIT_SUCCESS;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      parseResult != CommandLineProcessor::PARSE_SUCCESSFUL,
      std::invalid_argument, "Failed to parse command-line arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      numTrials < 0, std::invalid_argument, "numTrials must be nonnegative.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      numEltsPerProc < 0, std::invalid_argument,
      "numEltsPerProc must be nonnegative.");
  }

  // Derived benchmark parameters
  const int numMapCreateTrials = numTrials;
  const int numImportCreateTrials = numTrials;
  const int numVectorCreateTrials = numTrials;
  const int numImportExecTrials = numTrials;

  if (myRank == 0) {
    cout << endl << "---" << endl
         << "Command-line options:" << endl
         << "  numEltsPerProc: " << numEltsPerProc << endl
         << "  numTrials: " << numTrials << endl
         << "  runEpetra: " << (runEpetra ? "true" : "false") << endl
         << "  runTpetra: " << (runTpetra ? "true" : "false") << endl
         << endl;
  }

  // Run the benchmark
  Array<GO> srcGlobalElts, destGlobalElts;
  createGidLists (srcGlobalElts, destGlobalElts, numProcs, myRank,
                  numEltsPerProc, indexBase);
  if (runEpetra) {
    benchmarkEpetraImport (srcGlobalElts, destGlobalElts, indexBase, epetraComm,
                           numMapCreateTrials, numImportCreateTrials,
                           numVectorCreateTrials, numImportExecTrials);
  }
  if (runTpetra) {
    benchmarkTpetraImport (srcGlobalElts, destGlobalElts, indexBase, tpetraComm,
                           node, numMapCreateTrials, numImportCreateTrials,
                           numVectorCreateTrials, numImportExecTrials);
  }
  TimeMonitor::report (tpetraComm.ptr (), std::cout);
  return EXIT_SUCCESS;
}
