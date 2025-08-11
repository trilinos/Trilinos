// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_ConfigDefs.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_OrdinalTraits.hpp>
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

using ST = Tpetra::Vector<>::scalar_type;
using LO = Tpetra::Vector<>::local_ordinal_type;
using GO = int; // So Epetra and Tpetra can use the same GID lists

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
                       const int numMapCreateTrials,
                       const int numImportCreateTrials,
                       const int numVectorCreateTrials,
                       const int numImportExecTrials)
{
  using import_type = Tpetra::Import<LO, GO>;
  using map_type = Tpetra::Map<LO, GO>;
  using vector_type = Tpetra::Vector<ST, LO, GO>;
  const global_size_t INVALID =
    Teuchos::OrdinalTraits<global_size_t>::invalid ();

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
      srcMap = rcp (new map_type (INVALID, srcGlobalElts, indexBase, comm));
      destMap = rcp (new map_type (INVALID, destGlobalElts, indexBase, comm));
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

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    RCP<const Teuchos::Comm<int> > tpetraComm;
    tpetraComm = Tpetra::getDefaultComm ();

    const int numProcs = tpetraComm->getSize ();
    const int myRank = tpetraComm->getRank ();
    const int indexBase = 1; // of interest to Sierra

    // Benchmark parameters
    int numEltsPerProc = 100000;
    int numTrials = 100;
    bool runEpetra = false;
    bool runTpetra = true;

    CommandLineProcessor cmdp;
    cmdp.setOption ("numEltsPerProc", &numEltsPerProc,
                    "Number of global indices "
                    "owned by each process");
    cmdp.setOption ("numTrials", &numTrials,
                    "Number of times to repeat each "
                    "operation in a timing loop");
    cmdp.setOption ("runEpetra", "noEpetra", &runEpetra,
                    "Whether to run the Epetra benchmark");
    cmdp.setOption ("runTpetra", "noTpetra", &runTpetra,
                    "Whether to run the Tpetra benchmark");
    const CommandLineProcessor::EParseCommandLineReturn parseResult =
      cmdp.parse (argc, argv);
    if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED) {
      // The user specified --help at the command line to print help
      // with command-line arguments.  We printed help already, so
      // quit with a happy return code.
      return EXIT_SUCCESS;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION
        (parseResult != CommandLineProcessor::PARSE_SUCCESSFUL,
         std::invalid_argument,
         "Failed to parse command-line arguments.");
      TEUCHOS_TEST_FOR_EXCEPTION
        (numTrials < 0, std::invalid_argument,
         "numTrials must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION
        (numEltsPerProc < 0, std::invalid_argument,
         "numEltsPerProc must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION
        (runEpetra, std::invalid_argument, "Tpetra was not built with "
         "Epetra enable, so you cannot run the Epetra benchmark." );
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
    if (runTpetra) {
      benchmarkTpetraImport (srcGlobalElts, destGlobalElts, indexBase, tpetraComm,
                             numMapCreateTrials, numImportCreateTrials,
                             numVectorCreateTrials, numImportExecTrials);
    }
    TimeMonitor::report (tpetraComm.ptr (), std::cout);
  }
  return EXIT_SUCCESS;
}
