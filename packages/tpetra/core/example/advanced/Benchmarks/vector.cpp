// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This benchmark compares performance of common operations by
// Epetra_Vector and Teptra::Vector.  Both Epetra and Tpetra implement
// sparse matrix and dense vector data structures and computational
// kernels for users and other Trilinos data structures.  Both
// packages use MPI (Message Passing Interface) for distributed-memory
// parallelism.  Tpetra additionally uses Kokkos for shared-memory
// parallelism within an MPI process.

#include <Tpetra_ConfigDefs.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_TimeMonitor.hpp>

typedef Tpetra::global_size_t GST;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::CommandLineProcessor;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Time;
using Teuchos::TimeMonitor;

typedef Tpetra::Vector<>::scalar_type ST;
typedef Tpetra::Map<>::local_ordinal_type LO;
typedef Tpetra::Map<>::global_ordinal_type GO;
typedef Tpetra::Vector<>::node_type NT;

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
benchmarkTpetra (RCP<const Comm<int> > comm,
                 const int numIndPerProc,
                 const int indexBase,
                 const int numMapCreateTrials,
                 const int numVecCreateTrials,
                 const int numVecRandTrials,
                 const int numVecNormTrials,
                 const int numVecDotTrials)
{
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::Vector<ST, LO, GO, NT> vector_type;
  typedef typename vector_type::dot_type dot_type;
  typedef typename vector_type::mag_type mag_type;
  typedef Teuchos::ScalarTraits<ST> STS;
  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

  RCP<Time> mapCreateTimer = getTimer ("Tpetra: Map: Create");
  RCP<Time> vecCreateTimer = getTimer ("Tpetra: Vector: Create");
  RCP<Time> vecRandTimer = getTimer ("Tpetra: Vector: Rand");
  RCP<Time> vecNormTimer = getTimer ("Tpetra: Vector: Norm");
  RCP<Time> vecDotTimer = getTimer ("Tpetra: Vector: Dot");
  RCP<Time> vecAxpyTimer = getTimer ("Tpetra: Vector: Axpy");

  // Benchmark creation of a Map with a given number of indices per
  // process, telling the Map to compute the global number of indices.
  RCP<map_type> map;
  {
    TimeMonitor timeMon (*mapCreateTimer);
    for (int k = 0; k < numMapCreateTrials; ++k) {
      map = rcp (new map_type (INVALID, static_cast<size_t> (numIndPerProc),
                               static_cast<GO> (indexBase), comm));
    }
  }
  if (map.is_null ()) { // no Map create trials means no Map
    return;
  }

  // Benchmark creation of a Vector using the above Map.
  RCP<vector_type> x;
  {
    TimeMonitor timeMon (*vecCreateTimer);
    // This benchmarks both vector creation and vector destruction.
    for (int k = 0; k < numVecCreateTrials; ++k) {
      x = rcp (new vector_type (map));
    }
  }
  if (x.is_null ()) { // no Vector create trials means no Vector
    return;
  }

  // Benchmark filling a Vector with random data.
  {
    TimeMonitor timeMon (*vecRandTimer);
    for (int k = 0; k < numVecRandTrials; ++k) {
      x->randomize ();
    }
  }

  // Benchmark computing the 2-norm of a Vector.
  mag_type normResults[2];
  normResults[0] = 0.0;
  normResults[1] = 0.0;
  {
    TimeMonitor timeMon (*vecNormTimer);
    for (int k = 0; k < numVecNormTrials; ++k) {
      // "Confuse" the compiler so it doesn't optimize away the norm2() calls.
      normResults[k % 2] = x->norm2 ();
    }
  }

  // Benchmark computing the dot product of two Vectors.
  RCP<vector_type> y = rcp (new vector_type (map));
  y->randomize ();
  dot_type dotResults[2];
  dotResults[0] = 0.0;
  dotResults[1] = 0.0;
  {
    TimeMonitor timeMon (*vecDotTimer);
    for (int k = 0; k < numVecDotTrials; ++k) {
      // "Confuse" the compiler so it doesn't optimize away the dot() calls.
      dotResults[k % 2] = x->dot (*y);
    }
  }

  // Benchmark axpy (3-argument update).
  const ST HALF = STS::one () / (STS::one () + STS::one ());
  {
    TimeMonitor timeMon (*vecAxpyTimer);
    for (int k = 0; k < numVecDotTrials; ++k) {
      y->update (HALF, *x, STS::zero ());
    }
  }

  // Trick the compiler into not complaining that normResults and
  // dotResults never get used, by printing them to the equivalent of
  // /dev/null.
  Teuchos::oblackholestream blackHole;
  blackHole << "Norm results: " << normResults[0] << "," << normResults[1]
            << std::endl
            << "Dot results:  " << dotResults[0] << "," << dotResults[1]
            << std::endl;
}


int
main (int argc, char* argv[])
{
  using std::cout;
  using std::endl;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    RCP<const Teuchos::Comm<int> > tpetraComm;
    tpetraComm = Tpetra::getDefaultComm ();

    //const int numProcs = tpetraComm->getSize (); // unused
    const int myRank = tpetraComm->getRank ();
    const int indexBase = 0;

    // Benchmark parameters
    int numIndsPerProc = 100000;
    int numTrials = 1000;

    bool runEpetra = false;
    bool runTpetra = true;

    CommandLineProcessor cmdp;
    cmdp.setOption ("numIndsPerProc", &numIndsPerProc, "Number of global indices "
                    "owned by each process");
    cmdp.setOption ("numTrials", &numTrials, "Number of timing loop iterations for each event to time");
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
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION
        (parseResult != CommandLineProcessor::PARSE_SUCCESSFUL,
         std::invalid_argument, "Failed to parse command-line arguments.");
      TEUCHOS_TEST_FOR_EXCEPTION
        (numIndsPerProc < 0, std::invalid_argument,
         "numIndsPerProc must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION
        (runEpetra, std::invalid_argument, "Tpetra was not built with Epetra "
         "enabled, so you cannot run the Epetra benchmark." );
    }

    if (myRank == 0) {
      cout << endl << "---" << endl
           << "Command-line options:" << endl
           << "  numIndsPerProc: " << numIndsPerProc << endl
           << "  numTrials: " << numTrials << endl;
      cout << "  runEpetra: " << (runEpetra ? "true" : "false") << endl;
      cout << "  runTpetra: " << (runTpetra ? "true" : "false") << endl
           << endl;
    }

    const int numMapCreateTrials = numTrials;
    const int numVecCreateTrials = numTrials;
    const int numVecRandTrials = numTrials;
    const int numVecNormTrials = numTrials;
    const int numVecDotTrials = numTrials;

    // Run the benchmark
    if (runTpetra) {
      benchmarkTpetra (tpetraComm, numIndsPerProc, indexBase,
                       numMapCreateTrials,
                       numVecCreateTrials,
                       numVecRandTrials,
                       numVecNormTrials,
                       numVecDotTrials);
    }
    TimeMonitor::report (tpetraComm.ptr (), std::cout);
  }
  return EXIT_SUCCESS;
}
