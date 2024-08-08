//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#include <KokkosBlas.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Comm.hpp>
#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#else
#include <Teuchos_DefaultSerialComm.hpp>
#endif  // HAVE_MPI

using Teuchos::Comm;
using Teuchos::CommandLineProcessor;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Time;
using Teuchos::TimeMonitor;

// Create a new timer with the given name if it hasn't already been
// created, else get the previously created timer with that name.
RCP<Time> getTimer(const std::string& timerName) {
  RCP<Time> timer = TimeMonitor::lookupCounter(timerName);
  if (timer.is_null()) {
    timer = TimeMonitor::getNewCounter(timerName);
  }
  return timer;
}

bool benchmarkKokkos(std::ostream& out, const int lclNumRows, const int numTrials) {
  using std::endl;
  typedef Kokkos::View<double*, Kokkos::LayoutLeft> vector_type;

  RCP<Time> vecCreateTimer = getTimer("Kokkos: Vector: Create");
  RCP<Time> vecFillTimer   = getTimer("Kokkos: Vector: Fill");
  RCP<Time> vecDotTimer    = getTimer("Kokkos: Vector: Dot");

  // Benchmark creation of a Vector.
  vector_type x;
  {
    TimeMonitor timeMon(*vecCreateTimer);
    // This benchmarks both vector creation and vector destruction.
    for (int k = 0; k < numTrials; ++k) {
      x = vector_type("x", lclNumRows);
    }
  }

  // Benchmark filling a Vector.
  {
    TimeMonitor timeMon(*vecFillTimer);
    for (int k = 0; k < numTrials; ++k) {
      Kokkos::deep_copy(x, 1.0);
    }
  }

  vector_type y("y", lclNumRows);
  Kokkos::deep_copy(y, -1.0);

  // Benchmark computing the dot product of two Vectors.
  double dotResults[2];
  dotResults[0] = 0.0;
  dotResults[1] = 0.0;
  {
    TimeMonitor timeMon(*vecDotTimer);
    for (int k = 0; k < numTrials; ++k) {
      // "Confuse" the compiler so it doesn't optimize away the dot() calls.
      dotResults[k % 2] = KokkosBlas::dot(x, y);
    }
  }

  if (numTrials > 0) {
    const double expectedResult = static_cast<double>(lclNumRows) * -1.0;
    if (dotResults[0] != expectedResult) {
      out << "Kokkos dot product result is wrong!  Expected " << expectedResult << " but got " << dotResults[0]
          << " instead." << endl;
      return false;
    } else {
      return true;
    }
  }
  return true;
}

bool benchmarkRaw(std::ostream& out, const int lclNumRows, const int numTrials) {
  using std::endl;
  RCP<Time> vecCreateTimer = getTimer("Raw: Vector: Create");
  RCP<Time> vecFillTimer   = getTimer("Raw: Vector: Fill");
  RCP<Time> vecDotTimer    = getTimer("Raw: Vector: Dot");

  // Benchmark creation of a Vector.
  double* x = 0;
  {
    TimeMonitor timeMon(*vecCreateTimer);
    // This benchmarks both vector creation and vector destruction.
    for (int k = 0; k < numTrials; ++k) {
      x = new double[lclNumRows];
      memset(x, 0, lclNumRows * sizeof(double));
      if (k + 1 < numTrials) {
        delete[] x;
      }
    }
  }

  // Benchmark filling a Vector.
  {
    TimeMonitor timeMon(*vecFillTimer);
    for (int k = 0; k < numTrials; ++k) {
      for (int i = 0; i < lclNumRows; ++i) {
        x[i] = 1.0;
      }
    }
  }

  double* y = new double[lclNumRows];
  for (int i = 0; i < lclNumRows; ++i) {
    y[i] = -1.0;
  }

  // Benchmark computing the dot product of two Vectors.
  double dotResults[2];
  dotResults[0] = 0.0;
  dotResults[1] = 0.0;
  {
    TimeMonitor timeMon(*vecDotTimer);
    for (int k = 0; k < numTrials; ++k) {
      double sum = 0.0;
      for (int i = 0; i < lclNumRows; ++i) {
        sum += x[i] * y[i];
      }
      // "Confuse" the compiler so it doesn't optimize away the loops.
      dotResults[k % 2] = sum;
    }
  }

  if (x != NULL) {
    delete[] x;
    x = NULL;
  }
  if (y != NULL) {
    delete[] y;
    y = NULL;
  }

  if (numTrials == 0) {
    return true;  // trivially
  } else {        // numTrials > 0
    const double expectedResult = static_cast<double>(lclNumRows) * -1.0;
    if (dotResults[0] != expectedResult) {
      out << "Raw dot product result is wrong!  Expected " << expectedResult << " but got " << dotResults[0]
          << " instead." << endl;
      return false;
    } else {
      return true;
    }
  }
}

int main(int argc, char* argv[]) {
  using std::cout;
  using std::endl;
  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackHole);
  Kokkos::initialize(argc, argv);

#ifdef HAVE_MPI
  RCP<const Comm<int> > comm = rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
#else
  RCP<const Comm<int> > comm = rcp(new Teuchos::SerialComm<int>());
#endif  // HAVE_MPI

  // const int numProcs = comm->getSize (); // unused
  const int myRank = comm->getRank();

  // Benchmark parameters
  int lclNumRows = 100000;
  int numTrials  = 1000;

  bool runKokkos = true;
  bool runRaw    = true;

  CommandLineProcessor cmdp;
  cmdp.setOption("lclNumRows", &lclNumRows,
                 "Number of global indices "
                 "owned by each process");
  cmdp.setOption("numTrials", &numTrials, "Number of timing loop iterations for each event to time");
  cmdp.setOption("runKokkos", "noKokkos", &runKokkos, "Whether to run the Kokkos benchmark");
  cmdp.setOption("runRaw", "noRaw", &runRaw, "Whether to run the raw benchmark");
  const CommandLineProcessor::EParseCommandLineReturn parseResult = cmdp.parse(argc, argv);
  if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED) {
    // The user specified --help at the command line to print help
    // with command-line arguments.  We printed help already, so quit
    // with a happy return code.
    return EXIT_SUCCESS;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL, std::invalid_argument,
                               "Failed to parse command-line arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(lclNumRows < 0, std::invalid_argument, "lclNumRows must be nonnegative.");
  }

  if (myRank == 0) {
    cout << endl
         << "---" << endl
         << "Command-line options:" << endl
         << "  lclNumRows: " << lclNumRows << endl
         << "  numTrials: " << numTrials << endl
         << "  runKokkos: " << (runKokkos ? "true" : "false") << endl
         << "  runRaw: " << (runRaw ? "true" : "false") << endl
         << endl;
  }

  // Run the benchmark
  bool success = true;
  if (runKokkos) {
    const bool lclSuccess = benchmarkKokkos(cout, lclNumRows, numTrials);
    success               = success && lclSuccess;
  }
  if (runRaw) {
    const bool lclSuccess = benchmarkRaw(cout, lclNumRows, numTrials);
    success               = success && lclSuccess;
  }

  TimeMonitor::report(comm.ptr(), cout);
  if (success) {
    cout << "End Result: TEST PASSED" << endl;
  } else {
    cout << "End Result: TEST FAILED" << endl;
  }
  Kokkos::finalize();
  return EXIT_SUCCESS;
}
