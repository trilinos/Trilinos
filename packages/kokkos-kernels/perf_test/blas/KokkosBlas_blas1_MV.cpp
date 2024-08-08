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

bool benchmarkKokkos(std::ostream& out, const int numRows, const int numCols, const int numTrials) {
  using Kokkos::ALL;
  using Kokkos::subview;
  using std::endl;
#ifdef KOKKOS_ENABLE_SERIAL
  typedef Kokkos::Serial execution_space;
#else
  typedef Kokkos::View<double**, Kokkos::LayoutLeft>::execution_space execution_space;
#endif  // KOKKOS_ENABLE_SERIAL
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, execution_space> mv_type;
  bool success = true;

  RCP<Time> vecCreateTimer      = getTimer("Kokkos: MV: Create");
  RCP<Time> vecFillZeroTimer    = getTimer("Kokkos: MV: Fill zero");
  RCP<Time> vecFillNonZeroTimer = getTimer("Kokkos: MV: Fill nonzero");
  RCP<Time> vecNrm2Timer        = getTimer("Kokkos: MV: Nrm2 (contiguous)");
  RCP<Time> vecNrm2Timer2       = getTimer("Kokkos: MV: Nrm2 (noncontiguous)");
  RCP<Time> vecNrm1Timer        = getTimer("Kokkos: MV: Nrm1 (contiguous)");
  RCP<Time> vecNrm1Timer2       = getTimer("Kokkos: MV: Nrm1 (noncontiguous)");
  RCP<Time> vecDotTimer         = getTimer("Kokkos: MV: Dot (contiguous)");
  RCP<Time> vecDotTimer2        = getTimer("Kokkos: MV: Dot (noncontiguous)");
  RCP<Time> vecNrmInfTimer      = getTimer("Kokkos: MV: NrmInf (contiguous)");
  RCP<Time> vecNrmInfTimer2     = getTimer("Kokkos: MV: NrmInf (noncontiguous)");
  RCP<Time> vecAxpyTimer        = getTimer("Kokkos: MV: Axpy");
  RCP<Time> vecAxpbyTimer       = getTimer("Kokkos: MV: Axpby");
  RCP<Time> vecScalTimer        = getTimer("Kokkos: MV: Scal");

  // Benchmark creation of a MultiVector.
  mv_type x;
  {
    TimeMonitor timeMon(*vecCreateTimer);
    // This benchmarks both vector creation and vector destruction.
    for (int k = 0; k < numTrials; ++k) {
      x = mv_type("x", numRows, numCols);
    }
  }

  // Benchmark filling a Vector with zero.
  {
    TimeMonitor timeMon(*vecFillZeroTimer);
    for (int k = 0; k < numTrials; ++k) {
      KokkosBlas::fill(x, 0.0);
    }
  }

  // Benchmark filling a Vector with a nonzero value.
  {
    TimeMonitor timeMon(*vecFillNonZeroTimer);
    for (int k = 0; k < numTrials; ++k) {
      KokkosBlas::fill(x, 1.0);
    }
  }

  // Benchmark computing the (square of the) 2-norm of a MultiVector.
  typedef Kokkos::View<double*, Kokkos::DefaultHostExecutionSpace> norms_type;
  norms_type norms("norms", numCols);
  {
    TimeMonitor timeMon(*vecNrm2Timer);
    for (int k = 0; k < numTrials; ++k) {
      KokkosBlas::nrm2_squared(norms, x);
    }
  }

  if (numTrials > 0 && numCols > 0) {
    for (int j = 0; j < numCols; ++j) {
      const double expectedResult = static_cast<double>(numRows);
      if (norms(j) != expectedResult) {
        out << "Kokkos 2-norm (squared) result is wrong!  Expected " << expectedResult << " but got " << norms(j)
            << " instead." << endl;
        success = false;
      }
    }
  }

  // Benchmark computing the (square of the) 2-norm of a MultiVector,
  // using the 4-argument variant.
  {
    TimeMonitor timeMon(*vecNrm2Timer2);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        norms(j) = KokkosBlas::nrm2_squared(subview(x, ALL(), j));
      }
    }
  }

  if (numTrials > 0 && numCols > 0) {
    for (int j = 0; j < numCols; ++j) {
      const double expectedResult = static_cast<double>(numRows);
      if (norms(j) != expectedResult) {
        out << "Kokkos 2-norm (squared) result (3-arg variant) is wrong!  "
            << "Expected " << expectedResult << " but got " << norms(j) << " instead." << endl;
        success = false;
      }
    }
  }

  // Benchmark computing the 1-norm of a MultiVector.
  {
    TimeMonitor timeMon(*vecNrm1Timer);
    for (int k = 0; k < numTrials; ++k) {
      KokkosBlas::nrm1(norms, x);
    }
  }

  if (numTrials > 0 && numCols > 0) {
    for (int j = 0; j < numCols; ++j) {
      const double expectedResult = static_cast<double>(numRows);
      if (norms(j) != expectedResult) {
        out << "Kokkos 1-norm result is wrong!  Expected " << expectedResult << " but got " << norms(j) << " instead."
            << endl;
        success = false;
      }
    }
  }

  // Benchmark computing the 1-norm of a MultiVector, using the
  // 4-argument variant.
  {
    TimeMonitor timeMon(*vecNrm1Timer2);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        norms(j) = KokkosBlas::nrm1(subview(x, ALL(), j));
      }
    }
  }

  if (numTrials > 0 && numCols > 0) {
    for (int j = 0; j < numCols; ++j) {
      const double expectedResult = static_cast<double>(numRows);
      if (norms(j) != expectedResult) {
        out << "Kokkos 1-norm result (3-arg variant) is wrong!  "
            << "Expected " << expectedResult << " but got " << norms(j) << " instead." << endl;
        success = false;
      }
    }
  }

  // Benchmark computing the dot product of two MultiVectors.
  mv_type y("y", numRows, numCols);
  KokkosBlas::fill(y, -1.0);
  typedef Kokkos::View<double*, Kokkos::DefaultHostExecutionSpace> dots_type;
  dots_type dots("dots", numCols);
  {
    TimeMonitor timeMon(*vecDotTimer);
    for (int k = 0; k < numTrials; ++k) {
      KokkosBlas::dot(dots, x, y);
    }
  }

  if (numTrials > 0 && numCols > 0) {
    for (int j = 0; j < numCols; ++j) {
      const double expectedResult = static_cast<double>(numRows) * -1.0;
      if (dots(j) != expectedResult) {
        out << "Kokkos dot product result is wrong!  Expected " << expectedResult << " but got " << (j) << " instead."
            << endl;
        success = false;
      }
    }
  }

  // Benchmark computing the dot product of two MultiVectors,
  // using the variant that selects a column at a time.
  {
    TimeMonitor timeMon(*vecDotTimer2);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        dots(j) = KokkosBlas::dot(subview(x, ALL(), j), subview(y, ALL(), j));
      }
    }
  }

  if (numTrials > 0 && numCols > 0) {
    for (int j = 0; j < numCols; ++j) {
      const double expectedResult = static_cast<double>(numRows) * -1.0;
      if (dots(j) != expectedResult) {
        out << "Kokkos dot product result (5-arg variant) is wrong!  "
            << "Expected " << expectedResult << " but got " << (j) << " instead." << endl;
        success = false;
      }
    }
  }

  // Benchmark computing the inf-norm of a MultiVector.
  {
    TimeMonitor timeMon(*vecNrmInfTimer);
    for (int k = 0; k < numTrials; ++k) {
      KokkosBlas::nrminf(norms, x);
    }
  }

  if (numTrials > 0 && numCols > 0) {
    for (int j = 0; j < numCols; ++j) {
      const double expectedResult = 1.0;
      if (norms(j) != expectedResult) {
        out << "Kokkos inf-norm result is wrong!  Expected " << expectedResult << " but got " << norms(j) << " instead."
            << endl;
        success = false;
      }
    }
  }

  // Benchmark computing the inf-norm of a MultiVector, using the
  // 4-argument variant.
  {
    TimeMonitor timeMon(*vecNrmInfTimer2);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        norms(j) = KokkosBlas::nrminf(subview(x, ALL(), j));
      }
    }
  }

  if (numTrials > 0 && numCols > 0) {
    for (int j = 0; j < numCols; ++j) {
      const double expectedResult = 1.0;
      if (norms(j) != expectedResult) {
        out << "Kokkos inf-norm result (3-arg variant) is wrong!  "
            << "Expected " << expectedResult << " but got " << norms(j) << " instead." << endl;
        success = false;
      }
    }
  }

  // Benchmark y := alpha*x + beta*y for beta = 0 and beta != 0.
  {
    TimeMonitor timeMon(*vecAxpyTimer);
    const double alpha = 3.0;
    const double beta  = 0.0;

    for (int k = 0; k < numTrials; ++k) {
      KokkosBlas::axpby(alpha, x, beta, y);
    }
  }
  {
    TimeMonitor timeMon(*vecAxpbyTimer);
    const double alpha = 3.0;
    const double beta  = 4.0;

    for (int k = 0; k < numTrials; ++k) {
      KokkosBlas::axpby(alpha, x, beta, y);
    }
  }

  // Benchmark y := alpha*y.
  {
    TimeMonitor timeMon(*vecScalTimer);
    const double alpha = 0.5;

    for (int k = 0; k < numTrials; ++k) {
      KokkosBlas::scal(y, alpha, y);
    }
  }

  return success;
}

bool benchmarkRaw(std::ostream& out, const int numRows, const int numCols, const int numTrials) {
  using std::endl;
  RCP<Time> vecCreateTimer      = getTimer("Raw: MV: Create");
  RCP<Time> vecFillZeroTimer    = getTimer("Raw: MV: Fill zero");
  RCP<Time> vecFillNonzeroTimer = getTimer("Raw: MV: Fill nonzero");
  RCP<Time> vecNrm2Timer        = getTimer("Raw: MV: Nrm2");
  RCP<Time> vecNrm1Timer        = getTimer("Raw: MV: Nrm1");
  RCP<Time> vecDotTimer         = getTimer("Raw: MV: Dot");
  RCP<Time> vecNrmInfTimer      = getTimer("Raw: MV: NrmInf");
  RCP<Time> vecAxpyTimer        = getTimer("Raw: MV: Axpy");
  RCP<Time> vecAxpbyTimer       = getTimer("Raw: MV: Axpby");
  bool success                  = true;

  if (numTrials <= 0) {
    return success;  // trivial success
  }

  // Benchmark creation of a MultiVector.
  double* x = NULL;
  {
    TimeMonitor timeMon(*vecCreateTimer);
    // This benchmarks both vector creation and vector destruction.
    for (int k = 0; k < numTrials; ++k) {
      x = new double[numRows * numCols];
      memset(x, 0, numRows * numCols * sizeof(double));
      if (k + 1 < numTrials) {
        delete[] x;
      }
    }
  }

  // Benchmark filling a Vector with zeros.
  {
    TimeMonitor timeMon(*vecFillZeroTimer);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        double* x_j = x + numRows * j;
        for (int i = 0; i < numRows; ++i) {
          x_j[i] = 0.0;
        }
      }
    }
  }

  // Benchmark filling a Vector with a nonzero value.
  {
    TimeMonitor timeMon(*vecFillNonzeroTimer);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        double* x_j = x + numRows * j;
        for (int i = 0; i < numRows; ++i) {
          x_j[i] = 1.0;
        }
      }
    }
  }

  // Benchmark computing the (square of the) 2-norm of a MultiVector.
  double* norms = new double[numCols];
  {
    TimeMonitor timeMon(*vecNrm2Timer);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        double sum  = 0.0;
        double* x_j = x + numRows * j;
        for (int i = 0; i < numRows; ++i) {
          sum += x_j[i] * x_j[i];
        }
        norms[j] = sum;
      }
    }
  }

  // Benchmark computing the 1-norm of a MultiVector.
  {
    TimeMonitor timeMon(*vecNrm1Timer);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        double sum  = 0.0;
        double* x_j = x + numRows * j;
        for (int i = 0; i < numRows; ++i) {
          const double tmp = x_j[i] < 0.0 ? -x_j[i] : x_j[i];
          sum += tmp;
        }
        norms[j] = sum;
      }
    }
  }

  // Benchmark computing the inf-norm of a MultiVector.
  {
    TimeMonitor timeMon(*vecNrmInfTimer);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        double norm = 0.0;
        double* x_j = x + numRows * j;
        for (int i = 0; i < numRows; ++i) {
          const double tmp = x_j[i];
          if (norm < tmp) {
            norm = tmp;
          }
        }
        norms[j] = norm;
      }
    }
  }

  // Benchmark computing the dot product of two MultiVectors.
  double* y = new double[numRows * numCols];
  for (int j = 0; j < numCols; ++j) {
    double* y_j = y + numRows * j;
    for (int i = 0; i < numRows; ++i) {
      y_j[i] = -1.0;
    }
  }
  double* dots = new double[numCols];
  for (int j = 0; j < numCols; ++j) {
    dots[j] = 0.0;
  }
  {
    TimeMonitor timeMon(*vecDotTimer);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        double* const x_j = x + numRows * j;
        double* const y_j = y + numRows * j;
        double sum        = 0.0;
        for (int i = 0; i < numRows; ++i) {
          sum += x_j[i] * y_j[i];
        }
        dots[j] = sum;
      }
    }
  }

  // Benchmark y := alpha*x + beta*y for beta = 0 and beta != 0.
  {
    TimeMonitor timeMon(*vecAxpyTimer);
    const double alpha = 3.0;

    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        double* const x_j = x + numRows * j;
        double* const y_j = y + numRows * j;
        for (int i = 0; i < numRows; ++i) {
          y_j[i] += alpha * x_j[i];
        }
      }
    }
  }
  {
    TimeMonitor timeMon(*vecAxpbyTimer);
    const double alpha = 3.0;
    const double beta  = 4.0;

    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        double* const x_j = x + numRows * j;
        double* const y_j = y + numRows * j;
        for (int i = 0; i < numRows; ++i) {
          y_j[i] = alpha * x_j[i] + beta * y_j[i];
        }
      }
    }
  }

  if (norms != NULL) {
    delete[] norms;
    norms = NULL;
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
    return success;
  } else {  // numTrials > 0
    const double expectedResult = static_cast<double>(numRows) * -1.0;
    if (dots[0] != expectedResult) {
      out << "Raw dot product result is wrong!  Expected " << expectedResult << " but got " << dots[0] << " instead."
          << endl;
      if (dots != NULL) {
        delete[] dots;
      }
      success = false;
    } else {
      if (dots != NULL) {
        delete[] dots;
      }
    }
  }

  return success;
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
  int numRows    = 100000;
  int numCols    = 3;
  int numTrials  = 350;
  bool runKokkos = true;
  bool runRaw    = true;

  CommandLineProcessor cmdp;
  cmdp.setOption("numRows", &numRows, "Number of rows in the multivectors");
  cmdp.setOption("numCols", &numCols, "Number of columns in the multivectors");
  cmdp.setOption("numTrials", &numTrials,
                 "Number of timing loop iterations "
                 "for each event to time");
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
    TEUCHOS_TEST_FOR_EXCEPTION(numRows < 0, std::invalid_argument, "numRows must be nonnegative.");
    TEUCHOS_TEST_FOR_EXCEPTION(numCols < 0, std::invalid_argument, "numCols must be nonnegative.");
  }

  if (myRank == 0) {
    cout << endl
         << "---" << endl
         << "Command-line options:" << endl
         << "  numRows: " << numRows << endl
         << "  numCols: " << numCols << endl
         << "  numTrials: " << numTrials << endl
         << "  runKokkos: " << (runKokkos ? "true" : "false") << endl
         << "  runRaw: " << (runRaw ? "true" : "false") << endl
         << endl;
  }

  // Run the benchmark
  bool success = true;
  if (runKokkos) {
    const bool lclSuccess = benchmarkKokkos(cout, numRows, numCols, numTrials);
    success               = success && lclSuccess;
  }
  if (runRaw) {
    const bool lclSuccess = benchmarkRaw(cout, numRows, numCols, numTrials);
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
