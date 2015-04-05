// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#include <Kokkos_Blas1_MV.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Comm.hpp>
#ifdef HAVE_MPI
#  include <Teuchos_DefaultMpiComm.hpp>
#else
#  include <Teuchos_DefaultSerialComm.hpp>
#endif // HAVE_MPI

using Teuchos::Comm;
using Teuchos::CommandLineProcessor;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Time;
using Teuchos::TimeMonitor;

// Create a new timer with the given name if it hasn't already been
// created, else get the previously created timer with that name.
RCP<Time> getTimer (const std::string& timerName) {
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }
  return timer;
}

bool
benchmarkKokkos (std::ostream& out,
                 const int numRows,
                 const int numCols,
                 const int numTrials)
{
  using std::endl;
  typedef Kokkos::View<double**, Kokkos::LayoutLeft> mv_type;
  bool success = true;

  RCP<Time> vecCreateTimer = getTimer ("Kokkos: MV: Create");
  RCP<Time> vecFillTimer = getTimer ("Kokkos: MV: Fill");
  RCP<Time> vecNrm2Timer = getTimer ("Kokkos: MV: Nrm2 (2-arg)");
  RCP<Time> vecNrm2Timer2 = getTimer ("Kokkos: MV: Nrm2 (3-arg)");
  RCP<Time> vecNrm1Timer = getTimer ("Kokkos: MV: Nrm1 (2-arg)");
  RCP<Time> vecNrm1Timer2 = getTimer ("Kokkos: MV: Nrm1 (3-arg)");
  RCP<Time> vecDotTimer = getTimer ("Kokkos: MV: Dot (3-arg)");
  RCP<Time> vecDotTimer2 = getTimer ("Kokkos: MV: Dot (5-arg)");

  // Benchmark creation of a MultiVector.
  mv_type x;
  {
    TimeMonitor timeMon (*vecCreateTimer);
    // This benchmarks both vector creation and vector destruction.
    for (int k = 0; k < numTrials; ++k) {
      x = mv_type ("x", numRows, numCols);
    }
  }

  // Benchmark filling a Vector.
  {
    TimeMonitor timeMon (*vecFillTimer);
    for (int k = 0; k < numTrials; ++k) {
      KokkosBlas::fill (x, 1.0);
    }
  }

  // Benchmark computing the (square of the) 2-norm of a MultiVector.
  typedef Kokkos::View<double*, Kokkos::LayoutLeft> norms_type;
  norms_type norms ("norms", numCols);
  {
    TimeMonitor timeMon (*vecNrm2Timer);
    for (int k = 0; k < numTrials; ++k) {
      KokkosBlas::nrm2_squared (norms, x);
    }
  }

  if (numTrials > 0 && numCols > 0) {
    norms_type::HostMirror norms_h = Kokkos::create_mirror_view (norms);
    Kokkos::deep_copy (norms_h, norms);

    for (int j = 0; j < numCols; ++j) {
      const double expectedResult = static_cast<double> (numRows);
      if (norms_h(j) != expectedResult) {
        out << "Kokkos 2-norm (squared) result is wrong!  Expected "
            << expectedResult << " but got " << norms_h(j) << " instead."
            << endl;
        success = false;
      }
    }
  }

  // Benchmark computing the (square of the) 2-norm of a MultiVector,
  // using the 3-argument variant.
  {
    TimeMonitor timeMon (*vecNrm2Timer2);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        KokkosBlas::nrm2_squared (norms, x, j);
      }
    }
  }

  if (numTrials > 0 && numCols > 0) {
    norms_type::HostMirror norms_h = Kokkos::create_mirror_view (norms);
    Kokkos::deep_copy (norms_h, norms);

    for (int j = 0; j < numCols; ++j) {
      const double expectedResult = static_cast<double> (numRows);
      if (norms_h(j) != expectedResult) {
        out << "Kokkos 2-norm (squared) result (3-arg variant) is wrong!  "
            << "Expected " << expectedResult << " but got " << norms_h(j)
            << " instead." << endl;
        success = false;
      }
    }
  }

  // Benchmark computing the 1-norm of a MultiVector.
  {
    TimeMonitor timeMon (*vecNrm1Timer);
    for (int k = 0; k < numTrials; ++k) {
      KokkosBlas::nrm1 (norms, x);
    }
  }

  if (numTrials > 0 && numCols > 0) {
    norms_type::HostMirror norms_h = Kokkos::create_mirror_view (norms);
    Kokkos::deep_copy (norms_h, norms);

    for (int j = 0; j < numCols; ++j) {
      const double expectedResult = static_cast<double> (numRows);
      if (norms_h(j) != expectedResult) {
        out << "Kokkos 1-norm result is wrong!  Expected " << expectedResult
            << " but got " << norms_h(j) << " instead." << endl;
        success = false;
      }
    }
  }

  // Benchmark computing the 1-norm of a MultiVector, using the
  // 3-argument variant.
  {
    TimeMonitor timeMon (*vecNrm1Timer2);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        KokkosBlas::nrm1 (norms, x, j);
      }
    }
  }

  if (numTrials > 0 && numCols > 0) {
    norms_type::HostMirror norms_h = Kokkos::create_mirror_view (norms);
    Kokkos::deep_copy (norms_h, norms);

    for (int j = 0; j < numCols; ++j) {
      const double expectedResult = static_cast<double> (numRows);
      if (norms_h(j) != expectedResult) {
        out << "Kokkos 1-norm result (3-arg variant) is wrong!  "
            << "Expected " << expectedResult << " but got " << norms_h(j)
            << " instead." << endl;
        success = false;
      }
    }
  }


  mv_type y ("y", numRows, numCols);
  KokkosBlas::fill (y, -1.0);

  // Benchmark computing the dot product of two MultiVectors.
  typedef Kokkos::View<double*, Kokkos::LayoutLeft> dots_type;
  dots_type dots ("dots", numCols);
  {
    TimeMonitor timeMon (*vecDotTimer);
    for (int k = 0; k < numTrials; ++k) {
      KokkosBlas::dot (dots, x, y);
    }
  }

  if (numTrials > 0 && numCols > 0) {
    dots_type::HostMirror dots_h = Kokkos::create_mirror_view (dots);
    Kokkos::deep_copy (dots_h, dots);

    for (int j = 0; j < numCols; ++j) {
      const double expectedResult = static_cast<double> (numRows) * -1.0;
      if (dots_h(j) != expectedResult) {
        out << "Kokkos dot product result is wrong!  Expected " << expectedResult
            << " but got " << dots_h(j) << " instead." << endl;
        success = false;
      }
    }
  }

  // Benchmark computing the dot product of two MultiVectors,
  // using the variant that selects a column at a time.
  {
    TimeMonitor timeMon (*vecDotTimer2);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        KokkosBlas::dot (dots, x, j, y, j);
      }
    }
  }

  if (numTrials > 0 && numCols > 0) {
    dots_type::HostMirror dots_h = Kokkos::create_mirror_view (dots);
    Kokkos::deep_copy (dots_h, dots);

    for (int j = 0; j < numCols; ++j) {
      const double expectedResult = static_cast<double> (numRows) * -1.0;
      if (dots_h(j) != expectedResult) {
        out << "Kokkos dot product result (5-arg variant) is wrong!  "
            << "Expected " << expectedResult << " but got " << dots_h(j)
            << " instead." << endl;
        success = false;
      }
    }
  }

  return success;
}


bool
benchmarkRaw (std::ostream& out,
              const int numRows,
              const int numCols,
              const int numTrials)
{
  using std::endl;
  RCP<Time> vecCreateTimer = getTimer ("Raw: MV: Create");
  RCP<Time> vecFillTimer = getTimer ("Raw: MV: Fill");
  RCP<Time> vecNrm2Timer = getTimer ("Raw: MV: Nrm2");
  RCP<Time> vecNrm1Timer = getTimer ("Raw: MV: Nrm1");
  RCP<Time> vecDotTimer = getTimer ("Raw: MV: Dot");

  if (numTrials <= 0) {
    return true; // trivially
  }

  // Benchmark creation of a MultiVector.
  double* x = NULL;
  {
    TimeMonitor timeMon (*vecCreateTimer);
    // This benchmarks both vector creation and vector destruction.
    for (int k = 0; k < numTrials; ++k) {
      x = new double [numRows * numCols];
      memset (x, 0, numRows * numCols * sizeof (double));
      if (k + 1 < numTrials) {
        delete [] x;
      }
    }
  }

  // Benchmark filling a Vector.
  {
    TimeMonitor timeMon (*vecFillTimer);
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
  double* norms = new double [numCols];
  {
    TimeMonitor timeMon (*vecNrm2Timer);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        double sum = 0.0;
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
    TimeMonitor timeMon (*vecNrm1Timer);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        double sum = 0.0;
        double* x_j = x + numRows * j;
        for (int i = 0; i < numRows; ++i) {
          const double tmp = x_j[i] < 0.0 ? -x_j[i] : x_j[i];
          sum += tmp;
        }
        norms[j] = sum;
      }
    }
  }

  double* y = new double [numRows * numCols];
  for (int j = 0; j < numCols; ++j) {
    double* y_j = y + numRows * j;
    for (int i = 0; i < numRows; ++i) {
      y_j[i] = -1.0;
    }
  }

  // Benchmark computing the dot product of two MultiVectors.
  double* dots = new double [numCols];
  for (int j = 0; j < numCols; ++j) {
    dots[j] = 0.0;
  }
  {
    TimeMonitor timeMon (*vecDotTimer);
    for (int k = 0; k < numTrials; ++k) {
      for (int j = 0; j < numCols; ++j) {
        double* x_j = x + numRows * j;
        double* y_j = y + numRows * j;
        double sum = 0.0;
        for (int i = 0; i < numRows; ++i) {
          sum += x_j[i] * y_j[i];
        }
        dots[j] = sum;
      }
    }
  }
  if (norms != NULL) {
    delete [] norms;
    norms = NULL;
  }

  if (x != NULL) {
    delete [] x;
    x = NULL;
  }
  if (y != NULL) {
    delete [] y;
    y = NULL;
  }

  if (numTrials == 0) {
    return true; // trivially
  }
  else { // numTrials > 0
    const double expectedResult = static_cast<double> (numRows) * -1.0;
    if (dots[0] != expectedResult) {
      out << "Raw dot product result is wrong!  Expected " << expectedResult
          << " but got " << dots[0] << " instead." << endl;
      if (dots != NULL) {
        delete [] dots;
      }
      return false;
    } else {
      if (dots != NULL) {
        delete [] dots;
      }
      return true;
    }
  }
}


int
main (int argc, char* argv[])
{
  using std::cout;
  using std::endl;
  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
  Kokkos::initialize (argc, argv);

#ifdef HAVE_MPI
  RCP<const Comm<int> > comm = rcp (new Teuchos::MpiComm<int> (MPI_COMM_WORLD));
#else
  RCP<const Comm<int> > comm = rcp (new Teuchos::SerialComm<int> ());
#endif // HAVE_MPI

  //const int numProcs = comm->getSize (); // unused
  const int myRank = comm->getRank ();

  // Benchmark parameters
  int numRows = 100000;
  int numCols = 3;
  int numTrials = 350;
  bool runKokkos = true;
  bool runRaw = true;

  CommandLineProcessor cmdp;
  cmdp.setOption ("numRows", &numRows, "Number of rows in the multivectors");
  cmdp.setOption ("numCols", &numCols, "Number of columns in the multivectors");
  cmdp.setOption ("numTrials", &numTrials, "Number of timing loop iterations "
                  "for each event to time");
  cmdp.setOption ("runKokkos", "noKokkos", &runKokkos,
                  "Whether to run the Kokkos benchmark");
  cmdp.setOption ("runRaw", "noRaw", &runRaw,
                  "Whether to run the raw benchmark");
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
      numRows < 0, std::invalid_argument, "numRows must be nonnegative.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      numCols < 0, std::invalid_argument, "numCols must be nonnegative.");
  }

  if (myRank == 0) {
    cout << endl << "---" << endl
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
    const bool lclSuccess = benchmarkKokkos (cout, numRows, numCols, numTrials);
    success = success && lclSuccess;
  }
  if (runRaw) {
    const bool lclSuccess = benchmarkRaw (cout, numRows, numCols, numTrials);
    success = success && lclSuccess;
  }

  TimeMonitor::report (comm.ptr (), cout);
  if (success) {
    cout << "End Result: TEST PASSED" << endl;
  } else {
    cout << "End Result: TEST FAILED" << endl;
  }
  Kokkos::finalize ();
  return EXIT_SUCCESS;
}
