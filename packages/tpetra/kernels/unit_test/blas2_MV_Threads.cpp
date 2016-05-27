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

#include <TpetraKernels_Blas2_MV_GEMV.hpp>
#include <Teuchos_Comm.hpp>
#ifdef HAVE_MPI
#  include <Teuchos_DefaultMpiComm.hpp>
#else
#  include <Teuchos_DefaultSerialComm.hpp>
#endif // HAVE_MPI
#include <Teuchos_CommandLineProcessor.hpp>

using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::rcp;

namespace { // (anonymous)

// Test all BLAS 2 operations, for all Scalar types in {double, float,
// int, complex<double>, complex<float>}, for all Layout types in
// {LayoutLeft, LayoutRight}, and for Device<Threads, HostSpace>.
//
// The 'indent' string must be passed in by value, so that the caller
// doesn't see changes inside (for incrementing indent level).
bool
runThreadsTests (std::ostream& out,
                std::string indent,
                const bool testFloat,
                const bool testComplex,
                const bool printOnlyOnFailure,
                const bool debug)
{
  bool success = true;
  bool curSuccess = true;

  indent = indent + "  ";
  out << indent << "Test all BLAS 2 operations" << endl;

#ifdef KOKKOS_HAVE_PTHREAD
  {
    indent = indent + "  ";
    out << indent << "Test Threads, HostSpace" << endl;
    typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> device_type;
    curSuccess = testScalarsLayouts<device_type> (out, indent,
                                                  testFloat,
                                                  testComplex,
                                                  printOnlyOnFailure,
                                                  debug);
    success = success && curSuccess;
  }
#endif // KOKKOS_HAVE_PTHREAD

  return success;
}

} // namespace (anonymous)

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
  const int myRank = comm->getRank ();

  // mfh 10 May 2016: Currently, neither Scalar = float nor
  // Scalar = Kokkos::complex<float> work.  Skip them for now.
  bool testFloat = false;
  bool testComplex = true;
  bool printOnlyOnFailure = true;
  bool debug = false;

  const bool rejectUnrecognizedOptions = true;
  Teuchos::CommandLineProcessor cmdp (false, rejectUnrecognizedOptions);
  cmdp.setOption ("testFloat", "noTestFloat", &testFloat,
                  "Whether to test Scalar=float (and Scalar=complex<float>, "
                  "if testComplex is true)");
  cmdp.setOption ("testComplex", "noTestComplex", &testComplex,
                  "Whether to test complex arithmetic");
  cmdp.setOption ("printOnlyOnFailure", "noPrintOnlyOnFailure",
                  &printOnlyOnFailure,
                  "Whether to print per-test non-debug output only if that "
                  "test fails.");
  cmdp.setOption ("debug", "release", &debug,
                  "Whether to print extra debug output");

  if (cmdp.parse (argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    if (myRank == 0) {
      cout << "TEST FAILED to parse command-line arguments!" << endl;
    }
    return EXIT_FAILURE;
  }

  bool curSuccess = true;
  bool success = true;
  std::string indent = "";

  curSuccess = runThreadsTests (cout, indent + "  ",
                                testFloat,
                                testComplex,
                                printOnlyOnFailure,
                                debug);
  success = curSuccess && success;

  Kokkos::finalize ();

  if (success) {
    if (myRank == 0) {
      cout << "End Result: TEST PASSED" << endl;
    }
  }
  else {
    if (myRank == 0) {
      cout << "End Result: TEST FAILED" << endl;
    }
  }

  return EXIT_SUCCESS;
}
