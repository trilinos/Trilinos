
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

#include <TpetraKernels_Blas1_MV_UnitTests.hpp>
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

namespace KokkosBlas {
namespace Impl {

bool
testOverScalarsAndLayoutsAndDevices (std::ostream& out, const int numCols,
                                     const bool oneCol, const bool testComplex)
{
  using std::endl;
  bool curSuccess = true;
  bool success = true;

  out << endl << "Test with numCols=" << numCols << endl << endl;

  // testOverScalarsAndLayouts is responsible for checking whether the
  // execution space is initialized, and skipping the test if it is.

#ifdef KOKKOS_HAVE_SERIAL
  {
    typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> device_type;
    const char deviceName[] = "Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>";

    if (Kokkos::Serial::is_initialized ()) {
      curSuccess = testOverScalarsAndLayouts<device_type> (out, deviceName, numCols, oneCol, testComplex);
      success = success && curSuccess;
    } else {
      out << endl << "Serial NOT initialized; skipping test" << endl;
    }
  }
#endif // KOKKOS_HAVE_SERIAL
#ifdef KOKKOS_HAVE_OPENMP
  {
    typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> device_type;
    const char deviceName[] = "Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>";

    if (Kokkos::OpenMP::is_initialized ()) {
      curSuccess = testOverScalarsAndLayouts<device_type> (out, deviceName, numCols, oneCol, testComplex);
      success = success && curSuccess;
    } else {
      out << endl << "OpenMP NOT initialized; skipping test" << endl;
    }
  }
#endif // KOKKOS_HAVE_OPENMP
#ifdef KOKKOS_HAVE_PTHREAD
  {
    typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> device_type;
    const char deviceName[] = "Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>";

    if (Kokkos::Threads::is_initialized ()) {
      curSuccess = testOverScalarsAndLayouts<device_type> (out, deviceName, numCols, oneCol, testComplex);
      success = success && curSuccess;
    } else {
      out << endl << "Threads NOT initialized; skipping test" << endl;
    }
  }
#endif // KOKKOS_HAVE_PTHREAD
#ifdef KOKKOS_HAVE_CUDA
  {
    typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> device_type;
    const char deviceName[] = "Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>";

    curSuccess = testOverScalarsAndLayouts<device_type> (out, deviceName, numCols, oneCol, testComplex);
    success = success && curSuccess;
  }
  {
    typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace> device_type;
    const char deviceName[] = "Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>";

    curSuccess = testOverScalarsAndLayouts<device_type> (out, deviceName, numCols, oneCol, testComplex);
    out << endl << "Cuda NOT initialized; skipping test" << endl;
  }
#endif // KOKKOS_HAVE_CUDA

  return success;
}

} // namespace Impl
} // namespace KokkosBlas


int
main (int argc, char* argv[])
{
  using KokkosBlas::Impl::testOverScalarsAndLayoutsAndDevices;
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

  // Number of columns in the 2-D View(s) to test.
  int numCols = 3;
  bool oneCol = false;
  bool testComplex = true;

  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("numCols", &numCols,
                  "Number of columns in the 2-D View(s) to test");
  cmdp.setOption ("oneCol", "noOneCol", &oneCol, "Whether to test the 1-D View "
                  "(single-column) versions of the kernels");
  cmdp.setOption ("testComplex", "noTestComplex", &testComplex,
                  "Whether to test complex arithmetic");
  if (cmdp.parse (argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    if (myRank == 0) {
      cout << "TEST FAILED to parse command-line arguments!" << endl;
    }
    return EXIT_FAILURE;
  }

  bool curSuccess = true;
  bool success = true;

  // Always test with numCols=1 first.
  curSuccess = testOverScalarsAndLayoutsAndDevices (cout, 1, oneCol, testComplex);
  success = curSuccess && success;
  if (numCols != 1) {
    curSuccess = testOverScalarsAndLayoutsAndDevices (cout, numCols,
                                                      oneCol, testComplex);
    success = curSuccess && success;
  }
  if (success) {
    if (myRank == 0) {
      cout << "End Result: TEST PASSED" << endl;
    }
  } else {
    if (myRank == 0) {
      cout << "End Result: TEST FAILED" << endl;
    }
  }
  Kokkos::finalize ();
  return EXIT_SUCCESS;
}
