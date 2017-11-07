/*
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
*/

// This test tests whether the MPI implementation that Trilinos uses
// is CUDA aware.  See Trilinos GitHub issues #1571 and #1088 to learn
// what it means for an MPI implementation to be "CUDA aware," and why
// this matters for performance.
//
// The test will only build if CUDA is enabled.  If you want to
// exercise this test, you must build Tpetra with CUDA enabled, and
// set the environment variable TPETRA_ASSUME_CUDA_AWARE_MPI to some
// true value (e.g., "1" or "TRUE").  If you set the environment
// variable to some false value (e.g., "0" or "FALSE"), the test will
// run but will pass trivially.  This means that you may control the
// test's behavior at run time.

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Details_assumeMpiIsCudaAware.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Map.hpp" // creating a Map ensures Kokkos initialization
#include "Teuchos_CommHelpers.hpp"
#include "Kokkos_TeuchosCommAdapters.hpp"

namespace { // (anonymous)
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  typedef Tpetra::global_size_t GST;

#if ! defined(KOKKOS_ENABLE_CUDA) || ! defined(HAVE_TPETRA_INST_CUDA)
#  error "Building this test requires that Trilinos was built with CUDA enabled, and that Tpetra_INST_CUDA:BOOL=ON.  The latter should be true by default if the former is true.  Thus, if Trilinos was built with CUDA enabled, then you must have set some nondefault CMake option."
#endif // ! defined(KOKKOS_ENABLE_CUDA) && ! defined(HAVE_TPETRA_INST_CUDA)

  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> test_device_type;
  typedef Kokkos::Compat::KokkosCudaWrapperNode test_node_type;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<LO, GO, test_node_type> map_type;

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST( Comm, CudaAwareMpi )
  {
    constexpr bool debug = false;
    int lclSuccess = 1; // to be modified below
    int gblSuccess = 0; // output argument; to be modified below

    out << "Testing CUDA-awareness of the MPI implementation" << endl;
    Teuchos::OSTab tab1 (out);
    const bool assumeMpiIsCudaAware =
      Tpetra::Details::assumeMpiIsCudaAware (&out);
    {
      const bool fromBehavior =
        ::Tpetra::Details::Behavior::assumeMpiIsCudaAware ();
      TEST_EQUALITY( assumeMpiIsCudaAware, fromBehavior );
    }
    if (! assumeMpiIsCudaAware) {
      out << "Trilinos (or you, the user) asserts that MPI is NOT CUDA aware."
          << endl
          << "That's OK; we consider the test having passed in this case,"
          << endl
          << "but we won't run any of the test after this point."
          << endl;
      return;
    }

    RCP<const Comm<int> > comm = Tpetra::TestingUtilities::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();
    if (numProcs < 2) {
      out << "This test is more meaningful if run with at least 2 MPI "
        "processes.  You ran with only 1 MPI process." << endl;
    }
    const std::string prefix ([=] { // compare to Lisp's LET (yay const!)
        std::ostringstream os;
        os << "Proc " << myRank << ": ";
        return os.str ();
      } ());

    // Create a Map, just to initialize Kokkos.
    map_type map (numProcs, 1, 0, comm);

    // Create Views for sending and receiving data.
    Kokkos::View<int*, test_device_type> sendBuf ("sendBuf", 1);
    Kokkos::View<int*, test_device_type> recvBuf ("recvBuf", 1);
    auto recvBuf_h = Kokkos::create_mirror_view (recvBuf);

    out << "Test whether self messsages work correctly" << endl;
    {
      const int msgTag = 42; // just some nontrivial tag for MPI messages

      // Fill send buffer with some unique positive value.
      const int correctValue = myRank + 1;
      Kokkos::deep_copy (sendBuf, correctValue);
      // Fill receive buffer with a flag value, always negative.
      const int flagValue = -(myRank + 1);
      Kokkos::deep_copy (recvBuf, flagValue);
      if (debug) {
        std::ostringstream os;
        os << prefix << "Self message!" << endl;
        std::cerr << os.str ();
      }
      const int srcRank = myRank; // in this case
      const int tgtRank = myRank; // in this case
      auto recvReq = Teuchos::ireceive (recvBuf, srcRank, msgTag, *comm);
      auto sendReq = Teuchos::isend (sendBuf, tgtRank, msgTag, *comm);

      Teuchos::wait (*comm, outArg (sendReq));
      Teuchos::wait (*comm, outArg (recvReq));

      Kokkos::deep_copy (recvBuf_h, recvBuf);
      TEST_EQUALITY( recvBuf_h(0), correctValue );
    }

    // Make sure that everybody finished and got the right answer.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Self-message test failed on some process; aborting test!" << endl;
    }
    else {
      out << "Self-message test succeeded on all processes!" << endl;
    }

    if (numProcs > 1) {
      out << "Test whether messages between different processes work "
        "correctly" << endl;
      // Some tag, nontrivial, different than the tag used above (in
      // case something went wrong, we don't want messages from above
      // to get mixed up with messages below).
      const int msgTag = 43;
      const int correctValue = 1;
      const int flagValue = -1;

      // We only need to exercise the first two processes.  The rest
      // may stay idle.
      if (myRank == 0) { // sending process
        Kokkos::deep_copy (sendBuf, correctValue);
        const int tgtRank = 1;
        auto sendReq = Teuchos::isend (sendBuf, tgtRank, msgTag, *comm);
        Teuchos::wait (*comm, outArg (sendReq));
      }
      else if (myRank == 1) { // receiving process
        Kokkos::deep_copy (recvBuf, flagValue);
        const int srcRank = 0;
        auto recvReq = Teuchos::ireceive (recvBuf, srcRank, msgTag, *comm);
        Teuchos::wait (*comm, outArg (recvReq));
        Kokkos::deep_copy (recvBuf_h, recvBuf);
        TEST_EQUALITY( recvBuf_h(0), correctValue );
      }
    }

    // Make sure that everybody finished and got the right answer.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Other-message test failed on some process!" << endl;
    }
    else {
      out << "Other-message test succeeded on all processes!" << endl;
    }
  }

} // namespace (anonymous)


