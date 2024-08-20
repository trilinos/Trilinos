// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This test tests whether the MPI implementation that Trilinos uses
// is CUDA aware.  See Trilinos GitHub issues #1571 and #1088 to learn
// what it means for an MPI implementation to be "CUDA aware," and why
// this matters for performance.
//
// The test will only build if CUDA/HIP is enabled.  If you want to
// exercise this test, you must build Tpetra with CUDA or HIP enabled, and
// set the environment variable TPETRA_ASSUME_GPU_AWARE_MPI to some
// true value (e.g., "1" or "TRUE").  If you set the environment
// variable to some false value (e.g., "0" or "FALSE"), the test will
// run but will pass trivially.  This means that you may control the
// test's behavior at run time.

#include "Tpetra_TestingUtilities.hpp"
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

#if (! defined(KOKKOS_ENABLE_CUDA) || ! defined(HAVE_TPETRA_INST_CUDA) ) && (! defined(KOKKOS_ENABLE_HIP) || ! defined(HAVE_TPETRA_INST_HIP)) && (! defined(KOKKOS_ENABLE_SYCL) || ! defined(HAVE_TPETRA_INST_SYCL))
#  error "Building this test requires that Trilinos was built with CUDA, HIP or SYCL enabled, and that Tpetra_INST_CUDA:BOOL=ON, Tpetra_INST_HIP:BOOL=ON or Tpetra_INST_SYCL:BOOL=ON.  The latter should be true by default if the former is true.  Thus, if Trilinos was built with CUDA/HIP enabled, then you must have set some nondefault CMake option."
#endif

#if defined(KOKKOS_ENABLE_CUDA)
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> test_device_type;
  typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode test_node_type;
#elif defined(KOKKOS_ENABLE_HIP)
  typedef Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace> test_device_type;
  typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode test_node_type;
#elif defined(KOKKOS_ENABLE_SYCL)
  typedef Kokkos::Device<Kokkos::Experimental::SYCL, Kokkos::Experimental::SYCLDeviceUSMSpace> test_device_type;
  typedef Tpetra::KokkosCompat::KokkosSYCLWrapperNode test_node_type;
#endif
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<LO, GO, test_node_type> map_type;

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST( Comm, GPUAwareMpi )
  {
    constexpr bool debug = false;
    int lclSuccess = 1; // to be modified below
    int gblSuccess = 0; // output argument; to be modified below

    out << "Testing GPU-awareness of the MPI implementation" << endl;
    const bool assumeMpiIsGPUAware =
      ::Tpetra::Details::Behavior::assumeMpiIsGPUAware ();
    if (! assumeMpiIsGPUAware) {
      out << "Trilinos (or you, the user) asserts that MPI is NOT GPU aware."
          << endl
          << "That's OK; we consider the test having passed in this case,"
          << endl
          << "but we won't run any of the test after this point."
          << endl;
      return;
    }
    else {
      out << "Trilinos (or you, the user) asserts that MPI is GPU aware."
          << endl
          << "Beginning the test now."
          << endl;
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


