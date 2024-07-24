// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
// #include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_iallreduce.hpp"
#ifdef HAVE_TPETRACORE_MPI
#  include "Tpetra_Details_MpiTypeTraits.hpp"
#endif // HAVE_TPETRACORE_MPI
#include "Teuchos_CommHelpers.hpp"
#include "Tpetra_Map.hpp" // creating a Map ensures Kokkos initialization
#include "Kokkos_ArithTraits.hpp"

namespace {
  using Tpetra::TestingUtilities::getDefaultComm;
  // using Tpetra::Details::gathervPrint;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  typedef Tpetra::global_size_t GST;

  typedef Tpetra::Map<> map_type;
  typedef map_type::local_ordinal_type LO;
  typedef map_type::global_ordinal_type GO;

  //
  // UNIT TESTS
  //

  template<class Packet, class Device>
  void
  testIallreduce (bool& success,
                  Teuchos::FancyOStream& out,
                  const std::string& packetTypeName,
                  const std::string& deviceTypeName,
                  const LO lclNumPackets,
                  const Teuchos::Comm<int>& comm)
  {
    using Tpetra::Details::iallreduce;
    using Teuchos::reduceAll;
    typedef typename Kokkos::ArithTraits<Packet>::val_type val_type;
    typedef Kokkos::ArithTraits<val_type> STS;
    typedef typename STS::mag_type mag_type;
    typedef typename Device::device_type device_type;

    out << "Test iallreduce for Packet=" << packetTypeName
        << " and Device=" << deviceTypeName << endl;
    Teuchos::OSTab tab1 (out);

    TEST_ASSERT( lclNumPackets >= static_cast<LO> (1) );
    if (lclNumPackets < static_cast<LO> (1)) {
      out << "lclNumPackets must be a positive integer" << endl;
      return;
    }

    int lclSuccess = 1; // to be updated below
    int gblSuccess = 0; // output argument

    //const int myRank = comm.getRank ();
    const int numProcs = comm.getSize ();

    Kokkos::View<val_type*, device_type> sendbuf ("sendbuf", lclNumPackets);
    Kokkos::View<val_type*, device_type> recvbuf ("recvbuf", lclNumPackets);

    auto sendbuf_h = Kokkos::create_mirror_view (sendbuf); // save & reuse
    auto recvbuf_h = Kokkos::create_mirror_view (recvbuf); // save & reuse

    // Fill input buffer with values chosen so that we know what their
    // sum across processes should be.
    {
      val_type curVal = STS::one ();
      for (LO k = 0; k < lclNumPackets; ++k, curVal += STS::one ()) {
        sendbuf_h(k) = curVal;
      }
      Kokkos::deep_copy (sendbuf, sendbuf_h);
    }

    // Make a "back-up" of the send buffer with input values, just in
    // case iallreduce has a bug that corrupts it (hopefully not).
    Kokkos::View<val_type*, device_type> sendbuf_bak ("sendbuf_bak", lclNumPackets);
    Kokkos::deep_copy (sendbuf_bak, sendbuf);
    auto sendbuf_bak_h = Kokkos::create_mirror_view (sendbuf_bak);
    Kokkos::deep_copy (sendbuf_bak_h, sendbuf_bak);

    out << "Test version of iallreduce that takes a rank-1 Kokkos::View" << endl;

    Kokkos::View<const val_type*, device_type> sendbuf_const = sendbuf;

    auto request = iallreduce (sendbuf_const, recvbuf, Teuchos::REDUCE_SUM, comm);
    TEST_ASSERT( request.get () != NULL );
    if (request.get () != NULL) {
      TEST_NOTHROW( request->wait () );
    }

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0; // output argument
    reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "iallreduce call or wait failed on at least one process!" << endl;
      return;
    }

    out << "Make sure the input values were not changed" << endl;
    {
      Kokkos::deep_copy (sendbuf_h, sendbuf);
      for (LO k = 0; k < lclNumPackets; ++k) {
        TEST_EQUALITY( sendbuf_h(k), sendbuf_bak_h(k) );
      }
    }

    out << "Make sure the output values are correct" << endl;
    {
      // There's no direct cast from int to some Scalar types (e.g.,
      // std::complex<mag_type> or Kokkos::complex<mag_type>).  Thus,
      // we make an intermediate cast through (real-valued) mag_type.
      const val_type np = static_cast<val_type> (static_cast<mag_type> (numProcs));

      Kokkos::deep_copy (recvbuf_h, recvbuf);
      for (LO k = 0; k < lclNumPackets; ++k) {
        TEST_EQUALITY( recvbuf_h(k), sendbuf_bak_h(k) * np );
      }
    }

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0; // output argument
    reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Leaving test early, due to failure" << endl;
      return;
    }

#if MPI_VERSION >= 3
    out << endl << "Test #850 fix (can alias sendbuf and recvbuf if comm is an "
      "intracommunicator)" << endl;
    // Use recvbuf as both the send buffer and the receive buffer.
    {
      val_type curVal = STS::one ();
      for (LO k = 0; k < lclNumPackets; ++k, curVal += STS::one ()) {
        recvbuf_h(k) = curVal;
      }
      Kokkos::deep_copy (recvbuf, recvbuf_h);
    }

    Kokkos::View<const val_type*, device_type> recvbuf_const = recvbuf;
    request = iallreduce (recvbuf_const, recvbuf, Teuchos::REDUCE_SUM, comm);
    TEST_ASSERT( request.get () != NULL );
    if (request.get () != NULL) {
      TEST_NOTHROW( request->wait () );
    }

    out << "Make sure the output values are correct" << endl;
    {
      // There's no direct cast from int to some Scalar types (e.g.,
      // std::complex<mag_type> or Kokkos::complex<mag_type>).  Thus,
      // we make an intermediate cast through (real-valued) mag_type.
      const val_type np = static_cast<val_type> (static_cast<mag_type> (numProcs));

      Kokkos::deep_copy (recvbuf_h, recvbuf);
      for (LO k = 0; k < lclNumPackets; ++k) {
        TEST_EQUALITY( recvbuf_h(k), sendbuf_bak_h(k) * np );
      }
    }

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0; // output argument
    reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Leaving test early, due to failure" << endl;
      return;
    }
#endif // MPI_VERSION >= 3

    out << "Test rank-0 Kokkos::View output version of iallreduce" << endl;
    auto sendbuf_0d = Kokkos::subview (sendbuf, 0);
    auto recvbuf_0d = Kokkos::subview (recvbuf, 0);
    auto sendbuf_h_0d = Kokkos::subview (sendbuf_h, 0);
    auto recvbuf_h_0d = Kokkos::subview (recvbuf_h, 0);

    // Reset contents of output buffer (important, since it currently
    // contains the correct answer!).
    recvbuf_h_0d() = STS::zero ();
    Kokkos::deep_copy (recvbuf_0d, recvbuf_h_0d);

    // Fill input buffer with value chosen so that we know what its
    // sum across processes should be.
    sendbuf_h_0d() = STS::one ();
    Kokkos::deep_copy (sendbuf_0d, sendbuf_h_0d);

    // Make a "back-up" of the send buffer with input values, just in
    // case iallreduce has a bug that corrupts it (hopefully not).
    auto sendbuf_bak_0d = Kokkos::subview (sendbuf_bak, 0);
    auto sendbuf_bak_h_0d = Kokkos::subview (sendbuf_bak_h, 0);
    Kokkos::deep_copy (sendbuf_bak_h_0d, sendbuf_h_0d);
    Kokkos::deep_copy (sendbuf_bak_0d, sendbuf_bak_h_0d);

    Kokkos::View<const val_type, device_type> sendbuf_0d_const = sendbuf_0d;
    request = iallreduce (sendbuf_0d_const, recvbuf_0d, Teuchos::REDUCE_SUM, comm);
    TEST_ASSERT( request.get () != NULL );
    if (request.get () != NULL) {
      TEST_NOTHROW( request->wait () );
    }

    out << "Make sure the input value was not changed" << endl;
    {
      Kokkos::deep_copy (sendbuf_h_0d, sendbuf_0d);
      TEST_EQUALITY( sendbuf_h_0d(), sendbuf_bak_h_0d() );
    }

    out << "Make sure the output value is correct" << endl;
    {
      // There's no direct cast from int to some Scalar types (e.g.,
      // std::complex<mag_type> or Kokkos::complex<mag_type>).  Thus,
      // we make an intermediate cast through (real-valued) mag_type.
      const val_type np = static_cast<val_type> (static_cast<mag_type> (numProcs));

      Kokkos::deep_copy (recvbuf_h_0d, recvbuf_0d);
      TEST_EQUALITY( recvbuf_h_0d(), sendbuf_bak_h_0d() * np );
    }

    // TODO (mfh 14 Nov 2016) Add the following tests:
    //
    //   1. That we can launch two iallreduce calls without them
    //      getting mixed up (e.g., make sure the implementation does
    //      not use some kind of nonreentrant global state)
    //
    //   2. That we can launch an iallreduce, then launch a (blocking)
    //      all-reduce, and finally wait on the iallreduce, without
    //      the two calls getting mixed up, and without deadlock
    //
    //   1. That we can launch an iallreduce, then launch nonblocking
    //      point-to-point operations (irecv and isend), and not
    //      deadlock or get the calls' results mixed up, no matter in
    //      what (sensible) order we wait on the operations.
  }

  template<class Device>
  void
  testIallreducePackets (bool& success,
                         Teuchos::FancyOStream& out,
                         const std::string& deviceTypeName,
                         const LO lclNumPackets,
                         const Teuchos::Comm<int>& comm)
  {
    typedef typename Device::device_type device_type;

    testIallreduce<short, device_type> (success, out, "short", deviceTypeName, lclNumPackets, comm);
    testIallreduce<int, device_type> (success, out, "int", deviceTypeName, lclNumPackets, comm);
    testIallreduce<long, device_type> (success, out, "long", deviceTypeName, lclNumPackets, comm);
    testIallreduce<long long, device_type> (success, out, "long long", deviceTypeName, lclNumPackets, comm);

    testIallreduce<unsigned short, device_type> (success, out, "unsigned short", deviceTypeName, lclNumPackets, comm);
    testIallreduce<unsigned int, device_type> (success, out, "unsigned int", deviceTypeName, lclNumPackets, comm);
    testIallreduce<unsigned long, device_type> (success, out, "unsigned long", deviceTypeName, lclNumPackets, comm);
    testIallreduce<unsigned long long, device_type> (success, out, "unsigned long long", deviceTypeName, lclNumPackets, comm);

    testIallreduce<float, device_type> (success, out, "float", deviceTypeName, lclNumPackets, comm);
    testIallreduce<double, device_type> (success, out, "double", deviceTypeName, lclNumPackets, comm);

    // FIXME (mfh 17 Nov 2016) Currently, if we don't have built-in
    // MPI_Datatype for these two Kokkos::complex specializations,
    // then Teuchos::REDUCE_SUM does not work.  This is because the
    // conversion from Teuchos::EReductionType to MPI_Op only uses
    // built-in MPI_Op, but we would need a custom MPI_Op if we had to
    // construct a custom MPI_Datatype for Kokkos::complex.

#ifdef HAVE_TPETRACORE_MPI
    using Tpetra::Details::MpiTypeTraits;
    if (MpiTypeTraits<Kokkos::complex<float> >::isSpecialized &&
        ! MpiTypeTraits<Kokkos::complex<float> >::needsFree) {
#endif // HAVE_TPETRACORE_MPI
      testIallreduce<Kokkos::complex<float>, device_type> (success, out, "float", deviceTypeName, lclNumPackets, comm);
#ifdef HAVE_TPETRACORE_MPI
    }
#endif // HAVE_TPETRACORE_MPI

#ifdef HAVE_TPETRACORE_MPI
    if (MpiTypeTraits<Kokkos::complex<double> >::isSpecialized &&
        ! MpiTypeTraits<Kokkos::complex<double> >::needsFree) {
#endif // HAVE_TPETRACORE_MPI
      testIallreduce<Kokkos::complex<double>, device_type> (success, out, "double", deviceTypeName, lclNumPackets, comm);
#ifdef HAVE_TPETRACORE_MPI
    }
#endif // HAVE_TPETRACORE_MPI

    // FIXME (mfh 14 Nov 2016) Test for all enabled Scalar types, not
    // just for the subset above.
  }


  TEUCHOS_UNIT_TEST( iallreduce, basic )
  {
    typedef map_type::device_type device_type;

    out << "Testing Tpetra::Details::iallreduce" << endl;
    Teuchos::OSTab tab1 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    const LO lclNumPackets = 5;
    const GO gblNumPackets = static_cast<GO> (lclNumPackets) *
      static_cast<GO> (comm->getSize ());
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (static_cast<GST> (gblNumPackets),
                         static_cast<size_t> (lclNumPackets),
                         indexBase, comm));

    const std::string deviceTypeName ("default");
    testIallreducePackets<device_type> (success, out, deviceTypeName, lclNumPackets, *comm);
  }

  TEUCHOS_UNIT_TEST( iallreduce, single_int)
  {
    RCP<const Comm<int> > comm = getDefaultComm ();
    //Use iallreduce to compute sum(1 ... nranks)
    int input = 1 + comm->getRank();
    int correctOutput = 0;
    for(int i = 0; i < comm->getSize(); i++)
      correctOutput += (1 + i);
    int output = -1;
    auto req = Tpetra::iallreduce(input, output, Teuchos::REDUCE_SUM, *comm);
    TEST_ASSERT( req.get () != NULL );
    if (req.get () != NULL) {
      TEST_NOTHROW( req->wait () );
    }
    TEST_EQUALITY(output, correctOutput);
  }

} // namespace (anonymous)


