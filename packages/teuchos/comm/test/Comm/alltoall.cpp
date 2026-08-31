// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#ifdef HAVE_TEUCHOS_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#endif // HAVE_TEUCHOS_MPI
#include "Teuchos_UnitTestHarness.hpp"


template<class PacketType>
bool
testAlltoall (bool& success, std::ostream& out,
              const Teuchos::Comm<int>& comm)
{
  using Teuchos::alltoAll;
  using Teuchos::TypeNameTraits;
  using std::endl;
  typedef PacketType packet_type;

  const int myRank = comm.getRank();
  const int numProcs = comm.getSize();
  const int send_count = 2;
  const int recv_count = 2;
  const int total_count = send_count * numProcs;

  // Teuchos constructs the output stream such that it only prints on
  // Process 0 anyway.
  out << "Testing Teuchos::alltoAll<int, " << TypeNameTraits<packet_type>::name ()
      << "> with send count = " << send_count  << " to each MPI rank"
      << " and numProcs = " << numProcs << endl;

  const packet_type ZERO = Teuchos::ScalarTraits<packet_type>::zero ();
  const packet_type ONE = Teuchos::ScalarTraits<packet_type>::one ();
  packet_type sendBuf[total_count];
  packet_type recvBuf[total_count];
  for (int i = 0; i < total_count; ++i) {
    sendBuf[i] = myRank * total_count + i;
    recvBuf[i] = ZERO;
  }
  alltoAll<int, packet_type> (sendBuf, send_count, recvBuf, recv_count, comm);

  // Don't trust that any other Teuchos communicator wrapper functions
  // work here.  Instead, if building with MPI, use raw MPI.  If not
  // building with MPI, first test that comm has only one process,
  // then assume this in the test.

#ifdef HAVE_TEUCHOS_MPI
  using Teuchos::MpiComm;
  int err = MPI_SUCCESS;

  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  TEUCHOS_TEST_FOR_EXCEPTION
    (mpiComm == NULL, std::logic_error, "Building with MPI, but default "
     "communicator is not a Teuchos::MpiComm!");
  MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());

  // Use a barrier to make sure that every process got this far.
  err = MPI_Barrier (rawMpiComm);
  TEUCHOS_TEST_FOR_EXCEPTION
    (err != MPI_SUCCESS, std::logic_error, "MPI_Barrier failed!");

  // Recompute using MPI.  Use an all-reduce to simplify the test.
  packet_type sendBuf2[total_count];
  packet_type recvBuf2[total_count];
  for (int i = 0; i < total_count; ++i) {
    sendBuf2[i] = myRank * total_count + i;
    recvBuf2[i] = ZERO;
  }
  // FIXME (14 Jul 2015) Need to get the MPI_Datatype right for PacketType.
  MPI_Datatype rawMpiType;
  if (typeid (packet_type) == typeid (short)) {
    rawMpiType = MPI_SHORT;
  } else if (typeid (packet_type) == typeid (unsigned short)) {
    rawMpiType = MPI_UNSIGNED_SHORT;
  } else if (typeid (packet_type) == typeid (int)) {
    rawMpiType = MPI_INT;
  } else if (typeid (packet_type) == typeid (unsigned int)) {
    rawMpiType = MPI_UNSIGNED;
  } else if (typeid (packet_type) == typeid (long)) {
    rawMpiType = MPI_LONG;
  } else if (typeid (packet_type) == typeid (unsigned long)) {
    rawMpiType = MPI_UNSIGNED_LONG;
  } else if (typeid (packet_type) == typeid (long long)) {
    rawMpiType = MPI_LONG_LONG;
  } else if (typeid (packet_type) == typeid (unsigned long long)) {
    rawMpiType = MPI_UNSIGNED_LONG_LONG;
  } else if (typeid (packet_type) == typeid (float)) {
    rawMpiType = MPI_FLOAT;
  } else if (typeid (packet_type) == typeid (double)) {
    rawMpiType = MPI_DOUBLE;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::logic_error, "Unimplemented conversion from PacketType = "
       << TypeNameTraits<packet_type>::name () << " to MPI_Datatype.");
  }
  err = MPI_Alltoall (sendBuf2, send_count, rawMpiType,
                      recvBuf2, recv_count, rawMpiType,
                      rawMpiComm);
  TEUCHOS_TEST_FOR_EXCEPTION
    (err != MPI_SUCCESS, std::logic_error, "MPI_Alltoall failed!");
  for (int i = 0; i < total_count; ++i) {
    TEST_EQUALITY( recvBuf2[i], recvBuf[i] );
  }

  int lclSuccess = success ? 1 : 0;
  int gblSuccess = lclSuccess;
  err = MPI_Allreduce (&lclSuccess, &gblSuccess, 1, MPI_INT, MPI_MIN, rawMpiComm);
  TEUCHOS_TEST_FOR_EXCEPTION
    (err != MPI_SUCCESS, std::logic_error, "MPI_Alltoall failed!");
  success = gblSuccess == 1 ? true : false;
#else // HAVE_TEUCHOS_MPI
  TEUCHOS_TEST_FOR_EXCEPTION
    (comm.getSize () != 1, std::logic_error, "Not building with MPI, but "
     "communicator has size = " << comm.getSize () << " != 1.  We don't know "
     "how to test this case.");
  for (int i = 0; i < total_count; ++i) {
    TEST_EQUALITY( sendBuf[i], recvBuf[i] );
  }
#endif // HAVE_TEUCHOS_MPI

  return success;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Comm, AlltoAll, PacketType )
{
  using Teuchos::Comm;
  using Teuchos::RCP;
  typedef PacketType packet_type;

  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
  const bool curSuccess = testAlltoall<packet_type> (success, out, *comm);
  TEST_EQUALITY_CONST( curSuccess, true );
  success = success && curSuccess;
}

//
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Comm, AlltoAll, int )
