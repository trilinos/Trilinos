// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Details_PackTriples.hpp"
#ifdef HAVE_TPETRACORE_MPI
#  include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"
#endif // HAVE_TPETRACORE_MPI

namespace { // (anonymous)

#ifdef HAVE_TPETRACORE_MPI

typedef double scalar_type;
typedef int ordinal_type;

// Compute and send the total pack size in bytes.
void
sendPackSize (bool& success,
              Teuchos::FancyOStream& out,
              int& packSize, // output argument; pack size in bytes
              const int numEnt,
              const int destProc,
              const int tag,
              const Teuchos::Comm<int>& comm)
{
  using ::Tpetra::Details::countPackTriples;
  using ::Tpetra::Details::countPackTriplesCount;
  using ::Tpetra::Details::extractMpiCommFromTeuchos;
  int curPackSize = 0;
  int errCode = MPI_SUCCESS;
  packSize = 0;

  curPackSize = 0;
  errCode = countPackTriplesCount (comm, curPackSize, &out);
  TEST_ASSERT( errCode == MPI_SUCCESS );
  if (errCode != MPI_SUCCESS) {
    return;
  }
  packSize += curPackSize;

  curPackSize = 0;
  errCode =
    countPackTriples<scalar_type, ordinal_type> (numEnt, comm,
                                                 curPackSize, &out);
  TEST_ASSERT( errCode == MPI_SUCCESS );
  if (errCode != MPI_SUCCESS) {
    return;
  }
  packSize += curPackSize;

  MPI_Comm mpiComm = extractMpiCommFromTeuchos (comm);
  errCode = MPI_Send (&packSize, 1, MPI_INT, destProc, tag, mpiComm);
  TEST_ASSERT( errCode == MPI_SUCCESS );
}

void
recvPackSize (bool& success,
              Teuchos::FancyOStream& out,
              int& packSize,
              const int srcProc,
              const int tag,
              const Teuchos::Comm<int>& comm)
{
  using ::Tpetra::Details::extractMpiCommFromTeuchos;

  int errCode = MPI_SUCCESS;
  packSize = 0;
  MPI_Comm mpiComm = extractMpiCommFromTeuchos (comm);
  errCode = MPI_Recv (&packSize, 1, MPI_INT, srcProc, tag, mpiComm,
                      MPI_STATUS_IGNORE);
  TEST_ASSERT( errCode == MPI_SUCCESS );
}

void
packAndSendTriples (bool& success,
                    Teuchos::FancyOStream& out,
                    const ordinal_type gblRowInds[],
                    const ordinal_type gblColInds[],
                    const scalar_type vals[],
                    const int numEnt,
                    char outBuf[],
                    const int outBufSize,
                    const int destProc,
                    const int tag,
                    const Teuchos::Comm<int>& comm)
{
  using ::Tpetra::Details::packTriples;
  using ::Tpetra::Details::packTriplesCount;
  using ::Tpetra::Details::extractMpiCommFromTeuchos;

  int errCode = MPI_SUCCESS;
  int outBufCurPos = 0;
  const size_t minNumBytesForTriples =
    static_cast<size_t> (numEnt) *
    (2 * sizeof (ordinal_type) + sizeof (scalar_type));

  {
    std::ostringstream os;
    os << "packAndSendTriples: numEnt=" << numEnt << std::endl;
    std::cerr << os.str ();
  }
  errCode = packTriplesCount (numEnt, outBuf, outBufSize,
                              outBufCurPos, comm, &out);
  TEST_ASSERT( errCode == MPI_SUCCESS );
  if (errCode != MPI_SUCCESS) {
    return;
  }
  TEST_ASSERT( static_cast<size_t> (outBufCurPos) >= sizeof (int) );

  errCode = packTriples (gblRowInds, gblColInds, vals, numEnt, outBuf,
                         outBufSize, outBufCurPos, comm, &out);
  TEST_ASSERT( errCode == MPI_SUCCESS );
  if (errCode != MPI_SUCCESS) {
    return;
  }
  TEST_ASSERT( static_cast<size_t> (outBufCurPos) >=
               minNumBytesForTriples + sizeof (int) );

  MPI_Comm mpiComm = extractMpiCommFromTeuchos (comm);
  errCode = MPI_Send (outBuf, outBufSize, MPI_BYTE, destProc, tag, mpiComm);
  TEST_ASSERT( errCode == MPI_SUCCESS );
}

void
recvAndUnpackTriples (bool& success,
                      Teuchos::FancyOStream& out,
                      char inBuf[],
                      const int inBufSize,
                      ordinal_type gblRowInds[],
                      ordinal_type gblColInds[],
                      scalar_type vals[],
                      int& numEnt, // output argument!
                      const int srcProc,
                      const int tag,
                      const Teuchos::Comm<int>& comm)
{
  using ::Tpetra::Details::unpackTriples;
  using ::Tpetra::Details::unpackTriplesCount;
  using ::Tpetra::Details::extractMpiCommFromTeuchos;

  int errCode = MPI_SUCCESS;
  MPI_Comm mpiComm = extractMpiCommFromTeuchos (comm);
  MPI_Status status;
  errCode = MPI_Recv (inBuf, inBufSize, MPI_BYTE,
                      srcProc, tag, mpiComm, &status);
  TEST_ASSERT( errCode == MPI_SUCCESS );

  // Test that we received the expected number of bytes.
  int actualRecvCount = 0;
  errCode = MPI_Get_count (&status, MPI_BYTE, &actualRecvCount);
  TEST_ASSERT( errCode == MPI_SUCCESS );
  TEST_ASSERT( actualRecvCount == inBufSize );

  int inBufCurPos = 0;
  errCode = unpackTriplesCount (inBuf, inBufSize, inBufCurPos,
                                numEnt, comm, &out);
  TEST_ASSERT( errCode == MPI_SUCCESS );
  if (errCode != MPI_SUCCESS) {
    return;
  }
  TEST_ASSERT( inBufCurPos <= inBufSize );
  if (inBufCurPos > inBufSize) {
    return;
  }
  TEST_ASSERT( numEnt >= 0 );
  if (numEnt < 0) {
    return;
  }
  {
    std::ostringstream os;
    os << "recvAndUnpackTriples: numEnt=" << numEnt << std::endl;
    std::cerr << os.str ();
  }

  errCode = unpackTriples (inBuf, inBufSize, inBufCurPos, gblRowInds,
                           gblColInds, vals, numEnt, comm, &out);
  TEST_ASSERT( errCode == MPI_SUCCESS );
  TEST_ASSERT( inBufCurPos <= inBufSize );
}

void
testTriples (bool& success,
             Teuchos::FancyOStream& out,
             const Teuchos::Comm<int>& comm)
{
  constexpr int expectedNumEnt = 5;
  const ordinal_type expectedGblRowInds[] = {93, 418, 31, 666, 93};
  const ordinal_type expectedGblColInds[] = {7, 4, 6, 5, 3};
  const scalar_type expectedVals[] = {7.0, 4.0, 6.0, 5.0, 3.0};

  const int myRank = comm.getRank ();
  const int numProcs = comm.getSize ();

  // This test requires at least 2 MPI processes.
  TEST_ASSERT( numProcs >= 2 );
  if (numProcs < 2) {
    return;
  }

  const int tag = 101; // just something not 0

  if (myRank == 0) {
    const int srcProc = 1;
    int inBufSize = 0;
    recvPackSize (success, out, inBufSize, srcProc, tag, comm);
    // For the test, we just happen to know the correct lower bound.
    // In practice, the receiving process wouldn't know this unless told.
    TEST_ASSERT( static_cast<size_t> (inBufSize) >=
                 sizeof (int) +
                 expectedNumEnt * (2 * sizeof (ordinal_type) + sizeof (scalar_type)) );
    std::vector<char> inBuf (inBufSize);

    // For the test, we just happen to know the number of entries.  In
    // practice, we might need dynamic resizing, if the output arrays
    // are not long enough.
    std::vector<ordinal_type> gblRowInds (expectedNumEnt);
    std::vector<ordinal_type> gblColInds (expectedNumEnt);
    std::vector<scalar_type> vals (expectedNumEnt);
    int numEnt = 0; // output argument
    recvAndUnpackTriples (success, out, inBuf.data (), inBufSize,
                          gblRowInds.data (), gblColInds.data (), vals.data (),
                          numEnt, srcProc, tag, comm);
    TEST_ASSERT( numEnt == expectedNumEnt );
    if (numEnt == expectedNumEnt) {
      for (int k = 0; k < numEnt; ++k) {
        TEST_EQUALITY( expectedGblRowInds[k], gblRowInds[k] );
        TEST_EQUALITY( expectedGblColInds[k], gblColInds[k] );
        TEST_EQUALITY( expectedVals[k], vals[k] );
      }
    }
  }
  else if (myRank == 1) {
    const int destProc = 0;
    int outBufSize = 0;

    // I'm the sending process; I know how many entries I want to send.
    const int numEnt = expectedNumEnt;
    sendPackSize (success, out, outBufSize, numEnt, destProc, tag, comm);
    TEST_ASSERT( static_cast<size_t> (outBufSize) >=
                 sizeof (int) +
                 numEnt * (2 * sizeof (ordinal_type) + sizeof (scalar_type)) );

    std::vector<char> outBuf (outBufSize);
    packAndSendTriples (success, out, expectedGblRowInds, expectedGblColInds,
                        expectedVals, numEnt, outBuf.data (), outBufSize,
                        destProc, tag, comm);
  }
  //
  // This test only does anything significant on Processes 0 and 1.
  //
  {
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0; // output argument
    using ::Tpetra::Details::extractMpiCommFromTeuchos;
    MPI_Comm mpiComm = extractMpiCommFromTeuchos (comm);
    const int errCode = MPI_Allreduce (&lclSuccess, &gblSuccess, 1,
                                       MPI_INT, MPI_MIN, mpiComm);
    TEST_ASSERT( errCode == MPI_SUCCESS );
    TEST_ASSERT( gblSuccess == 1 );
  }
}

TEUCHOS_UNIT_TEST( PackTriples, intDouble )
{
  using Tpetra::TestingUtilities::getDefaultComm;
  using std::endl;

  out << "Testing packTriples and unpackTriples" << endl;
  Teuchos::OSTab tab1 (out);

  testTriples (success, out, * (getDefaultComm ()));
}

#endif // HAVE_TPETRACORE_MPI

} // namespace (anonymous)

