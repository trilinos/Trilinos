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

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Details_PackTriples.hpp"
#ifdef HAVE_TPETRACORE_MPI
#  include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"
#endif // HAVE_TPETRACORE_MPI

namespace { // (anonymous)

#ifdef HAVE_TPETRACORE_MPI

typedef double scalar_type;
typedef int ordinal_type;

void
sendPackSize (bool& success,
              Teuchos::FancyOStream& out,
              int& packSize,
              const int numEnt,
              const int destProc,
              const int tag,
              const Teuchos::Comm<int>& comm)
{
  using ::Tpetra::Details::countPackTriples;
  using ::Tpetra::Details::extractMpiCommFromTeuchos;

  int errCode = MPI_SUCCESS;
  packSize = 0;
  errCode =
    countPackTriples<scalar_type, ordinal_type> (numEnt, comm, packSize, &out);
  TEST_ASSERT( errCode == MPI_SUCCESS );

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
  using ::Tpetra::Details::extractMpiCommFromTeuchos;

  int errCode = MPI_SUCCESS;
  int outBufCurPos = 0;
  errCode = packTriples (gblRowInds, gblColInds, vals, numEnt, outBuf,
                         outBufSize, outBufCurPos, comm, &out);
  TEST_ASSERT( errCode == MPI_SUCCESS );
  TEST_ASSERT( static_cast<size_t> (outBufCurPos) >=
               static_cast<size_t> (numEnt) * (2 * sizeof (ordinal_type) + sizeof (scalar_type)));

  MPI_Comm mpiComm = extractMpiCommFromTeuchos (comm);
  errCode = MPI_Send (outBuf, outBufSize, MPI_BYTE, destProc, tag, mpiComm);
}

void
recvAndUnpackTriples (bool& success,
                      Teuchos::FancyOStream& out,
                      char inBuf[],
                      const int inBufSize,
                      ordinal_type gblRowInds[],
                      ordinal_type gblColInds[],
                      scalar_type vals[],
                      const int numEnt, // we would know this from DistObject
                      const int srcProc,
                      const int tag,
                      const Teuchos::Comm<int>& comm)
{
  using ::Tpetra::Details::unpackTriples;
  using ::Tpetra::Details::extractMpiCommFromTeuchos;

  int errCode = MPI_SUCCESS;
  MPI_Comm mpiComm = extractMpiCommFromTeuchos (comm);
  MPI_Status status;
  errCode = MPI_Recv (inBuf, inBufSize, MPI_BYTE,
                      srcProc, tag, mpiComm, &status);

  TEST_ASSERT( errCode == MPI_SUCCESS );
  int count = 0;
  errCode = MPI_Get_count (&status, MPI_BYTE, &count);
  TEST_ASSERT( errCode == MPI_SUCCESS );
  TEST_ASSERT( count <= inBufSize );
  // NOTE (mfh 16 Jan 2017) This assumes that sizeof is the right
  // way to measure how much space an ordinal_type or scalar_type
  // value occupies.
  TEST_ASSERT( static_cast<size_t> (count) >=
               static_cast<size_t> (numEnt * (2 * sizeof (ordinal_type) + sizeof (scalar_type))) );

  int inBufCurPos = 0;
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
  constexpr int numEnt = 5;
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
                 numEnt * (2 * sizeof (ordinal_type) + sizeof (scalar_type)) );
    std::vector<char> inBuf (inBufSize);

    // For the test, we just happen to know numEnt.
    // Tpetra::DistObject would normally get this from Distributor's
    // nonconstant-number-of-packets path.
    std::vector<ordinal_type> gblRowInds (numEnt);
    std::vector<ordinal_type> gblColInds (numEnt);
    std::vector<scalar_type> vals (numEnt);
    recvAndUnpackTriples (success, out, inBuf.data (), inBufSize,
                          gblRowInds.data (), gblColInds.data (), vals.data (),
                          numEnt, srcProc, tag, comm);
    for (int k = 0; k < numEnt; ++k) {
      TEST_EQUALITY( expectedGblRowInds[k], gblRowInds[k] );
      TEST_EQUALITY( expectedGblColInds[k], gblColInds[k] );
      TEST_EQUALITY( expectedVals[k], vals[k] );
    }
  }
  else if (myRank == 1) {
    const int destProc = 0;
    int outBufSize = 0;
    sendPackSize (success, out, outBufSize, numEnt, destProc, tag, comm);
    TEST_ASSERT( static_cast<size_t> (outBufSize) >=
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

