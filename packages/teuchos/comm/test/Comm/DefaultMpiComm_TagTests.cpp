/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
// @HEADER
*/

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"

namespace std { 

template <typename Packet>
ostream & operator<< (ostream& os, const pair<Packet, Packet>& arg)
{
  os << "(" << arg.first << "," << arg.second << ")";
  return os;
}

} // namespace std

namespace {
using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::CommRequest;
using Teuchos::CommStatus;
using Teuchos::ireceive;
using Teuchos::send;
using Teuchos::MpiComm;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::waitAll;
using std::endl;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP ();
  clp.addOutputSetupOptions (true);
}

TEUCHOS_UNIT_TEST( DefaultMpiComm, Tag )
{
  typedef ArrayRCP<int>::size_type size_type;

  RCP<const Comm<int> > comm = rcp (new MpiComm<int> (MPI_COMM_WORLD));
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  // If there is only one process, then the left and right neighbors
  // are the same, namely the calling process.  MPI allows a process
  // to send to and receive from itself.  On the other hand, we want
  // to test blocking sends, and we want to post two sends, so we only
  // allow > 1 processes.
  if (numProcs == 1) {
    out << "numProcs == 1, test passes trivially." << endl;
    return;
  }
  out << "Test setup" << endl;

  // If there are only two processes, the left neighbor and the right
  // neighbor are the same, namely the other process.
  int leftNeighbor = (myRank - 1) % numProcs;
  int rightNeighbor = (myRank + 1) % numProcs;
  // C doesn't guarantee nonnegativity of the result of the % operator.
  if (leftNeighbor < 0) {
    leftNeighbor += numProcs;
  }
  if (rightNeighbor < 0) {
    rightNeighbor += numProcs;
  }
  Array<int> expectedSourceRanks (2); // expected source ranks for receives
  expectedSourceRanks[0] = leftNeighbor;
  expectedSourceRanks[1] = rightNeighbor;
  std::sort (expectedSourceRanks.begin (), expectedSourceRanks.end ());

  // Receive buffer, with subbuffers for each neighbor.
  ArrayRCP<int> recvBuf (2);
  ArrayRCP<int> leftRecvBuf = recvBuf.persistingView (0, 1);
  ArrayRCP<int> rightRecvBuf = recvBuf.persistingView (1, 1);
  
  // Send buffer, with subbuffers for each neighbor.
  Array<int> sendBuf (2);
  ArrayView<int> leftSendBuf = sendBuf.view (0, 1);
  ArrayView<int> rightSendBuf = sendBuf.view (1, 1);

  ////////////////////////////////////////////////////////////////////
  // First round of messages
  ////////////////////////////////////////////////////////////////////
  out << "Round 1 of messages" << endl;

  // Tag to use for the first set of messages.
  const int tag1 = 42;

  // Fill receive buffer with error flags.
  for (size_type k = 0; k < recvBuf.size (); ++k) {
    recvBuf[k] = -1;
  }
  
  // Send my process rank plus the current tag to all neighbors.
  for (size_type k = 0; k < sendBuf.size (); ++k) {
    sendBuf[k] = myRank + tag1;
  }
  
  // Post receives from left and right neighbors.
  Array<RCP<CommRequest<int> > > requests (2);
  requests[0] = ireceive<int, int> (leftRecvBuf, leftNeighbor, tag1, *comm);
  requests[1] = ireceive<int, int> (rightRecvBuf, rightNeighbor, tag1, *comm);

  // Post sends to left and right neighbors.
  send<int, int> (leftSendBuf.getRawPtr (), as<int> (leftSendBuf.size ()), 
		  leftNeighbor, tag1, *comm);
  send<int, int> (rightSendBuf.getRawPtr (), as<int> (rightSendBuf.size ()),
		  rightNeighbor, tag1, *comm);

  // Wait for the receives to complete.
  Array<RCP<CommStatus<int> > > statuses (2);
  waitAll (*comm, requests (), statuses ());

  // Make sure the source ranks are correct.
  Array<int> sourceRanks (2);
  for (size_type k = 0; k < 2; ++k) {
    sourceRanks[k] = statuses[k]->getSourceRank ();
  }
  std::sort (sourceRanks.begin (), sourceRanks.end ());
  TEST_EQUALITY( sourceRanks.size (), expectedSourceRanks.size () );
  for (size_type k = 0; k < sourceRanks.size (); ++k) {
    TEST_EQUALITY( sourceRanks[k], expectedSourceRanks[k] );
  }

  // Make sure the source tags are correct.
  for (size_type k = 0; k < statuses.size (); ++k) {
    TEST_EQUALITY( statuses[k]->getTag (), tag1 );
  }

  // Make sure the message contents are correct.
  TEST_EQUALITY( leftRecvBuf[0], leftNeighbor + tag1 );
  TEST_EQUALITY( rightRecvBuf[0], rightNeighbor + tag1 );

  ////////////////////////////////////////////////////////////////////
  // Second round of messages
  ////////////////////////////////////////////////////////////////////
  out << "Round 2 of messages" << endl;

  // Tag to use for the second set of messages.
  const int tag2 = 100;

  // Fill receive buffer with error flags.
  for (size_type k = 0; k < recvBuf.size (); ++k) {
    recvBuf[k] = -1;
  }
  
  // Send my process rank plus the current tag to all neighbors.
  for (size_type k = 0; k < sendBuf.size (); ++k) {
    sendBuf[k] = myRank + tag2;
  }
  
  // Post receives from left and right neighbors.
  requests[0] = ireceive<int, int> (leftRecvBuf, leftNeighbor, tag2, *comm);
  requests[1] = ireceive<int, int> (rightRecvBuf, rightNeighbor, tag2, *comm);

  // Post sends to left and right neighbors.
  send<int, int> (leftSendBuf.getRawPtr (), as<int> (leftSendBuf.size ()), 
		  leftNeighbor, tag2, *comm);
  send<int, int> (rightSendBuf.getRawPtr (), as<int> (rightSendBuf.size ()),
		  rightNeighbor, tag2, *comm);

  // Wait for the receives to complete.
  waitAll (*comm, requests (), statuses ());

  // Make sure the source ranks are correct.
  for (size_type k = 0; k < 2; ++k) {
    sourceRanks[k] = statuses[k]->getSourceRank ();
  }
  std::sort (sourceRanks.begin (), sourceRanks.end ());
  TEST_EQUALITY( sourceRanks.size (), expectedSourceRanks.size () );
  for (size_type k = 0; k < sourceRanks.size (); ++k) {
    TEST_EQUALITY( sourceRanks[k], expectedSourceRanks[k] );
  }

  // Make sure the source tags are correct.
  for (size_type k = 0; k < statuses.size (); ++k) {
    TEST_EQUALITY( statuses[k]->getTag (), tag2 );
  }

  // Make sure the message contents are correct.
  TEST_EQUALITY( leftRecvBuf[0], leftNeighbor + tag2 );
  TEST_EQUALITY( rightRecvBuf[0], rightNeighbor + tag2 );

  ////////////////////////////////////////////////////////////////////
  // Third round of messages
  ////////////////////////////////////////////////////////////////////
  out << "Round 3 of messages" << endl;

  // In this round, we try again with the first tag.  This will tell
  // us if any first-round messages got mixed up with second-round
  // messages.
  const int tag3 = tag1;

  // Fill receive buffer with error flags.
  for (size_type k = 0; k < recvBuf.size (); ++k) {
    recvBuf[k] = -1;
  }
  
  // Send my process rank plus the current tag to all neighbors.
  for (size_type k = 0; k < sendBuf.size (); ++k) {
    sendBuf[k] = myRank + tag3;
  }
  
  // Post receives from left and right neighbors.
  requests[0] = ireceive<int, int> (leftRecvBuf, leftNeighbor, tag3, *comm);
  requests[1] = ireceive<int, int> (rightRecvBuf, rightNeighbor, tag3, *comm);

  // Post sends to left and right neighbors.
  send<int, int> (leftSendBuf.getRawPtr (), as<int> (leftSendBuf.size ()), 
		  leftNeighbor, tag3, *comm);
  send<int, int> (rightSendBuf.getRawPtr (), as<int> (rightSendBuf.size ()),
		  rightNeighbor, tag3, *comm);

  // Wait for the receives to complete.
  waitAll (*comm, requests (), statuses ());

  // Make sure the source ranks are correct.
  for (size_type k = 0; k < 2; ++k) {
    sourceRanks[k] = statuses[k]->getSourceRank ();
  }
  std::sort (sourceRanks.begin (), sourceRanks.end ());
  TEST_EQUALITY( sourceRanks.size (), expectedSourceRanks.size () );
  for (size_type k = 0; k < sourceRanks.size (); ++k) {
    TEST_EQUALITY( sourceRanks[k], expectedSourceRanks[k] );
  }

  // Make sure the source tags are correct.
  for (size_type k = 0; k < statuses.size (); ++k) {
    TEST_EQUALITY( statuses[k]->getTag (), tag3 );
  }

  // Make sure the message contents are correct.
  TEST_EQUALITY( leftRecvBuf[0], leftNeighbor + tag3 );
  TEST_EQUALITY( rightRecvBuf[0], rightNeighbor + tag3 );

  ////////////////////////////////////////////////////////////////////
  // Final check
  ////////////////////////////////////////////////////////////////////
  out << "Final check" << endl;

  // At this point, if we do a barrier, all the processes should reach
  // it.  None should hang.  If this test times out, it probably means
  // that not all the processes reached this point.
  comm->barrier ();
  out << "All processes successfully completed this test." << endl;
}

} // namespace
