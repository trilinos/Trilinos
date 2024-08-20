// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
using Teuchos::isend;
using Teuchos::MpiComm;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::send;
using Teuchos::waitAll;
using std::endl;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP ();
  clp.addOutputSetupOptions (true);
}

TEUCHOS_UNIT_TEST( MpiCommTag, IrecvSend )
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


TEUCHOS_UNIT_TEST( MpiCommTag, IrecvIsend )
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
  ArrayRCP<int> sendBuf (2);
  ArrayRCP<int> leftSendBuf = sendBuf.persistingView (0, 1);
  ArrayRCP<int> rightSendBuf = sendBuf.persistingView (1, 1);

  // Requests for both nonblocking receives and nonblocking sends.
  Array<RCP<CommRequest<int> > > requests (4);
  // Statuses for both nonblocking receives and nonblocking sends.  We
  // only use these to test that the ranks of the received messages
  // were correct.
  Array<RCP<CommStatus<int> > > statuses (4);

  ////////////////////////////////////////////////////////////////////
  // First round of messages
  ////////////////////////////////////////////////////////////////////
  out << "Round 1 of messages" << endl;

  // Tag to use for the first set of messages.
  const int tag1 = 101;

  // Fill receive buffer with error flags.
  for (size_type k = 0; k < recvBuf.size (); ++k) {
    recvBuf[k] = -1;
  }

  // Send my process rank plus the current tag to all neighbors.
  for (size_type k = 0; k < sendBuf.size (); ++k) {
    sendBuf[k] = myRank + tag1;
  }

  // Post receives from left and right neighbors.
  requests[0] = ireceive<int, int> (leftRecvBuf, leftNeighbor, tag1, *comm);
  requests[1] = ireceive<int, int> (rightRecvBuf, rightNeighbor, tag1, *comm);

  // Post sends to left and right neighbors.
  requests[2] = isend<int, int> (leftSendBuf, leftNeighbor, tag1, *comm);
  requests[3] = isend<int, int> (rightSendBuf, rightNeighbor, tag1, *comm);

  // Wait for the receives to complete.
  waitAll (*comm, requests (), statuses ());

  // Make sure the source tags are correct.
  for (size_type k = 0; k < 2; ++k) {
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
  const int tag2 = 202;

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
  requests[2] = isend<int, int> (leftSendBuf, leftNeighbor, tag2, *comm);
  requests[3] = isend<int, int> (rightSendBuf, rightNeighbor, tag2, *comm);

  // Wait for the receives to complete.
  waitAll (*comm, requests (), statuses ());

  // Make sure the source tags are correct.
  for (size_type k = 0; k < 2; ++k) {
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
  requests[2] = isend<int, int> (leftSendBuf, leftNeighbor, tag3, *comm);
  requests[3] = isend<int, int> (rightSendBuf, rightNeighbor, tag3, *comm);

  // Wait for the receives to complete.
  waitAll (*comm, requests (), statuses ());

  // Make sure the source tags are correct.
  for (size_type k = 0; k < 2; ++k) {
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
