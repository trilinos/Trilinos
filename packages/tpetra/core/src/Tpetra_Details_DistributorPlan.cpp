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
// ************************************************************************
// @HEADER

#include "Tpetra_Details_DistributorPlan.hpp"

#include "Tpetra_Util.hpp"
#include <numeric>

namespace Tpetra {
namespace Details {

std::string
DistributorSendTypeEnumToString (EDistributorSendType sendType)
{
  if (sendType == DISTRIBUTOR_ISEND) {
    return "Isend";
  }
  else if (sendType == DISTRIBUTOR_RSEND) {
    return "Rsend";
  }
  else if (sendType == DISTRIBUTOR_SEND) {
    return "Send";
  }
  else if (sendType == DISTRIBUTOR_SSEND) {
    return "Ssend";
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid "
      "EDistributorSendType enum value " << sendType << ".");
  }
}

std::string
DistributorHowInitializedEnumToString (EDistributorHowInitialized how)
{
  switch (how) {
  case Details::DISTRIBUTOR_NOT_INITIALIZED:
    return "Not initialized yet";
  case Details::DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_SENDS:
    return "By createFromSends";
  case Details::DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_RECVS:
    return "By createFromRecvs";
  case Details::DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_SENDS_N_RECVS:
    return "By createFromSendsAndRecvs";
  case Details::DISTRIBUTOR_INITIALIZED_BY_REVERSE:
    return "By createReverseDistributor";
  case Details::DISTRIBUTOR_INITIALIZED_BY_COPY:
    return "By copy constructor";
  default:
    return "INVALID";
  }
}

DistributorPlan::DistributorPlan(Teuchos::RCP<const Teuchos::Comm<int>> comm)
  : comm_(comm),
    howInitialized_(DISTRIBUTOR_NOT_INITIALIZED),
    sendType_(DISTRIBUTOR_SEND),
    barrierBetweenRecvSend_(barrierBetween_default),
    useDistinctTags_(useDistinctTags_default),
    sendMessageToSelf_(false),
    numSendsToOtherProcs_(0),
    maxSendLength_(0),
    numReceives_(0),
    totalReceiveLength_(0)
{ }

DistributorPlan::DistributorPlan(const DistributorPlan& otherPlan)
  : comm_(otherPlan.comm_),
    howInitialized_(DISTRIBUTOR_INITIALIZED_BY_COPY),
    sendType_(otherPlan.sendType_),
    barrierBetweenRecvSend_(otherPlan.barrierBetweenRecvSend_),
    useDistinctTags_(otherPlan.useDistinctTags_),
    sendMessageToSelf_(otherPlan.sendMessageToSelf_),
    numSendsToOtherProcs_(otherPlan.numSendsToOtherProcs_),
    procIdsToSendTo_(otherPlan.procIdsToSendTo_),
    startsTo_(otherPlan.startsTo_),
    lengthsTo_(otherPlan.lengthsTo_),
    maxSendLength_(otherPlan.maxSendLength_),
    indicesTo_(otherPlan.indicesTo_),
    numReceives_(otherPlan.numReceives_),
    totalReceiveLength_(otherPlan.totalReceiveLength_),
    lengthsFrom_(otherPlan.lengthsFrom_),
    procsFrom_(otherPlan.procsFrom_),
    startsFrom_(otherPlan.startsFrom_),
    indicesFrom_(otherPlan.indicesFrom_)
{ }

int DistributorPlan::getTag(const int pathTag) const {
  return useDistinctTags_ ? pathTag : comm_->getTag();
}

void DistributorPlan::computeReceives()
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::as;
  using Teuchos::CommStatus;
  using Teuchos::CommRequest;
  using Teuchos::ireceive;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_SUM;
  using Teuchos::receive;
  using Teuchos::reduce;
  using Teuchos::scatter;
  using Teuchos::send;
  using Teuchos::waitAll;

  const int myRank = comm_->getRank();
  const int numProcs = comm_->getSize();

  // MPI tag for nonblocking receives and blocking sends in this method.
  const int pathTag = 2;
  const int tag = getTag(pathTag);

  // toProcsFromMe[i] == the number of messages sent by this process
  // to process i.  The data in numSendsToOtherProcs_, procIdsToSendTo_, and lengthsTo_
  // concern the contiguous sends.  Therefore, each process will be
  // listed in procIdsToSendTo_ at most once, and so toProcsFromMe[i] will
  // either be 0 or 1.
  {
    Array<int> toProcsFromMe (numProcs, 0);
#ifdef HAVE_TEUCHOS_DEBUG
    bool counting_error = false;
#endif // HAVE_TEUCHOS_DEBUG
    for (size_t i = 0; i < (numSendsToOtherProcs_ + (sendMessageToSelf_ ? 1 : 0)); ++i) {
#ifdef HAVE_TEUCHOS_DEBUG
      if (toProcsFromMe[procIdsToSendTo_[i]] != 0) {
        counting_error = true;
      }
#endif // HAVE_TEUCHOS_DEBUG
      toProcsFromMe[procIdsToSendTo_[i]] = 1;
    }
#ifdef HAVE_TEUCHOS_DEBUG
    SHARED_TEST_FOR_EXCEPTION(counting_error, std::logic_error,
        "Tpetra::Distributor::computeReceives: There was an error on at least "
        "one process in counting the number of messages send by that process to "
        "the other processs.  Please report this bug to the Tpetra developers.",
        *comm_);
#endif // HAVE_TEUCHOS_DEBUG

    // Compute the number of receives that this process needs to
    // post.  The number of receives includes any self sends (i.e.,
    // messages sent by this process to itself).
    //
    // (We will use numReceives_ this below to post exactly that
    // number of receives, with MPI_ANY_SOURCE as the sending rank.
    // This will tell us from which processes this process expects
    // to receive, and how many packets of data we expect to receive
    // from each process.)
    //
    // toProcsFromMe[i] is the number of messages sent by this
    // process to process i.  Compute the sum (elementwise) of all
    // the toProcsFromMe arrays on all processes in the
    // communicator.  If the array x is that sum, then if this
    // process has rank j, x[j] is the number of messages sent
    // to process j, that is, the number of receives on process j
    // (including any messages sent by process j to itself).
    //
    // Yes, this requires storing and operating on an array of
    // length P, where P is the number of processes in the
    // communicator.  Epetra does this too.  Avoiding this O(P)
    // memory bottleneck would require some research.
    //
    // mfh 09 Jan 2012, 15 Jul 2015: There are three ways to
    // implement this O(P) memory algorithm.
    //
    //   1. Use MPI_Reduce and MPI_Scatter: reduce on the root
    //      process (0) from toProcsFromMe, to numRecvsOnEachProc.
    //      Then, scatter the latter, so that each process p gets
    //      numRecvsOnEachProc[p].
    //
    //   2. Like #1, but use MPI_Reduce_scatter instead of
    //      MPI_Reduce and MPI_Scatter.  MPI_Reduce_scatter might be
    //      optimized to reduce the number of messages, but
    //      MPI_Reduce_scatter is more general than we need (it
    //      allows the equivalent of MPI_Scatterv).  See Bug 6336.
    //
    //   3. Do an all-reduce on toProcsFromMe, and let my process
    //      (with rank myRank) get numReceives_ from
    //      toProcsFromMe[myRank].  The HPCCG miniapp uses the
    //      all-reduce method.
    //
    // Approaches 1 and 3 have the same critical path length.
    // However, #3 moves more data.  This is because the final
    // result is just one integer, but #3 moves a whole array of
    // results to all the processes.  This is why we use Approach 1
    // here.
    //
    // mfh 12 Apr 2013: See discussion in createFromSends() about
    // how we could use this communication to propagate an error
    // flag for "free" in a release build.

    const int root = 0; // rank of root process of the reduction
    Array<int> numRecvsOnEachProc; // temp; only needed on root
    if (myRank == root) {
      numRecvsOnEachProc.resize (numProcs);
    }
    int numReceivesAsInt = 0; // output
    reduce<int, int> (toProcsFromMe.getRawPtr (),
        numRecvsOnEachProc.getRawPtr (),
        numProcs, REDUCE_SUM, root, *comm_);
    scatter<int, int> (numRecvsOnEachProc.getRawPtr (), 1,
        &numReceivesAsInt, 1, root, *comm_);
    numReceives_ = static_cast<size_t> (numReceivesAsInt);
  }

  // Now we know numReceives_, which is this process' number of
  // receives.  Allocate the lengthsFrom_ and procsFrom_ arrays
  // with this number of entries.
  lengthsFrom_.assign (numReceives_, 0);
  procsFrom_.assign (numReceives_, 0);

  //
  // Ask (via nonblocking receive) each process from which we are
  // receiving how many packets we should expect from it in the
  // communication pattern.
  //

  // At this point, numReceives_ includes any self message that
  // there may be.  At the end of this routine, we'll subtract off
  // the self message (if there is one) from numReceives_.  In this
  // routine, we don't need to receive a message from ourselves in
  // order to figure out our lengthsFrom_ and source process ID; we
  // can just ask ourselves directly.  Thus, the actual number of
  // nonblocking receives we post here does not include the self
  // message.
  const size_t actualNumReceives = numReceives_ - (sendMessageToSelf_ ? 1 : 0);

  // Teuchos' wrapper for nonblocking receives requires receive
  // buffers that it knows won't go away.  This is why we use RCPs,
  // one RCP per nonblocking receive request.  They get allocated in
  // the loop below.
  Array<RCP<CommRequest<int> > > requests (actualNumReceives);
  Array<ArrayRCP<size_t> > lengthsFromBuffers (actualNumReceives);
  Array<RCP<CommStatus<int> > > statuses (actualNumReceives);

  // Teuchos::Comm treats a negative process ID as MPI_ANY_SOURCE
  // (receive data from any process).
#ifdef HAVE_MPI
  const int anySourceProc = MPI_ANY_SOURCE;
#else
  const int anySourceProc = -1;
#endif

  // Post the (nonblocking) receives.
  for (size_t i = 0; i < actualNumReceives; ++i) {
    // Once the receive completes, we can ask the corresponding
    // CommStatus object (output by wait()) for the sending process'
    // ID (which we'll assign to procsFrom_[i] -- don't forget to
    // do that!).
    lengthsFromBuffers[i].resize (1);
    lengthsFromBuffers[i][0] = as<size_t> (0);
    requests[i] = ireceive<int, size_t> (lengthsFromBuffers[i], anySourceProc,
        tag, *comm_);
  }

  // Post the sends: Tell each process to which we are sending how
  // many packets it should expect from us in the communication
  // pattern.  We could use nonblocking sends here, as long as we do
  // a waitAll() on all the sends and receives at once.
  //
  // We assume that numSendsToOtherProcs_ and sendMessageToSelf_ have already been
  // set.  The value of numSendsToOtherProcs_ (my process' number of sends) does
  // not include any message that it might send to itself.
  for (size_t i = 0; i < numSendsToOtherProcs_ + (sendMessageToSelf_ ? 1 : 0); ++i) {
    if (procIdsToSendTo_[i] != myRank) {
      // Send a message to procIdsToSendTo_[i], telling that process that
      // this communication pattern will send that process
      // lengthsTo_[i] blocks of packets.
      const size_t* const lengthsTo_i = &lengthsTo_[i];
      send<int, size_t> (lengthsTo_i, 1, as<int> (procIdsToSendTo_[i]), tag, *comm_);
    }
    else {
      // We don't need a send in the self-message case.  If this
      // process will send a message to itself in the communication
      // pattern, then the last element of lengthsFrom_ and
      // procsFrom_ corresponds to the self-message.  Of course
      // this process knows how long the message is, and the process
      // ID is its own process ID.
      lengthsFrom_[numReceives_-1] = lengthsTo_[i];
      procsFrom_[numReceives_-1] = myRank;
    }
  }

  //
  // Wait on all the receives.  When they arrive, check the status
  // output of wait() for the receiving process ID, unpack the
  // request buffers into lengthsFrom_, and set procsFrom_ from the
  // status.
  //
  waitAll (*comm_, requests (), statuses ());
  for (size_t i = 0; i < actualNumReceives; ++i) {
    lengthsFrom_[i] = *lengthsFromBuffers[i];
    procsFrom_[i] = statuses[i]->getSourceRank ();
  }

  // Sort the procsFrom_ array, and apply the same permutation to
  // lengthsFrom_.  This ensures that procsFrom_[i] and
  // lengthsFrom_[i] refers to the same thing.
  sort2 (procsFrom_.begin(), procsFrom_.end(), lengthsFrom_.begin());

  // Compute indicesFrom_
  totalReceiveLength_ =
    std::accumulate (lengthsFrom_.begin (), lengthsFrom_.end (), 0);
  indicesFrom_.clear ();

  startsFrom_.clear ();
  startsFrom_.reserve (numReceives_);
  for (size_t i = 0, j = 0; i < numReceives_; ++i) {
    startsFrom_.push_back(j);
    j += lengthsFrom_[i];
  }

  if (sendMessageToSelf_) {
    --numReceives_;
  }
}

}
}
