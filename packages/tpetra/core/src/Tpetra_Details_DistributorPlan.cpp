// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_DistributorPlan.hpp"

#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include <numeric>

namespace Tpetra {
namespace Details {

std::string
DistributorSendTypeEnumToString (EDistributorSendType sendType)
{
  if (sendType == DISTRIBUTOR_ISEND) {
    return "Isend";
  }
  else if (sendType == DISTRIBUTOR_SEND) {
    return "Send";
  }
  else if (sendType == DISTRIBUTOR_ALLTOALL) {
    return "Alltoall";
  }
#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  else if (sendType == DISTRIBUTOR_MPIADVANCE_ALLTOALL) {
    return "MpiAdvanceAlltoall";
  }
  else if (sendType == DISTRIBUTOR_MPIADVANCE_NBRALLTOALLV) {
    return "MpiAdvanceNbralltoallv";
  }
#endif
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
#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
    mpixComm_(Teuchos::null),
#endif
    howInitialized_(DISTRIBUTOR_NOT_INITIALIZED),
    reversePlan_(Teuchos::null),
    sendType_(DISTRIBUTOR_SEND),
    sendMessageToSelf_(false),
    numSendsToOtherProcs_(0),
    maxSendLength_(0),
    numReceives_(0),
    totalReceiveLength_(0)
{ }

DistributorPlan::DistributorPlan(const DistributorPlan& otherPlan)
  : comm_(otherPlan.comm_),
#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
    mpixComm_(otherPlan.mpixComm_),
#endif
    howInitialized_(DISTRIBUTOR_INITIALIZED_BY_COPY),
    reversePlan_(otherPlan.reversePlan_),
    sendType_(otherPlan.sendType_),
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

size_t DistributorPlan::createFromSends(const Teuchos::ArrayView<const int>& exportProcIDs) {
  using Teuchos::outArg;
  using Teuchos::REDUCE_MAX;
  using Teuchos::reduceAll;
  using std::endl;
  const char rawPrefix[] = "Tpetra::DistributorPlan::createFromSends";

  const size_t numExports = exportProcIDs.size();
  const int myProcID = comm_->getRank();
  const int numProcs = comm_->getSize();
  const bool debug = Details::Behavior::debug("Distributor");

  // exportProcIDs tells us the communication pattern for this
  // distributor.  It dictates the way that the export data will be
  // interpreted in doPosts().  We want to perform at most one
  // send per process in doPosts; this is for two reasons:
  //   * minimize latency / overhead in the comm routines (nice)
  //   * match the number of receives and sends between processes
  //     (necessary)
  //
  // Teuchos::Comm requires that the data for a send are contiguous
  // in a send buffer.  Therefore, if the data in the send buffer
  // for doPosts() are not contiguous, they will need to be copied
  // into a contiguous buffer.  The user has specified this
  // noncontiguous pattern and we can't do anything about it.
  // However, if they do not provide an efficient pattern, we will
  // warn them if one of the following compile-time options has been
  // set:
  //   * HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS
  //
  // If the data are contiguous, then we can post the sends in situ
  // (i.e., without needing to copy them into a send buffer).
  //
  // Determine contiguity. There are a number of ways to do this:
  // * If the export IDs are sorted, then all exports to a
  //   particular proc must be contiguous. This is what Epetra does.
  // * If the export ID of the current export already has been
  //   listed, then the previous listing should correspond to the
  //   same export.  This tests contiguity, but not sortedness.
  //
  // Both of these tests require O(n), where n is the number of
  // exports. However, the latter will positively identify a greater
  // portion of contiguous patterns.  We use the latter method.
  //
  // Check to see if values are grouped by procs without gaps
  // If so, indices_to -> 0.

  if (debug) {
    // Test whether any process in the communicator got an invalid
    // process ID.  If badID != -1 on this process, then it equals
    // this process' rank.  The max of all badID over all processes
    // is the max rank which has an invalid process ID.
    int badID = -1;
    for (size_t i = 0; i < numExports; ++i) {
      const int exportID = exportProcIDs[i];
      if (exportID >= numProcs || exportID < 0) {
        badID = myProcID;
        break;
      }
    }
    int gbl_badID;
    reduceAll<int, int> (*comm_, REDUCE_MAX, badID, outArg (gbl_badID));
    TEUCHOS_TEST_FOR_EXCEPTION
      (gbl_badID >= 0, std::runtime_error, rawPrefix << "Proc "
        << gbl_badID << ", perhaps among other processes, got a bad "
        "send process ID.");
  }

  // Set up data structures for quick traversal of arrays.
  // This contains the number of sends for each process ID.
  //
  // FIXME (mfh 20 Mar 2014) This is one of a few places in Tpetra
  // that create an array of length the number of processes in the
  // communicator (plus one).  Given how this code uses this array,
  // it should be straightforward to replace it with a hash table or
  // some other more space-efficient data structure.  In practice,
  // most of the entries of starts should be zero for a sufficiently
  // large process count, unless the communication pattern is dense.
  // Note that it's important to be able to iterate through keys (i
  // for which starts[i] is nonzero) in increasing order.
  Teuchos::Array<size_t> starts (numProcs + 1, 0);

  // numActive is the number of sends that are not Null
  size_t numActive = 0;
  int needSendBuff = 0; // Boolean
  
  for (size_t i = 0; i < numExports; ++i) {
    const int exportID = exportProcIDs[i];
    if (exportID >= 0) {
      // exportID is a valid process ID.  Increment the number of
      // messages this process will send to that process.
      ++starts[exportID];

      // If we're sending more than one message to process exportID,
      // then it is possible that the data are not contiguous.
      // Check by seeing if the previous process ID in the list
      // (exportProcIDs[i-1]) is the same.  It's safe to use i-1,
      // because if starts[exportID] > 1, then i must be > 1 (since
      // the starts array was filled with zeros initially).

      // null entries break continuity.
      // e.g.,  [ 0, 0, 0, 1, -99, 1, 2, 2, 2] is not contiguous
      if (needSendBuff == 0 && starts[exportID] > 1 &&
          exportID != exportProcIDs[i-1]) {
        needSendBuff = 1;
      }
      ++numActive;
    }
  }

#if defined(HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS)
  {
    int global_needSendBuff;
    reduceAll<int, int> (*comm_, REDUCE_MAX, needSendBuff,
        outArg (global_needSendBuff));
    TPETRA_EFFICIENCY_WARNING(
        global_needSendBuff != 0,
        "::createFromSends: Grouping export IDs together by process rank often "
        "improves performance.");
  }
#endif

  // Determine from the caller's data whether or not the current
  // process should send (a) message(s) to itself.
  if (starts[myProcID] != 0) {
    sendMessageToSelf_ = true;
  }
  else {
    sendMessageToSelf_ = false;
  }

  if (! needSendBuff) {
    // grouped by proc, no send buffer or indicesTo_ needed
    numSendsToOtherProcs_ = 0;
    // Count total number of sends, i.e., total number of procs to
    // which we are sending.  This includes myself, if applicable.
    for (int i = 0; i < numProcs; ++i) {
      if (starts[i]) {
        ++numSendsToOtherProcs_;
      }
    }

    // Not only do we not need these, but we must clear them, as
    // empty status of indicesTo is a flag used later.
    indicesTo_.resize(0);
    // Size these to numSendsToOtherProcs_; note, at the moment, numSendsToOtherProcs_
    // includes self sends.  Set their values to zeros.
    procIdsToSendTo_.assign(numSendsToOtherProcs_,0);
    startsTo_.assign(numSendsToOtherProcs_,0);
    lengthsTo_.assign(numSendsToOtherProcs_,0);

    // set startsTo to the offset for each send (i.e., each proc ID)
    // set procsTo to the proc ID for each send
    // in interpreting this code, remember that we are assuming contiguity
    // that is why index skips through the ranks
    {
      size_t procIndex = 0;
      for (size_t i = 0; i < numSendsToOtherProcs_; ++i) {
        while (exportProcIDs[procIndex] < 0) {
          ++procIndex; // skip all negative proc IDs
        }
        startsTo_[i] = procIndex;
        int procID = exportProcIDs[procIndex];
        procIdsToSendTo_[i] = procID;
        procIndex += starts[procID];
      }
    }
    // sort the startsTo and proc IDs together, in ascending order, according
    // to proc IDs
    if (numSendsToOtherProcs_ > 0) {
      sort2(procIdsToSendTo_.begin(), procIdsToSendTo_.end(), startsTo_.begin());
    }
    // compute the maximum send length
    maxSendLength_ = 0;
    for (size_t i = 0; i < numSendsToOtherProcs_; ++i) {
      int procID = procIdsToSendTo_[i];
      lengthsTo_[i] = starts[procID];
      if ((procID != myProcID) && (lengthsTo_[i] > maxSendLength_)) {
        maxSendLength_ = lengthsTo_[i];
      }
    }
  }
  else {
    // not grouped by proc, need send buffer and indicesTo_

    // starts[i] is the number of sends to proc i
    // numActive equals number of sends total, \sum_i starts[i]

    // this loop starts at starts[1], so explicitly check starts[0]
    if (starts[0] == 0 ) {
      numSendsToOtherProcs_ = 0;
    }
    else {
      numSendsToOtherProcs_ = 1;
    }
    for (Teuchos::Array<size_t>::iterator i=starts.begin()+1,
        im1=starts.begin();
        i != starts.end(); ++i)
    {
      if (*i != 0) ++numSendsToOtherProcs_;
      *i += *im1;
      im1 = i;
    }
    // starts[i] now contains the number of exports to procs 0 through i

    for (Teuchos::Array<size_t>::reverse_iterator ip1=starts.rbegin(),
        i=starts.rbegin()+1;
        i != starts.rend(); ++i)
    {
      *ip1 = *i;
      ip1 = i;
    }
    starts[0] = 0;
    // starts[i] now contains the number of exports to procs 0 through
    // i-1, i.e., all procs before proc i

    indicesTo_.resize(numActive);

    for (size_t i = 0; i < numExports; ++i) {
      if (exportProcIDs[i] >= 0) {
        // record the offset to the sendBuffer for this export
        indicesTo_[starts[exportProcIDs[i]]] = i;
        // now increment the offset for this proc
        ++starts[exportProcIDs[i]];
      }
    }
    // our send buffer will contain the export data for each of the procs
    // we communicate with, in order by proc id
    // sendBuffer = {proc_0_data, proc_1_data, ..., proc_np-1_data}
    // indicesTo now maps each export to the location in our send buffer
    // associated with the export
    // data for export i located at sendBuffer[indicesTo[i]]
    //
    // starts[i] once again contains the number of exports to
    // procs 0 through i
    for (int proc = numProcs-1; proc != 0; --proc) {
      starts[proc] = starts[proc-1];
    }
    starts.front() = 0;
    starts[numProcs] = numActive;
    //
    // starts[proc] once again contains the number of exports to
    // procs 0 through proc-1
    // i.e., the start of my data in the sendBuffer

    // this contains invalid data at procs we don't care about, that is okay
    procIdsToSendTo_.resize(numSendsToOtherProcs_);
    startsTo_.resize(numSendsToOtherProcs_);
    lengthsTo_.resize(numSendsToOtherProcs_);

    // for each group of sends/exports, record the destination proc,
    // the length, and the offset for this send into the
    // send buffer (startsTo_)
    maxSendLength_ = 0;
    size_t snd = 0;
    for (int proc = 0; proc < numProcs; ++proc ) {
      if (starts[proc+1] != starts[proc]) {
        lengthsTo_[snd] = starts[proc+1] - starts[proc];
        startsTo_[snd] = starts[proc];
        // record max length for all off-proc sends
        if ((proc != myProcID) && (lengthsTo_[snd] > maxSendLength_)) {
          maxSendLength_ = lengthsTo_[snd];
        }
        procIdsToSendTo_[snd] = proc;
        ++snd;
      }
    }
  }

  if (sendMessageToSelf_) {
    --numSendsToOtherProcs_;
  }

  // Invert map to see what msgs are received and what length
  computeReceives();

#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  initializeMpiAdvance();
#endif

  // createFromRecvs() calls createFromSends(), but will set
  // howInitialized_ again after calling createFromSends().
  howInitialized_ = Details::DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_SENDS;

  return totalReceiveLength_;
}

void DistributorPlan::createFromRecvs(const Teuchos::ArrayView<const int>& remoteProcIDs)
{
  createFromSends(remoteProcIDs);
  *this = *getReversePlan();
  howInitialized_ = Details::DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_RECVS;
}

void DistributorPlan::createFromSendsAndRecvs(const Teuchos::ArrayView<const int>& exportProcIDs,
                                              const Teuchos::ArrayView<const int>& remoteProcIDs)
{
  // note the exportProcIDs and remoteProcIDs _must_ be a list that has
  // an entry for each GID. If the export/remoteProcIDs is taken from
  // the getProcs{From|To} lists that are extracted from a previous distributor,
  // it will generate a wrong answer, because those lists have a unique entry
  // for each processor id. A version of this with lengthsTo and lengthsFrom
  // should be made.

  howInitialized_ = Tpetra::Details::DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_SENDS_N_RECVS;

  int myProcID = comm_->getRank ();
  int numProcs = comm_->getSize();

  const size_t numExportIDs = exportProcIDs.size();
  Teuchos::Array<size_t> starts (numProcs + 1, 0);

  size_t numActive = 0;
  int needSendBuff = 0; // Boolean

  for(size_t i = 0; i < numExportIDs; i++ )
  {
    if( needSendBuff==0 && i && (exportProcIDs[i] < exportProcIDs[i-1]) )
      needSendBuff = 1;
    if( exportProcIDs[i] >= 0 )
    {
      ++starts[ exportProcIDs[i] ];
      ++numActive;
    }
  }

  sendMessageToSelf_ = ( starts[myProcID] != 0 ) ? 1 : 0;

  numSendsToOtherProcs_ = 0;

  if( needSendBuff ) //grouped by processor, no send buffer or indicesTo_ needed
  {
    if (starts[0] == 0 ) {
      numSendsToOtherProcs_ = 0;
    }
    else {
      numSendsToOtherProcs_ = 1;
    }
    for (Teuchos::Array<size_t>::iterator i=starts.begin()+1,
        im1=starts.begin();
        i != starts.end(); ++i)
    {
      if (*i != 0) ++numSendsToOtherProcs_;
      *i += *im1;
      im1 = i;
    }
    // starts[i] now contains the number of exports to procs 0 through i

    for (Teuchos::Array<size_t>::reverse_iterator ip1=starts.rbegin(),
        i=starts.rbegin()+1;
        i != starts.rend(); ++i)
    {
      *ip1 = *i;
      ip1 = i;
    }
    starts[0] = 0;
    // starts[i] now contains the number of exports to procs 0 through
    // i-1, i.e., all procs before proc i

    indicesTo_.resize(numActive);

    for (size_t i = 0; i < numExportIDs; ++i) {
      if (exportProcIDs[i] >= 0) {
        // record the offset to the sendBuffer for this export
        indicesTo_[starts[exportProcIDs[i]]] = i;
        // now increment the offset for this proc
        ++starts[exportProcIDs[i]];
      }
    }
    for (int proc = numProcs-1; proc != 0; --proc) {
      starts[proc] = starts[proc-1];
    }
    starts.front() = 0;
    starts[numProcs] = numActive;
    procIdsToSendTo_.resize(numSendsToOtherProcs_);
    startsTo_.resize(numSendsToOtherProcs_);
    lengthsTo_.resize(numSendsToOtherProcs_);
    maxSendLength_ = 0;
    size_t snd = 0;
    for (int proc = 0; proc < numProcs; ++proc ) {
      if (starts[proc+1] != starts[proc]) {
        lengthsTo_[snd] = starts[proc+1] - starts[proc];
        startsTo_[snd] = starts[proc];
        // record max length for all off-proc sends
        if ((proc != myProcID) && (lengthsTo_[snd] > maxSendLength_)) {
          maxSendLength_ = lengthsTo_[snd];
        }
        procIdsToSendTo_[snd] = proc;
        ++snd;
      }
    }
  }
  else {
    // grouped by proc, no send buffer or indicesTo_ needed
    numSendsToOtherProcs_ = 0;
    // Count total number of sends, i.e., total number of procs to
    // which we are sending.  This includes myself, if applicable.
    for (int i = 0; i < numProcs; ++i) {
      if (starts[i]) {
        ++numSendsToOtherProcs_;
      }
    }

    // Not only do we not need these, but we must clear them, as
    // empty status of indicesTo is a flag used later.
    indicesTo_.resize(0);
    // Size these to numSendsToOtherProcs_; note, at the moment, numSendsToOtherProcs_
    // includes self sends.  Set their values to zeros.
    procIdsToSendTo_.assign(numSendsToOtherProcs_,0);
    startsTo_.assign(numSendsToOtherProcs_,0);
    lengthsTo_.assign(numSendsToOtherProcs_,0);

    // set startsTo to the offset for each send (i.e., each proc ID)
    // set procsTo to the proc ID for each send
    // in interpreting this code, remember that we are assuming contiguity
    // that is why index skips through the ranks
    {
      size_t procIndex = 0;
      for (size_t i = 0; i < numSendsToOtherProcs_; ++i) {
        while (exportProcIDs[procIndex] < 0) {
          ++procIndex; // skip all negative proc IDs
        }
        startsTo_[i] = procIndex;
        int procID = exportProcIDs[procIndex];
        procIdsToSendTo_[i] = procID;
        procIndex += starts[procID];
      }
    }
    // sort the startsTo and proc IDs together, in ascending order, according
    // to proc IDs
    if (numSendsToOtherProcs_ > 0) {
      sort2(procIdsToSendTo_.begin(), procIdsToSendTo_.end(), startsTo_.begin());
    }
    // compute the maximum send length
    maxSendLength_ = 0;
    for (size_t i = 0; i < numSendsToOtherProcs_; ++i) {
      int procID = procIdsToSendTo_[i];
      lengthsTo_[i] = starts[procID];
      if ((procID != myProcID) && (lengthsTo_[i] > maxSendLength_)) {
        maxSendLength_ = lengthsTo_[i];
      }
    }
  }


  numSendsToOtherProcs_ -= sendMessageToSelf_;
  std::vector<int> recv_list;
  recv_list.reserve(numSendsToOtherProcs_); //reserve an initial guess for size needed

  int last_pid=-2;
  for(int i=0; i<remoteProcIDs.size(); i++) {
    if(remoteProcIDs[i]>last_pid) {
      recv_list.push_back(remoteProcIDs[i]);
      last_pid = remoteProcIDs[i];
    }
    else if (remoteProcIDs[i]<last_pid)
      throw std::runtime_error("Tpetra::Distributor:::createFromSendsAndRecvs expected RemotePIDs to be in sorted order");
  }
  numReceives_ = recv_list.size();
  if(numReceives_) {
    procsFrom_.assign(numReceives_,0);
    lengthsFrom_.assign(numReceives_,0);
    indicesFrom_.assign(numReceives_,0);
    startsFrom_.assign(numReceives_,0);
  }
  for(size_t i=0,j=0; i<numReceives_; ++i) {
    int jlast=j;
    procsFrom_[i]  = recv_list[i];
    startsFrom_[i] = j;
    for( ; j<(size_t)remoteProcIDs.size() &&
        remoteProcIDs[jlast]==remoteProcIDs[j]  ; j++){;}
    lengthsFrom_[i] = j-jlast;
  }
  totalReceiveLength_ = remoteProcIDs.size();
  indicesFrom_.clear ();
  numReceives_-=sendMessageToSelf_;
  
#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  initializeMpiAdvance();
#endif
}

Teuchos::RCP<DistributorPlan> DistributorPlan::getReversePlan() const {
  if (reversePlan_.is_null()) createReversePlan();
  return reversePlan_;
}

void DistributorPlan::createReversePlan() const
{
  reversePlan_ = Teuchos::rcp(new DistributorPlan(comm_));
  reversePlan_->howInitialized_ = Details::DISTRIBUTOR_INITIALIZED_BY_REVERSE;
  reversePlan_->sendType_ = sendType_;

  // The total length of all the sends of this DistributorPlan.  We
  // calculate it because it's the total length of all the receives
  // of the reverse DistributorPlan.
  size_t totalSendLength =
    std::accumulate(lengthsTo_.begin(), lengthsTo_.end(), 0);

  // The maximum length of any of the receives of this DistributorPlan.
  // We calculate it because it's the maximum length of any of the
  // sends of the reverse DistributorPlan.
  size_t maxReceiveLength = 0;
  const int myProcID = comm_->getRank();
  for (size_t i=0; i < numReceives_; ++i) {
    if (procsFrom_[i] != myProcID) {
      // Don't count receives for messages sent by myself to myself.
      if (lengthsFrom_[i] > maxReceiveLength) {
        maxReceiveLength = lengthsFrom_[i];
      }
    }
  }

  reversePlan_->sendMessageToSelf_ = sendMessageToSelf_;
  reversePlan_->numSendsToOtherProcs_ = numReceives_;
  reversePlan_->procIdsToSendTo_ = procsFrom_;
  reversePlan_->startsTo_ = startsFrom_;
  reversePlan_->lengthsTo_ = lengthsFrom_;
  reversePlan_->maxSendLength_ = maxReceiveLength;
  reversePlan_->indicesTo_ = indicesFrom_;
  reversePlan_->numReceives_ = numSendsToOtherProcs_;
  reversePlan_->totalReceiveLength_ = totalSendLength;
  reversePlan_->lengthsFrom_ = lengthsTo_;
  reversePlan_->procsFrom_ = procIdsToSendTo_;
  reversePlan_->startsFrom_ = startsTo_;
  reversePlan_->indicesFrom_ = indicesTo_;

#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  // is there a smarter way to do this
  reversePlan_->initializeMpiAdvance();
#endif
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

  const int mpiTag = DEFAULT_MPI_TAG;

  // toProcsFromMe[i] == the number of messages sent by this process
  // to process i.  The data in numSendsToOtherProcs_, procIdsToSendTo_, and lengthsTo_
  // concern the contiguous sends.  Therefore, each process will be
  // listed in procIdsToSendTo_ at most once, and so toProcsFromMe[i] will
  // either be 0 or 1.
  {
    Array<int> toProcsFromMe (numProcs, 0);
#ifdef HAVE_TPETRA_DEBUG
    bool counting_error = false;
#endif // HAVE_TPETRA_DEBUG
    for (size_t i = 0; i < (numSendsToOtherProcs_ + (sendMessageToSelf_ ? 1 : 0)); ++i) {
#ifdef HAVE_TPETRA_DEBUG
      if (toProcsFromMe[procIdsToSendTo_[i]] != 0) {
        counting_error = true;
      }
#endif // HAVE_TPETRA_DEBUG
      toProcsFromMe[procIdsToSendTo_[i]] = 1;
    }
#ifdef HAVE_TPETRA_DEBUG
    // Note that SHARED_TEST_FOR_EXCEPTION does a global reduction
    SHARED_TEST_FOR_EXCEPTION(counting_error, std::logic_error,
        "Tpetra::Distributor::computeReceives: There was an error on at least "
        "one process in counting the number of messages send by that process to "
        "the other processs.  Please report this bug to the Tpetra developers.",
        *comm_);
#endif // HAVE_TPETRA_DEBUG

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
        mpiTag, *comm_);
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
      send<int, size_t> (lengthsTo_i, 1, as<int> (procIdsToSendTo_[i]), mpiTag, *comm_);
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

void DistributorPlan::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& plist)
{
  using Teuchos::FancyOStream;
  using Teuchos::getIntegralValue;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using std::endl;

  if (! plist.is_null()) {
    RCP<const ParameterList> validParams = getValidParameters ();
    plist->validateParametersAndSetDefaults (*validParams);

    const Details::EDistributorSendType sendType =
      getIntegralValue<Details::EDistributorSendType> (*plist, "Send type");

    // Now that we've validated the input list, save the results.
    sendType_ = sendType;

    // ParameterListAcceptor semantics require pointer identity of the
    // sublist passed to setParameterList(), so we save the pointer.
    this->setMyParamList (plist);
  }
}

Teuchos::Array<std::string> distributorSendTypes()
{
  Teuchos::Array<std::string> sendTypes;
  sendTypes.push_back ("Isend");
  sendTypes.push_back ("Send");
  sendTypes.push_back ("Alltoall");
#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  sendTypes.push_back ("MpiAdvanceAlltoall");
  sendTypes.push_back ("MpiAdvanceNbralltoallv");
#endif
  return sendTypes;
}

Teuchos::RCP<const Teuchos::ParameterList>
DistributorPlan::getValidParameters() const
{
  using Teuchos::Array;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::setStringToIntegralParameter;

  Array<std::string> sendTypes = distributorSendTypes ();
  const std::string defaultSendType ("Send");
  Array<Details::EDistributorSendType> sendTypeEnums;
  sendTypeEnums.push_back (Details::DISTRIBUTOR_ISEND);
  sendTypeEnums.push_back (Details::DISTRIBUTOR_SEND);
  sendTypeEnums.push_back (Details::DISTRIBUTOR_ALLTOALL);
#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  sendTypeEnums.push_back (Details::DISTRIBUTOR_MPIADVANCE_ALLTOALL);
  sendTypeEnums.push_back (Details::DISTRIBUTOR_MPIADVANCE_NBRALLTOALLV);
#endif

  RCP<ParameterList> plist = parameterList ("Tpetra::Distributor");

  setStringToIntegralParameter<Details::EDistributorSendType> ("Send type",
      defaultSendType, "When using MPI, the variant of send to use in "
      "do[Reverse]Posts()", sendTypes(), sendTypeEnums(), plist.getRawPtr());
  plist->set ("Timer Label","","Label for Time Monitor output");

  return Teuchos::rcp_const_cast<const ParameterList> (plist);
}

#if defined(HAVE_TPETRACORE_MPI_ADVANCE)

// Used by Teuchos::RCP to clean up an owned MPIX_Comm*
struct MpixCommDeallocator {
  void free(MPIX_Comm **comm) const {
    MPIX_Comm_free(*comm);
  }
};

void DistributorPlan::initializeMpiAdvance() {

  // assert the mpix communicator is null. if this is not the case we will figure out why
  TEUCHOS_ASSERT(mpixComm_.is_null());

  // use the members to initialize the graph for neightborhood mode, or just the MPIX communicator for non-neighborhood mode
  Teuchos::RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm_);
  Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawComm = mpiComm->getRawMpiComm();
  int err = 0;
  if (sendType_ == DISTRIBUTOR_MPIADVANCE_ALLTOALL) {
    MPIX_Comm **mpixComm = new(MPIX_Comm*);
    err = MPIX_Comm_init(mpixComm, (*rawComm)());
    mpixComm_ = Teuchos::RCP(mpixComm,
      MpixCommDeallocator(),
      true /*take ownership*/
    );
  }
  else if (sendType_ == DISTRIBUTOR_MPIADVANCE_NBRALLTOALLV) {
    int numRecvs = (int)(numReceives_ + (sendMessageToSelf_ ? 1 : 0));
    int *sourceRanks = procsFrom_.data();
    
    // int *sourceWeights = static_cast<int*>(lengthsFrom_.data());// lengthsFrom_ may not be int
    const int *sourceWeights = MPI_UNWEIGHTED;
    int numSends = (int)(numSendsToOtherProcs_ + (sendMessageToSelf_ ? 1 : 0));
    int *destRanks = procIdsToSendTo_.data();

    // int *destWeights = static_cast<int*>(lengthsTo_.data()); // lengthsTo_ may not be int
    const int *destWeights = MPI_UNWEIGHTED; // lengthsTo_ may not be int

    MPIX_Comm **mpixComm = new(MPIX_Comm*);
    err = MPIX_Dist_graph_create_adjacent((*rawComm)(), numRecvs, sourceRanks, sourceWeights, numSends, destRanks, destWeights, MPI_INFO_NULL, false, mpixComm);
    mpixComm_ = Teuchos::RCP(mpixComm,
      MpixCommDeallocator(),
      true /*take ownership*/
    );
  }

  TEUCHOS_ASSERT(err == 0);
}
#endif

}
}
