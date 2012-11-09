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

#include "Tpetra_Distributor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

namespace Tpetra {

  Array<std::string>
  distributorSendTypes ()
  {
    Array<std::string> sendTypes;
    sendTypes.push_back ("Isend");
    sendTypes.push_back ("Rsend");
    sendTypes.push_back ("Send");
    sendTypes.push_back ("Ssend");
    return sendTypes;
  }

  Distributor::Distributor (const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
    : comm_(comm)
    , sendType_ (DISTRIBUTOR_SEND)
    , barrierBetween_ (true)
    , numExports_(0)
    , selfMessage_(false)
    , numSends_(0)
    , maxSendLength_(0)
    , numReceives_(0)
    , totalReceiveLength_(0)
  {
    using Teuchos::getFancyOStream;
    using Teuchos::oblackholestream;
    using Teuchos::rcp;

    // Always start by making sure the Distributor won't print anything.
    this->setVerbLevel (Teuchos::VERB_NONE);
    this->setOStream (getFancyOStream (rcp (new oblackholestream)));
  }

  Distributor::Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                            const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : comm_(comm)
    , sendType_ (DISTRIBUTOR_SEND)
    , barrierBetween_ (true)
    , numExports_(0)
    , selfMessage_(false)
    , numSends_(0)
    , maxSendLength_(0)
    , numReceives_(0)
    , totalReceiveLength_(0)
  {
    TEUCHOS_TEST_FOR_EXCEPTION (plist.is_null(), std::invalid_argument, "The "
      "two-argument Distributor constructor requires that the input "
      "RCP<ParameterList> be nonnull.  If you don't know what parameters to "
      "set, you can either call the one-argument constructor, or supply a "
      "nonnull but empty ParameterList.  Both of these options will set default "
      "parameters.");
    using Teuchos::getFancyOStream;
    using Teuchos::oblackholestream;
    using Teuchos::rcp;

    // Always start by making sure the Distributor won't print anything.
    this->setVerbLevel (Teuchos::VERB_NONE);
    this->setOStream (getFancyOStream (rcp (new oblackholestream)));
    // Setting parameters may override these, if there is a
    // "VerboseObject" sublist.
    this->setParameterList (plist);
  }

  Distributor::Distributor (const Distributor & distributor)
    : comm_(distributor.comm_)
    , sendType_ (distributor.sendType_)
    , barrierBetween_ (distributor.barrierBetween_)
    , numExports_(distributor.numExports_)
    , selfMessage_(distributor.selfMessage_)
    , numSends_(distributor.numSends_)
    , maxSendLength_(distributor.maxSendLength_)
    , numReceives_(distributor.numReceives_)
    , totalReceiveLength_(distributor.totalReceiveLength_)
    , reverseDistributor_(distributor.reverseDistributor_)
  {
    using Teuchos::getFancyOStream;
    using Teuchos::oblackholestream;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

    this->setVerbLevel (Teuchos::VERB_NONE);
    this->setOStream (getFancyOStream (rcp (new oblackholestream)));

    // Clone the right-hand side's ParameterList, so that this' list
    // is decoupled from the right-hand side's list.  We don't need to
    // do validation, since the right-hand side already has validated
    // its parameters, so just call setMyParamList().  Note that this
    // won't work if the right-hand side doesn't have a list set yet,
    // so we first check for null.
    RCP<const ParameterList> rhsList = distributor.getParameterList ();
    if (! rhsList.is_null ()) {
      this->setMyParamList (parameterList (* rhsList));
    }
  }

  Distributor::~Distributor()
  {
    // We shouldn't have any outstanding communication requests at
    // this point.
    TEUCHOS_TEST_FOR_EXCEPTION(requests_.size() != 0, std::runtime_error,
      "Tpetra::Distributor: Destructor called with " << requests_.size()
      << " outstanding posts (unfulfilled communication requests).  There "
      "should be none at this point.  Please report this bug to the Tpetra "
      "developers.");
  }

  void
  Distributor::setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    using Teuchos::FancyOStream;
    using Teuchos::getIntegralValue;
    using Teuchos::includesVerbLevel;
    using Teuchos::OSTab;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using std::endl;

    RCP<const ParameterList> validParams = getValidParameters ();
    plist->validateParametersAndSetDefaults (*validParams);

    const bool barrierBetween =
      plist->get<bool> ("Barrier between receives and sends");
    const EDistributorSendType sendType =
      getIntegralValue<EDistributorSendType> (*plist, "Send type");

    // We check this property explicitly, since we haven't yet learned
    // how to make a validator that can cross-check properties.
    // Later, turn this into a validator so that it can be embedded in
    // the valid ParameterList and used in Optika.
    TEUCHOS_TEST_FOR_EXCEPTION(! barrierBetween && sendType == DISTRIBUTOR_RSEND,
      std::invalid_argument, "If you use ready sends, you must include a "
      "barrier between receives and sends.  Ready sends require that their "
      "corresponding receives have already been posted, and the only way to "
      "guarantee that in general is with a barrier.");

    if (plist->isSublist ("VerboseObject")) {
      // Read the "VerboseObject" sublist for (optional) verbosity
      // settings.  We've set defaults already in Distributor's
      // constructor, so we don't need this sublist to exist.
      Teuchos::readVerboseObjectSublist (&*plist, this);
    }

    // Now that we've validated the input list, save the results.
    sendType_ = sendType;
    barrierBetween_ = barrierBetween;

#ifdef HAVE_TEUCHOS_DEBUG
    // Prepare for verbose output, if applicable.
    Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel ();
    RCP<FancyOStream> out = this->getOStream ();
    const int myRank = comm_->getRank ();
    // We only want one process to print verbose output here.
    const bool doPrint = out.get () && (myRank == 0) &&
      includesVerbLevel (verbLevel, Teuchos::VERB_EXTREME, true);

    if (doPrint) {
      *out << "Distributor::setParameterList" << endl;
      OSTab tab = this->getOSTab(); // Add one tab level
      *out << "sendType_=" << DistributorSendTypeEnumToString (sendType_)
           << ", barrierBetween_=" << barrierBetween_ << endl;
    }
#endif // HAVE_TEUCHOS_DEBUG

    // ParameterListAcceptor semantics require pointer identity of the
    // sublist passed to setParameterList(), so we save the pointer.
    this->setMyParamList (plist);
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  Distributor::getValidParameters () const
  {
    using Teuchos::Array;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::setStringToIntegralParameter;

    const bool barrierBetween = false;

    Array<std::string> sendTypes = distributorSendTypes ();
    const std::string defaultSendType ("Send");
    Array<EDistributorSendType> sendTypeEnums;
    sendTypeEnums.push_back (DISTRIBUTOR_ISEND);
    sendTypeEnums.push_back (DISTRIBUTOR_RSEND);
    sendTypeEnums.push_back (DISTRIBUTOR_SEND);
    sendTypeEnums.push_back (DISTRIBUTOR_SSEND);

    RCP<ParameterList> plist = parameterList ("Tpetra::Distributor");
    plist->set ("Barrier between receives and sends", barrierBetween,
                "Whether to execute a barrier between receives and sends in do"
                "[Reverse]Posts().  Required for correctness when \"Send type\""
                "=\"Rsend\", otherwise correct but not recommended.");
    setStringToIntegralParameter<EDistributorSendType> ("Send type",
      defaultSendType, "When using MPI, the variant of MPI_Send to use in "
      "do[Reverse]Posts()", sendTypes(), sendTypeEnums(), plist.getRawPtr());

    Teuchos::setupVerboseObjectSublist (&*plist);
    return Teuchos::rcp_const_cast<const ParameterList> (plist);
  }


  size_t Distributor::getTotalReceiveLength() const
  { return totalReceiveLength_; }

  size_t Distributor::getNumReceives() const
  { return numReceives_; }

  bool Distributor::hasSelfMessage() const
  { return selfMessage_; }

  size_t Distributor::getNumSends() const
  { return numSends_; }

  size_t Distributor::getMaxSendLength() const
  { return maxSendLength_; }

  Teuchos::ArrayView<const int> Distributor::getImagesFrom() const
  { return imagesFrom_; }

  Teuchos::ArrayView<const size_t> Distributor::getLengthsFrom() const
  { return lengthsFrom_; }

  Teuchos::ArrayView<const int> Distributor::getImagesTo() const
  { return imagesTo_; }

  Teuchos::ArrayView<const size_t> Distributor::getLengthsTo() const
  { return lengthsTo_; }

  const Teuchos::RCP<Distributor> &
  Distributor::getReverse() const {
    if (reverseDistributor_ == Teuchos::null) {
      // need to create reverse distributor
      createReverseDistributor();
    }
    return reverseDistributor_;
  }


  void Distributor::createReverseDistributor() const {

    reverseDistributor_ = Teuchos::rcp(new Distributor(comm_));

    // The total length of all the sends of this Distributor.  We
    // calculate it because it's the total length of all the receives
    // of the reverse Distributor.
    size_t totalSendLength = std::accumulate(lengthsTo_.begin(),lengthsTo_.end(),0);

    // The maximum length of any of the receives of this Distributor.
    // We calculate it because it's the maximum length of any of the
    // sends of the reverse Distributor.
    size_t maxReceiveLength = 0;
    const int myImageID = comm_->getRank();
    for (size_t i=0; i < numReceives_; ++i) {
      if (imagesFrom_[i] != myImageID) {
        // Don't count receives for messages sent by myself to myself.
        if (lengthsFrom_[i] > maxReceiveLength) {
          maxReceiveLength = lengthsFrom_[i];
        }
      }
    }

    // Initialize all of reverseDistributor's data members.  This
    // mainly just involves flipping "send" and "receive," or the
    // equivalent "to" and "from."
    reverseDistributor_->lengthsTo_ = lengthsFrom_;
    reverseDistributor_->imagesTo_ = imagesFrom_;
    reverseDistributor_->indicesTo_ = indicesFrom_;
    reverseDistributor_->startsTo_ = startsFrom_;
    reverseDistributor_->lengthsFrom_ = lengthsTo_;
    reverseDistributor_->imagesFrom_ = imagesTo_;
    reverseDistributor_->indicesFrom_ = indicesTo_;
    reverseDistributor_->startsFrom_ = startsTo_;
    reverseDistributor_->numSends_ = numReceives_;
    reverseDistributor_->numReceives_ = numSends_;
    reverseDistributor_->selfMessage_ = selfMessage_;
    reverseDistributor_->maxSendLength_ = maxReceiveLength;
    reverseDistributor_->totalReceiveLength_ = totalSendLength;
    // Note: technically, I am my reverse distributor's reverse distributor, but
    //       we will not set this up, as it gives us an opportunity to test
    //       that reverseDistributor is an inverse operation w.r.t. value semantics of distributors
    // Note: numExports_ was not copied
  }


  void Distributor::doWaits() {
    using Teuchos::FancyOStream;
    using Teuchos::includesVerbLevel;
    using Teuchos::is_null;
    using Teuchos::OSTab;
    using Teuchos::RCP;
    using Teuchos::waitAll;
    using std::endl;

#ifdef HAVE_TEUCHOS_DEBUG
    // Prepare for verbose output, if applicable.
    Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel ();
    RCP<FancyOStream> out = this->getOStream ();
    const bool doPrint = out.get () &&
      includesVerbLevel (verbLevel, Teuchos::VERB_EXTREME, true);
    const int myRank = comm_->getRank ();

    if (doPrint) {
      *out << "Distributor::doWaits (Proc " << myRank << "):" << endl;
    }
    OSTab tab = this->getOSTab(); // Add one tab level
#endif // HAVE_TEUCHOS_DEBUG

    if (requests_.size() > 0) {
      waitAll (*comm_, requests_());

#ifdef HAVE_TEUCHOS_DEBUG
      if (doPrint) {
        *out << "Proc " << myRank << ": waitAll completed " << requests_.size()
             << " requests" << endl;
      }

      // Make sure that waitAll() nulled out all the requests.
      using Teuchos::Array;
      using Teuchos::CommRequest;
      using Teuchos::RCP;
      for (Array<RCP<CommRequest<int> > >::const_iterator it = requests_.begin();
           it != requests_.end(); ++it)
      {
        TEUCHOS_TEST_FOR_EXCEPTION( ! is_null (*it), std::runtime_error,
          Teuchos::typeName(*this) << "::doWaits(): Communication requests "
          "should all be null aftr calling Teuchos::waitAll() on them, but "
          "at least one request is not null.");
      }
#endif // HAVE_TEUCHOS_DEBUG
      // Restore the invariant that requests_.size() is the number of
      // outstanding nonblocking communication requests.
      requests_.resize (0);
    }
  }


  void Distributor::doReverseWaits()
  {
    // call doWaits() on the reverse Distributor, if it exists
    if (! reverseDistributor_.is_null()) {
      reverseDistributor_->doWaits();
    }
  }

  std::string Distributor::description() const
  {
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    return oss.str();
  }

  void Distributor::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
  {
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;
    const int myImageID = comm_->getRank();
    const int numImages = comm_->getSize();
    Teuchos::OSTab tab(out);

    if (vl == VERB_NONE) {
      return;
    } else {
      if (myImageID == 0) {
        // VERB_LOW and higher prints description() (on Proc 0 only).
        out << this->description() << endl;
      }
      if (vl == VERB_LOW) {
        return;
      } else {
        // vl > VERB_LOW lets each image print its data.  We assume
        // that all images can print to the given output stream, and
        // execute barriers to make it more likely that the output
        // will be in the right order.
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << "[Node " << myImageID << " of " << numImages << "]" << endl;
            out << " selfMessage: " << hasSelfMessage() << endl;
            out << " numSends: " << getNumSends() << endl;
            if (vl == VERB_HIGH || vl == VERB_EXTREME) {
              out << " imagesTo: " << toString(imagesTo_) << endl;
              out << " lengthsTo: " << toString(lengthsTo_) << endl;
              out << " maxSendLength: " << getMaxSendLength() << endl;
            }
            if (vl == VERB_EXTREME) {
              out << " startsTo: " << toString(startsTo_) << endl;
              out << " indicesTo: " << toString(indicesTo_) << endl;
            }
            if (vl == VERB_HIGH || vl == VERB_EXTREME) {
              out << " numReceives: " << getNumReceives() << endl;
              out << " totalReceiveLength: " << getTotalReceiveLength() << endl;
              out << " lengthsFrom: " << toString(lengthsFrom_) << endl;
              out << " startsFrom: " << toString(startsFrom_) << endl;
              out << " imagesFrom: " << toString(imagesFrom_) << endl;
            }
            // Last output is a flush; it leaves a space and also
            // helps synchronize output.
            out << std::flush;
          } // if it's my image's turn to print
          // Execute barriers to give output time to synchronize.
          // One barrier generally isn't enough.
          comm_->barrier();
          comm_->barrier();
          comm_->barrier();
        } // for each image
      }
    }
  }

  void
  Distributor::computeReceives ()
  {
    using Teuchos::Array;
    using Teuchos::CommStatus;
    using Teuchos::CommRequest;
    using Teuchos::ireceive;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_SUM;
    using Teuchos::receive;
    using Teuchos::reduceAllAndScatter;
    using Teuchos::send;
    using Teuchos::waitAll;

    const int myRank = comm_->getRank();
    const int numProcs = comm_->getSize();

    // toNodesFromMe[i] == the number of messages sent by this process
    // to process i.  The data in numSends_, imagesTo_, lengthsTo_
    // concern the contiguous sends.  Therefore, each process will be
    // listed in imagesTo_ at most once.
    {
      Array<size_t> toNodesFromMe (numProcs,0);
#ifdef HAVE_TEUCHOS_DEBUG
      bool counting_error = false;
#endif // HAVE_TEUCHOS_DEBUG
      for (size_t i = 0; i < (numSends_ + (selfMessage_ ? 1 : 0)); ++i) {
#ifdef HAVE_TEUCHOS_DEBUG
        if (toNodesFromMe[imagesTo_[i]] != 0) {
          counting_error = true;
        }
#endif // HAVE_TEUCHOS_DEBUG
        toNodesFromMe[imagesTo_[i]] = 1;
      }
#ifdef HAVE_TEUCHOS_DEBUG
      SHARED_TEST_FOR_EXCEPTION(counting_error, std::logic_error,
        "Tpetra::Distributor::computeReceives: There was an error on at least "
        "one process in counting the number of messages send by that process to "
        "the other processs.  Please report this bug to the Tpetra developers.",
        *comm_);
#endif // HAVE_TEUCHOS_DEBUG

      // Each process will get back only one item (hence, counts =
      // ones) from the array of global sums, namely that entry
      // corresponding to the process, and detailing how many receives
      // it has.  This total includes self sends.
      //
      // mfh 09 Jan 2012: The reduceAllAndScatter really isn't
      // necessary here.  Since counts is just all ones, we could
      // replace this with an all-reduce on toNodesFromMe, and let my
      // process (with rank myRank) get numReceives_ from
      // toNodesFromMe[myRank].  The HPCCG miniapp uses the all-reduce
      // method.  It could be possible that reduceAllAndScatter is
      // faster, but it also makes the code more complicated, and it
      // can't be _asymptotically_ faster (MPI_Allreduce has twice the
      // critical path length of MPI_Reduce, so reduceAllAndScatter
      // can't be more than twice as fast as the all-reduce, even if
      // the scatter is free).
      Array<int> counts (numProcs, 1);
      reduceAllAndScatter (*comm_, REDUCE_SUM, numProcs, &toNodesFromMe[0],
                           &counts[0], &numReceives_);
    }

    // Now we know numReceives_, which is this process' number of
    // receives.  Allocate the lengthsFrom_ and imagesFrom_ arrays
    // with this number of entries.
    lengthsFrom_.assign (numReceives_, 0);
    imagesFrom_.assign (numReceives_, 0);

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
    const size_t actualNumReceives = numReceives_ - (selfMessage_ ? 1 : 0);

    // Teuchos' wrapper for nonblocking receives requires receive
    // buffers that it knows won't go away.  This is why we use RCPs,
    // one RCP per nonblocking receive request.  They get allocated in
    // the loop below.
    Array<RCP<CommRequest<int> > > requests (actualNumReceives);
    Array<RCP<size_t> > lengthsFromBuffers (actualNumReceives);
    Array<RCP<CommStatus<int> > > statuses (actualNumReceives);

    // Teuchos::Comm treats a negative process ID as MPI_ANY_SOURCE
    // (receive data from any process).
    const int anySourceProc = -1;

    for (size_t i = 0; i < actualNumReceives; ++i) {
      lengthsFromBuffers[i] = rcp (new size_t (0));
      // Once the receive completes, we can ask the corresponding
      // CommStatus object (output by wait()) for the sending process'
      // ID (which we'll assign to imagesFrom_[i] -- don't forget to
      // do that!).
      requests[i] = ireceive (*comm_, lengthsFromBuffers[i], anySourceProc);
    }

    // Tell each process to which we are sending how many packets it
    // should expect from us in the communication pattern.
    //
    // We assume that numSends_ and selfMessage_ have already been
    // set.  The value of numSends_ (my process' number of sends) does
    // not include any message that it might send to itself.
    for (size_t i = 0; i < numSends_ + (selfMessage_ ? 1 : 0); ++i) {
      if (imagesTo_[i] != myRank ) {
        // Send a message to imagesTo_[i], telling that process that
        // this communication pattern will send that process
        // lengthsTo_[i] blocks of packets.
        send (*comm_, lengthsTo_[i], imagesTo_[i]);
      }
      else {
        // We don't need a send in the self-message case.  If this
        // process will send a message to itself in the communication
        // pattern, then the last element of lengthsFrom_ and
        // imagesFrom_ corresponds to the self-message.  Of course
        // this process knows how long the message is, and the process
        // ID is its own process ID.
        lengthsFrom_[numReceives_-1] = lengthsTo_[i];
        imagesFrom_[numReceives_-1] = myRank;
      }
    }

    //
    // Wait on all the receives.  When they arrive, check the status
    // output of wait() for the receiving process ID, unpack the
    // request buffers into lengthsFrom_, and set imagesFrom_ from the
    // status.
    //
    waitAll (*comm_, requests (), statuses ());
    for (size_t i = 0; i < actualNumReceives; ++i) {
      lengthsFrom_[i] = *lengthsFromBuffers[i];
      imagesFrom_[i] = statuses[i]->getSourceRank ();
    }

#ifdef HAVE_TEUCHOS_DEBUG
    comm_->barrier();
#endif // HAVE_TEUCHOS_DEBUG

    // Sort the imagesFrom_ array, and apply the same permutation to
    // lengthsFrom_.  This ensures that imagesFrom_[i] and
    // lengthsFrom_[i] refers to the same thing.
    sort2 (imagesFrom_.begin(), imagesFrom_.end(), lengthsFrom_.begin());

    // Compute indicesFrom_
    totalReceiveLength_ = std::accumulate (lengthsFrom_.begin(), lengthsFrom_.end(), 0);
    indicesFrom_.clear ();
    indicesFrom_.reserve (totalReceiveLength_);
    for (size_t i = 0; i < totalReceiveLength_; ++i) {
      indicesFrom_.push_back(i);
    }

    startsFrom_.clear ();
    startsFrom_.reserve (numReceives_);
    for (size_t i = 0, j = 0; i < numReceives_; ++i) {
      startsFrom_.push_back(j);
      j += lengthsFrom_[i];
    }

    if (selfMessage_) {
      --numReceives_;
    }

#ifdef HAVE_TEUCHOS_DEBUG
    comm_->barrier();
#endif // HAVE_TEUCHOS_DEBUG
  }

  size_t
  Distributor::createFromSends (const Teuchos::ArrayView<const int> &exportNodeIDs)
  {
    using Teuchos::outArg;
    numExports_ = exportNodeIDs.size();

    const int myImageID = comm_->getRank();
    const int numImages = comm_->getSize();

    // exportNodeIDs tells us the communication pattern for this
    // distributor.  It dictates the way that the export data will be
    // interpreted in doPosts().  We want to perform at most one
    // communication per node; this is for two reasons:
    //   * minimize latency/overhead in the comm routines (nice)
    //   * match the number of receives and sends between nodes (necessary)
    //
    // Teuchos::Comm requires that the data for a send is contiguous
    // in a send buffer.  Therefore, if the data in the send buffer
    // for doPosts() is not contiguous, it will need to be copied into
    // a contiguous buffer.  The user has specified this noncontiguous
    // pattern and we can't do anything about it.  However, if they do
    // not provide an efficient pattern, we will warn them if one of
    // the following compile-time options has been set:
    //   * HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS
    //   * HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS
    //
    // If the data is contiguous, then we can post the sends in situ
    // (i.e., without needing to copy them into a send buffer).
    //
    // Determine contiguity. There are a number of ways to do this:
    // * If the export IDs are sorted, then all exports to a
    //   particular node must be contiguous. This is what Epetra does.
    // * If the export ID of the current export already has been
    //   listed, then the previous listing should correspond to the
    //   same export.  This tests contiguity, but not sortedness.
    //
    // Both of these tests require O(n), where n is the number of
    // exports. However, the latter will positively identify a greater
    // portion of contiguous patterns. We will use the latter method.
    //
    // Check to see if values are grouped by images without gaps
    // If so, indices_to -> 0.

    // Set up data structures for quick traversal of arrays.
    // This contains the number of sends for each image id.
    Teuchos::Array<size_t> starts (numImages + 1, 0);

    // numActive is the number of sends that are not Null
    size_t numActive = 0;
    char needSendBuff = 0;

    int badID = -1;
    for (size_t i = 0; i < numExports_; ++i) {
      int exportID = exportNodeIDs[i];
      if (exportID >= numImages) {
        badID = myImageID;
        break;
      }
      else if (exportID >= 0) {
        // increment starts[exportID]
        ++starts[exportID];
        // if after incrementing it is greater than one, check that the
        // previous export went to this node
        // this is a safe comparison, because starts[exportID] > 1
        // implies that i > 1.
        // null entries break continuity.
        // e.g.,  [ 0, 0, 0, 1, -99, 1, 2, 2, 2] is not contiguous
        if (needSendBuff==0 && starts[exportID] > 1 && exportID != exportNodeIDs[i-1]) {
          needSendBuff = 1;
        }
        ++numActive;
      }
    }
    {
      int gbl_badID;
      Teuchos::reduceAll(*comm_,Teuchos::REDUCE_MAX,badID,outArg(gbl_badID));
      TEUCHOS_TEST_FOR_EXCEPTION(gbl_badID >= 0, std::runtime_error,
          Teuchos::typeName(*this) << "::createFromSends(): bad node id listed on node " << gbl_badID << ".");
    }

#   if defined(HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS) || defined(HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS)
    {
      char global_needSendBuff;
      Teuchos::reduceAll (*comm_, Teuchos::REDUCE_MAX, needSendBuff, Teuchos::ptr (&global_needSendBuff));
      TPETRA_EFFICIENCY_WARNING(global_needSendBuff,std::runtime_error,
          "::createFromSends(): Grouping export IDs together leads to improved performance.");
    }
#   endif

    // Determine from the caller's data whether or not the current
    // process should send (a) message(s) to itself.
    if (starts[myImageID] != 0) {
      selfMessage_ = true;
    }
    else {
      selfMessage_ = false;
    }


#ifdef HAVE_TEUCHOS_DEBUG
    bool index_neq_numActive = false;
    bool send_neq_numSends = false;
#endif
    if (!needSendBuff) {
      // grouped by image, no send buffer or indicesTo_ needed
      numSends_ = 0;
      // Count total number of sends, i.e., total number of images to
      // which we are sending.  This includes myself, if applicable.
      for (int i=0; i < numImages; ++i) {
        if (starts[i]) ++numSends_;
      }

      // Not only do we not need these, but we must clear them, as
      // empty status of indicesTo is a flag used later.
      indicesTo_.resize(0);
      // Size these to numSends_; note, at the moment, numSends_
      // includes self sends.  Set their values to zeros.
      imagesTo_.assign(numSends_,0);
      startsTo_.assign(numSends_,0);
      lengthsTo_.assign(numSends_,0);

      // set startsTo to the offset for each send (i.e., each image ID)
      // set imagesTo to the image ID for each send
      // in interpreting this code, remember that we are assuming contiguity
      // that is why index skips through the ranks
      {
        size_t index = 0, nodeIndex = 0;
        for (size_t i = 0; i < numSends_; ++i) {
          while (exportNodeIDs[nodeIndex] < 0) {
            ++nodeIndex; // skip all negative node IDs
          }
          startsTo_[i] = nodeIndex;
          int imageID = exportNodeIDs[nodeIndex];
          imagesTo_[i] = imageID;
          index     += starts[imageID];
          nodeIndex += starts[imageID];
        }
#ifdef HAVE_TEUCHOS_DEBUG
        if (index != numActive) {
          index_neq_numActive = true;
        }
#endif
      }
      // sort the startsTo and image IDs together, in ascending order, according
      // to image IDs
      if (numSends_ > 0) {
        sort2(imagesTo_.begin(), imagesTo_.end(), startsTo_.begin());
      }
      // compute the maximum send length
      maxSendLength_ = 0;
      for (size_t i = 0; i < numSends_; ++i) {
        int imageID = imagesTo_[i];
        lengthsTo_[i] = starts[imageID];
        if ((imageID != myImageID) && (lengthsTo_[i] > maxSendLength_)) {
          maxSendLength_ = lengthsTo_[i];
        }
      }
    }
    else {
      // not grouped by image, need send buffer and indicesTo_

      // starts[i] is the number of sends to node i
      // numActive equals number of sends total, \sum_i starts[i]

      // this loop starts at starts[1], so explicitly check starts[0]
      if (starts[0] == 0 ) {
        numSends_ = 0;
      }
      else {
        numSends_ = 1;
      }
      for (Teuchos::Array<size_t>::iterator i=starts.begin()+1,
                                            im1=starts.begin();
           i != starts.end(); ++i)
      {
        if (*i != 0) ++numSends_;
        *i += *im1;
        im1 = i;
      }
      // starts[i] now contains the number of exports to nodes 0 through i

      for (Teuchos::Array<size_t>::reverse_iterator ip1=starts.rbegin(),
                                                      i=starts.rbegin()+1;
           i != starts.rend(); ++i)
      {
        *ip1 = *i;
        ip1 = i;
      }
      starts[0] = 0;
      // starts[i] now contains the number of exports to nodes 0 through
      // i-1, i.e., all nodes before node i

      indicesTo_.resize(numActive);

      for (size_t i = 0; i < numExports_; ++i) {
        if (exportNodeIDs[i] >= 0) {
          // record the offset to the sendBuffer for this export
          indicesTo_[starts[exportNodeIDs[i]]] = i;
          // now increment the offset for this node
          ++starts[exportNodeIDs[i]];
        }
      }
      // our send buffer will contain the export data for each of the nodes
      // we communicate with, in order by node id
      // sendBuffer = {node_0_data, node_1_data, ..., node_np-1_data}
      // indicesTo now maps each export to the location in our send buffer
      // associated with the export
      // data for export i located at sendBuffer[indicesTo[i]]
      //
      // starts[i] once again contains the number of exports to
      // nodes 0 through i
      for (int node = numImages-1; node != 0; --node) {
        starts[node] = starts[node-1];
      }
      starts.front() = 0;
      starts[numImages] = numActive;
      //
      // starts[node] once again contains the number of exports to
      // nodes 0 through node-1
      // i.e., the start of my data in the sendBuffer

      // this contains invalid data at nodes we don't care about, that is okay
      imagesTo_.resize(numSends_);
      startsTo_.resize(numSends_);
      lengthsTo_.resize(numSends_);

      // for each group of sends/exports, record the destination node,
      // the length, and the offset for this send into the
      // send buffer (startsTo_)
      maxSendLength_ = 0;
      size_t snd = 0;
      for (int node = 0; node < numImages; ++node ) {
        if (starts[node+1] != starts[node]) {
          lengthsTo_[snd] = starts[node+1] - starts[node];
          startsTo_[snd] = starts[node];
          // record max length for all off-node sends
          if ((node != myImageID) && (lengthsTo_[snd] > maxSendLength_)) {
            maxSendLength_ = lengthsTo_[snd];
          }
          imagesTo_[snd] = node;
          ++snd;
        }
      }
#ifdef HAVE_TEUCHOS_DEBUG
      if (snd != numSends_) {
        send_neq_numSends = true;
      }
#endif
    }
#ifdef HAVE_TEUCHOS_DEBUG
        SHARED_TEST_FOR_EXCEPTION(index_neq_numActive, std::logic_error,
            "Tpetra::Distributor::createFromSends: logic error. Please notify the Tpetra team.",*comm_);
        SHARED_TEST_FOR_EXCEPTION(send_neq_numSends, std::logic_error,
            "Tpetra::Distributor::createFromSends: logic error. Please notify the Tpetra team.",*comm_);
#endif

    if (selfMessage_) --numSends_;

    // Invert map to see what msgs are received and what length
    computeReceives();

    return totalReceiveLength_;
  }

} // namespace Tpetra
