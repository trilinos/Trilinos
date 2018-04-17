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
#include "Tpetra_Details_gathervPrint.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
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
  } // namespace Details

  Teuchos::Array<std::string>
  distributorSendTypes ()
  {
    Teuchos::Array<std::string> sendTypes;
    sendTypes.push_back ("Isend");
    sendTypes.push_back ("Rsend");
    sendTypes.push_back ("Send");
    sendTypes.push_back ("Ssend");
    return sendTypes;
  }

  // We set default values of Distributor's Boolean parameters here,
  // in this one place.  That way, if we want to change the default
  // value of a parameter, we don't have to search the whole file to
  // ensure a consistent setting.
  namespace {
    // Default value of the "Debug" parameter.
    const bool tpetraDistributorDebugDefault = false;
    // Default value of the "Barrier between receives and sends" parameter.
    const bool barrierBetween_default = false;
    // Default value of the "Use distinct tags" parameter.
    const bool useDistinctTags_default = true;
  } // namespace (anonymous)

  int Distributor::getTag (const int pathTag) const {
    return useDistinctTags_ ? pathTag : comm_->getTag ();
  }


#ifdef TPETRA_DISTRIBUTOR_TIMERS
  void Distributor::makeTimers () {
    const std::string name_doPosts3 = "Tpetra::Distributor: doPosts(3)";
    const std::string name_doPosts4 = "Tpetra::Distributor: doPosts(4)";
    const std::string name_doWaits = "Tpetra::Distributor: doWaits";
    const std::string name_doPosts3_recvs = "Tpetra::Distributor: doPosts(3): recvs";
    const std::string name_doPosts4_recvs = "Tpetra::Distributor: doPosts(4): recvs";
    const std::string name_doPosts3_barrier = "Tpetra::Distributor: doPosts(3): barrier";
    const std::string name_doPosts4_barrier = "Tpetra::Distributor: doPosts(4): barrier";
    const std::string name_doPosts3_sends = "Tpetra::Distributor: doPosts(3): sends";
    const std::string name_doPosts4_sends = "Tpetra::Distributor: doPosts(4): sends";

    timer_doPosts3_ = Teuchos::TimeMonitor::getNewTimer (name_doPosts3);
    timer_doPosts4_ = Teuchos::TimeMonitor::getNewTimer (name_doPosts4);
    timer_doWaits_ = Teuchos::TimeMonitor::getNewTimer (name_doWaits);
    timer_doPosts3_recvs_ = Teuchos::TimeMonitor::getNewTimer (name_doPosts3_recvs);
    timer_doPosts4_recvs_ = Teuchos::TimeMonitor::getNewTimer (name_doPosts4_recvs);
    timer_doPosts3_barrier_ = Teuchos::TimeMonitor::getNewTimer (name_doPosts3_barrier);
    timer_doPosts4_barrier_ = Teuchos::TimeMonitor::getNewTimer (name_doPosts4_barrier);
    timer_doPosts3_sends_ = Teuchos::TimeMonitor::getNewTimer (name_doPosts3_sends);
    timer_doPosts4_sends_ = Teuchos::TimeMonitor::getNewTimer (name_doPosts4_sends);
  }
#endif // TPETRA_DISTRIBUTOR_TIMERS

  void
  Distributor::init (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                     const Teuchos::RCP<Teuchos::FancyOStream>& out,
                     const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    const bool verbose = Tpetra::Details::Behavior::verbose("Distributor");
    this->out_ = out.is_null () ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) : out;
    if (! plist.is_null ()) {
      this->setParameterList (plist);
    }

#ifdef TPETRA_DISTRIBUTOR_TIMERS
    makeTimers ();
#endif // TPETRA_DISTRIBUTOR_TIMERS

    if (verbose || debug_) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (out_.is_null (), std::logic_error, "Tpetra::Distributor::init: "
         "verbose and/or debug_ are true but out_ (pointer to the output "
         "stream) is NULL.  Please report this bug to the Tpetra developers.");
    }
    if (verbose) {
      Teuchos::OSTab tab (out_);
      std::ostringstream os;
      os << comm_->getRank ()
         << ": Distributor ctor done" << std::endl;
      *out_ << os.str ();
    }
  }

  Distributor::Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
    : comm_ (comm)
    , howInitialized_ (Details::DISTRIBUTOR_NOT_INITIALIZED)
    , sendType_ (Details::DISTRIBUTOR_SEND)
    , barrierBetween_ (barrierBetween_default)
    , debug_ (tpetraDistributorDebugDefault)
    , selfMessage_ (false)
    , numSends_ (0)
    , maxSendLength_ (0)
    , numReceives_ (0)
    , totalReceiveLength_ (0)
    , lastRoundBytesSend_ (0)
    , lastRoundBytesRecv_ (0)
    , useDistinctTags_ (useDistinctTags_default)
  {
    init (comm, Teuchos::null, Teuchos::null);
  }

  Distributor::Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                            const Teuchos::RCP<Teuchos::FancyOStream>& out)
    : comm_ (comm)
    , howInitialized_ (Details::DISTRIBUTOR_NOT_INITIALIZED)
    , sendType_ (Details::DISTRIBUTOR_SEND)
    , barrierBetween_ (barrierBetween_default)
    , debug_ (tpetraDistributorDebugDefault)
    , selfMessage_ (false)
    , numSends_ (0)
    , maxSendLength_ (0)
    , numReceives_ (0)
    , totalReceiveLength_ (0)
    , lastRoundBytesSend_ (0)
    , lastRoundBytesRecv_ (0)
    , useDistinctTags_ (useDistinctTags_default)
  {
    init (comm, out, Teuchos::null);
  }

  Distributor::Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                            const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : comm_ (comm)
    , howInitialized_ (Details::DISTRIBUTOR_NOT_INITIALIZED)
    , sendType_ (Details::DISTRIBUTOR_SEND)
    , barrierBetween_ (barrierBetween_default)
    , debug_ (tpetraDistributorDebugDefault)
    , selfMessage_ (false)
    , numSends_ (0)
    , maxSendLength_ (0)
    , numReceives_ (0)
    , totalReceiveLength_ (0)
    , lastRoundBytesSend_ (0)
    , lastRoundBytesRecv_ (0)
    , useDistinctTags_ (useDistinctTags_default)
  {
    init (comm, Teuchos::null, plist);
  }

  Distributor::Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                            const Teuchos::RCP<Teuchos::FancyOStream>& out,
                            const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : comm_ (comm)
    , howInitialized_ (Details::DISTRIBUTOR_NOT_INITIALIZED)
    , sendType_ (Details::DISTRIBUTOR_SEND)
    , barrierBetween_ (barrierBetween_default)
    , debug_ (tpetraDistributorDebugDefault)
    , selfMessage_ (false)
    , numSends_ (0)
    , maxSendLength_ (0)
    , numReceives_ (0)
    , totalReceiveLength_ (0)
    , lastRoundBytesSend_ (0)
    , lastRoundBytesRecv_ (0)
    , useDistinctTags_ (useDistinctTags_default)
  {
    init (comm, out, plist);
  }

  Distributor::Distributor (const Distributor & distributor)
    : comm_ (distributor.comm_)
    , out_ (distributor.out_)
    , howInitialized_ (Details::DISTRIBUTOR_INITIALIZED_BY_COPY)
    , sendType_ (distributor.sendType_)
    , barrierBetween_ (distributor.barrierBetween_)
    , debug_ (distributor.debug_)
    , selfMessage_ (distributor.selfMessage_)
    , numSends_ (distributor.numSends_)
    , procsTo_ (distributor.procsTo_)
    , startsTo_ (distributor.startsTo_)
    , lengthsTo_ (distributor.lengthsTo_)
    , maxSendLength_ (distributor.maxSendLength_)
    , indicesTo_ (distributor.indicesTo_)
    , numReceives_ (distributor.numReceives_)
    , totalReceiveLength_ (distributor.totalReceiveLength_)
    , lengthsFrom_ (distributor.lengthsFrom_)
    , procsFrom_ (distributor.procsFrom_)
    , startsFrom_ (distributor.startsFrom_)
    , indicesFrom_ (distributor.indicesFrom_)
    , reverseDistributor_ (distributor.reverseDistributor_)
    , lastRoundBytesSend_ (distributor.lastRoundBytesSend_)
    , lastRoundBytesRecv_ (distributor.lastRoundBytesRecv_)
    , useDistinctTags_ (distributor.useDistinctTags_)
  {
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

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

#ifdef TPETRA_DISTRIBUTOR_TIMERS
    makeTimers ();
#endif // TPETRA_DISTRIBUTOR_TIMERS

    const bool verbose = Tpetra::Details::Behavior::verbose("Distributor");
    if (verbose || debug_) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (out_.is_null (), std::logic_error, "Tpetra::Distributor::init: "
         "verbose and/or debug_ are true but out_ (pointer to the output "
         "stream) is NULL.  Please report this bug to the Tpetra developers.");
      Teuchos::OSTab tab (out_);
    }
    if (verbose) {
      std::ostringstream os;
      os << comm_->getRank ()
         << ": Distributor copy ctor done" << std::endl;
      *out_ << os.str ();
    }
  }

  void Distributor::swap (Distributor& rhs) {
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;

    std::swap (comm_, rhs.comm_);
    std::swap (out_, rhs.out_);
    std::swap (howInitialized_, rhs.howInitialized_);
    std::swap (sendType_, rhs.sendType_);
    std::swap (barrierBetween_, rhs.barrierBetween_);
    std::swap (debug_, rhs.debug_);
    std::swap (selfMessage_, rhs.selfMessage_);
    std::swap (numSends_, rhs.numSends_);
    std::swap (procsTo_, rhs.procsTo_);
    std::swap (startsTo_, rhs.startsTo_);
    std::swap (lengthsTo_, rhs.lengthsTo_);
    std::swap (maxSendLength_, rhs.maxSendLength_);
    std::swap (indicesTo_, rhs.indicesTo_);
    std::swap (numReceives_, rhs.numReceives_);
    std::swap (totalReceiveLength_, rhs.totalReceiveLength_);
    std::swap (lengthsFrom_, rhs.lengthsFrom_);
    std::swap (procsFrom_, rhs.procsFrom_);
    std::swap (startsFrom_, rhs.startsFrom_);
    std::swap (indicesFrom_, rhs.indicesFrom_);
    std::swap (reverseDistributor_, rhs.reverseDistributor_);
    std::swap (lastRoundBytesSend_, rhs.lastRoundBytesSend_);
    std::swap (lastRoundBytesRecv_, rhs.lastRoundBytesRecv_);
    std::swap (useDistinctTags_, rhs.useDistinctTags_);

    // Swap parameter lists.  If they are the same object, make a deep
    // copy first, so that modifying one won't modify the other one.
    RCP<ParameterList> lhsList = this->getNonconstParameterList ();
    RCP<ParameterList> rhsList = rhs.getNonconstParameterList ();
    if (lhsList.getRawPtr () == rhsList.getRawPtr () && ! rhsList.is_null ()) {
      rhsList = parameterList (*rhsList);
    }
    if (! rhsList.is_null ()) {
      this->setMyParamList (rhsList);
    }
    if (! lhsList.is_null ()) {
      rhs.setMyParamList (lhsList);
    }

    // We don't need to swap timers, because all instances of
    // Distributor use the same timers.
  }

  Distributor::~Distributor()
  {
    // mfh 10 May 2017: We shouldn't have any outstanding
    // communication requests at this point.  It would be legitimate
    // to check here, and report an error if requests_.size() != 0.
    // However, throwing in a destructor is bad form.  See #1303:
    //
    // https://github.com/trilinos/Trilinos/issues/1303
    //
    // If someone wants to restore the error check, please don't
    // throw; instead, use MPI_Abort (or exit() in a non-MPI build).

    // TEUCHOS_TEST_FOR_EXCEPTION(requests_.size() != 0, std::runtime_error,
    //   "Tpetra::Distributor: Destructor called with " << requests_.size()
    //   << " outstanding posts (unfulfilled communication requests).  There "
    //   "should be none at this point.  Please report this bug to the Tpetra "
    //   "developers.");
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
    const Details::EDistributorSendType sendType =
      getIntegralValue<Details::EDistributorSendType> (*plist, "Send type");
    const bool useDistinctTags = plist->get<bool> ("Use distinct tags");
    const bool debug = plist->get<bool> ("Debug");
    {
      // mfh 03 May 2016: We keep this option only for backwards
      // compatibility, but it must always be true.  See discussion of
      // Github Issue #227.
      const bool enable_cuda_rdma =
        plist->get<bool> ("Enable MPI CUDA RDMA support");
      TEUCHOS_TEST_FOR_EXCEPTION
        (! enable_cuda_rdma, std::invalid_argument, "Tpetra::Distributor::"
         "setParameterList: " << "You specified \"Enable MPI CUDA RDMA "
         "support\" = false.  This is no longer valid.  You don't need to "
         "specify this option any more; Tpetra assumes it is always true.  "
         "This is a very light assumption on the MPI implementation, and in "
         "fact does not actually involve hardware or system RDMA support.  "
         "Tpetra just assumes that the MPI implementation can tell whether a "
         "pointer points to host memory or CUDA device memory.");
    }

    // We check this property explicitly, since we haven't yet learned
    // how to make a validator that can cross-check properties.
    // Later, turn this into a validator so that it can be embedded in
    // the valid ParameterList and used in Optika.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! barrierBetween && sendType == Details::DISTRIBUTOR_RSEND,
      std::invalid_argument, "Tpetra::Distributor::setParameterList: " << endl
      << "You specified \"Send type\"=\"Rsend\", but turned off the barrier "
      "between receives and sends." << endl << "This is invalid; you must "
      "include the barrier if you use ready sends." << endl << "Ready sends "
      "require that their corresponding receives have already been posted, "
      "and the only way to guarantee that in general is with a barrier.");

    // Now that we've validated the input list, save the results.
    sendType_ = sendType;
    barrierBetween_ = barrierBetween;
    useDistinctTags_ = useDistinctTags;
    debug_ = debug;

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

    const bool barrierBetween = barrierBetween_default;
    const bool useDistinctTags = useDistinctTags_default;
    const bool debug = tpetraDistributorDebugDefault;

    Array<std::string> sendTypes = distributorSendTypes ();
    const std::string defaultSendType ("Send");
    Array<Details::EDistributorSendType> sendTypeEnums;
    sendTypeEnums.push_back (Details::DISTRIBUTOR_ISEND);
    sendTypeEnums.push_back (Details::DISTRIBUTOR_RSEND);
    sendTypeEnums.push_back (Details::DISTRIBUTOR_SEND);
    sendTypeEnums.push_back (Details::DISTRIBUTOR_SSEND);

    RCP<ParameterList> plist = parameterList ("Tpetra::Distributor");
    plist->set ("Barrier between receives and sends", barrierBetween,
                "Whether to execute a barrier between receives and sends in do"
                "[Reverse]Posts().  Required for correctness when \"Send type\""
                "=\"Rsend\", otherwise correct but not recommended.");
    setStringToIntegralParameter<Details::EDistributorSendType> ("Send type",
      defaultSendType, "When using MPI, the variant of send to use in "
      "do[Reverse]Posts()", sendTypes(), sendTypeEnums(), plist.getRawPtr());
    plist->set ("Use distinct tags", useDistinctTags, "Whether to use distinct "
                "MPI message tags for different code paths.  Highly recommended"
                " to avoid message collisions.");
    plist->set ("Debug", debug, "Whether to print copious debugging output on "
                "all processes.");
    plist->set ("Timer Label","","Label for Time Monitor output");
    plist->set ("Enable MPI CUDA RDMA support", true, "Assume that MPI can "
                "tell whether a pointer points to host memory or CUDA device "
                "memory.  You don't need to specify this option any more; "
                "Tpetra assumes it is always true.  This is a very light "
                "assumption on the MPI implementation, and in fact does not "
                "actually involve hardware or system RDMA support.");

    // mfh 24 Dec 2015: Tpetra no longer inherits from
    // Teuchos::VerboseObject, so it doesn't need the "VerboseObject"
    // sublist.  However, we retain the "VerboseObject" sublist
    // anyway, for backwards compatibility (otherwise the above
    // validation would fail with an invalid parameter name, should
    // the user still want to provide this list).
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

  Teuchos::ArrayView<const int> Distributor::getProcsFrom() const
  { return procsFrom_; }

  Teuchos::ArrayView<const size_t> Distributor::getLengthsFrom() const
  { return lengthsFrom_; }

  Teuchos::ArrayView<const int> Distributor::getProcsTo() const
  { return procsTo_; }

  Teuchos::ArrayView<const size_t> Distributor::getLengthsTo() const
  { return lengthsTo_; }

  Teuchos::RCP<Distributor>
  Distributor::getReverse() const {
    if (reverseDistributor_.is_null ()) {
      createReverseDistributor ();
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (reverseDistributor_.is_null (), std::logic_error, "The reverse "
       "Distributor is null after createReverseDistributor returned.  "
       "Please report this bug to the Tpetra developers.");
    return reverseDistributor_;
  }


  void
  Distributor::createReverseDistributor() const
  {
    reverseDistributor_ = Teuchos::rcp (new Distributor (comm_, out_));
    reverseDistributor_->howInitialized_ = Details::DISTRIBUTOR_INITIALIZED_BY_REVERSE;
    reverseDistributor_->sendType_ = sendType_;
    reverseDistributor_->barrierBetween_ = barrierBetween_;
    reverseDistributor_->debug_ = debug_;

    // The total length of all the sends of this Distributor.  We
    // calculate it because it's the total length of all the receives
    // of the reverse Distributor.
    size_t totalSendLength =
      std::accumulate (lengthsTo_.begin(), lengthsTo_.end(), 0);

    // The maximum length of any of the receives of this Distributor.
    // We calculate it because it's the maximum length of any of the
    // sends of the reverse Distributor.
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

    // Initialize all of reverseDistributor's data members.  This
    // mainly just involves flipping "send" and "receive," or the
    // equivalent "to" and "from."

    reverseDistributor_->selfMessage_ = selfMessage_;
    reverseDistributor_->numSends_ = numReceives_;
    reverseDistributor_->procsTo_ = procsFrom_;
    reverseDistributor_->startsTo_ = startsFrom_;
    reverseDistributor_->lengthsTo_ = lengthsFrom_;
    reverseDistributor_->maxSendLength_ = maxReceiveLength;
    reverseDistributor_->indicesTo_ = indicesFrom_;
    reverseDistributor_->numReceives_ = numSends_;
    reverseDistributor_->totalReceiveLength_ = totalSendLength;
    reverseDistributor_->lengthsFrom_ = lengthsTo_;
    reverseDistributor_->procsFrom_ = procsTo_;
    reverseDistributor_->startsFrom_ = startsTo_;
    reverseDistributor_->indicesFrom_ = indicesTo_;

    // requests_: Allocated on demand.
    // reverseDistributor_: See note below

    // mfh 31 Mar 2016: These are statistics, kept on calls to
    // doPostsAndWaits or doReversePostsAndWaits.  They weren't here
    // when I started, and I didn't add them, so I don't know if they
    // are accurate.
    reverseDistributor_->lastRoundBytesSend_ = 0;
    reverseDistributor_->lastRoundBytesRecv_ = 0;

    reverseDistributor_->useDistinctTags_ = useDistinctTags_;

    // I am my reverse Distributor's reverse Distributor.
    // Thus, it would be legit to do the following:
    //
    // reverseDistributor_->reverseDistributor_ = Teuchos::rcp (this, false);
    //
    // (Note use of a "weak reference" to avoid a circular RCP
    // dependency.)  The only issue is that if users hold on to the
    // reverse Distributor but let go of the forward one, this
    // reference won't be valid anymore.  However, the reverse
    // Distributor is really an implementation detail of Distributor
    // and not meant to be used directly, so we don't need to do this.
    reverseDistributor_->reverseDistributor_ = Teuchos::null;
  }


  void Distributor::doWaits() {
    using Teuchos::Array;
    using Teuchos::CommRequest;
    using Teuchos::FancyOStream;
    using Teuchos::includesVerbLevel;
    using Teuchos::is_null;
    using Teuchos::OSTab;
    using Teuchos::RCP;
    using Teuchos::waitAll;
    using std::endl;

    Teuchos::OSTab tab (out_);

#ifdef TPETRA_DISTRIBUTOR_TIMERS
    Teuchos::TimeMonitor timeMon (*timer_doWaits_);
#endif // TPETRA_DISTRIBUTOR_TIMERS

    const int myRank = comm_->getRank ();
    const bool verbose = Tpetra::Details::Behavior::verbose("Distributor");

    if (verbose) {
      std::ostringstream os;
      os << myRank << ": doWaits: # reqs = "
         << requests_.size () << endl;
      *out_ << os.str ();
    }

    if (requests_.size() > 0) {
      waitAll (*comm_, requests_());

#ifdef HAVE_TEUCHOS_DEBUG
      // Make sure that waitAll() nulled out all the requests.
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

#ifdef HAVE_TEUCHOS_DEBUG
    {
      const int localSizeNonzero = (requests_.size () != 0) ? 1 : 0;
      int globalSizeNonzero = 0;
      Teuchos::reduceAll<int, int> (*comm_, Teuchos::REDUCE_MAX,
                                    localSizeNonzero,
                                    Teuchos::outArg (globalSizeNonzero));
      TEUCHOS_TEST_FOR_EXCEPTION(
        globalSizeNonzero != 0, std::runtime_error,
        "Tpetra::Distributor::doWaits: After waitAll, at least one process has "
        "a nonzero number of outstanding posts.  There should be none at this "
        "point.  Please report this bug to the Tpetra developers.");
    }
#endif // HAVE_TEUCHOS_DEBUG

    if (verbose) {
      std::ostringstream os;
      os << myRank << ": doWaits done" << endl;
      *out_ << os.str ();
    }
  }

  void Distributor::doReverseWaits() {
    // call doWaits() on the reverse Distributor, if it exists
    if (! reverseDistributor_.is_null()) {
      reverseDistributor_->doWaits();
    }
  }

  std::string Distributor::description () const {
    std::ostringstream out;

    out << "\"Tpetra::Distributor\": {";
    const std::string label = this->getObjectLabel ();
    if (label != "") {
      out << "Label: " << label << ", ";
    }
    out << "How initialized: "
        << Details::DistributorHowInitializedEnumToString (howInitialized_)
        << ", Parameters: {"
        << "Send type: "
        << DistributorSendTypeEnumToString (sendType_)
        << ", Barrier between receives and sends: "
        << (barrierBetween_ ? "true" : "false")
        << ", Use distinct tags: "
        << (useDistinctTags_ ? "true" : "false")
        << ", Debug: " << (debug_ ? "true" : "false")
        << "}}";
    return out.str ();
  }

  std::string
  Distributor::
  localDescribeToString (const Teuchos::EVerbosityLevel vl) const
  {
    using Teuchos::toString;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    using std::endl;

    // This preserves current behavior of Distributor.
    if (vl <= Teuchos::VERB_LOW || comm_.is_null ()) {
      return std::string ();
    }

    auto outStringP = Teuchos::rcp (new std::ostringstream ());
    auto outp = Teuchos::getFancyOStream (outStringP); // returns RCP
    Teuchos::FancyOStream& out = *outp;

    const int myRank = comm_->getRank ();
    const int numProcs = comm_->getSize ();
    out << "Process " << myRank << " of " << numProcs << ":" << endl;
    Teuchos::OSTab tab1 (out);

    out << "selfMessage: " << hasSelfMessage () << endl;
    out << "numSends: " << getNumSends () << endl;
    if (vl == VERB_HIGH || vl == VERB_EXTREME) {
      out << "procsTo: " << toString (procsTo_) << endl;
      out << "lengthsTo: " << toString (lengthsTo_) << endl;
      out << "maxSendLength: " << getMaxSendLength () << endl;
    }
    if (vl == VERB_EXTREME) {
      out << "startsTo: " << toString (startsTo_) << endl;
      out << "indicesTo: " << toString (indicesTo_) << endl;
    }
    if (vl == VERB_HIGH || vl == VERB_EXTREME) {
      out << "numReceives: " << getNumReceives () << endl;
      out << "totalReceiveLength: " << getTotalReceiveLength () << endl;
      out << "lengthsFrom: " << toString (lengthsFrom_) << endl;
      out << "startsFrom: " << toString (startsFrom_) << endl;
      out << "procsFrom: " << toString (procsFrom_) << endl;
    }

    out.flush (); // make sure the ostringstream got everything
    return outStringP->str ();
  }

  void
  Distributor::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using std::endl;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    const Teuchos::EVerbosityLevel vl =
      (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

    if (vl == VERB_NONE) {
      return; // don't print anything
    }
    // If this Distributor's Comm is null, then the the calling
    // process does not participate in Distributor-related collective
    // operations with the other processes.  In that case, it is not
    // even legal to call this method.  The reasonable thing to do in
    // that case is nothing.
    if (comm_.is_null ()) {
      return;
    }
    const int myRank = comm_->getRank ();
    const int numProcs = comm_->getSize ();

    // Only Process 0 should touch the output stream, but this method
    // in general may need to do communication.  Thus, we may need to
    // preserve the current tab level across multiple "if (myRank ==
    // 0) { ... }" inner scopes.  This is why we sometimes create
    // OSTab instances by pointer, instead of by value.  We only need
    // to create them by pointer if the tab level must persist through
    // multiple inner scopes.
    Teuchos::RCP<Teuchos::OSTab> tab0, tab1;

    if (myRank == 0) {
      // At every verbosity level but VERB_NONE, Process 0 prints.
      // By convention, describe() always begins with a tab before
      // printing.
      tab0 = Teuchos::rcp (new Teuchos::OSTab (out));
      // We quote the class name because it contains colons.
      // This makes the output valid YAML.
      out << "\"Tpetra::Distributor\":" << endl;
      tab1 = Teuchos::rcp (new Teuchos::OSTab (out));

      const std::string label = this->getObjectLabel ();
      if (label != "") {
        out << "Label: " << label << endl;
      }
      out << "Number of processes: " << numProcs << endl
          << "How initialized: "
          << Details::DistributorHowInitializedEnumToString (howInitialized_)
          << endl;
      {
        out << "Parameters: " << endl;
        Teuchos::OSTab tab2 (out);
        out << "\"Send type\": "
            << DistributorSendTypeEnumToString (sendType_) << endl
            << "\"Barrier between receives and sends\": "
            << (barrierBetween_ ? "true" : "false") << endl
            << "\"Use distinct tags\": "
            << (useDistinctTags_ ? "true" : "false") << endl
            << "\"Debug\": " << (debug_ ? "true" : "false") << endl;
      }
    } // if myRank == 0

    // This is collective over the Map's communicator.
    if (vl > VERB_LOW) {
      const std::string lclStr = this->localDescribeToString (vl);
      Tpetra::Details::gathervPrint (out, lclStr, *comm_);
    }

    out << "Reverse Distributor:";
    if (reverseDistributor_.is_null ()) {
      out << " null" << endl;
    }
    else {
      out << endl;
      reverseDistributor_->describe (out, vl);
    }
  }

  void
  Distributor::computeReceives ()
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
    using std::endl;

    Teuchos::OSTab tab (out_);
    const int myRank = comm_->getRank();
    const int numProcs = comm_->getSize();

    // MPI tag for nonblocking receives and blocking sends in this method.
    const int pathTag = 2;
    const int tag = this->getTag (pathTag);
    const bool verbose = Tpetra::Details::Behavior::verbose("Distributor");

    if (verbose) {
      std::ostringstream os;
      os << myRank << ": computeReceives: "
        "{selfMessage_: " << (selfMessage_ ? "true" : "false")
         << ", tag: " << tag << "}" << endl;
      *out_ << os.str ();
    }

    // toProcsFromMe[i] == the number of messages sent by this process
    // to process i.  The data in numSends_, procsTo_, and lengthsTo_
    // concern the contiguous sends.  Therefore, each process will be
    // listed in procsTo_ at most once, and so toProcsFromMe[i] will
    // either be 0 or 1.
    {
      Array<int> toProcsFromMe (numProcs, 0);
#ifdef HAVE_TEUCHOS_DEBUG
      bool counting_error = false;
#endif // HAVE_TEUCHOS_DEBUG
      for (size_t i = 0; i < (numSends_ + (selfMessage_ ? 1 : 0)); ++i) {
#ifdef HAVE_TEUCHOS_DEBUG
        if (toProcsFromMe[procsTo_[i]] != 0) {
          counting_error = true;
        }
#endif // HAVE_TEUCHOS_DEBUG
        toProcsFromMe[procsTo_[i]] = 1;
      }
#ifdef HAVE_TEUCHOS_DEBUG
      SHARED_TEST_FOR_EXCEPTION(counting_error, std::logic_error,
        "Tpetra::Distributor::computeReceives: There was an error on at least "
        "one process in counting the number of messages send by that process to "
        "the other processs.  Please report this bug to the Tpetra developers.",
        *comm_);
#endif // HAVE_TEUCHOS_DEBUG

      if (verbose) {
        std::ostringstream os;
        os << myRank << ": computeReceives: Calling reduce and scatter" << endl;
        *out_ << os.str ();
      }

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
    const size_t actualNumReceives = numReceives_ - (selfMessage_ ? 1 : 0);

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

    if (verbose) {
      std::ostringstream os;
      os << myRank << ": computeReceives: Posting "
         << actualNumReceives << " irecvs" << endl;
      *out_ << os.str ();
    }

    // Post the (nonblocking) receives.
    for (size_t i = 0; i < actualNumReceives; ++i) {
      // Once the receive completes, we can ask the corresponding
      // CommStatus object (output by wait()) for the sending process'
      // ID (which we'll assign to procsFrom_[i] -- don't forget to
      // do that!).
      lengthsFromBuffers[i].resize (1);
      lengthsFromBuffers[i][0] = as<size_t> (0);
      requests[i] = ireceive<int, size_t> (lengthsFromBuffers[i], anySourceProc, tag, *comm_);
      if (verbose) {
        std::ostringstream os;
        os << myRank << ": computeReceives: "
          "Posted any-proc irecv w/ specified tag " << tag << endl;
        *out_ << os.str ();
      }
    }

    if (verbose) {
      std::ostringstream os;
      os << myRank << ": computeReceives: "
        "posting " << numSends_ << " sends" << endl;
      *out_ << os.str ();
    }
    // Post the sends: Tell each process to which we are sending how
    // many packets it should expect from us in the communication
    // pattern.  We could use nonblocking sends here, as long as we do
    // a waitAll() on all the sends and receives at once.
    //
    // We assume that numSends_ and selfMessage_ have already been
    // set.  The value of numSends_ (my process' number of sends) does
    // not include any message that it might send to itself.
    for (size_t i = 0; i < numSends_ + (selfMessage_ ? 1 : 0); ++i) {
      if (procsTo_[i] != myRank) {
        // Send a message to procsTo_[i], telling that process that
        // this communication pattern will send that process
        // lengthsTo_[i] blocks of packets.
        const size_t* const lengthsTo_i = &lengthsTo_[i];
        send<int, size_t> (lengthsTo_i, 1, as<int> (procsTo_[i]), tag, *comm_);
        if (verbose) {
          std::ostringstream os;
          os << myRank << ": computeReceives: "
            "Posted send to Proc " << procsTo_[i] << " w/ specified tag "
             << tag << endl;
          *out_ << os.str ();
        }
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

    if (verbose) {
      std::ostringstream os;
      os << myRank << ": computeReceives: waitAll on "
         << requests.size () << " requests" << endl;
      *out_ << os.str ();
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

    if (verbose) {
      std::ostringstream os;
      os << myRank << ": computeReceives: done" << endl;
      *out_ << os.str ();
    }
  }

  size_t
  Distributor::createFromSends (const Teuchos::ArrayView<const int> &exportProcIDs)
  {
    using Teuchos::outArg;
    using Teuchos::REDUCE_MAX;
    using Teuchos::reduceAll;
    using std::endl;

    Teuchos::OSTab tab (out_);
    const size_t numExports = exportProcIDs.size();
    const int myProcID = comm_->getRank();
    const int numProcs = comm_->getSize();
    const bool verbose = Tpetra::Details::Behavior::verbose("Distributor");
    if (verbose) {
      std::ostringstream os;
      os << myProcID << ": createFromSends" << endl;
      *out_ << os.str ();
    }

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
    //   * HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS
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

#ifdef HAVE_TPETRA_DEBUG
    int badID = -1; // only used in a debug build
#endif // HAVE_TPETRA_DEBUG
    for (size_t i = 0; i < numExports; ++i) {
      const int exportID = exportProcIDs[i];
      if (exportID >= numProcs) {
#ifdef HAVE_TPETRA_DEBUG
        badID = myProcID;
#endif // HAVE_TPETRA_DEBUG
        break;
      }
      else if (exportID >= 0) {
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
        if (needSendBuff==0 && starts[exportID] > 1 && exportID != exportProcIDs[i-1]) {
          needSendBuff = 1;
        }
        ++numActive;
      }
    }

#ifdef HAVE_TPETRA_DEBUG
    // Test whether any process in the communicator got an invalid
    // process ID.  If badID != -1 on this process, then it equals
    // this process' rank.  The max of all badID over all processes is
    // the max rank which has an invalid process ID.
    {
      int gbl_badID;
      reduceAll<int, int> (*comm_, REDUCE_MAX, badID, outArg (gbl_badID));
      TEUCHOS_TEST_FOR_EXCEPTION(gbl_badID >= 0, std::runtime_error,
        Teuchos::typeName(*this) << "::createFromSends(): Process  " << gbl_badID
        << ", perhaps among other processes, got a bad send process ID.");
    }
#else
    // FIXME (mfh 12 Apr 2013, 15 Jul 2015) Rather than simply
    // ignoring this information, we should think about how to pass it
    // along so that all the processes find out about it.  In a
    // release build with efficiency warnings turned off, the next
    // collective communication happens in computeReceives().  We
    // could figure out how to encode the error flag in that
    // operation, for example by adding an extra entry to the
    // collective's output array that encodes the error condition (0
    // on all processes if no error, else 1 on any process with the
    // error, so that the sum will produce a nonzero value if any
    // process had an error).  I'll defer this change for now and
    // recommend instead that people with troubles try a debug build.
#endif // HAVE_TPETRA_DEBUG

#if defined(HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS) || defined(HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS)
    {
      int global_needSendBuff;
      reduceAll<int, int> (*comm_, REDUCE_MAX, needSendBuff,
                            outArg (global_needSendBuff));
      TPETRA_EFFICIENCY_WARNING(
        global_needSendBuff != 0, std::runtime_error,
        "::createFromSends: Grouping export IDs together by process rank often "
        "improves performance.");
    }
#endif

    // Determine from the caller's data whether or not the current
    // process should send (a) message(s) to itself.
    if (starts[myProcID] != 0) {
      selfMessage_ = true;
    }
    else {
      selfMessage_ = false;
    }

#ifdef HAVE_TEUCHOS_DEBUG
    bool index_neq_numActive = false;
    bool send_neq_numSends = false;
#endif
    if (! needSendBuff) {
      // grouped by proc, no send buffer or indicesTo_ needed
      numSends_ = 0;
      // Count total number of sends, i.e., total number of procs to
      // which we are sending.  This includes myself, if applicable.
      for (int i = 0; i < numProcs; ++i) {
        if (starts[i]) {
          ++numSends_;
        }
      }

      // Not only do we not need these, but we must clear them, as
      // empty status of indicesTo is a flag used later.
      indicesTo_.resize(0);
      // Size these to numSends_; note, at the moment, numSends_
      // includes self sends.  Set their values to zeros.
      procsTo_.assign(numSends_,0);
      startsTo_.assign(numSends_,0);
      lengthsTo_.assign(numSends_,0);

      // set startsTo to the offset for each send (i.e., each proc ID)
      // set procsTo to the proc ID for each send
      // in interpreting this code, remember that we are assuming contiguity
      // that is why index skips through the ranks
      {
        size_t index = 0, procIndex = 0;
        for (size_t i = 0; i < numSends_; ++i) {
          while (exportProcIDs[procIndex] < 0) {
            ++procIndex; // skip all negative proc IDs
          }
          startsTo_[i] = procIndex;
          int procID = exportProcIDs[procIndex];
          procsTo_[i] = procID;
          index     += starts[procID];
          procIndex += starts[procID];
        }
#ifdef HAVE_TEUCHOS_DEBUG
        if (index != numActive) {
          index_neq_numActive = true;
        }
#endif
      }
      // sort the startsTo and proc IDs together, in ascending order, according
      // to proc IDs
      if (numSends_ > 0) {
        sort2(procsTo_.begin(), procsTo_.end(), startsTo_.begin());
      }
      // compute the maximum send length
      maxSendLength_ = 0;
      for (size_t i = 0; i < numSends_; ++i) {
        int procID = procsTo_[i];
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
      procsTo_.resize(numSends_);
      startsTo_.resize(numSends_);
      lengthsTo_.resize(numSends_);

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
          procsTo_[snd] = proc;
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

    if (verbose) {
      std::ostringstream os;
      os << myProcID << ": createFromSends: done" << endl;
      *out_ << os.str ();
    }

    // createFromRecvs() calls createFromSends(), but will set
    // howInitialized_ again after calling createFromSends().
    howInitialized_ = Details::DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_SENDS;

    return totalReceiveLength_;
  }

  void
  Distributor::
  createFromSendsAndRecvs (const Teuchos::ArrayView<const int>& exportProcIDs,
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

    selfMessage_ = ( starts[myProcID] != 0 ) ? 1 : 0;

    numSends_ = 0;

    if( needSendBuff ) //grouped by processor, no send buffer or indicesTo_ needed
      {
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
        procsTo_.resize(numSends_);
        startsTo_.resize(numSends_);
        lengthsTo_.resize(numSends_);
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
            procsTo_[snd] = proc;
            ++snd;
          }
        }
      }
    else {
      // grouped by proc, no send buffer or indicesTo_ needed
      numSends_ = 0;
      // Count total number of sends, i.e., total number of procs to
      // which we are sending.  This includes myself, if applicable.
      for (int i = 0; i < numProcs; ++i) {
        if (starts[i]) {
          ++numSends_;
        }
      }

      // Not only do we not need these, but we must clear them, as
      // empty status of indicesTo is a flag used later.
      indicesTo_.resize(0);
      // Size these to numSends_; note, at the moment, numSends_
      // includes self sends.  Set their values to zeros.
      procsTo_.assign(numSends_,0);
      startsTo_.assign(numSends_,0);
      lengthsTo_.assign(numSends_,0);

      // set startsTo to the offset for each send (i.e., each proc ID)
      // set procsTo to the proc ID for each send
      // in interpreting this code, remember that we are assuming contiguity
      // that is why index skips through the ranks
      {
        size_t index = 0, procIndex = 0;
        for (size_t i = 0; i < numSends_; ++i) {
          while (exportProcIDs[procIndex] < 0) {
            ++procIndex; // skip all negative proc IDs
          }
          startsTo_[i] = procIndex;
          int procID = exportProcIDs[procIndex];
          procsTo_[i] = procID;
          index     += starts[procID];
          procIndex += starts[procID];
        }
      }
      // sort the startsTo and proc IDs together, in ascending order, according
      // to proc IDs
      if (numSends_ > 0) {
        sort2(procsTo_.begin(), procsTo_.end(), startsTo_.begin());
      }
      // compute the maximum send length
      maxSendLength_ = 0;
      for (size_t i = 0; i < numSends_; ++i) {
        int procID = procsTo_[i];
        lengthsTo_[i] = starts[procID];
        if ((procID != myProcID) && (lengthsTo_[i] > maxSendLength_)) {
          maxSendLength_ = lengthsTo_[i];
        }
      }
    }


    numSends_ -= selfMessage_;
    std::vector<int> recv_list;
    recv_list.reserve(numSends_); //reserve an initial guess for size needed

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
    indicesFrom_.reserve (totalReceiveLength_);
    for (size_t i = 0; i < totalReceiveLength_; ++i) {
      indicesFrom_.push_back(i);
    }

    numReceives_-=selfMessage_;
  }

} // namespace Tpetra
