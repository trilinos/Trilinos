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

#include "Tpetra_Distributor.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Details_makeValidVerboseStream.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include <numeric>

namespace Tpetra {
  // We set default values of Distributor's Boolean parameters here,
  // in this one place.  That way, if we want to change the default
  // value of a parameter, we don't have to search the whole file to
  // ensure a consistent setting.
  namespace {
    // Default value of the "Debug" parameter.
    const bool tpetraDistributorDebugDefault = false;
  } // namespace (anonymous)

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

  Distributor::
  Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
               const Teuchos::RCP<Teuchos::FancyOStream>& /* out */,
               const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : plan_(comm)
    , lastRoundBytesSend_ (0)
    , lastRoundBytesRecv_ (0)
  {
    this->setParameterList(plist);
  }

  Distributor::
  Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
    : Distributor (comm, Teuchos::null, Teuchos::null)
  {}

  Distributor::
  Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
               const Teuchos::RCP<Teuchos::FancyOStream>& out)
    : Distributor (comm, out, Teuchos::null)
  {}

  Distributor::
  Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
               const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : Distributor (comm, Teuchos::null, plist)
  {}

  Distributor::
  Distributor (const Distributor& distributor)
    : plan_(distributor.plan_)
    , actor_(distributor.actor_)
    , verbose_ (distributor.verbose_)
    , reverseDistributor_ (distributor.reverseDistributor_)
    , lastRoundBytesSend_ (distributor.lastRoundBytesSend_)
    , lastRoundBytesRecv_ (distributor.lastRoundBytesRecv_)
  {
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<const ParameterList> rhsList = distributor.getParameterList ();
    RCP<ParameterList> newList = rhsList.is_null () ? Teuchos::null :
      Teuchos::parameterList (*rhsList);
    this->setParameterList (newList);
  }

  void Distributor::swap (Distributor& rhs) {
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;

    std::swap (plan_, rhs.plan_);
    std::swap (actor_, rhs.actor_);
    std::swap (verbose_, rhs.verbose_);
    std::swap (reverseDistributor_, rhs.reverseDistributor_);
    std::swap (lastRoundBytesSend_, rhs.lastRoundBytesSend_);
    std::swap (lastRoundBytesRecv_, rhs.lastRoundBytesRecv_);

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

  bool
  Distributor::getVerbose()
  {
    return Details::Behavior::verbose("Distributor") ||
      Details::Behavior::verbose("Tpetra::Distributor");
  }

  std::unique_ptr<std::string>
  Distributor::
  createPrefix(const char methodName[]) const
  {
    return Details::createPrefix(
      plan_.comm_.getRawPtr(), "Distributor", methodName);
  }

  void
  Distributor::
  setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    using ::Tpetra::Details::Behavior;
    using Teuchos::FancyOStream;
    using Teuchos::getIntegralValue;
    using Teuchos::includesVerbLevel;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using std::endl;

    if (! plist.is_null()) {
      RCP<const ParameterList> validParams = getValidParameters ();
      plist->validateParametersAndSetDefaults (*validParams);

      const bool barrierBetween =
        plist->get<bool> ("Barrier between receives and sends");
      const Details::EDistributorSendType sendType =
        getIntegralValue<Details::EDistributorSendType> (*plist, "Send type");
      const bool useDistinctTags = plist->get<bool> ("Use distinct tags");
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
      TEUCHOS_TEST_FOR_EXCEPTION
        (! barrierBetween && sendType == Details::DISTRIBUTOR_RSEND,
         std::invalid_argument, "Tpetra::Distributor::setParameterList: " << endl
         << "You specified \"Send type\"=\"Rsend\", but turned off the barrier "
         "between receives and sends." << endl << "This is invalid; you must "
         "include the barrier if you use ready sends." << endl << "Ready sends "
         "require that their corresponding receives have already been posted, "
         "and the only way to guarantee that in general is with a barrier.");

      // Now that we've validated the input list, save the results.
      plan_.sendType_ = sendType;
      plan_.barrierBetweenRecvSend_ = barrierBetween;
      plan_.useDistinctTags_ = useDistinctTags;

      // ParameterListAcceptor semantics require pointer identity of the
      // sublist passed to setParameterList(), so we save the pointer.
      this->setMyParamList (plist);
    }
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  Distributor::getValidParameters () const
  {
    using Teuchos::Array;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::setStringToIntegralParameter;

    const bool barrierBetween = Details::barrierBetween_default;
    const bool useDistinctTags = Details::useDistinctTags_default;
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
  { return plan_.totalReceiveLength_; }

  size_t Distributor::getNumReceives() const
  { return plan_.numReceives_; }

  bool Distributor::hasSelfMessage() const
  { return plan_.sendMessageToSelf_; }

  size_t Distributor::getNumSends() const
  { return plan_.numSendsToOtherProcs_; }

  size_t Distributor::getMaxSendLength() const
  { return plan_.maxSendLength_; }

  Teuchos::ArrayView<const int> Distributor::getProcsFrom() const
  { return plan_.procsFrom_; }

  Teuchos::ArrayView<const size_t> Distributor::getLengthsFrom() const
  { return plan_.lengthsFrom_; }

  Teuchos::ArrayView<const int> Distributor::getProcsTo() const
  { return plan_.procIdsToSendTo_; }

  Teuchos::ArrayView<const size_t> Distributor::getLengthsTo() const
  { return plan_.lengthsTo_; }

  Teuchos::RCP<Distributor>
  Distributor::getReverse(bool create) const {
    if (reverseDistributor_.is_null () && create) {
      createReverseDistributor ();
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (reverseDistributor_.is_null () && create, std::logic_error, "The reverse "
       "Distributor is null after createReverseDistributor returned.  "
       "Please report this bug to the Tpetra developers.");
    return reverseDistributor_;
  }


  void
  Distributor::createReverseDistributor() const
  {
    reverseDistributor_ = Teuchos::rcp(new Distributor(plan_.comm_));
    reverseDistributor_->plan_.howInitialized_ = Details::DISTRIBUTOR_INITIALIZED_BY_REVERSE;
    reverseDistributor_->plan_.sendType_ = plan_.sendType_;
    reverseDistributor_->plan_.barrierBetweenRecvSend_ = plan_.barrierBetweenRecvSend_;
    reverseDistributor_->verbose_ = verbose_;

    // The total length of all the sends of this Distributor.  We
    // calculate it because it's the total length of all the receives
    // of the reverse Distributor.
    size_t totalSendLength =
      std::accumulate (plan_.lengthsTo_.begin(), plan_.lengthsTo_.end(), 0);

    // The maximum length of any of the receives of this Distributor.
    // We calculate it because it's the maximum length of any of the
    // sends of the reverse Distributor.
    size_t maxReceiveLength = 0;
    const int myProcID = plan_.comm_->getRank();
    for (size_t i=0; i < plan_.numReceives_; ++i) {
      if (plan_.procsFrom_[i] != myProcID) {
        // Don't count receives for messages sent by myself to myself.
        if (plan_.lengthsFrom_[i] > maxReceiveLength) {
          maxReceiveLength = plan_.lengthsFrom_[i];
        }
      }
    }

    // Initialize all of reverseDistributor's data members.  This
    // mainly just involves flipping "send" and "receive," or the
    // equivalent "to" and "from."

    reverseDistributor_->plan_.sendMessageToSelf_ = plan_.sendMessageToSelf_;
    reverseDistributor_->plan_.numSendsToOtherProcs_ = plan_.numReceives_;
    reverseDistributor_->plan_.procIdsToSendTo_ = plan_.procsFrom_;
    reverseDistributor_->plan_.startsTo_ = plan_.startsFrom_;
    reverseDistributor_->plan_.lengthsTo_ = plan_.lengthsFrom_;
    reverseDistributor_->plan_.maxSendLength_ = maxReceiveLength;
    reverseDistributor_->plan_.indicesTo_ = plan_.indicesFrom_;
    reverseDistributor_->plan_.numReceives_ = plan_.numSendsToOtherProcs_;
    reverseDistributor_->plan_.totalReceiveLength_ = totalSendLength;
    reverseDistributor_->plan_.lengthsFrom_ = plan_.lengthsTo_;
    reverseDistributor_->plan_.procsFrom_ = plan_.procIdsToSendTo_;
    reverseDistributor_->plan_.startsFrom_ = plan_.startsTo_;
    reverseDistributor_->plan_.indicesFrom_ = plan_.indicesTo_;

    // requests_: Allocated on demand.
    // reverseDistributor_: See note below

    // mfh 31 Mar 2016: These are statistics, kept on calls to
    // doPostsAndWaits or doReversePostsAndWaits.  They weren't here
    // when I started, and I didn't add them, so I don't know if they
    // are accurate.
    reverseDistributor_->lastRoundBytesSend_ = 0;
    reverseDistributor_->lastRoundBytesRecv_ = 0;

    reverseDistributor_->plan_.useDistinctTags_ = plan_.useDistinctTags_;

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

  void
  Distributor::doWaits()
  {
    using Teuchos::Array;
    using Teuchos::CommRequest;
    using Teuchos::FancyOStream;
    using Teuchos::includesVerbLevel;
    using Teuchos::is_null;
    using Teuchos::RCP;
    using Teuchos::waitAll;
    using std::endl;

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor timeMon (*actor_.timer_doWaits_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

    const bool debug = Details::Behavior::debug("Distributor");

    std::unique_ptr<std::string> prefix;
    if (verbose_) {
      prefix = createPrefix("doWaits");
      std::ostringstream os;
      os << *prefix << "Start: requests_.size(): "
         << actor_.requests_.size() << endl;
      std::cerr << os.str();
    }

    if (actor_.requests_.size() > 0) {
      waitAll(*plan_.comm_, actor_.requests_());

      if (debug) {
        // Make sure that waitAll() nulled out all the requests.
        for (auto it = actor_.requests_.begin(); it != actor_.requests_.end(); ++it) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (! is_null(*it), std::runtime_error,
             "Tpetra::Distributor::doWaits: Communication requests "
             "should all be null aftr calling Teuchos::waitAll on "
             "them, but at least one request is not null.");
        }
      }
      // Restore the invariant that requests_.size() is the number of
      // outstanding nonblocking communication requests.
      actor_.requests_.resize (0);
    }

    if (debug) {
      const int localSizeNonzero = (actor_.requests_.size () != 0) ? 1 : 0;
      int globalSizeNonzero = 0;
      Teuchos::reduceAll<int, int> (*plan_.comm_, Teuchos::REDUCE_MAX,
                                    localSizeNonzero,
                                    Teuchos::outArg (globalSizeNonzero));
      TEUCHOS_TEST_FOR_EXCEPTION(
        globalSizeNonzero != 0, std::runtime_error,
        "Tpetra::Distributor::doWaits: After waitAll, at least one process has "
        "a nonzero number of outstanding posts.  There should be none at this "
        "point.  Please report this bug to the Tpetra developers.");
    }

    if (verbose_) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
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
        << Details::DistributorHowInitializedEnumToString (plan_.howInitialized_)
        << ", Parameters: {"
        << "Send type: "
        << DistributorSendTypeEnumToString (plan_.sendType_)
        << ", Barrier between receives and sends: "
        << (plan_.barrierBetweenRecvSend_ ? "true" : "false")
        << ", Use distinct tags: "
        << (plan_.useDistinctTags_ ? "true" : "false")
        << ", Debug: " << (verbose_ ? "true" : "false")
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
    if (vl <= Teuchos::VERB_LOW || plan_.comm_.is_null ()) {
      return std::string ();
    }

    auto outStringP = Teuchos::rcp (new std::ostringstream ());
    auto outp = Teuchos::getFancyOStream (outStringP); // returns RCP
    Teuchos::FancyOStream& out = *outp;

    const int myRank = plan_.comm_->getRank ();
    const int numProcs = plan_.comm_->getSize ();
    out << "Process " << myRank << " of " << numProcs << ":" << endl;
    Teuchos::OSTab tab1 (out);

    out << "selfMessage: " << hasSelfMessage () << endl;
    out << "numSends: " << getNumSends () << endl;
    if (vl == VERB_HIGH || vl == VERB_EXTREME) {
      out << "procsTo: " << toString (plan_.procIdsToSendTo_) << endl;
      out << "lengthsTo: " << toString (plan_.lengthsTo_) << endl;
      out << "maxSendLength: " << getMaxSendLength () << endl;
    }
    if (vl == VERB_EXTREME) {
      out << "startsTo: " << toString (plan_.startsTo_) << endl;
      out << "indicesTo: " << toString (plan_.indicesTo_) << endl;
    }
    if (vl == VERB_HIGH || vl == VERB_EXTREME) {
      out << "numReceives: " << getNumReceives () << endl;
      out << "totalReceiveLength: " << getTotalReceiveLength () << endl;
      out << "lengthsFrom: " << toString (plan_.lengthsFrom_) << endl;
      out << "startsFrom: " << toString (plan_.startsFrom_) << endl;
      out << "procsFrom: " << toString (plan_.procsFrom_) << endl;
    }

    out.flush (); // make sure the ostringstream got everything
    return outStringP->str ();
  }

  void
  Distributor::
  describe (Teuchos::FancyOStream& out,
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
    if (plan_.comm_.is_null ()) {
      return;
    }
    const int myRank = plan_.comm_->getRank ();
    const int numProcs = plan_.comm_->getSize ();

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
          << Details::DistributorHowInitializedEnumToString (plan_.howInitialized_)
          << endl;
      {
        out << "Parameters: " << endl;
        Teuchos::OSTab tab2 (out);
        out << "\"Send type\": "
            << DistributorSendTypeEnumToString (plan_.sendType_) << endl
            << "\"Barrier between receives and sends\": "
            << (plan_.barrierBetweenRecvSend_ ? "true" : "false") << endl
            << "\"Use distinct tags\": "
            << (plan_.useDistinctTags_ ? "true" : "false") << endl
            << "\"Debug\": " << (verbose_ ? "true" : "false") << endl;
      }
    } // if myRank == 0

    // This is collective over the Map's communicator.
    if (vl > VERB_LOW) {
      const std::string lclStr = this->localDescribeToString (vl);
      Tpetra::Details::gathervPrint (out, lclStr, *plan_.comm_);
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

  size_t
  Distributor::
  createFromSends(const Teuchos::ArrayView<const int>& exportProcIDs)
  {
    plan_.createFromSends(exportProcIDs);
  }

  void
  Distributor::
  createFromSendsAndRecvs (const Teuchos::ArrayView<const int>& exportProcIDs,
                           const Teuchos::ArrayView<const int>& remoteProcIDs)
  {
    std::unique_ptr<std::string> prefix;
    if (verbose_) {
      prefix = createPrefix("createFromSendsAndRecvs");
      std::ostringstream os;
      os << *prefix << "Start" << std::endl;
      std::cerr << os.str();
    }

    // note the exportProcIDs and remoteProcIDs _must_ be a list that has
    // an entry for each GID. If the export/remoteProcIDs is taken from
    // the getProcs{From|To} lists that are extracted from a previous distributor,
    // it will generate a wrong answer, because those lists have a unique entry
    // for each processor id. A version of this with lengthsTo and lengthsFrom
    // should be made.

    plan_.howInitialized_ = Tpetra::Details::DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_SENDS_N_RECVS;


    int myProcID = plan_.comm_->getRank ();
    int numProcs = plan_.comm_->getSize();

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

    plan_.sendMessageToSelf_ = ( starts[myProcID] != 0 ) ? 1 : 0;

    plan_.numSendsToOtherProcs_ = 0;

    if( needSendBuff ) //grouped by processor, no send buffer or indicesTo_ needed
      {
        if (starts[0] == 0 ) {
          plan_.numSendsToOtherProcs_ = 0;
        }
        else {
          plan_.numSendsToOtherProcs_ = 1;
        }
        for (Teuchos::Array<size_t>::iterator i=starts.begin()+1,
               im1=starts.begin();
             i != starts.end(); ++i)
          {
            if (*i != 0) ++plan_.numSendsToOtherProcs_;
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

        plan_.indicesTo_.resize(numActive);

        for (size_t i = 0; i < numExportIDs; ++i) {
          if (exportProcIDs[i] >= 0) {
            // record the offset to the sendBuffer for this export
            plan_.indicesTo_[starts[exportProcIDs[i]]] = i;
            // now increment the offset for this proc
            ++starts[exportProcIDs[i]];
          }
        }
        for (int proc = numProcs-1; proc != 0; --proc) {
          starts[proc] = starts[proc-1];
        }
        starts.front() = 0;
        starts[numProcs] = numActive;
        plan_.procIdsToSendTo_.resize(plan_.numSendsToOtherProcs_);
        plan_.startsTo_.resize(plan_.numSendsToOtherProcs_);
        plan_.lengthsTo_.resize(plan_.numSendsToOtherProcs_);
        plan_.maxSendLength_ = 0;
        size_t snd = 0;
        for (int proc = 0; proc < numProcs; ++proc ) {
          if (starts[proc+1] != starts[proc]) {
            plan_.lengthsTo_[snd] = starts[proc+1] - starts[proc];
            plan_.startsTo_[snd] = starts[proc];
            // record max length for all off-proc sends
            if ((proc != myProcID) && (plan_.lengthsTo_[snd] > plan_.maxSendLength_)) {
              plan_.maxSendLength_ = plan_.lengthsTo_[snd];
            }
            plan_.procIdsToSendTo_[snd] = proc;
            ++snd;
          }
        }
      }
    else {
      // grouped by proc, no send buffer or indicesTo_ needed
      plan_.numSendsToOtherProcs_ = 0;
      // Count total number of sends, i.e., total number of procs to
      // which we are sending.  This includes myself, if applicable.
      for (int i = 0; i < numProcs; ++i) {
        if (starts[i]) {
          ++plan_.numSendsToOtherProcs_;
        }
      }

      // Not only do we not need these, but we must clear them, as
      // empty status of indicesTo is a flag used later.
      plan_.indicesTo_.resize(0);
      // Size these to numSendsToOtherProcs_; note, at the moment, numSendsToOtherProcs_
      // includes self sends.  Set their values to zeros.
      plan_.procIdsToSendTo_.assign(plan_.numSendsToOtherProcs_,0);
      plan_.startsTo_.assign(plan_.numSendsToOtherProcs_,0);
      plan_.lengthsTo_.assign(plan_.numSendsToOtherProcs_,0);

      // set startsTo to the offset for each send (i.e., each proc ID)
      // set procsTo to the proc ID for each send
      // in interpreting this code, remember that we are assuming contiguity
      // that is why index skips through the ranks
      {
        size_t index = 0, procIndex = 0;
        for (size_t i = 0; i < plan_.numSendsToOtherProcs_; ++i) {
          while (exportProcIDs[procIndex] < 0) {
            ++procIndex; // skip all negative proc IDs
          }
          plan_.startsTo_[i] = procIndex;
          int procID = exportProcIDs[procIndex];
          plan_.procIdsToSendTo_[i] = procID;
          index     += starts[procID];
          procIndex += starts[procID];
        }
      }
      // sort the startsTo and proc IDs together, in ascending order, according
      // to proc IDs
      if (plan_.numSendsToOtherProcs_ > 0) {
        sort2(plan_.procIdsToSendTo_.begin(), plan_.procIdsToSendTo_.end(), plan_.startsTo_.begin());
      }
      // compute the maximum send length
      plan_.maxSendLength_ = 0;
      for (size_t i = 0; i < plan_.numSendsToOtherProcs_; ++i) {
        int procID = plan_.procIdsToSendTo_[i];
        plan_.lengthsTo_[i] = starts[procID];
        if ((procID != myProcID) && (plan_.lengthsTo_[i] > plan_.maxSendLength_)) {
          plan_.maxSendLength_ = plan_.lengthsTo_[i];
        }
      }
    }


    plan_.numSendsToOtherProcs_ -= plan_.sendMessageToSelf_;
    std::vector<int> recv_list;
    recv_list.reserve(plan_.numSendsToOtherProcs_); //reserve an initial guess for size needed

    int last_pid=-2;
    for(int i=0; i<remoteProcIDs.size(); i++) {
    if(remoteProcIDs[i]>last_pid) {
      recv_list.push_back(remoteProcIDs[i]);
      last_pid = remoteProcIDs[i];
    }
    else if (remoteProcIDs[i]<last_pid)
      throw std::runtime_error("Tpetra::Distributor:::createFromSendsAndRecvs expected RemotePIDs to be in sorted order");
    }
    plan_.numReceives_ = recv_list.size();
    if(plan_.numReceives_) {
      plan_.procsFrom_.assign(plan_.numReceives_,0);
      plan_.lengthsFrom_.assign(plan_.numReceives_,0);
      plan_.indicesFrom_.assign(plan_.numReceives_,0);
      plan_.startsFrom_.assign(plan_.numReceives_,0);
    }
    for(size_t i=0,j=0; i<plan_.numReceives_; ++i) {
      int jlast=j;
      plan_.procsFrom_[i]  = recv_list[i];
      plan_.startsFrom_[i] = j;
      for( ; j<(size_t)remoteProcIDs.size() &&
             remoteProcIDs[jlast]==remoteProcIDs[j]  ; j++){;}
      plan_.lengthsFrom_[i] = j-jlast;
    }
    plan_.totalReceiveLength_ = remoteProcIDs.size();
    plan_.indicesFrom_.clear ();
    plan_.numReceives_-=plan_.sendMessageToSelf_;

    if (verbose_) {
      std::ostringstream os;
      os << *prefix << "Done" << std::endl;
      std::cerr << os.str();
    }
  }

} // namespace Tpetra
