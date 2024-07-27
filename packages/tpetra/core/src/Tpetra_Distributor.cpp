// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    sendTypes.push_back ("Send");
    sendTypes.push_back ("Alltoall");
#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
    sendTypes.push_back ("MpiAdvanceAlltoall");
    sendTypes.push_back ("MpiAdvanceNbralltoallv");
#endif
    return sendTypes;
  }

  Distributor::
  Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
               const Teuchos::RCP<Teuchos::FancyOStream>& /* out */,
               const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : plan_(comm)
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
      plan_.getComm().getRawPtr(), "Distributor", methodName);
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

      // ParameterListAcceptor semantics require pointer identity of the
      // sublist passed to setParameterList(), so we save the pointer.
      this->setMyParamList (plist);

      RCP<ParameterList> planParams(plist);
      planParams->remove("Debug", false);
      planParams->remove("VerboseObject", false);
      plan_.setParameterList(planParams);
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

    const bool debug = tpetraDistributorDebugDefault;

    Array<std::string> sendTypes = distributorSendTypes ();
    const std::string defaultSendType ("Send");
    Array<Details::EDistributorSendType> sendTypeEnums;
    sendTypeEnums.push_back (Details::DISTRIBUTOR_ISEND);
    sendTypeEnums.push_back (Details::DISTRIBUTOR_SEND);
    sendTypeEnums.push_back (Details::DISTRIBUTOR_ALLTOALL);
#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
    sendTypeEnums.push_back(Details::DISTRIBUTOR_MPIADVANCE_ALLTOALL);
    sendTypeEnums.push_back(Details::DISTRIBUTOR_MPIADVANCE_NBRALLTOALLV);
#endif

    RCP<ParameterList> plist = parameterList ("Tpetra::Distributor");
    setStringToIntegralParameter<Details::EDistributorSendType> ("Send type",
      defaultSendType, "When using MPI, the variant of send to use in "
      "do[Reverse]Posts()", sendTypes(), sendTypeEnums(), plist.getRawPtr());
    plist->set ("Debug", debug, "Whether to print copious debugging output on "
                "all processes.");
    plist->set ("Timer Label","","Label for Time Monitor output");

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
  { return plan_.getTotalReceiveLength(); }

  size_t Distributor::getNumReceives() const
  { return plan_.getNumReceives(); }

  bool Distributor::hasSelfMessage() const
  { return plan_.hasSelfMessage(); }

  size_t Distributor::getNumSends() const
  { return plan_.getNumSends(); }

  size_t Distributor::getMaxSendLength() const
  { return plan_.getMaxSendLength(); }

  Teuchos::ArrayView<const int> Distributor::getProcsFrom() const
  { return plan_.getProcsFrom(); }

  Teuchos::ArrayView<const size_t> Distributor::getLengthsFrom() const
  { return plan_.getLengthsFrom(); }

  Teuchos::ArrayView<const int> Distributor::getProcsTo() const
  { return plan_.getProcsTo(); }

  Teuchos::ArrayView<const size_t> Distributor::getLengthsTo() const
  { return plan_.getLengthsTo(); }

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
    reverseDistributor_ = Teuchos::rcp(new Distributor(plan_.getComm()));
    reverseDistributor_->plan_ = *plan_.getReversePlan();
    reverseDistributor_->verbose_ = verbose_;

    // requests_: Allocated on demand.
    // reverseDistributor_: See note below

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
    actor_.doWaits(plan_);
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
        << Details::DistributorHowInitializedEnumToString (plan_.howInitialized())
        << ", Parameters: {"
        << "Send type: "
        << DistributorSendTypeEnumToString (plan_.getSendType())
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
    if (vl <= Teuchos::VERB_LOW || plan_.getComm().is_null ()) {
      return std::string ();
    }

    auto outStringP = Teuchos::rcp (new std::ostringstream ());
    auto outp = Teuchos::getFancyOStream (outStringP); // returns RCP
    Teuchos::FancyOStream& out = *outp;

    const int myRank = plan_.getComm()->getRank ();
    const int numProcs = plan_.getComm()->getSize ();
    out << "Process " << myRank << " of " << numProcs << ":" << endl;
    Teuchos::OSTab tab1 (out);

    out << "selfMessage: " << hasSelfMessage() << endl;
    out << "numSends: " << getNumSends() << endl;
    if (vl == VERB_HIGH || vl == VERB_EXTREME) {
      out << "procsTo: " << toString (plan_.getProcsTo()) << endl;
      out << "lengthsTo: " << toString (plan_.getLengthsTo()) << endl;
      out << "maxSendLength: " << getMaxSendLength() << endl;
    }
    if (vl == VERB_EXTREME) {
      out << "startsTo: " << toString (plan_.getStartsTo()) << endl;
      out << "indicesTo: " << toString (plan_.getIndicesTo()) << endl;
    }
    if (vl == VERB_HIGH || vl == VERB_EXTREME) {
      out << "numReceives: " << getNumReceives() << endl;
      out << "totalReceiveLength: " << getTotalReceiveLength() << endl;
      out << "lengthsFrom: " << toString (plan_.getLengthsFrom()) << endl;
      out << "procsFrom: " << toString (plan_.getProcsFrom()) << endl;
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
    if (plan_.getComm().is_null ()) {
      return;
    }
    const int myRank = plan_.getComm()->getRank ();
    const int numProcs = plan_.getComm()->getSize ();

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
          << Details::DistributorHowInitializedEnumToString (plan_.howInitialized())
          << endl;
      {
        out << "Parameters: " << endl;
        Teuchos::OSTab tab2 (out);
        out << "\"Send type\": "
            << DistributorSendTypeEnumToString (plan_.getSendType()) << endl
            << "\"Debug\": " << (verbose_ ? "true" : "false") << endl;
      }
    } // if myRank == 0

    // This is collective over the Map's communicator.
    if (vl > VERB_LOW) {
      const std::string lclStr = this->localDescribeToString (vl);
      Tpetra::Details::gathervPrint (out, lclStr, *plan_.getComm());
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
    return plan_.createFromSends(exportProcIDs);
  }

  void
  Distributor::
  createFromSendsAndRecvs (const Teuchos::ArrayView<const int>& exportProcIDs,
                           const Teuchos::ArrayView<const int>& remoteProcIDs)
  {
    plan_.createFromSendsAndRecvs(exportProcIDs, remoteProcIDs);
  }

} // namespace Tpetra
