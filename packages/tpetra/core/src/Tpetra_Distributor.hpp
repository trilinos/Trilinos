// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DISTRIBUTOR_HPP
#define TPETRA_DISTRIBUTOR_HPP

#include "Tpetra_Details_DistributorActor.hpp"
#include "Tpetra_Details_DistributorPlan.hpp"

#include "Tpetra_Util.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Tpetra_Details_Behavior.hpp"

#include "KokkosCompat_View.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_TeuchosCommAdapters.hpp"
#include <memory>
#include <sstream>
#include <type_traits>

namespace Tpetra {

  /// \brief Valid values for Distributor's "Send type" parameter.
  ///
  /// This is mainly useful as an implementation detail of
  /// Distributor.  You may use it if you would like a programmatic
  /// way to get all possible values of the "Send type" parameter of
  /// Distributor.
  Teuchos::Array<std::string> distributorSendTypes ();

  /// \class Distributor
  /// \brief Sets up and executes a communication plan for a Tpetra DistObject.
  ///
  /// \note Most Tpetra users do not need to know about this class.
  ///
  /// This class encapsulates the general information and
  /// communication services needed for subclasses of \c DistObject
  /// (such as CrsMatrix and MultiVector) to do data redistribution
  /// (Import and Export) operations.  It is an implementation detail
  /// of Import and Export; in particular; it actually does the
  /// communication.
  ///
  /// Here is the typical way to use this class:
  /// 1. Create a Distributor.  (The constructor is inexpensive.)
  /// 2. Set up the Distributor once, using one of the two "plan
  ///    creation" methods: either createFromSends(), or
  ///    createFromRecvs().  This may be more expensive and
  ///    communication-intensive than Step 3.
  /// 3. Communicate the data by calling doPostsAndWaits() (forward
  ///    mode), or doReversePostsAndWaits() (reverse mode).  You may
  ///    do this multiple times with the same Distributor instance.
  ///
  /// Step 2 is expensive, but you can amortize its cost over multiple
  /// uses of the Distributor for communication (Step 3).  You may
  /// also separate out "posts" (invoking nonblocking communication)
  /// and "waits" (waiting for that communication to complete), by
  /// calling doPosts() (resp. doReversePosts()), then doWaits()
  /// (resp. doReverseWaits()).  This is useful if you have local work
  /// to do between the posts and waits, because it may overlap
  /// communication with computation.  Whether it actually <i>does</i>
  /// overlap, depends on both the MPI implementation and your choice
  /// of parameters for the Distributor.
  ///
  /// Instances of Distributor take the following parameters that
  /// control communication and debug output:
  /// - "Send type" (<tt>std::string</tt>): When using MPI, the
  ///   variant of MPI_Send to use in do[Reverse]Posts().  Valid
  ///   values include "Isend", 
  ///   and "Send".  The
  ///   default is "Send".  (The receive type is always MPI_Irecv, a
  ///   nonblocking receive.  Since we post receives first before
  ///   sends, this prevents deadlock, even if MPI_Send blocks and
  ///   does not buffer.)
  /// - "Debug" (\c bool): If true, print copious debugging output on
  ///   all processes in the Distributor's communicator.  This is
  ///   useful only for debugging Distributor and other Tpetra classes
  ///   that use it (like Import and Export).  If the Distributor was
  ///   created using one of the constructors that takes a
  ///   Teuchos::FancyOStream, it will write debugging output to that
  ///   stream.  Otherwise, it will write debugging output to stderr.
  ///   Currently, the "Debug" parameter overrides "VerboseObject"
  ///   (see below).
  /// - "VerboseObject" (sublist): Optional sublist for controlling
  ///   behavior of Distributor as a Teuchos::VerboseObject.  This is
  ///   currently useful only for debugging.  This sublist takes
  ///   optional parameters "Verbosity Level" (std::string) and
  ///   "Output File" (std::string).  "Verbosity Level" has six valid
  ///   values: "VERB_DEFAULT", "VERB_NONE", "VERB_LOW",
  ///   "VERB_MEDIUM", "VERB_HIGH", and "VERB_EXTREME", with
  ///   increasing verbosity starting with "VERB_NONE".  "Output File"
  ///   is the name of a file to use for output; "none" means don't
  ///   open a file, but write to the default output stream.
  class Distributor :
    public Teuchos::Describable,
    public Teuchos::ParameterListAcceptorDefaultBase {
  public:
    //! @name Constructors and destructor
    //@{

    /// \brief Construct using the specified communicator and default parameters.
    ///
    /// \param comm [in] Communicator used by the Distributor.  MUST
    ///   be nonnull.
    ///
    /// The constructor doesn't actually set up the distribution
    /// pattern.  You need to call one of the "gather / scatter
    /// 'constructors'" to do that.
    explicit Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm);

    /// \brief Construct using the specified communicator and default
    ///   parameters, with an output stream
    ///
    /// \param comm [in] Communicator used by the Distributor.  MUST
    ///   be nonnull.
    /// \param out [in/out] Output stream (for debugging output).  If
    ///   null, the default is \c std::cerr.
    ///
    /// The constructor doesn't actually set up the distribution
    /// pattern.  You need to call one of the "gather / scatter
    /// 'constructors'" to do that.
    Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                 const Teuchos::RCP<Teuchos::FancyOStream>& out);

    /// \brief Construct using the specified communicator and ParameterList.
    ///
    /// \param comm [in] Communicator used by the Distributor.  MUST
    ///   be nonnull.
    /// \param plist [in/out] List of parameters controlling how the
    ///   Distributor performs communication, and specifying debug
    ///   options.  If null, all parameters take their default values.
    ///   Please see the class documentation for a list of all
    ///   accepted parameters and their default values.
    ///
    /// The constructor doesn't actually set up the distribution
    /// pattern.  You need to call one of the "gather / scatter
    /// 'constructors'" to do that.
    Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                 const Teuchos::RCP<Teuchos::ParameterList>& plist);

    /// \brief Construct using the specified communicator and
    ///   ParameterList, with an output stream
    ///
    /// \param comm [in] Communicator used by the Distributor.  MUST
    ///   be nonnull.
    /// \param out [in/out] Output stream (for debugging output).  If
    ///   null, the default is \c std::cerr.
    /// \param plist [in/out] List of parameters controlling how the
    ///   Distributor performs communication, and specifying debug
    ///   options.  If null, all parameters take their default values.
    ///   Please see the class documentation for a list of all
    ///   accepted parameters and their default values.
    ///
    /// The constructor doesn't actually set up the distribution
    /// pattern.  You need to call one of the "gather / scatter
    /// 'constructors'" to do that.
    Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                 const Teuchos::RCP<Teuchos::FancyOStream>& out,
                 const Teuchos::RCP<Teuchos::ParameterList>& plist);

    //! Copy constructor.
    Distributor (const Distributor& distributor);

    /// \brief Destructor (virtual for memory safety).
    ///
    /// \pre No outstanding communication requests.
    ///   (We could check, but see GitHub Issue #1303.)
    virtual ~Distributor () = default;

    /// \brief Swap the contents of rhs with those of *this.
    ///
    /// This is useful in Import's setUnion() method.  It avoids the
    /// overhead of copying arrays, since it can use std::swap on the
    /// arrays.
    void swap (Distributor& rhs);

    //@}
    //! @name Implementation of ParameterListAcceptorDefaultBase
    //@{

    /// \brief Set Distributor parameters.
    ///
    /// Please see the class documentation for a list of all accepted
    /// parameters and their default values.
    void setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist);

    /// \brief List of valid Distributor parameters.
    ///
    /// Please see the class documentation for a list of all accepted
    /// parameters and their default values.
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters () const;

    //@}
    //! \name Gather / scatter "constructors"
    //@{

    /// \brief Set up Distributor using list of process ranks to which
    ///   this process will send.
    ///
    /// Take a list of process ranks and construct a plan for
    /// efficiently scattering to those processes.  Return the number
    /// of processes which will send me (the calling process) data.
    ///
    /// \param exportProcIDs [in] List of ranks of the processes that
    ///   will get the exported data.  If there is a process rank
    ///   greater than or equal to the number of processes, all
    ///   processes will throw an <tt>std::runtime_error</tt>
    ///   exception.  Process ranks less than zero are ignored; their
    ///   placement corresponds to null sends in any future
    ///   exports. That is, if <tt>exportProcIDs[0] == -1</tt>, then
    ///   the corresponding position in the export array is ignored
    ///   during a call to doPosts() or doPostsAndWaits().  For this
    ///   reason, a negative entry is sufficient to break contiguity.
    ///
    /// \return Number of imports this process will be receiving.
    size_t createFromSends (const Teuchos::ArrayView<const int>& exportProcIDs);

    /// \brief Set up Distributor using list of process ranks from which to receive.
    ///
    /// Take a list of process ranks and construct a plan for
    /// efficiently scattering to those processes.  Return the number
    /// and list of IDs being sent by me (the calling process).
    ///
    /// Import invokes this method in order to create a Distributor
    /// from a list of receive neighbors and IDs.  A common use case
    /// for this process is setting up sends and receives for the
    /// remote entries of the source vector in a distributed sparse
    /// matrix-vector multiply.  The Mantevo HPCCG miniapp shows an
    /// annotated and simplified version of this process for that
    /// special case.
    ///
    /// \param remoteIDs [in] List of remote IDs wanted.
    ///
    /// \param remoteProcIDs [in] The ranks of the process that will
    ///   send the remote IDs listed in \c remoteIDs. Process ranks
    ///   less than zero are ignored; their placement corresponds to
    ///   null sends in any future exports.  If there is a process
    ///   rank greater than or equal to the number of processes, all
    ///   processes will throw an <tt>std::runtime_error</tt>
    ///   exception.
    ///
    /// \param exportIDs [out] List of IDs that need to be sent from
    ///   this process.
    ///
    /// \param exportProcIDs [out] The ranks of the processes that
    ///   will get the exported IDs in \c exportIDs.
    ///
    /// The \c exportGIDs and \c exportProcIDs arrays are resized by
    /// the Distributor, which is why they are passed in as a nonconst
    /// Teuchos::Array reference.
    template <class Ordinal>
    void
    createFromRecvs (const Teuchos::ArrayView<const Ordinal>& remoteIDs,
                     const Teuchos::ArrayView<const int>& remoteProcIDs,
                     Teuchos::Array<Ordinal>& exportIDs,
                     Teuchos::Array<int>& exportProcIDs);

    /// \brief Set up Distributor using list of process ranks to which
    ///   to send, and list of process ranks from which to receive.
    ///
    /// \param exportProcIDs [in] List of process ranks to which this
    ///   process must send a message.
    /// \param remoteProcIDs [in] List of process ranks from which
    ///   this process must receive a message.
    void
    createFromSendsAndRecvs (const Teuchos::ArrayView<const int>& exportProcIDs,
                             const Teuchos::ArrayView<const int>& remoteProcIDs);

    //@}
    //! @name Attribute accessor methods
    //@{

    /// \brief The number of processes from which we will receive data.
    ///
    /// The count does <i>not</i> include the calling process.
    size_t getNumReceives() const;

    /// \brief The number of processes to which we will send data.
    ///
    /// The count does <i>not</i> include the calling process.
    size_t getNumSends() const;

    //! Whether the calling process will send or receive messages to itself.
    bool hasSelfMessage() const;

    //! Maximum number of values this process will send to another single process.
    size_t getMaxSendLength() const;

    //! Total number of values this process will receive from other processes.
    size_t getTotalReceiveLength() const;

    /// \brief Ranks of the processes sending values to this process.
    ///
    /// This is a nonpersisting view.  It will last only as long as
    /// this Distributor instance does.
    Teuchos::ArrayView<const int> getProcsFrom() const;

    /// \brief Ranks of the processes to which this process will send values.
    ///
    /// This is a nonpersisting view.  It will last only as long as
    /// this Distributor instance does.
    Teuchos::ArrayView<const int> getProcsTo() const;

    /// \brief Number of values this process will receive from each process.
    ///
    /// This process will receive <tt>getLengthsFrom[i]</tt> values
    /// from process <tt>getProcsFrom[i]</tt>.
    ///
    /// This is a nonpersisting view.  It will last only as long as
    /// this Distributor instance does.
    Teuchos::ArrayView<const size_t> getLengthsFrom() const;

    /// \brief Number of values this process will send to each process.
    ///
    /// This process will send <tt>getLengthsTo[i]</tt> values to
    /// process <tt>getProcsTo[i]</tt>.
    ///
    /// This is a nonpersisting view.  It will last only as long as
    /// this Distributor instance does.
    Teuchos::ArrayView<const size_t> getLengthsTo() const;

    /// \brief Return an enum indicating whether and how a Distributor was initialized.
    ///
    /// This is an implementation detail of Tpetra.  Please do not
    /// call this method or rely on it existing in your code.
    Details::EDistributorHowInitialized howInitialized() const {
      return plan_.howInitialized();
    }

    //@}
    //! @name Reverse communication methods
    //@{

    /// \brief A reverse communication plan Distributor.
    ///
    /// The first time this method is called, it creates a Distributor
    /// with the reverse communication plan of <tt>*this</tt>.  On
    /// subsequent calls, it returns the cached reverse Distributor.
    ///
    /// Most users do not need to call this method.  If you invoke
    /// doReversePosts() or doReversePostsAndWaits(), the reverse
    /// Distributor will be created automatically if it does not yet
    /// exist.
    Teuchos::RCP<Distributor> getReverse(bool create=true) const;

    //@}
    //! @name Methods for executing a communication plan
    //@{

    /// Wait on any outstanding nonblocking message requests to complete.
    ///
    /// This method is for forward mode communication only, that is,
    /// after calling doPosts().  For reverse mode communication
    /// (after calling doReversePosts()), call doReverseWaits()
    /// instead.
    void doWaits ();

    /// Wait on any outstanding nonblocking message requests to complete.
    ///
    /// This method is for reverse mode communication only, that is,
    /// after calling doReversePosts().  For forward mode
    /// communication (after calling doPosts()), call doWaits()
    /// instead.
    void doReverseWaits ();

    /// \brief Execute the (forward) communication plan.
    ///
    /// Call this version of the method when you have the same number
    /// of Packets for each LID (local ID) to send or receive.
    ///
    /// \tparam Packet The type of data to send and receive.
    ///
    /// \param exports [in] Contains the values to be sent by this
    ///   process.  On exit from this method, it's OK to modify the
    ///   entries of this buffer.
    ///
    /// \param numPackets [in] The number of Packets per export /
    ///   import.  This version of the routine assumes that each LID
    ///   has the same number of Packets associated with it.  (\c
    ///   MultiVector is an example of a DistObject subclass
    ///   satisfying this property.)
    ///
    /// \param imports [out] On entry, buffer must be large enough to
    ///   accomodate the data exported (sent) to us.  On exit,
    ///   contains the values exported to us.
    template <class ExpView, class ImpView>
    typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
    doPostsAndWaits (
      const ExpView &exports,
      size_t numPackets,
      const ImpView &imports);

    /// \brief Execute the (forward) communication plan.
    ///
    /// Call this version of the method when you have possibly
    /// different numbers of Packets for each LID (local ID) to send
    /// or receive.
    ///
    /// \tparam Packet The type of data to send and receive.
    ///
    /// \param exports [in] Contains the values to be sent by this
    ///   process.  On exit from this method, it's OK to modify the
    ///   entries of this buffer.
    ///
    /// \param numExportPacketsPerLID [in] The number of packets for
    ///   each export LID (i.e., each LID to be sent).
    ///
    /// \param imports [out] On entry, buffer must be large enough to
    ///   accomodate the data exported (sent) to us.  On exit,
    ///   contains the values exported to us.
    ///
    /// \param numImportPacketsPerLID [in] The number of packets for
    ///   each import LID (i.e., each LID to be received).
    template <class ExpView, class ImpView>
    typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
    doPostsAndWaits (const ExpView &exports,
                     const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                     const ImpView &imports,
                     const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

    /// \brief Post the data for a forward plan, but do not execute the waits yet.
    ///
    /// Call this overload when you have the same number of Packets
    /// for each LID to send or receive.
    ///
    /// \tparam Packet The type of data to send and receive.
    ///
    /// \param exports [in] Contains the values to be sent by this
    ///   process.  This is an ArrayRCP and not an ArrayView so that
    ///   we have the freedom to use nonblocking sends if we wish.  Do
    ///   not modify the data in this array until \c doWaits() has
    ///   completed.
    ///
    /// \param numPackets [in] The number of Packets per export /
    ///   import.  (Same as the three-argument version of
    ///   doPostsAndWaits().)
    ///
    /// \param imports [out] On entry, buffer must be large enough to
    ///   accomodate the data exported (sent) to us.  This is an
    ///   ArrayRCP and not an ArrayView so that we have the freedom to
    ///   use nonblocking sends if we wish.  Do not modify the data in
    ///   this array until \c doWaits() has completed.  Upon
    ///   completion of \c doWaits(), this buffer contains the values
    ///   exported to us.
    template <class ExpView, class ImpView>
    typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
    doPosts (const ExpView &exports,
             size_t numPackets,
             const ImpView &imports);

    /// \brief Post the data for a forward plan, but do not execute the waits yet.
    ///
    /// Call this overload when you have possibly different numbers of
    /// Packets for each LID to send or receive.
    ///
    /// \tparam Packet The type of data to send and receive.
    ///
    /// \param exports [in] Same as in the three-argument version of
    ///   doPosts().
    ///
    /// \param numExportPacketsPerLID [in] Same as in the
    ///   four-argument version of doPostsAndWaits().
    ///
    /// \param imports [out] Same as in the three-argument version of
    ///   doPosts().
    ///
    /// \param numImportPacketsPerLID [in] Same as in the
    ///   four-argument version of doPostsAndWaits().
    template <class ExpView, class ImpView>
    typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
    doPosts (const ExpView &exports,
             const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
             const ImpView &imports,
             const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

    /// \brief Execute the reverse communication plan.
    ///
    /// This method takes the same arguments as the three-argument
    /// version of doPostsAndWaits().
    template <class ExpView, class ImpView>
    typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
    doReversePostsAndWaits (const ExpView &exports,
                            size_t numPackets,
                            const ImpView &imports);

    /// \brief Execute the reverse communication plan.
    ///
    /// This method takes the same arguments as the four-argument
    /// version of doPostsAndWaits().
    template <class ExpView, class ImpView>
    typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
    doReversePostsAndWaits (const ExpView &exports,
                            const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                            const ImpView &imports,
                            const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

    /// \brief Post the data for a reverse plan, but do not execute the waits yet.
    ///
    /// This method takes the same arguments as the three-argument
    /// version of doPosts().
    template <class ExpView, class ImpView>
    typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
    doReversePosts (const ExpView &exports,
                    size_t numPackets,
                    const ImpView &imports);

    /// \brief Post the data for a reverse plan, but do not execute the waits yet.
    ///
    /// This method takes the same arguments as the four-argument
    /// version of doPosts().
    template <class ExpView, class ImpView>
    typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
    doReversePosts (const ExpView &exports,
                    const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                    const ImpView &imports,
                    const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

    //@}
    //! @name Implementation of Teuchos::Describable
    //@{

    //! Return a one-line description of this object.
    std::string description() const;

    /// \brief Describe this object in a human-readable way to the
    ///   given output stream.
    ///
    /// You must call this method as a collective over all processes
    /// in this object's communicator.
    ///
    /// \param out [out] Output stream to which to write.  Only
    ///   Process 0 in this object's communicator may write to the
    ///   output stream.
    ///
    /// \param verbLevel [in] Verbosity level.  This also controls
    ///   whether this method does any communication.  At verbosity
    ///   levels higher (greater) than Teuchos::VERB_LOW, this method
    ///   behaves as a collective over the object's communicator.
    ///
    /// Teuchos::FancyOStream wraps std::ostream.  It adds features
    /// like tab levels.  If you just want to wrap std::cout, try
    /// this:
    /// \code
    /// auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::out));
    /// \endcode
    void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel =
                Teuchos::Describable::verbLevel_default) const;
    //@}

    /// \brief Get this Distributor's DistributorPlan
    ///
    /// FIXME: Delete this method when it's no longer needed for non-blocking
    ///        communication in the DistObject
    const Details::DistributorPlan& getPlan() const { return plan_; }
  private:
    Details::DistributorPlan plan_;
    Details::DistributorActor actor_;

    //! @name Parameters read in from the Teuchos::ParameterList
    //@{

    //! Get default value of verbose_ (see below).
    static bool getVerbose();

    /// \brief Get prefix for verbose debug output.
    ///
    /// \brief methodName [in] Name of the method in which you want to
    ///   print verbose debug output.
    std::unique_ptr<std::string>
    createPrefix(const char methodName[]) const;

    //! Whether to print copious debug output to stderr on all processes.
    bool verbose_ = getVerbose();
    //@}

    /// \brief The reverse distributor.
    ///
    /// This is created on demand in getReverse() and cached for
    /// later reuse.  This is why it is declared "mutable".
    mutable Teuchos::RCP<Distributor> reverseDistributor_;

    /// \brief Compute send (GID,PID) pairs from receive (GID,PID) pairs.
    ///
    /// GID means "global ID" and PID means "process ID" (rank, in MPI
    /// terms).
    ///
    /// \param importGIDs [in] GIDs to receive by my process.
    /// \param importProcIDs [in] Process IDs from which to receive by
    ///   my process.
    /// \param exportGIDs [out] GIDs to send by my process.  Resized if
    ///   necessary.
    /// \param exportProcIDs [out] Process IDs to which to send by my
    ///   process.  Resized if necessary.
    template <class Ordinal>
    void computeSends (const Teuchos::ArrayView<const Ordinal> &remoteGIDs,
                       const Teuchos::ArrayView<const int> &remoteProcIDs,
                       Teuchos::Array<Ordinal> &exportGIDs,
                       Teuchos::Array<int> &exportProcIDs);

    //! Create a distributor for the reverse communication pattern.
    void createReverseDistributor() const;


    /// \brief Print the calling process' verbose describe()
    ///   information to the given output string.
    ///
    /// This is an implementation detail of describe().
    std::string
    localDescribeToString (const Teuchos::EVerbosityLevel vl) const;
  }; // class Distributor

  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
  Distributor::
  doPostsAndWaits (const ExpView& exports,
                   size_t numPackets,
                   const ImpView& imports)
  {
    actor_.doPostsAndWaits(plan_, exports, numPackets, imports);
  }

  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
  Distributor::
  doPostsAndWaits(const ExpView& exports,
                  const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                  const ImpView& imports,
                  const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID)
  {
    actor_.doPostsAndWaits(plan_, exports, numExportPacketsPerLID, imports, numImportPacketsPerLID);
  }


  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
  Distributor::
  doPosts (const ExpView &exports,
           size_t numPackets,
           const ImpView &imports)
  {
    actor_.doPosts(plan_, exports, numPackets, imports);
  }

  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
  Distributor::
  doPosts (const ExpView &exports,
           const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
           const ImpView &imports,
           const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID)
  {
    actor_.doPosts(plan_, exports, numExportPacketsPerLID, imports, numImportPacketsPerLID);
  }

  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
  Distributor::
  doReversePostsAndWaits (const ExpView& exports,
                          size_t numPackets,
                          const ImpView& imports)
  {
    doReversePosts (exports, numPackets, imports);
    doReverseWaits ();
  }

  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
  Distributor::
  doReversePostsAndWaits (const ExpView& exports,
                          const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                          const ImpView& imports,
                          const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID)
  {
    doReversePosts (exports, numExportPacketsPerLID, imports,
                    numImportPacketsPerLID);
    doReverseWaits ();
  }

  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
  Distributor::
  doReversePosts (const ExpView &exports,
                  size_t numPackets,
                  const  ImpView &imports)
  {
    // FIXME (mfh 29 Mar 2012) WHY?
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! plan_.getIndicesTo().is_null(), std::runtime_error,
      "Tpetra::Distributor::doReversePosts(3 args): Can only do "
      "reverse communication when original data are blocked by process.");
    if (reverseDistributor_.is_null ()) {
      createReverseDistributor ();
    }
    reverseDistributor_->doPosts (exports, numPackets, imports);
  }

  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::is_view<ExpView>::value && Kokkos::is_view<ImpView>::value)>::type
  Distributor::
  doReversePosts (const ExpView &exports,
                  const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                  const ImpView &imports,
                  const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID)
  {
    // FIXME (mfh 29 Mar 2012) WHY?
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! plan_.getIndicesTo().is_null(), std::runtime_error,
      "Tpetra::Distributor::doReversePosts(3 args): Can only do "
      "reverse communication when original data are blocked by process.");
    if (reverseDistributor_.is_null ()) {
      createReverseDistributor ();
    }
    reverseDistributor_->doPosts (exports, numExportPacketsPerLID,
                                  imports, numImportPacketsPerLID);
  }

  template <class OrdinalType>
  void Distributor::
  computeSends(const Teuchos::ArrayView<const OrdinalType>& importGIDs,
               const Teuchos::ArrayView<const int>& importProcIDs,
               Teuchos::Array<OrdinalType>& exportGIDs,
               Teuchos::Array<int>& exportProcIDs)
  {
    // NOTE (mfh 19 Apr 2012): There was a note on this code saying:
    // "assumes that size_t >= Ordinal".  The code certainly does
    // assume that sizeof(size_t) >= sizeof(OrdinalType) as well as
    // sizeof(size_t) >= sizeof(int).  This is because it casts the
    // OrdinalType elements of importGIDs (along with their
    // corresponding process IDs, as int) to size_t, and does a
    // doPostsAndWaits<size_t>() to send the packed data.
    using Teuchos::ArrayView;
    using std::endl;
    using size_type = typename ArrayView<const OrdinalType>::size_type;
    const char errPrefix[] = "Tpetra::Distributor::computeSends: ";
    const char suffix[] =
      "  Please report this bug to the Tpetra developers.";

    const int myRank = plan_.getComm()->getRank ();

    TEUCHOS_TEST_FOR_EXCEPTION
      (importGIDs.size () != importProcIDs.size (),
       std::invalid_argument, errPrefix << "On Process " << myRank
       << ": importProcIDs.size()=" << importProcIDs.size()
       << " != importGIDs.size()=" << importGIDs.size() << ".");

    const size_type numImports = importProcIDs.size();
    Kokkos::View<size_t*, Kokkos::HostSpace> importObjs("importObjs", 2*numImports);
    // Pack pairs (importGIDs[i], my process ID) to send into importObjs.
    for (size_type i = 0; i < numImports; ++i) {
      importObjs[2*i]   = static_cast<size_t>(importGIDs[i]);
      importObjs[2*i+1] = static_cast<size_t>(myRank);
    }
    //
    // Use a temporary Distributor to send the (importGIDs[i], myRank)
    // pairs to importProcIDs[i].
    //
    Distributor tempPlan(plan_.getComm());
    // mfh 20 Mar 2014: An extra-cautious cast from unsigned to
    // signed, in order to forestall any possible causes for Bug 6069.
    const size_t numExportsAsSizeT =
      tempPlan.createFromSends(importProcIDs);
    const size_type numExports =
      static_cast<size_type>(numExportsAsSizeT);
    TEUCHOS_TEST_FOR_EXCEPTION
      (numExports < 0, std::logic_error, errPrefix <<
       "tempPlan.createFromSends() returned numExports="
       << numExportsAsSizeT << " as a size_t, which overflows to "
       << numExports << " when cast to " <<
       Teuchos::TypeNameTraits<size_type>::name () << "." << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION
      (size_type(tempPlan.getTotalReceiveLength()) != numExports,
       std::logic_error, errPrefix << "tempPlan.getTotalReceiveLength()="
       << tempPlan.getTotalReceiveLength () << " != numExports="
       << numExports  << "." << suffix);

    if (numExports > 0) {
      exportGIDs.resize(numExports);
      exportProcIDs.resize(numExports);
    }

    // exportObjs: Packed receive buffer.  (exportObjs[2*i],
    // exportObjs[2*i+1]) will give the (GID, PID) pair for export i,
    // after tempPlan.doPostsAndWaits(...) finishes below.
    //
    // FIXME (mfh 19 Mar 2014) This only works if OrdinalType fits in
    // size_t.  This issue might come up, for example, on a 32-bit
    // machine using 64-bit global indices.  I will add a check here
    // for that case.
    static_assert(sizeof(size_t) >= sizeof(OrdinalType),
      "Tpetra::Distributor::computeSends: "
      "sizeof(size_t) < sizeof(OrdinalType).");

    TEUCHOS_TEST_FOR_EXCEPTION
      (tempPlan.getTotalReceiveLength () < size_t(numExports),
       std::logic_error,
       errPrefix << "tempPlan.getTotalReceiveLength()="
       << tempPlan.getTotalReceiveLength() << " < numExports="
       << numExports << "." << suffix);

    Kokkos::View<size_t*, Kokkos::HostSpace> exportObjs("exportObjs", tempPlan.getTotalReceiveLength() * 2);
    tempPlan.doPostsAndWaits(importObjs, 2, exportObjs);

    // Unpack received (GID, PID) pairs into exportIDs resp. exportProcIDs.
    for (size_type i = 0; i < numExports; ++i) {
      exportGIDs[i] = static_cast<OrdinalType> (exportObjs[2*i]);
      exportProcIDs[i] = static_cast<int> (exportObjs[2*i+1]);
    }
  }

  template <class OrdinalType>
  void Distributor::
  createFromRecvs (const Teuchos::ArrayView<const OrdinalType> &remoteGIDs,
                   const Teuchos::ArrayView<const int> &remoteProcIDs,
                   Teuchos::Array<OrdinalType> &exportGIDs,
                   Teuchos::Array<int> &exportProcIDs)
  {
    using std::endl;
    const char errPrefix[] = "Tpetra::Distributor::createFromRecvs: ";
    const int myRank = plan_.getComm()->getRank();

    std::unique_ptr<std::string> prefix;
    if (verbose_) {
      prefix = createPrefix("createFromRecvs");
      std::ostringstream os;
      os << *prefix << "Start" << endl;
      std::cerr << os.str();
    }

    const bool debug = Details::Behavior::debug("Distributor");
    if (debug) {
      using Teuchos::outArg;
      using Teuchos::REDUCE_MAX;
      using Teuchos::reduceAll;
      // In debug mode, first test locally, then do an all-reduce to
      // make sure that all processes passed.
      const int errProc =
        (remoteGIDs.size () != remoteProcIDs.size ()) ? myRank : -1;
      int maxErrProc = -1;
      reduceAll(*plan_.getComm(), REDUCE_MAX, errProc, outArg(maxErrProc));
      TEUCHOS_TEST_FOR_EXCEPTION
        (maxErrProc != -1, std::runtime_error, errPrefix << "Lists "
         "of remote IDs and remote process IDs must have the same "
         "size on all participating processes.  Maximum process ID "
         "with error: " << maxErrProc << ".");
    }
    else { // in non-debug mode, just test locally
      // NOTE (mfh 13 Feb 2020) This needs to throw std::runtime_error
      // in order to make an existing Distributor unit test pass.
      TEUCHOS_TEST_FOR_EXCEPTION
        (remoteGIDs.size() != remoteProcIDs.size(), std::runtime_error,
         errPrefix << "On Process " << myRank << ": "
         "remoteGIDs.size()=" << remoteGIDs.size() <<
         " != remoteProcIDs.size()=" << remoteProcIDs.size() << ".");
    }

    computeSends(remoteGIDs, remoteProcIDs, exportGIDs, exportProcIDs);

    plan_.createFromRecvs(remoteProcIDs);

    if (verbose_) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
  }

} // namespace Tpetra

#endif // TPETRA_DISTRIBUTOR_HPP
