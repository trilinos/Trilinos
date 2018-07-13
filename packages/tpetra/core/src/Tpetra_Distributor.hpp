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

#ifndef TPETRA_DISTRIBUTOR_HPP
#define TPETRA_DISTRIBUTOR_HPP

#include "Tpetra_Util.hpp"
#include <Teuchos_as.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Tpetra_Details_Behavior.hpp>

// If TPETRA_DISTRIBUTOR_TIMERS is defined, Distributor will time
// doPosts (both versions) and doWaits, and register those timers with
// Teuchos::TimeMonitor so that summarize() or report() will show
// results.

// #ifndef TPETRA_DISTRIBUTOR_TIMERS
// #  define TPETRA_DISTRIBUTOR_TIMERS 1
// #endif // TPETRA_DISTRIBUTOR_TIMERS

#ifdef TPETRA_DISTRIBUTOR_TIMERS
#  undef TPETRA_DISTRIBUTOR_TIMERS
#endif // TPETRA_DISTRIBUTOR_TIMERS

#include "KokkosCompat_View.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_TeuchosCommAdapters.hpp"
#include <type_traits>

namespace Tpetra {

  namespace Details {
    /// \brief The type of MPI send that Distributor should use.
    ///
    /// This is an implementation detail of Distributor.  Please do
    /// not rely on these values in your code.
    enum EDistributorSendType {
      DISTRIBUTOR_ISEND, // Use MPI_Isend (Teuchos::isend)
      DISTRIBUTOR_RSEND, // Use MPI_Rsend (Teuchos::readySend)
      DISTRIBUTOR_SEND,  // Use MPI_Send (Teuchos::send)
      DISTRIBUTOR_SSEND  // Use MPI_Ssend (Teuchos::ssend)
    };

    /// \brief Convert an EDistributorSendType enum value to a string.
    ///
    /// This is an implementation detail of Distributor.  Please do
    /// not rely on this function in your code.
    std::string
    DistributorSendTypeEnumToString (EDistributorSendType sendType);

    /// \brief Enum indicating how and whether a Distributor was initialized.
    ///
    /// This is an implementation detail of Distributor.  Please do
    /// not rely on these values in your code.
    enum EDistributorHowInitialized {
      DISTRIBUTOR_NOT_INITIALIZED, // Not initialized yet
      DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_SENDS, // By createFromSends
      DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_RECVS, // By createFromRecvs
      DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_SENDS_N_RECVS, // By createFromSendsAndRecvs
      DISTRIBUTOR_INITIALIZED_BY_REVERSE, // By createReverseDistributor
      DISTRIBUTOR_INITIALIZED_BY_COPY, // By copy constructor
    };

    /// \brief Convert an EDistributorHowInitialized enum value to a string.
    ///
    /// This is an implementation detail of Distributor.  Please do
    /// not rely on this function in your code.
    std::string
    DistributorHowInitializedEnumToString (EDistributorHowInitialized how);

  } // namespace Details

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
  /// - "Barrier between receives and sends" (<tt>bool</tt>):
  ///   Whether to execute a barrier between receives and sends in
  ///   do[Reverse]Posts().  A barrier is required for correctness
  ///   when the "Send type" parameter is "Rsend".  Otherwise, a
  ///   barrier is correct and may be useful for debugging, but not
  ///   recommended, since it introduces useless synchronization.
  /// - "Send type" (<tt>std::string</tt>): When using MPI, the
  ///   variant of MPI_Send to use in do[Reverse]Posts().  Valid
  ///   values include "Isend", "Rsend", "Send", and "Ssend".  The
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

    //! Destructor (virtual for memory safety).
    virtual ~Distributor ();

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
    Details::EDistributorHowInitialized howInitialized () const {
      return howInitialized_;
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
    Teuchos::RCP<Distributor> getReverse() const;

    //@}
    //! @name Methods for executing a communication plan
    //@{

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
    template <class Packet>
    void
    doPostsAndWaits (const Teuchos::ArrayView<const Packet> &exports,
                     size_t numPackets,
                     const Teuchos::ArrayView<Packet> &imports);

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
    template <class Packet>
    void
    doPostsAndWaits (const Teuchos::ArrayView<const Packet> &exports,
                     const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                     const Teuchos::ArrayView<Packet> &imports,
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
    template <class Packet>
    void
    doPosts (const Teuchos::ArrayRCP<const Packet> &exports,
             size_t numPackets,
             const Teuchos::ArrayRCP<Packet> &imports);

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
    template <class Packet>
    void
    doPosts (const Teuchos::ArrayRCP<const Packet> &exports,
             const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
             const Teuchos::ArrayRCP<Packet> &imports,
             const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

    /// Wait on any outstanding nonblocking message requests to complete.
    ///
    /// This method is for forward mode communication only, that is,
    /// after calling doPosts().  For reverse mode communication
    /// (after calling doReversePosts()), call doReverseWaits()
    /// instead.
    void doWaits ();

    /// \brief Execute the reverse communication plan.
    ///
    /// This method takes the same arguments as the three-argument
    /// version of doPostsAndWaits().
    template <class Packet>
    void
    doReversePostsAndWaits (const Teuchos::ArrayView<const Packet> &exports,
                            size_t numPackets,
                            const Teuchos::ArrayView<Packet> &imports);

    /// \brief Execute the reverse communication plan.
    ///
    /// This method takes the same arguments as the four-argument
    /// version of doPostsAndWaits().
    template <class Packet>
    void
    doReversePostsAndWaits (const Teuchos::ArrayView<const Packet> &exports,
                            const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                            const Teuchos::ArrayView<Packet> &imports,
                            const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

    /// \brief Post the data for a reverse plan, but do not execute the waits yet.
    ///
    /// This method takes the same arguments as the three-argument
    /// version of doPosts().
    template <class Packet>
    void
    doReversePosts (const Teuchos::ArrayRCP<const Packet> &exports,
                    size_t numPackets,
                    const Teuchos::ArrayRCP<Packet> &imports);

    /// \brief Post the data for a reverse plan, but do not execute the waits yet.
    ///
    /// This method takes the same arguments as the four-argument
    /// version of doPosts().
    template <class Packet>
    void
    doReversePosts (const Teuchos::ArrayRCP<const Packet> &exports,
                    const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                    const Teuchos::ArrayRCP<Packet> &imports,
                    const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

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
    typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
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
    typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
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
    typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
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
    typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
    doPosts (const ExpView &exports,
             const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
             const ImpView &imports,
             const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

    /// \brief Execute the reverse communication plan.
    ///
    /// This method takes the same arguments as the three-argument
    /// version of doPostsAndWaits().
    template <class ExpView, class ImpView>
    typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
    doReversePostsAndWaits (const ExpView &exports,
                            size_t numPackets,
                            const ImpView &imports);

    /// \brief Execute the reverse communication plan.
    ///
    /// This method takes the same arguments as the four-argument
    /// version of doPostsAndWaits().
    template <class ExpView, class ImpView>
    typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
    doReversePostsAndWaits (const ExpView &exports,
                            const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                            const ImpView &imports,
                            const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

    /// \brief Post the data for a reverse plan, but do not execute the waits yet.
    ///
    /// This method takes the same arguments as the three-argument
    /// version of doPosts().
    template <class ExpView, class ImpView>
    typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
    doReversePosts (const ExpView &exports,
                    size_t numPackets,
                    const ImpView &imports);

    /// \brief Post the data for a reverse plan, but do not execute the waits yet.
    ///
    /// This method takes the same arguments as the four-argument
    /// version of doPosts().
    template <class ExpView, class ImpView>
    typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
    doReversePosts (const ExpView &exports,
                    const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                    const ImpView &imports,
                    const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

    /// \brief Information on the last call to do/doReverse
    ///
    /// Returns the amount of data sent & recv'd on this processor in bytes
    void getLastDoStatistics(size_t & bytes_sent,  size_t & bytes_recvd) const{
      bytes_sent  = lastRoundBytesSend_;
      bytes_recvd = lastRoundBytesRecv_;
    }

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

  private:
    //! The communicator over which to perform distributions.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;

    //! Output stream for debug output.
    Teuchos::RCP<Teuchos::FancyOStream> out_;

    //! How the Distributor was initialized (if it was).
    Details::EDistributorHowInitialized howInitialized_;

    //! @name Parameters read in from the Teuchos::ParameterList
    //@{

    //! The variant of send to use in do[Reverse]Posts().
    Details::EDistributorSendType sendType_;

    //! Whether to do a barrier between receives and sends in do[Reverse]Posts().
    bool barrierBetween_;

    //! Whether to print copious debug output to stderr on all processes.
    bool debug_;
    //@}

    /// \brief Whether I am supposed to send a message to myself.
    ///
    /// This is set in createFromSends or createReverseDistributor.
    bool selfMessage_;

    /// \brief The number of sends from this process to other process.
    ///
    /// This is always less than or equal to the number of processes.
    /// It does <i>not</i> count self receives (that is, messages from
    /// the calling process to itself).
    ///
    /// This value is computed by the createFromSends() method.  That
    /// method first includes self receives in the count, but at the
    /// end subtracts one if selfMessage_ is true.
    size_t numSends_;

    /// \brief List of process IDs to which to send.
    ///
    /// This array has length numSends_ + selfMessage_ (that is, it
    /// includes the self message, if there is one).
    Teuchos::Array<int> procsTo_;

    /// \brief Starting index of the block of Packets to send to each process.
    ///
    /// Given an export buffer that contains all of the data being
    /// sent by this process, the block of Packets to send to process
    /// p will start at position startsTo_[p].
    ///
    /// This array has length numSends_ + selfMessage_ (that is, it
    /// includes the self message, if there is one).
    Teuchos::Array<size_t> startsTo_;

    /// \brief Length (in number of Packets) of my process' send to each process.
    ///
    /// lengthsTo_[p] is the length of my process' send to process p.
    /// This array has length numSends_ + selfMessage_ (that is, it
    /// includes the self message, if there is one).
    Teuchos::Array<size_t> lengthsTo_;

    /// \brief The maximum send length (in number of Packets) to another process.
    ///
    /// maxSendLength_ = max(lengthsTo_[p]) for p != my process rank.
    size_t maxSendLength_;

    /// \brief Offset (by message, not by number of Packets) into exports array.
    ///
    /// This array is used by both versions of doPosts().  In that
    /// method, <tt>indicesTo_[j]*numPackets</tt> is the offset into
    /// the <tt>exports</tt> array, where <tt>j = startsTo_[p]</tt>
    /// and p is an index iterating through the sends in reverse order
    /// (starting with the process rank right before the self message,
    /// if there is a self message, else the largest process rank to
    /// which this process sends).
    ///
    /// This array is only used if export data are not blocked (laid
    /// out) by process rank, that is, if we need to use a send
    /// buffer.  Otherwise, this array has no entries.  (In fact,
    /// Distributor currently uses this in both overloads of doPosts()
    /// to test whether data are laid out by process.)
    Teuchos::Array<size_t> indicesTo_;

    /// \brief The number of messages received by my process from other processes.
    ///
    /// This does <i>not</i> count self receives.  If selfMessage_ is
    /// true, the actual number of receives is one more (we assume
    /// that we only receive zero or one messages from ourself).
    ///
    /// This value is computed by the \c computeReceives() method.
    /// That method first includes self receives in the count, but at
    /// the end subtracts one if selfMessage_ is true.
    size_t numReceives_;

    /// \brief sum(lengthsFrom_)
    ///
    /// This is computed by \c createFromSends() and is used to
    /// allocate the receive buffer.  The reverse communicator's total
    /// receive length is the total send length of the forward
    /// communicator.
    size_t totalReceiveLength_;

    /// \brief Array of lengths of incoming messages.
    ///
    /// This array has length numReceives_ + selfMessage_.  Incoming
    /// message i from process procsFrom_[i] has length
    /// lengthsFrom_[i].
    Teuchos::Array<size_t> lengthsFrom_;

    /// \brief Array of ranks of the process from which the calling
    ///   process will receive a message.
    ///
    /// This array has length numReceives_ + selfMessage_.  Incoming
    /// message i was sent by process procsFrom_[i].
    Teuchos::Array<int> procsFrom_;

    /// \brief Array of offsets of incoming messages.
    ///
    /// This array has length numReceives_ + selfMessage_.  It is an
    /// exclusive prefix sum of lengthsFrom_.  It is only used for
    /// constructing the reverse Distributor.
    Teuchos::Array<size_t> startsFrom_;

    /// \brief List that becomes the reverse communicator's indicesTo_.
    ///
    /// Array of length totalReceiveLength_.  Allocated and filled in
    /// computeReceives() as [0, 1, ..., totalReceiveLength_-1].  When
    /// creating the reverse Distributor, this is assigned to the
    /// reverse Distributor's indicesTo_.
    Teuchos::Array<size_t> indicesFrom_;

    /// \brief Communication requests associated with nonblocking receives and sends.
    ///
    /// \note To implementers: Distributor uses requests_.size() as
    ///   the number of outstanding nonblocking receives and sends.
    ///   This means you should always resize to zero after completing
    ///   receive and send requests.
    Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest<int> > > requests_;

    /// \brief The reverse distributor.
    ///
    /// This is created on demand in getReverse() and cached for
    /// later reuse.  This is why it is declared "mutable".
    mutable Teuchos::RCP<Distributor> reverseDistributor_;

    /// \brief The number of bytes sent by this proc in the last call to do/doReverse
    size_t lastRoundBytesSend_;

    /// \brief The number of bytes received by this proc in the last call to do/doReverse
    size_t lastRoundBytesRecv_;

#ifdef TPETRA_DISTRIBUTOR_TIMERS
    Teuchos::RCP<Teuchos::Time> timer_doPosts3_;
    Teuchos::RCP<Teuchos::Time> timer_doPosts4_;
    Teuchos::RCP<Teuchos::Time> timer_doWaits_;
    Teuchos::RCP<Teuchos::Time> timer_doPosts3_recvs_;
    Teuchos::RCP<Teuchos::Time> timer_doPosts4_recvs_;
    Teuchos::RCP<Teuchos::Time> timer_doPosts3_barrier_;
    Teuchos::RCP<Teuchos::Time> timer_doPosts4_barrier_;
    Teuchos::RCP<Teuchos::Time> timer_doPosts3_sends_;
    Teuchos::RCP<Teuchos::Time> timer_doPosts4_sends_;

    //! Make the instance's timers.  (Call only in constructor.)
    void makeTimers ();
#endif // TPETRA_DISTRIBUTOR_TIMERS

    /// \brief Whether to use different tags for different code paths.
    ///
    /// There are currently three code paths in Distributor that post
    /// receives and sends:
    ///
    /// 1. Three-argument variant of doPosts()
    /// 2. Four-argument variant of doPosts()
    /// 3. computeReceives()
    ///
    /// If this option is true, Distributor will use a distinct
    /// message tag for each of these paths.
    bool useDistinctTags_;

    //! Get the tag to use for receives and sends.
    ///
    /// See useDistinctTags_.  This is called in doPosts() (both
    /// variants) and computeReceives().
    int getTag (const int pathTag) const;

    /// \brief Initialize using the specified communicator and ParameterList.
    ///
    /// This method is only meant to be called by the constructor.
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
    /// This method doesn't actually set up the distribution pattern.
    /// You need to call one of the "gather / scatter 'constructors'"
    /// to do that.
    void
    init (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
          const Teuchos::RCP<Teuchos::FancyOStream>& out,
          const Teuchos::RCP<Teuchos::ParameterList>& plist);

    /// \brief Compute receive info from sends.
    ///
    /// This method computes numReceives_, lengthsFrom_, procsFrom_,
    /// totalReceiveLength_, indicesFrom_, and startsFrom_.
    ///
    /// \note This method currently ignores the sendType_ and
    ///   barrierBetween_ parameters, and always uses ireceive() /
    ///   send() for communication of the process IDs from which our
    ///   process is receiving and their corresponding receive packet
    ///   counts.
    void computeReceives ();

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


  template <class Packet>
  void Distributor::
  doPostsAndWaits (const Teuchos::ArrayView<const Packet>& exports,
                   size_t numPackets,
                   const Teuchos::ArrayView<Packet>& imports)
  {
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;
    typedef typename ArrayRCP<const Packet>::size_type size_type;

    TEUCHOS_TEST_FOR_EXCEPTION(
      requests_.size () != 0, std::runtime_error, "Tpetra::Distributor::"
      "doPostsAndWaits(3 args): There are " << requests_.size () <<
      " outstanding nonblocking messages pending.  It is incorrect to call "
      "this method with posts outstanding.");

    // doPosts() accepts the exports and imports arrays as ArrayRCPs,
    // requiring that the memory location is persisting (as is
    // necessary for nonblocking receives).  However, it need only
    // persist until doWaits() completes, so it is safe for us to use
    // a nonpersisting reference in this case.  The use of a
    // nonpersisting reference is purely a performance optimization.

    //const Packet* exportsPtr = exports.getRawPtr();
    //ArrayRCP<const Packet> exportsArcp (exportsPtr, static_cast<size_type> (0),
    //                                    exports.size(), false);
    ArrayRCP<const Packet> exportsArcp (exports.getRawPtr (),
                                        static_cast<size_type> (0),
                                        exports.size(), false);

    // For some reason, neither of the options below (that use arcp)
    // compile for Packet=std::complex<double> with GCC 4.5.1.  The
    // issue only arises with the exports array.  This is why we
    // construct a separate nonowning ArrayRCP.

    // doPosts (arcp<const Packet> (exports.getRawPtr(), 0, exports.size(), false),
    //              numPackets,
    //              arcp<Packet> (imports.getRawPtr(), 0, imports.size(), false));
    // doPosts (arcp<const Packet> (exportsPtr, 0, exports.size(), false),
    //              numPackets,
    //              arcp<Packet> (imports.getRawPtr(), 0, imports.size(), false));
    doPosts (exportsArcp,
             numPackets,
             arcp<Packet> (imports.getRawPtr (), 0, imports.size (), false));
    doWaits ();

    lastRoundBytesSend_ = exports.size () * sizeof (Packet);
    lastRoundBytesRecv_ = imports.size () * sizeof (Packet);
  }

  template <class Packet>
  void Distributor::
  doPostsAndWaits (const Teuchos::ArrayView<const Packet>& exports,
                   const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                   const Teuchos::ArrayView<Packet> &imports,
                   const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID)
  {
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;

    TEUCHOS_TEST_FOR_EXCEPTION(
      requests_.size () != 0, std::runtime_error, "Tpetra::Distributor::"
      "doPostsAndWaits: There are " << requests_.size () << " outstanding "
      "nonblocking messages pending.  It is incorrect to call doPostsAndWaits "
      "with posts outstanding.");

    // doPosts() accepts the exports and imports arrays as ArrayRCPs,
    // requiring that the memory location is persisting (as is
    // necessary for nonblocking receives).  However, it need only
    // persist until doWaits() completes, so it is safe for us to use
    // a nonpersisting reference in this case.

    // mfh 04 Apr 2012: For some reason, calling arcp<const Packet>
    // for Packet=std::complex<T> (e.g., T=float) fails to compile
    // with some versions of GCC.  The issue only arises with the
    // exports array.  This is why we construct a separate nonowning
    // ArrayRCP.
    typedef typename ArrayRCP<const Packet>::size_type size_type;
    ArrayRCP<const Packet> exportsArcp (exports.getRawPtr (),
                                        static_cast<size_type> (0),
                                        exports.size (), false);
    // mfh 04 Apr 2012: This is the offending code.  This statement
    // would normally be in place of "exportsArcp" in the
    // doPosts() call below.
    //arcp<const Packet> (exports.getRawPtr(), 0, exports.size(), false),
    doPosts (exportsArcp,
             numExportPacketsPerLID,
             arcp<Packet> (imports.getRawPtr (), 0, imports.size (), false),
             numImportPacketsPerLID);
    doWaits ();

    lastRoundBytesSend_ = exports.size () * sizeof (Packet);
    lastRoundBytesRecv_ = imports.size () * sizeof (Packet);
  }


  template <class Packet>
  void Distributor::
  doPosts (const Teuchos::ArrayRCP<const Packet>& exports,
           size_t numPackets,
           const Teuchos::ArrayRCP<Packet>& imports)
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::FancyOStream;
    using Teuchos::includesVerbLevel;
    using Teuchos::ireceive;
    using Teuchos::isend;
    using Teuchos::OSTab;
    using Teuchos::readySend;
    using Teuchos::send;
    using Teuchos::ssend;
    using Teuchos::TypeNameTraits;
    using Teuchos::typeName;
    using std::endl;
    typedef Array<size_t>::size_type size_type;

    const bool verbose = Tpetra::Details::Behavior::verbose("Distributor");

#ifdef TPETRA_DISTRIBUTOR_TIMERS
    Teuchos::TimeMonitor timeMon (*timer_doPosts3_);
#endif // TPETRA_DISTRIBUTOR_TIMERS

    const int myRank = comm_->getRank ();
    // Run-time configurable parameters that come from the input
    // ParameterList set by setParameterList().
    const Details::EDistributorSendType sendType = sendType_;
    const bool doBarrier = barrierBetween_;

    Teuchos::OSTab tab0 (out_);
    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myRank
         << ": Distributor::doPosts(3 args, Teuchos::ArrayRCP)" << endl;
      *out_ << os.str ();
    }
    Teuchos::OSTab tab1 (out_);

    TEUCHOS_TEST_FOR_EXCEPTION(
      sendType == Details::DISTRIBUTOR_RSEND && ! doBarrier, std::logic_error,
      "Tpetra::Distributor::doPosts(3 args, Teuchos::ArrayRCP): Ready-send "
      "version requires a barrier between posting receives and posting ready "
      "sends.  This should have been checked before.  "
      "Please report this bug to the Tpetra developers.");

    size_t selfReceiveOffset = 0;

    // mfh 30 Mar 2016: See Github Issue #227 to see why we need to
    // check whether we're doing reverse mode before checking the
    // length of the imports array.
    if (howInitialized_ != Details::DISTRIBUTOR_INITIALIZED_BY_REVERSE) {
      // Each message has the same number of packets.
      //
      // FIXME (mfh 18 Jul 2014): Relaxing this test from strict
      // inequality to a less-than seems to have fixed Bug 6170.  It's
      // OK for the 'imports' array to be longer than it needs to be;
      // I'm just curious why it would be.
      const size_t totalNumImportPackets = totalReceiveLength_ * numPackets;
      TEUCHOS_TEST_FOR_EXCEPTION
        (static_cast<size_t> (imports.size ()) < totalNumImportPackets,
         std::invalid_argument,
         "Tpetra::Distributor::doPosts(3 args, Teuchos::ArrayRCP): "
         "The 'imports' array must have enough entries to hold the expected number "
         "of import packets.  imports.size() = " << imports.size () << " < "
         "totalNumImportPackets = " << totalNumImportPackets << ".");
    }

    // MPI tag for nonblocking receives and blocking sends in this
    // method.  Some processes might take the "fast" path
    // (indicesTo_.empty()) and others might take the "slow" path for
    // the same doPosts() call, so the path tag must be the same for
    // both.
    const int pathTag = 0;
    const int tag = this->getTag (pathTag);

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION
      (requests_.size () != 0,
       std::logic_error,
       "Tpetra::Distributor::doPosts(3 args, Teuchos::ArrayRCP): Process "
       << myRank << ": requests_.size() = " << requests_.size () << " != 0.");
#endif // HAVE_TPETRA_DEBUG

    // Distributor uses requests_.size() as the number of outstanding
    // nonblocking message requests, so we resize to zero to maintain
    // this invariant.
    //
    // numReceives_ does _not_ include the self message, if there is
    // one.  Here, we do actually send a message to ourselves, so we
    // include any self message in the "actual" number of receives to
    // post.
    //
    // NOTE (mfh 19 Mar 2012): Epetra_MpiDistributor::DoPosts()
    // doesn't (re)allocate its array of requests.  That happens in
    // CreateFromSends(), ComputeRecvs_(), DoReversePosts() (on
    // demand), or Resize_().
    const size_type actualNumReceives = as<size_type> (numReceives_) +
      as<size_type> (selfMessage_ ? 1 : 0);
    requests_.resize (0);

    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myRank << ": doPosts(3 args, Teuchos::ArrayRCP, "
         << (indicesTo_.empty () ? "fast" : "slow") << "): Post receives"
         << endl;
      *out_ << os.str ();
    }

    // Post the nonblocking receives.  It's common MPI wisdom to post
    // receives before sends.  In MPI terms, this means favoring
    // adding to the "posted queue" (of receive requests) over adding
    // to the "unexpected queue" (of arrived messages not yet matched
    // with a receive).
    {
#ifdef TPETRA_DISTRIBUTOR_TIMERS
      Teuchos::TimeMonitor timeMonRecvs (*timer_doPosts3_recvs_);
#endif // TPETRA_DISTRIBUTOR_TIMERS

      size_t curBufOffset = 0;
      for (size_type i = 0; i < actualNumReceives; ++i) {
        const size_t curBufLen = lengthsFrom_[i] * numPackets;
        if (procsFrom_[i] != myRank) {
          if (verbose) {
            std::ostringstream os;
            os << "Proc " << myRank << ": doPosts(3 args, Teuchos::ArrayRCP, "
               << (indicesTo_.empty () ? "fast" : "slow") << "): "
               << "Post irecv: {source: " << procsFrom_[i]
               << ", tag: " << tag << "}" << endl;
            *out_ << os.str ();
          }
          // If my process is receiving these packet(s) from another
          // process (not a self-receive):
          //
          // 1. Set up the persisting view (recvBuf) of the imports
          //    array, given the offset and size (total number of
          //    packets from process procsFrom_[i]).
          // 2. Start the Irecv and save the resulting request.
          TEUCHOS_TEST_FOR_EXCEPTION(
            curBufOffset + curBufLen > static_cast<size_t> (imports.size ()),
            std::logic_error,
            "Tpetra::Distributor::doPosts(3 args, Teuchos::ArrayRCP): "
            "Exceeded size of 'imports' array in packing loop on Process " <<
            myRank << ".  imports.size() = " << imports.size () << " < "
            "curBufOffset(" << curBufOffset << ") + curBufLen(" << curBufLen
            << ").");
          ArrayRCP<Packet> recvBuf =
            imports.persistingView (curBufOffset, curBufLen);
          requests_.push_back (ireceive<int, Packet> (recvBuf, procsFrom_[i],
                                                      tag, *comm_));
        }
        else { // Receiving from myself
          selfReceiveOffset = curBufOffset; // Remember the self-recv offset
        }
        curBufOffset += curBufLen;
      }
    }

    if (doBarrier) {
#ifdef TPETRA_DISTRIBUTOR_TIMERS
      Teuchos::TimeMonitor timeMonBarrier (*timer_doPosts3_barrier_);
#endif // TPETRA_DISTRIBUTOR_TIMERS

      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myRank << ": doPosts(3 args, Teuchos::ArrayRCP, "
           << (indicesTo_.empty () ? "fast" : "slow") << "): Barrier" << endl;
        *out_ << os.str ();
      }
      // If we are using ready sends (MPI_Rsend) below, we need to do
      // a barrier before we post the ready sends.  This is because a
      // ready send requires that its matching receive has already
      // been posted before the send has been posted.  The only way to
      // guarantee that in this case is to use a barrier.
      comm_->barrier ();
    }

#ifdef TPETRA_DISTRIBUTOR_TIMERS
    Teuchos::TimeMonitor timeMonSends (*timer_doPosts3_sends_);
#endif // TPETRA_DISTRIBUTOR_TIMERS

    // setup scan through procsTo_ list starting with higher numbered procs
    // (should help balance message traffic)
    //
    // FIXME (mfh 20 Feb 2013) Why haven't we precomputed this?
    // It doesn't depend on the input at all.
    size_t numBlocks = numSends_ + selfMessage_;
    size_t procIndex = 0;
    while ((procIndex < numBlocks) && (procsTo_[procIndex] < myRank)) {
      ++procIndex;
    }
    if (procIndex == numBlocks) {
      procIndex = 0;
    }

    size_t selfNum = 0;
    size_t selfIndex = 0;

    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myRank
         << ": doPosts(3 args, Teuchos::ArrayRCP, "
         << (indicesTo_.empty () ? "fast" : "slow") << "): Post sends" << endl;
      *out_ << os.str ();
    }

    if (indicesTo_.empty()) {
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myRank
           << ": doPosts(3 args, Teuchos::ArrayRCP, fast): posting sends" << endl;
        *out_ << os.str ();
      }

      // Data are already blocked (laid out) by process, so we don't
      // need a separate send buffer (besides the exports array).
      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + procIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (procsTo_[p] != myRank) {
          if (verbose) {
            std::ostringstream os;
            os << "Proc " << myRank
               << ": doPosts(3 args, Teuchos::ArrayRCP, fast): Post send: "
              "{target: " << procsTo_[p] << ", tag: " << tag << "}" << endl;
            *out_ << os.str ();
          }

          ArrayView<const Packet> tmpSend =
            exports.view (startsTo_[p]*numPackets, lengthsTo_[p]*numPackets);

          if (sendType == Details::DISTRIBUTOR_SEND) {
            send<int, Packet> (tmpSend.getRawPtr (),
                               as<int> (tmpSend.size ()),
                               procsTo_[p], tag, *comm_);
          }
          else if (sendType == Details::DISTRIBUTOR_ISEND) {
            ArrayRCP<const Packet> tmpSendBuf =
              exports.persistingView (startsTo_[p] * numPackets,
                                      lengthsTo_[p] * numPackets);
            requests_.push_back (isend<int, Packet> (tmpSendBuf, procsTo_[p],
                                                     tag, *comm_));
          }
          else if (sendType == Details::DISTRIBUTOR_RSEND) {
            readySend<int, Packet> (tmpSend.getRawPtr (),
                                    as<int> (tmpSend.size ()),
                                    procsTo_[p], tag, *comm_);
          }
          else if (sendType == Details::DISTRIBUTOR_SSEND) {
            ssend<int, Packet> (tmpSend.getRawPtr (),
                                as<int> (tmpSend.size ()),
                                procsTo_[p], tag, *comm_);
          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::logic_error,
              "Tpetra::Distributor::doPosts(3 args, Teuchos::ArrayRCP): "
              "Invalid send type.  We should never get here.  "
              "Please report this bug to the Tpetra developers.");
          }
        }
        else { // "Sending" the message to myself
          selfNum = p;
        }
      }

      if (selfMessage_) {
        if (verbose) {
          std::ostringstream os;
          os << "Proc " << myRank
             << ": doPosts(3 args, Teuchos::ArrayRCP, fast): Self-send" << endl;
          *out_ << os.str ();
        }
        // This is how we "send a message to ourself": we copy from
        // the export buffer to the import buffer.  That saves
        // Teuchos::Comm implementations other than MpiComm (in
        // particular, SerialComm) the trouble of implementing self
        // messages correctly.  (To do this right, SerialComm would
        // need internal buffer space for messages, keyed on the
        // message's tag.)
        std::copy (exports.begin()+startsTo_[selfNum]*numPackets,
                   exports.begin()+startsTo_[selfNum]*numPackets+lengthsTo_[selfNum]*numPackets,
                   imports.begin()+selfReceiveOffset);
      }
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myRank
           << ": doPosts(3 args, Teuchos::ArrayRCP, fast) done" << endl;
        *out_ << os.str ();
      }
    }
    else { // data are not blocked by proc, use send buffer
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myRank
           << ": doPosts(3 args, Teuchos::ArrayRCP, slow): posting sends" << endl;
        *out_ << os.str ();
      }

      // FIXME (mfh 05 Mar 2013) This is broken for Isend (nonblocking
      // sends), because the buffer is only long enough for one send.
      ArrayRCP<Packet> sendArray (maxSendLength_ * numPackets); // send buffer

      TEUCHOS_TEST_FOR_EXCEPTION(
        sendType == Details::DISTRIBUTOR_ISEND, std::logic_error,
        "Tpetra::Distributor::doPosts(3 args, Teuchos::ArrayRCP): "
        "The \"send buffer\" code path doesn't currently work with nonblocking sends.");

      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + procIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (procsTo_[p] != myRank) {
          if (verbose) {
            std::ostringstream os;
            os << "Proc " << myRank
               << ": doPosts(3 args, Teuchos::ArrayRCP, slow): Post send: "
              "{target: " << procsTo_[p] << ", tag: " << tag << "}" << endl;
            *out_ << os.str ();
          }

          typename ArrayView<const Packet>::iterator srcBegin, srcEnd;
          size_t sendArrayOffset = 0;
          size_t j = startsTo_[p];
          for (size_t k = 0; k < lengthsTo_[p]; ++k, ++j) {
            srcBegin = exports.begin() + indicesTo_[j]*numPackets;
            srcEnd   = srcBegin + numPackets;
            std::copy (srcBegin, srcEnd, sendArray.begin()+sendArrayOffset);
            sendArrayOffset += numPackets;
          }
          ArrayView<const Packet> tmpSend =
            sendArray.view (0, lengthsTo_[p]*numPackets);

          if (sendType == Details::DISTRIBUTOR_SEND) {
            send<int, Packet> (tmpSend.getRawPtr (),
                               as<int> (tmpSend.size ()),
                               procsTo_[p], tag, *comm_);
          }
          else if (sendType == Details::DISTRIBUTOR_ISEND) {
            ArrayRCP<const Packet> tmpSendBuf =
              sendArray.persistingView (0, lengthsTo_[p] * numPackets);
            requests_.push_back (isend<int, Packet> (tmpSendBuf, procsTo_[p],
                                                     tag, *comm_));
          }
          else if (sendType == Details::DISTRIBUTOR_RSEND) {
            readySend<int, Packet> (tmpSend.getRawPtr (),
                                    as<int> (tmpSend.size ()),
                                    procsTo_[p], tag, *comm_);
          }
          else if (sendType == Details::DISTRIBUTOR_SSEND) {
            ssend<int, Packet> (tmpSend.getRawPtr (),
                                as<int> (tmpSend.size ()),
                                procsTo_[p], tag, *comm_);
          }
          else {
            TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::logic_error,
              "Tpetra::Distributor::doPosts(3 args, Teuchos::ArrayRCP): "
              "Invalid send type.  We should never get here.  "
              "Please report this bug to the Tpetra developers.");
          }
        }
        else { // "Sending" the message to myself
          selfNum = p;
          selfIndex = startsTo_[p];
        }
      }

      if (selfMessage_) {
        if (verbose) {
          std::ostringstream os;
          os << "Proc " << myRank
             << ": doPosts(3 args, Teuchos::ArrayRCP, slow): Self-send" << endl;
          *out_ << os.str ();
        }
        for (size_t k = 0; k < lengthsTo_[selfNum]; ++k) {
          std::copy (exports.begin()+indicesTo_[selfIndex]*numPackets,
                     exports.begin()+indicesTo_[selfIndex]*numPackets + numPackets,
                     imports.begin() + selfReceiveOffset);
          ++selfIndex;
          selfReceiveOffset += numPackets;
        }
      }
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myRank
           << ": doPosts(3 args, Teuchos::ArrayRCP, slow) done" << endl;
        *out_ << os.str ();
      }
    }

    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myRank << ": doPosts done" << endl;
      *out_ << os.str ();
    }
  }

  template <class Packet>
  void Distributor::
  doPosts (const Teuchos::ArrayRCP<const Packet>& exports,
           const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
           const Teuchos::ArrayRCP<Packet>& imports,
           const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID)
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::ireceive;
    using Teuchos::isend;
    using Teuchos::readySend;
    using Teuchos::send;
    using Teuchos::ssend;
    using Teuchos::TypeNameTraits;
#ifdef HAVE_TEUCHOS_DEBUG
    using Teuchos::OSTab;
#endif // HAVE_TEUCHOS_DEBUG
    using std::endl;
    typedef Array<size_t>::size_type size_type;

    Teuchos::OSTab tab (out_);
    const bool verbose = Tpetra::Details::Behavior::verbose("Distributor");

#ifdef TPETRA_DISTRIBUTOR_TIMERS
    Teuchos::TimeMonitor timeMon (*timer_doPosts4_);
#endif // TPETRA_DISTRIBUTOR_TIMERS

    // Run-time configurable parameters that come from the input
    // ParameterList set by setParameterList().
    const Details::EDistributorSendType sendType = sendType_;
    const bool doBarrier = barrierBetween_;

// #ifdef HAVE_TEUCHOS_DEBUG
//     // Prepare for verbose output, if applicable.
//     Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel ();
//     Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream ();
//     const bool doPrint = out.get () && (comm_->getRank () == 0) &&
//       includesVerbLevel (verbLevel, Teuchos::VERB_EXTREME, true);

//     if (doPrint) {
//       // Only need one process to print out parameters.
//       *out << "Distributor::doPosts (4 args)" << endl;
//     }
//     // Add one tab level.  We declare this outside the doPrint scopes
//     // so that the tab persists until the end of this method.
//     Teuchos::OSTab tab = this->getOSTab ();
//     if (doPrint) {
//       *out << "Parameters:" << endl;
//       {
//         OSTab tab2 (out);
//         *out << "sendType: " << DistributorSendTypeEnumToString (sendType)
//              << endl << "barrierBetween: " << doBarrier << endl;
//       }
//     }
// #endif // HAVE_TEUCHOS_DEBUG

    TEUCHOS_TEST_FOR_EXCEPTION(
      sendType == Details::DISTRIBUTOR_RSEND && ! doBarrier,
      std::logic_error,
      "Tpetra::Distributor::doPosts(4 args, Teuchos::ArrayRCP): Ready-send "
      "version requires a barrier between posting receives and posting ready "
      "ends.  This should have been checked before.  "
      "Please report this bug to the Tpetra developers.");

    const int myProcID = comm_->getRank ();
    size_t selfReceiveOffset = 0;

#ifdef HAVE_TEUCHOS_DEBUG
    // Different messages may have different numbers of packets.
    size_t totalNumImportPackets = 0;
    for (size_t ii = 0; ii < static_cast<size_t> (numImportPacketsPerLID.size ()); ++ii) {
      totalNumImportPackets += numImportPacketsPerLID[ii];
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (imports.size ()) < totalNumImportPackets,
      std::runtime_error,
      "Tpetra::Distributor::doPosts(4 args, Teuchos::ArrayRCP): The 'imports' "
      "array must have enough entries to hold the expected number of import "
      "packets.  imports.size() = " << imports.size() << " < "
      "totalNumImportPackets = " << totalNumImportPackets << ".");
#endif // HAVE_TEUCHOS_DEBUG

    // MPI tag for nonblocking receives and blocking sends in this
    // method.  Some processes might take the "fast" path
    // (indicesTo_.empty()) and others might take the "slow" path for
    // the same doPosts() call, so the path tag must be the same for
    // both.
    const int pathTag = 1;
    const int tag = this->getTag (pathTag);

    if (debug_) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        requests_.size () != 0,
        std::logic_error,
        "Tpetra::Distributor::doPosts(4 args, Teuchos::ArrayRCP): Process "
        << myProcID << ": requests_.size() = " << requests_.size () << " != 0.");
    }
    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myProcID << ": doPosts(4 args, Teuchos::ArrayRCP, "
         << (indicesTo_.empty () ? "fast" : "slow") << ")" << endl;
      *out_ << os.str ();
    }

    // Distributor uses requests_.size() as the number of outstanding
    // nonblocking message requests, so we resize to zero to maintain
    // this invariant.
    //
    // numReceives_ does _not_ include the self message, if there is
    // one.  Here, we do actually send a message to ourselves, so we
    // include any self message in the "actual" number of receives to
    // post.
    //
    // NOTE (mfh 19 Mar 2012): Epetra_MpiDistributor::DoPosts()
    // doesn't (re)allocate its array of requests.  That happens in
    // CreateFromSends(), ComputeRecvs_(), DoReversePosts() (on
    // demand), or Resize_().
    const size_type actualNumReceives = as<size_type> (numReceives_) +
      as<size_type> (selfMessage_ ? 1 : 0);
    requests_.resize (0);

    // Post the nonblocking receives.  It's common MPI wisdom to post
    // receives before sends.  In MPI terms, this means favoring
    // adding to the "posted queue" (of receive requests) over adding
    // to the "unexpected queue" (of arrived messages not yet matched
    // with a receive).
    {
#ifdef TPETRA_DISTRIBUTOR_TIMERS
      Teuchos::TimeMonitor timeMonRecvs (*timer_doPosts4_recvs_);
#endif // TPETRA_DISTRIBUTOR_TIMERS

      size_t curBufferOffset = 0;
      size_t curLIDoffset = 0;
      for (size_type i = 0; i < actualNumReceives; ++i) {
        size_t totalPacketsFrom_i = 0;
        for (size_t j = 0; j < lengthsFrom_[i]; ++j) {
          totalPacketsFrom_i += numImportPacketsPerLID[curLIDoffset+j];
        }
        curLIDoffset += lengthsFrom_[i];
        if (procsFrom_[i] != myProcID && totalPacketsFrom_i) {
          // If my process is receiving these packet(s) from another
          // process (not a self-receive), and if there is at least
          // one packet to receive:
          //
          // 1. Set up the persisting view (recvBuf) into the imports
          //    array, given the offset and size (total number of
          //    packets from process procsFrom_[i]).
          // 2. Start the Irecv and save the resulting request.
          ArrayRCP<Packet> recvBuf =
            imports.persistingView (curBufferOffset, totalPacketsFrom_i);
          requests_.push_back (ireceive<int, Packet> (recvBuf, procsFrom_[i],
                                                      tag, *comm_));
        }
        else { // Receiving these packet(s) from myself
          selfReceiveOffset = curBufferOffset; // Remember the offset
        }
        curBufferOffset += totalPacketsFrom_i;
      }
    }

    if (doBarrier) {
#ifdef TPETRA_DISTRIBUTOR_TIMERS
      Teuchos::TimeMonitor timeMonBarrier (*timer_doPosts4_barrier_);
#endif // TPETRA_DISTRIBUTOR_TIMERS
      // If we are using ready sends (MPI_Rsend) below, we need to do
      // a barrier before we post the ready sends.  This is because a
      // ready send requires that its matching receive has already
      // been posted before the send has been posted.  The only way to
      // guarantee that in this case is to use a barrier.
      comm_->barrier ();
    }

#ifdef TPETRA_DISTRIBUTOR_TIMERS
    Teuchos::TimeMonitor timeMonSends (*timer_doPosts4_sends_);
#endif // TPETRA_DISTRIBUTOR_TIMERS

    // setup arrays containing starting-offsets into exports for each send,
    // and num-packets-to-send for each send.
    Array<size_t> sendPacketOffsets(numSends_,0), packetsPerSend(numSends_,0);
    size_t maxNumPackets = 0;
    size_t curPKToffset = 0;
    for (size_t pp=0; pp<numSends_; ++pp) {
      sendPacketOffsets[pp] = curPKToffset;
      size_t numPackets = 0;
      for (size_t j=startsTo_[pp]; j<startsTo_[pp]+lengthsTo_[pp]; ++j) {
        numPackets += numExportPacketsPerLID[j];
      }
      if (numPackets > maxNumPackets) maxNumPackets = numPackets;
      packetsPerSend[pp] = numPackets;
      curPKToffset += numPackets;
    }

    // setup scan through procsTo_ list starting with higher numbered procs
    // (should help balance message traffic)
    size_t numBlocks = numSends_+ selfMessage_;
    size_t procIndex = 0;
    while ((procIndex < numBlocks) && (procsTo_[procIndex] < myProcID)) {
      ++procIndex;
    }
    if (procIndex == numBlocks) {
      procIndex = 0;
    }

    size_t selfNum = 0;
    size_t selfIndex = 0;

    if (indicesTo_.empty()) {
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myProcID
           << ": doPosts(4 args, Teuchos::ArrayRCP, fast): posting sends" << endl;
        *out_ << os.str ();
      }

      // Data are already blocked (laid out) by process, so we don't
      // need a separate send buffer (besides the exports array).
      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + procIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (procsTo_[p] != myProcID && packetsPerSend[p] > 0) {
          ArrayView<const Packet> tmpSend =
            exports.view (sendPacketOffsets[p], packetsPerSend[p]);

          if (sendType == Details::DISTRIBUTOR_SEND) { // the default, so put it first
            send<int, Packet> (tmpSend.getRawPtr (),
                               as<int> (tmpSend.size ()),
                               procsTo_[p], tag, *comm_);
          }
          else if (sendType == Details::DISTRIBUTOR_RSEND) {
            readySend<int, Packet> (tmpSend.getRawPtr (),
                                    as<int> (tmpSend.size ()),
                                    procsTo_[p], tag, *comm_);
          }
          else if (sendType == Details::DISTRIBUTOR_ISEND) {
            ArrayRCP<const Packet> tmpSendBuf =
              exports.persistingView (sendPacketOffsets[p], packetsPerSend[p]);
            requests_.push_back (isend<int, Packet> (tmpSendBuf, procsTo_[p],
                                                     tag, *comm_));
          }
          else if (sendType == Details::DISTRIBUTOR_SSEND) {
            ssend<int, Packet> (tmpSend.getRawPtr (),
                                as<int> (tmpSend.size ()),
                                procsTo_[p], tag, *comm_);
          }
          else {
            TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::logic_error,
              "Tpetra::Distributor::doPosts(4 args, Teuchos::ArrayRCP): "
              "Invalid send type.  We should never get here.  Please report "
              "this bug to the Tpetra developers.");
          }
        }
        else { // "Sending" the message to myself
          selfNum = p;
        }
      }

      if (selfMessage_) {
        std::copy (exports.begin()+sendPacketOffsets[selfNum],
                   exports.begin()+sendPacketOffsets[selfNum]+packetsPerSend[selfNum],
                   imports.begin()+selfReceiveOffset);
      }
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myProcID
           << ": doPosts(4 args, Teuchos::ArrayRCP, fast) done" << endl;
        *out_ << os.str ();
      }
    }
    else { // data are not blocked by proc, use send buffer
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myProcID
           << ": doPosts(4 args, Teuchos::ArrayRCP, slow): posting sends" << endl;
        *out_ << os.str ();
      }

      // FIXME (mfh 05 Mar 2013) This may be broken for Isend.
      ArrayRCP<Packet> sendArray (maxNumPackets); // send buffer

      TEUCHOS_TEST_FOR_EXCEPTION(
        sendType == Details::DISTRIBUTOR_ISEND,
        std::logic_error,
        "Tpetra::Distributor::doPosts(4 args, Teuchos::ArrayRCP): "
        "The \"send buffer\" code path may not necessarily work with nonblocking sends.");

      Array<size_t> indicesOffsets (numExportPacketsPerLID.size(), 0);
      size_t ioffset = 0;
      for (int j=0; j<numExportPacketsPerLID.size(); ++j) {
        indicesOffsets[j] = ioffset;
        ioffset += numExportPacketsPerLID[j];
      }

      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + procIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (procsTo_[p] != myProcID) {
          typename ArrayView<const Packet>::iterator srcBegin, srcEnd;
          size_t sendArrayOffset = 0;
          size_t j = startsTo_[p];
          size_t numPacketsTo_p = 0;
          for (size_t k = 0; k < lengthsTo_[p]; ++k, ++j) {
            srcBegin = exports.begin() + indicesOffsets[j];
            srcEnd   = srcBegin + numExportPacketsPerLID[j];
            numPacketsTo_p += numExportPacketsPerLID[j];
            std::copy (srcBegin, srcEnd, sendArray.begin()+sendArrayOffset);
            sendArrayOffset += numExportPacketsPerLID[j];
          }
          if (numPacketsTo_p > 0) {
            ArrayView<const Packet> tmpSend =
              sendArray.view (0, numPacketsTo_p);

            if (sendType == Details::DISTRIBUTOR_RSEND) {
              readySend<int, Packet> (tmpSend.getRawPtr (),
                                      as<int> (tmpSend.size ()),
                                      procsTo_[p], tag, *comm_);
            }
            else if (sendType == Details::DISTRIBUTOR_ISEND) {
              ArrayRCP<const Packet> tmpSendBuf =
                sendArray.persistingView (0, numPacketsTo_p);
              requests_.push_back (isend<int, Packet> (tmpSendBuf, procsTo_[p],
                                                       tag, *comm_));
            }
            else if (sendType == Details::DISTRIBUTOR_SSEND) {
              ssend<int, Packet> (tmpSend.getRawPtr (),
                                  as<int> (tmpSend.size ()),
                                  procsTo_[p], tag, *comm_);
            }
            else { // if (sendType == Details::DISTRIBUTOR_SSEND)
              send<int, Packet> (tmpSend.getRawPtr (),
                                 as<int> (tmpSend.size ()),
                                 procsTo_[p], tag, *comm_);
            }
          }
        }
        else { // "Sending" the message to myself
          selfNum = p;
          selfIndex = startsTo_[p];
        }
      }

      if (selfMessage_) {
        for (size_t k = 0; k < lengthsTo_[selfNum]; ++k) {
          std::copy (exports.begin()+indicesOffsets[selfIndex],
                     exports.begin()+indicesOffsets[selfIndex]+numExportPacketsPerLID[selfIndex],
                     imports.begin() + selfReceiveOffset);
          selfReceiveOffset += numExportPacketsPerLID[selfIndex];
          ++selfIndex;
        }
      }
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myProcID
           << ": doPosts(4 args, Teuchos::ArrayRCP, slow) done" << endl;
        *out_ << os.str ();
      }
    }
  }

  template <class Packet>
  void Distributor::
  doReversePostsAndWaits (const Teuchos::ArrayView<const Packet>& exports,
                          size_t numPackets,
                          const Teuchos::ArrayView<Packet>& imports)
  {
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

    // doReversePosts() takes exports and imports as ArrayRCPs,
    // requiring that the memory locations are persisting.  However,
    // they need only persist within the scope of that routine, so it
    // is safe for us to use nonpersisting references in this case.

    // mfh 04 Apr 2012: For some reason, calling arcp<const Packet>
    // for Packet=std::complex<T> (e.g., T=float) fails to compile
    // with some versions of GCC.  The issue only arises with the
    // exports array.  This is why we construct a separate nonowning
    // ArrayRCP.
    typedef typename ArrayRCP<const Packet>::size_type size_type;
    ArrayRCP<const Packet> exportsArcp (exports.getRawPtr(), as<size_type> (0),
                                        exports.size(), false);
    // mfh 04 Apr 2012: This is the offending code.  This statement
    // would normally be in place of "exportsArcp" in the
    // doReversePosts() call below.
    //arcp<const Packet> (exports.getRawPtr(), 0, exports.size(), false)
    doReversePosts (exportsArcp,
                    numPackets,
                    arcp<Packet> (imports.getRawPtr (), 0, imports.size (), false));
    doReverseWaits ();

    lastRoundBytesSend_ = exports.size() * sizeof(Packet);
    lastRoundBytesRecv_ = imports.size() * sizeof(Packet);
  }

  template <class Packet>
  void Distributor::
  doReversePostsAndWaits (const Teuchos::ArrayView<const Packet>& exports,
                          const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                          const Teuchos::ArrayView<Packet> &imports,
                          const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID)
  {
    using Teuchos::as;
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;

    TEUCHOS_TEST_FOR_EXCEPTION(
      requests_.size () != 0, std::runtime_error, "Tpetra::Distributor::"
      "doReversePostsAndWaits(4 args): There are " << requests_.size ()
      << " outstanding nonblocking messages pending.  It is incorrect to call "
      "this method with posts outstanding.");

    // doReversePosts() accepts the exports and imports arrays as
    // ArrayRCPs, requiring that the memory location is persisting (as
    // is necessary for nonblocking receives).  However, it need only
    // persist until doReverseWaits() completes, so it is safe for us
    // to use a nonpersisting reference in this case.  The use of a
    // nonpersisting reference is purely a performance optimization.

    // mfh 02 Apr 2012: For some reason, calling arcp<const Packet>
    // for Packet=std::complex<double> fails to compile with some
    // versions of GCC.  The issue only arises with the exports array.
    // This is why we construct a separate nonowning ArrayRCP.
    typedef typename ArrayRCP<const Packet>::size_type size_type;
    ArrayRCP<const Packet> exportsArcp (exports.getRawPtr (), as<size_type> (0),
                                        exports.size (), false);
    doReversePosts (exportsArcp,
                    numExportPacketsPerLID,
                    arcp<Packet> (imports.getRawPtr (), 0, imports.size (), false),
                    numImportPacketsPerLID);
    doReverseWaits ();

    lastRoundBytesSend_ = exports.size() * sizeof(Packet);
    lastRoundBytesRecv_ = imports.size() * sizeof(Packet);
  }

  template <class Packet>
  void Distributor::
  doReversePosts (const Teuchos::ArrayRCP<const Packet>& exports,
                  size_t numPackets,
                  const Teuchos::ArrayRCP<Packet>& imports)
  {
    // FIXME (mfh 29 Mar 2012) WHY?
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! indicesTo_.empty (), std::runtime_error,
      "Tpetra::Distributor::doReversePosts(3 args): Can only do reverse "
      "communication when original data are blocked by process.");
    if (reverseDistributor_.is_null ()) {
      createReverseDistributor ();
    }
    reverseDistributor_->doPosts (exports, numPackets, imports);
  }

  template <class Packet>
  void Distributor::
  doReversePosts (const Teuchos::ArrayRCP<const Packet>& exports,
                  const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                  const Teuchos::ArrayRCP<Packet>& imports,
                  const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID)
  {
    // FIXME (mfh 29 Mar 2012) WHY?
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! indicesTo_.empty (), std::runtime_error,
      "Tpetra::Distributor::doReversePosts(3 args): Can only do reverse "
      "communication when original data are blocked by process.");
    if (reverseDistributor_.is_null ()) {
      createReverseDistributor ();
    }
    reverseDistributor_->doPosts (exports, numExportPacketsPerLID,
                                  imports, numImportPacketsPerLID);
  }

  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
  Distributor::
  doPostsAndWaits (const ExpView& exports,
                   size_t numPackets,
                   const ImpView& imports)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    const bool verbose = Tpetra::Details::Behavior::verbose("Distributor");

    RCP<Teuchos::OSTab> tab0, tab1;
    if (verbose) {
      tab0 = rcp (new Teuchos::OSTab (out_));
      const int myRank = comm_->getRank ();
      std::ostringstream os;
      os << "Proc " << myRank
         << ": Distributor::doPostsAndWaits(3 args, Kokkos): "
         << "{sendType: " << DistributorSendTypeEnumToString (sendType_)
         << ", barrierBetween: " << barrierBetween_ << "}" << endl;
      *out_ << os.str ();
      tab1 = rcp (new Teuchos::OSTab (out_));
    }

    TEUCHOS_TEST_FOR_EXCEPTION(
      requests_.size () != 0, std::runtime_error, "Tpetra::Distributor::"
      "doPostsAndWaits(3 args): There are " << requests_.size () <<
      " outstanding nonblocking messages pending.  It is incorrect to call "
      "this method with posts outstanding.");

    if (verbose) {
      const int myRank = comm_->getRank ();
      std::ostringstream os;
      os << "Proc " << myRank
         << ": Distributor::doPostsAndWaits: Call doPosts" << endl;
      *out_ << os.str ();
    }
    doPosts (exports, numPackets, imports);
    if (verbose) {
      const int myRank = comm_->getRank ();
      std::ostringstream os;
      os << "Proc " << myRank
         << ": Distributor::doPostsAndWaits: Call doWaits" << endl;
      *out_ << os.str ();
    }
    doWaits ();
  }

  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
  Distributor::
  doPostsAndWaits (const ExpView& exports,
                   const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                   const ImpView& imports,
                   const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      requests_.size () != 0, std::runtime_error,
      "Tpetra::Distributor::doPostsAndWaits(4 args): There are "
      << requests_.size () << " outstanding nonblocking messages pending.  "
      "It is incorrect to call this method with posts outstanding.");

    doPosts (exports, numExportPacketsPerLID, imports, numImportPacketsPerLID);
    doWaits ();
  }


  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
  Distributor::
  doPosts (const ExpView &exports,
           size_t numPackets,
           const ImpView &imports)
  {
    using Teuchos::Array;
    using Teuchos::as;
    using Teuchos::FancyOStream;
    using Teuchos::includesVerbLevel;
    using Teuchos::ireceive;
    using Teuchos::isend;
    using Teuchos::OSTab;
    using Teuchos::readySend;
    using Teuchos::send;
    using Teuchos::ssend;
    using Teuchos::TypeNameTraits;
    using Teuchos::typeName;
    using std::endl;
    using Kokkos::Compat::create_const_view;
    using Kokkos::Compat::create_view;
    using Kokkos::Compat::subview_offset;
    using Kokkos::Compat::deep_copy_offset;
    typedef Array<size_t>::size_type size_type;
    typedef ExpView exports_view_type;
    typedef ImpView imports_view_type;

    const bool verbose = Tpetra::Details::Behavior::verbose("Distributor");
#ifdef KOKKOS_ENABLE_CUDA
    static_assert (! std::is_same<typename ExpView::memory_space, Kokkos::CudaUVMSpace>::value &&
                   ! std::is_same<typename ImpView::memory_space, Kokkos::CudaUVMSpace>::value,
                   "Please do not use Tpetra::Distributor with UVM "
                   "allocations.  See GitHub issue #1088.");
#endif // KOKKOS_ENABLE_CUDA

#ifdef TPETRA_DISTRIBUTOR_TIMERS
    Teuchos::TimeMonitor timeMon (*timer_doPosts3_);
#endif // TPETRA_DISTRIBUTOR_TIMERS

    const int myRank = comm_->getRank ();
    // Run-time configurable parameters that come from the input
    // ParameterList set by setParameterList().
    const Details::EDistributorSendType sendType = sendType_;
    const bool doBarrier = barrierBetween_;

    Teuchos::OSTab tab0 (out_);
    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myRank
         << ": Distributor::doPosts(3 args, Kokkos)" << endl;
      *out_ << os.str ();
    }
    Teuchos::OSTab tab1 (out_);

    TEUCHOS_TEST_FOR_EXCEPTION(
      sendType == Details::DISTRIBUTOR_RSEND && ! doBarrier,
      std::logic_error,
      "Tpetra::Distributor::doPosts(3 args, Kokkos): Ready-send version "
      "requires a barrier between posting receives and posting ready sends.  "
      "This should have been checked before.  "
      "Please report this bug to the Tpetra developers.");

    size_t selfReceiveOffset = 0;

    // mfh 30 Mar 2016: See Github Issue #227 to see why we need to
    // check whether we're doing reverse mode before checking the
    // length of the imports array.
    if (false /* howInitialized_ != Details::DISTRIBUTOR_INITIALIZED_BY_REVERSE */) {
      // Each message has the same number of packets.
      const size_t totalNumImportPackets = totalReceiveLength_ * numPackets;

      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myRank << ": doPosts: totalNumImportPackets = " <<
          totalNumImportPackets << " = " << totalReceiveLength_ << " * " <<
          numPackets << "; imports.extent(0) = " << imports.extent (0)
           << endl;
        *out_ << os.str ();
      }

#ifdef HAVE_TPETRA_DEBUG
      // mfh 31 Mar 2016: Extra special all-reduce check to help diagnose #227.
      {
        const size_t importBufSize = static_cast<size_t> (imports.extent (0));
        const int lclBad = (importBufSize < totalNumImportPackets) ? 1 : 0;
        int gblBad = 0;
        using Teuchos::reduceAll;
        using Teuchos::REDUCE_MAX;
        using Teuchos::outArg;
        reduceAll (*comm_, REDUCE_MAX, lclBad, outArg (gblBad));
        TEUCHOS_TEST_FOR_EXCEPTION
          (gblBad != 0,
           std::runtime_error,
           "Tpetra::Distributor::doPosts(3 args, Kokkos): "
           "On one or more MPI processes, the 'imports' array "
           "does not have enough entries to hold the expected number of "
           "import packets.  ");
      }
#else
      TEUCHOS_TEST_FOR_EXCEPTION
        (static_cast<size_t> (imports.extent (0)) < totalNumImportPackets,
         std::runtime_error,
         "Tpetra::Distributor::doPosts(3 args, Kokkos): The 'imports' "
         "array must have enough entries to hold the expected number of import "
         "packets.  imports.extent(0) = " << imports.extent (0) << " < "
         "totalNumImportPackets = " << totalNumImportPackets << " = "
         "totalReceiveLength_ (" << totalReceiveLength_ << ") * numPackets ("
         << numPackets << ").");
#endif // HAVE_TPETRA_DEBUG
    }

    // MPI tag for nonblocking receives and blocking sends in this
    // method.  Some processes might take the "fast" path
    // (indicesTo_.empty()) and others might take the "slow" path for
    // the same doPosts() call, so the path tag must be the same for
    // both.
    const int pathTag = 0;
    const int tag = this->getTag (pathTag);

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION
      (requests_.size () != 0,
       std::logic_error,
       "Tpetra::Distributor::doPosts(3 args, Kokkos): Process "
       << myRank << ": requests_.size() = " << requests_.size () << " != 0.");
#endif // HAVE_TPETRA_DEBUG

    // Distributor uses requests_.size() as the number of outstanding
    // nonblocking message requests, so we resize to zero to maintain
    // this invariant.
    //
    // numReceives_ does _not_ include the self message, if there is
    // one.  Here, we do actually send a message to ourselves, so we
    // include any self message in the "actual" number of receives to
    // post.
    //
    // NOTE (mfh 19 Mar 2012): Epetra_MpiDistributor::DoPosts()
    // doesn't (re)allocate its array of requests.  That happens in
    // CreateFromSends(), ComputeRecvs_(), DoReversePosts() (on
    // demand), or Resize_().
    const size_type actualNumReceives = as<size_type> (numReceives_) +
      as<size_type> (selfMessage_ ? 1 : 0);
    requests_.resize (0);

    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myRank << ": doPosts(3 args, Kokkos, "
         << (indicesTo_.empty () ? "fast" : "slow") << "): Post receives"
         << endl;
      *out_ << os.str ();
    }

    // Post the nonblocking receives.  It's common MPI wisdom to post
    // receives before sends.  In MPI terms, this means favoring
    // adding to the "posted queue" (of receive requests) over adding
    // to the "unexpected queue" (of arrived messages not yet matched
    // with a receive).
    {
#ifdef TPETRA_DISTRIBUTOR_TIMERS
      Teuchos::TimeMonitor timeMonRecvs (*timer_doPosts3_recvs_);
#endif // TPETRA_DISTRIBUTOR_TIMERS

      size_t curBufferOffset = 0;
      for (size_type i = 0; i < actualNumReceives; ++i) {
        const size_t curBufLen = lengthsFrom_[i] * numPackets;
        if (procsFrom_[i] != myRank) {
          if (verbose) {
            std::ostringstream os;
            os << "Proc " << myRank << ": doPosts(3 args, Kokkos, "
               << (indicesTo_.empty () ? "fast" : "slow") << "): "
               << "Post irecv: {source: " << procsFrom_[i]
               << ", tag: " << tag << "}" << endl;
            *out_ << os.str ();
          }
          // If my process is receiving these packet(s) from another
          // process (not a self-receive):
          //
          // 1. Set up the persisting view (recvBuf) of the imports
          //    array, given the offset and size (total number of
          //    packets from process procsFrom_[i]).
          // 2. Start the Irecv and save the resulting request.
          TEUCHOS_TEST_FOR_EXCEPTION(
            curBufferOffset + curBufLen > static_cast<size_t> (imports.size ()),
            std::logic_error, "Tpetra::Distributor::doPosts(3 args, Kokkos): "
            "Exceeded size of 'imports' array in packing loop on Process " <<
            myRank << ".  imports.size() = " << imports.size () << " < "
            "curBufferOffset(" << curBufferOffset << ") + curBufLen(" <<
            curBufLen << ").");
          imports_view_type recvBuf =
            subview_offset (imports, curBufferOffset, curBufLen);
          requests_.push_back (ireceive<int> (recvBuf, procsFrom_[i],
                                              tag, *comm_));
        }
        else { // Receiving from myself
          selfReceiveOffset = curBufferOffset; // Remember the self-recv offset
        }
        curBufferOffset += curBufLen;
      }
    }

    if (doBarrier) {
#ifdef TPETRA_DISTRIBUTOR_TIMERS
      Teuchos::TimeMonitor timeMonBarrier (*timer_doPosts3_barrier_);
#endif // TPETRA_DISTRIBUTOR_TIMERS

      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myRank << ": doPosts(3 args, Kokkos, "
           << (indicesTo_.empty () ? "fast" : "slow") << "): Barrier" << endl;
        *out_ << os.str ();
      }
      // If we are using ready sends (MPI_Rsend) below, we need to do
      // a barrier before we post the ready sends.  This is because a
      // ready send requires that its matching receive has already
      // been posted before the send has been posted.  The only way to
      // guarantee that in this case is to use a barrier.
      comm_->barrier ();
    }

#ifdef TPETRA_DISTRIBUTOR_TIMERS
    Teuchos::TimeMonitor timeMonSends (*timer_doPosts3_sends_);
#endif // TPETRA_DISTRIBUTOR_TIMERS

    // setup scan through procsTo_ list starting with higher numbered procs
    // (should help balance message traffic)
    //
    // FIXME (mfh 20 Feb 2013) Why haven't we precomputed this?
    // It doesn't depend on the input at all.
    size_t numBlocks = numSends_ + selfMessage_;
    size_t procIndex = 0;
    while ((procIndex < numBlocks) && (procsTo_[procIndex] < myRank)) {
      ++procIndex;
    }
    if (procIndex == numBlocks) {
      procIndex = 0;
    }

    size_t selfNum = 0;
    size_t selfIndex = 0;

    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myRank << ": doPosts(3 args, Kokkos, "
         << (indicesTo_.empty () ? "fast" : "slow") << "): Post sends" << endl;
      *out_ << os.str ();
    }

    if (indicesTo_.empty()) {
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myRank
           << ": doPosts(3 args, Kokkos, fast): posting sends" << endl;
        *out_ << os.str ();
      }

      // Data are already blocked (laid out) by process, so we don't
      // need a separate send buffer (besides the exports array).
      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + procIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (procsTo_[p] != myRank) {
          if (verbose) {
            std::ostringstream os;
            os << "Proc " << myRank << ": doPosts(3 args, Kokkos, fast): Post send: "
              "{target: " << procsTo_[p] << ", tag: " << tag << "}" << endl;
            *out_ << os.str ();
          }
          // if (debug_) {
          //   const size_t off = startsTo_[p] * numPackets;
          //   const size_t len = lengthsTo_[p] * numPackets;
          //   TEUCHOS_TEST_FOR_EXCEPTION
          //     (static_cast<size_t> (off + len) >
          //      static_cast<size_t> (exports.size ()), std::logic_error,
          //      "doPosts: off=" << off << " + len=" << len << " > "
          //      "exports.size()=" << exports.size () << ".");
          // }

          exports_view_type tmpSend = subview_offset(
            exports, startsTo_[p]*numPackets, lengthsTo_[p]*numPackets);

          if (sendType == Details::DISTRIBUTOR_SEND) {
            send<int> (tmpSend,
                       as<int> (tmpSend.size ()),
                       procsTo_[p], tag, *comm_);
          }
          else if (sendType == Details::DISTRIBUTOR_ISEND) {
            exports_view_type tmpSendBuf =
              subview_offset (exports, startsTo_[p] * numPackets,
                              lengthsTo_[p] * numPackets);
            requests_.push_back (isend<int> (tmpSendBuf, procsTo_[p],
                                             tag, *comm_));
          }
          else if (sendType == Details::DISTRIBUTOR_RSEND) {
            readySend<int> (tmpSend,
                            as<int> (tmpSend.size ()),
                            procsTo_[p], tag, *comm_);
          }
          else if (sendType == Details::DISTRIBUTOR_SSEND) {
            ssend<int> (tmpSend,
                        as<int> (tmpSend.size ()),
                        procsTo_[p], tag, *comm_);
          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(
              true,
              std::logic_error,
              "Tpetra::Distributor::doPosts(3 args, Kokkos): "
              "Invalid send type.  We should never get here.  "
              "Please report this bug to the Tpetra developers.");
          }
        }
        else { // "Sending" the message to myself
          selfNum = p;
        }
      }

      if (selfMessage_) {
        if (verbose) {
          std::ostringstream os;
          os << "Proc " << myRank
             << ": doPosts(3 args, Kokkos, fast): Self-send" << endl;
          *out_ << os.str ();
        }
        // This is how we "send a message to ourself": we copy from
        // the export buffer to the import buffer.  That saves
        // Teuchos::Comm implementations other than MpiComm (in
        // particular, SerialComm) the trouble of implementing self
        // messages correctly.  (To do this right, SerialComm would
        // need internal buffer space for messages, keyed on the
        // message's tag.)
        deep_copy_offset(imports, exports, selfReceiveOffset,
                         startsTo_[selfNum]*numPackets,
                         lengthsTo_[selfNum]*numPackets);
      }
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myRank << ": doPosts(3 args, Kokkos, fast) done" << endl;
        *out_ << os.str ();
      }
    }
    else { // data are not blocked by proc, use send buffer
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myRank
           << ": doPosts(3 args, Kokkos, slow): posting sends" << endl;
        *out_ << os.str ();
      }

      typedef typename ExpView::non_const_value_type Packet;
      typedef typename ExpView::array_layout Layout;
      typedef typename ExpView::device_type Device;
      typedef typename ExpView::memory_traits Mem;
      Kokkos::View<Packet*,Layout,Device,Mem> sendArray ("sendArray",
                                               maxSendLength_ * numPackets);

      // FIXME (mfh 05 Mar 2013) This is broken for Isend (nonblocking
      // sends), because the buffer is only long enough for one send.
      TEUCHOS_TEST_FOR_EXCEPTION(
        sendType == Details::DISTRIBUTOR_ISEND,
        std::logic_error,
        "Tpetra::Distributor::doPosts(3 args, Kokkos): The \"send buffer\" code path "
        "doesn't currently work with nonblocking sends.");

      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + procIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (procsTo_[p] != myRank) {
          if (verbose) {
            std::ostringstream os;
            os << "Proc " << myRank
               << ": doPosts(3 args, Kokkos, slow): Post send: {target: "
               << procsTo_[p] << ", tag: " << tag << "}" << endl;
            *out_ << os.str ();
          }

          size_t sendArrayOffset = 0;
          size_t j = startsTo_[p];
          for (size_t k = 0; k < lengthsTo_[p]; ++k, ++j) {
            deep_copy_offset(sendArray, exports, sendArrayOffset,
                             indicesTo_[j]*numPackets, numPackets);
            sendArrayOffset += numPackets;
          }
          ImpView tmpSend =
            subview_offset(sendArray, size_t(0), lengthsTo_[p]*numPackets);

          if (sendType == Details::DISTRIBUTOR_SEND) {
            send<int> (tmpSend,
                       as<int> (tmpSend.size ()),
                       procsTo_[p], tag, *comm_);
          }
          else if (sendType == Details::DISTRIBUTOR_ISEND) {
            exports_view_type tmpSendBuf =
              subview_offset (sendArray, size_t(0), lengthsTo_[p] * numPackets);
            requests_.push_back (isend<int> (tmpSendBuf, procsTo_[p],
                                             tag, *comm_));
          }
          else if (sendType == Details::DISTRIBUTOR_RSEND) {
            readySend<int> (tmpSend,
                            as<int> (tmpSend.size ()),
                            procsTo_[p], tag, *comm_);
          }
          else if (sendType == Details::DISTRIBUTOR_SSEND) {
            ssend<int> (tmpSend,
                        as<int> (tmpSend.size ()),
                        procsTo_[p], tag, *comm_);
          }
          else {
            TEUCHOS_TEST_FOR_EXCEPTION(
              true,
              std::logic_error,
              "Tpetra::Distributor::doPosts(3 args, Kokkos): "
              "Invalid send type.  We should never get here.  "
              "Please report this bug to the Tpetra developers.");
          }
        }
        else { // "Sending" the message to myself
          selfNum = p;
          selfIndex = startsTo_[p];
        }
      }

      if (selfMessage_) {
        if (verbose) {
          std::ostringstream os;
          os << "Proc " << myRank
             << ": doPosts(3 args, Kokkos, slow): Self-send" << endl;
          *out_ << os.str ();
        }
        for (size_t k = 0; k < lengthsTo_[selfNum]; ++k) {
          deep_copy_offset(imports, exports, selfReceiveOffset,
                           indicesTo_[selfIndex]*numPackets, numPackets);
          ++selfIndex;
          selfReceiveOffset += numPackets;
        }
      }
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myRank
           << ": doPosts(3 args, Kokkos, slow) done" << endl;
        *out_ << os.str ();
      }
    }

    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myRank << ": doPosts done" << endl;
      *out_ << os.str ();
    }
  }

  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
  Distributor::
  doPosts (const ExpView &exports,
           const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
           const ImpView &imports,
           const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID)
  {
    using Teuchos::Array;
    using Teuchos::as;
    using Teuchos::ireceive;
    using Teuchos::isend;
    using Teuchos::readySend;
    using Teuchos::send;
    using Teuchos::ssend;
    using Teuchos::TypeNameTraits;
#ifdef HAVE_TEUCHOS_DEBUG
    using Teuchos::OSTab;
#endif // HAVE_TEUCHOS_DEBUG
    using std::endl;
    using Kokkos::Compat::create_const_view;
    using Kokkos::Compat::create_view;
    using Kokkos::Compat::subview_offset;
    using Kokkos::Compat::deep_copy_offset;
    typedef Array<size_t>::size_type size_type;
    typedef ExpView exports_view_type;
    typedef ImpView imports_view_type;
    const bool verbose = Tpetra::Details::Behavior::verbose("Distributor");

#ifdef KOKKOS_ENABLE_CUDA
    static_assert (! std::is_same<typename ExpView::memory_space, Kokkos::CudaUVMSpace>::value &&
                   ! std::is_same<typename ImpView::memory_space, Kokkos::CudaUVMSpace>::value,
                   "Please do not use Tpetra::Distributor with UVM "
                   "allocations.  See GitHub issue #1088.");
#endif // KOKKOS_ENABLE_CUDA

    Teuchos::OSTab tab (out_);

#ifdef TPETRA_DISTRIBUTOR_TIMERS
    Teuchos::TimeMonitor timeMon (*timer_doPosts4_);
#endif // TPETRA_DISTRIBUTOR_TIMERS

    // Run-time configurable parameters that come from the input
    // ParameterList set by setParameterList().
    const Details::EDistributorSendType sendType = sendType_;
    const bool doBarrier = barrierBetween_;

// #ifdef HAVE_TEUCHOS_DEBUG
//     // Prepare for verbose output, if applicable.
//     Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel ();
//     RCP<Teuchos::FancyOStream> out = this->getOStream ();
//     const bool doPrint = out.get () && (comm_->getRank () == 0) &&
//       includesVerbLevel (verbLevel, Teuchos::VERB_EXTREME, true);

//     if (doPrint) {
//       // Only need one process to print out parameters.
//       *out << "Distributor::doPosts (4 args)" << endl;
//     }
//     // Add one tab level.  We declare this outside the doPrint scopes
//     // so that the tab persists until the end of this method.
//     Teuchos::OSTab tab = this->getOSTab ();
//     if (doPrint) {
//       *out << "Parameters:" << endl;
//       {
//         OSTab tab2 (out);
//         *out << "sendType: " << DistributorSendTypeEnumToString (sendType)
//              << endl << "barrierBetween: " << doBarrier << endl;
//       }
//     }
// #endif // HAVE_TEUCHOS_DEBUG

    TEUCHOS_TEST_FOR_EXCEPTION(
      sendType == Details::DISTRIBUTOR_RSEND && ! doBarrier,
      std::logic_error, "Tpetra::Distributor::doPosts(4 args, Kokkos): Ready-send "
      "version requires a barrier between posting receives and posting ready "
      "sends.  This should have been checked before.  "
      "Please report this bug to the Tpetra developers.");

    const int myProcID = comm_->getRank ();
    size_t selfReceiveOffset = 0;

#ifdef HAVE_TEUCHOS_DEBUG
    // Different messages may have different numbers of packets.
    size_t totalNumImportPackets = 0;
    for (size_type ii = 0; ii < numImportPacketsPerLID.size (); ++ii) {
      totalNumImportPackets += numImportPacketsPerLID[ii];
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      imports.extent (0) < totalNumImportPackets, std::runtime_error,
      "Tpetra::Distributor::doPosts(4 args, Kokkos): The 'imports' array must have "
      "enough entries to hold the expected number of import packets.  "
      "imports.extent(0) = " << imports.extent (0) << " < "
      "totalNumImportPackets = " << totalNumImportPackets << ".");
#endif // HAVE_TEUCHOS_DEBUG

    // MPI tag for nonblocking receives and blocking sends in this
    // method.  Some processes might take the "fast" path
    // (indicesTo_.empty()) and others might take the "slow" path for
    // the same doPosts() call, so the path tag must be the same for
    // both.
    const int pathTag = 1;
    const int tag = this->getTag (pathTag);

    if (debug_) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        requests_.size () != 0, std::logic_error, "Tpetra::Distributor::"
        "doPosts(4 args, Kokkos): Process " << myProcID << ": requests_.size () = "
        << requests_.size () << " != 0.");
    }
    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myProcID << ": doPosts(4 args, Kokkos, "
         << (indicesTo_.empty () ? "fast" : "slow") << ")" << endl;
      *out_ << os.str ();
    }

    // Distributor uses requests_.size() as the number of outstanding
    // nonblocking message requests, so we resize to zero to maintain
    // this invariant.
    //
    // numReceives_ does _not_ include the self message, if there is
    // one.  Here, we do actually send a message to ourselves, so we
    // include any self message in the "actual" number of receives to
    // post.
    //
    // NOTE (mfh 19 Mar 2012): Epetra_MpiDistributor::DoPosts()
    // doesn't (re)allocate its array of requests.  That happens in
    // CreateFromSends(), ComputeRecvs_(), DoReversePosts() (on
    // demand), or Resize_().
    const size_type actualNumReceives = as<size_type> (numReceives_) +
      as<size_type> (selfMessage_ ? 1 : 0);
    requests_.resize (0);

    // Post the nonblocking receives.  It's common MPI wisdom to post
    // receives before sends.  In MPI terms, this means favoring
    // adding to the "posted queue" (of receive requests) over adding
    // to the "unexpected queue" (of arrived messages not yet matched
    // with a receive).
    {
#ifdef TPETRA_DISTRIBUTOR_TIMERS
      Teuchos::TimeMonitor timeMonRecvs (*timer_doPosts4_recvs_);
#endif // TPETRA_DISTRIBUTOR_TIMERS

      size_t curBufferOffset = 0;
      size_t curLIDoffset = 0;
      for (size_type i = 0; i < actualNumReceives; ++i) {
        size_t totalPacketsFrom_i = 0;
        for (size_t j = 0; j < lengthsFrom_[i]; ++j) {
          totalPacketsFrom_i += numImportPacketsPerLID[curLIDoffset+j];
        }
        curLIDoffset += lengthsFrom_[i];
        if (procsFrom_[i] != myProcID && totalPacketsFrom_i) {
          // If my process is receiving these packet(s) from another
          // process (not a self-receive), and if there is at least
          // one packet to receive:
          //
          // 1. Set up the persisting view (recvBuf) into the imports
          //    array, given the offset and size (total number of
          //    packets from process procsFrom_[i]).
          // 2. Start the Irecv and save the resulting request.
          imports_view_type recvBuf =
            subview_offset (imports, curBufferOffset, totalPacketsFrom_i);
          requests_.push_back (ireceive<int> (recvBuf, procsFrom_[i],
                                              tag, *comm_));
        }
        else { // Receiving these packet(s) from myself
          selfReceiveOffset = curBufferOffset; // Remember the offset
        }
        curBufferOffset += totalPacketsFrom_i;
      }
    }

    if (doBarrier) {
#ifdef TPETRA_DISTRIBUTOR_TIMERS
      Teuchos::TimeMonitor timeMonBarrier (*timer_doPosts4_barrier_);
#endif // TPETRA_DISTRIBUTOR_TIMERS
      // If we are using ready sends (MPI_Rsend) below, we need to do
      // a barrier before we post the ready sends.  This is because a
      // ready send requires that its matching receive has already
      // been posted before the send has been posted.  The only way to
      // guarantee that in this case is to use a barrier.
      comm_->barrier ();
    }

#ifdef TPETRA_DISTRIBUTOR_TIMERS
    Teuchos::TimeMonitor timeMonSends (*timer_doPosts4_sends_);
#endif // TPETRA_DISTRIBUTOR_TIMERS

    // setup arrays containing starting-offsets into exports for each send,
    // and num-packets-to-send for each send.
    Array<size_t> sendPacketOffsets(numSends_,0), packetsPerSend(numSends_,0);
    size_t maxNumPackets = 0;
    size_t curPKToffset = 0;
    for (size_t pp=0; pp<numSends_; ++pp) {
      sendPacketOffsets[pp] = curPKToffset;
      size_t numPackets = 0;
      for (size_t j=startsTo_[pp]; j<startsTo_[pp]+lengthsTo_[pp]; ++j) {
        numPackets += numExportPacketsPerLID[j];
      }
      if (numPackets > maxNumPackets) maxNumPackets = numPackets;
      packetsPerSend[pp] = numPackets;
      curPKToffset += numPackets;
    }

    // setup scan through procsTo_ list starting with higher numbered procs
    // (should help balance message traffic)
    size_t numBlocks = numSends_+ selfMessage_;
    size_t procIndex = 0;
    while ((procIndex < numBlocks) && (procsTo_[procIndex] < myProcID)) {
      ++procIndex;
    }
    if (procIndex == numBlocks) {
      procIndex = 0;
    }

    size_t selfNum = 0;
    size_t selfIndex = 0;
    if (indicesTo_.empty()) {
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myProcID
           << ": doPosts(4 args, Kokkos, fast): posting sends" << endl;
        *out_ << os.str ();
      }

      // Data are already blocked (laid out) by process, so we don't
      // need a separate send buffer (besides the exports array).
      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + procIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (procsTo_[p] != myProcID && packetsPerSend[p] > 0) {
          exports_view_type tmpSend =
            subview_offset(exports, sendPacketOffsets[p], packetsPerSend[p]);

          if (sendType == Details::DISTRIBUTOR_SEND) { // the default, so put it first
            send<int> (tmpSend,
                       as<int> (tmpSend.size ()),
                       procsTo_[p], tag, *comm_);
          }
          else if (sendType == Details::DISTRIBUTOR_RSEND) {
            readySend<int> (tmpSend,
                            as<int> (tmpSend.size ()),
                            procsTo_[p], tag, *comm_);
          }
          else if (sendType == Details::DISTRIBUTOR_ISEND) {
            exports_view_type tmpSendBuf =
              subview_offset (exports, sendPacketOffsets[p], packetsPerSend[p]);
            requests_.push_back (isend<int> (tmpSendBuf, procsTo_[p],
                                             tag, *comm_));
          }
          else if (sendType == Details::DISTRIBUTOR_SSEND) {
            ssend<int> (tmpSend,
                        as<int> (tmpSend.size ()),
                        procsTo_[p], tag, *comm_);
          }
          else {
            TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::logic_error,
              "Tpetra::Distributor::doPosts(4 args, Kokkos): "
              "Invalid send type.  We should never get here.  "
              "Please report this bug to the Tpetra developers.");
          }
        }
        else { // "Sending" the message to myself
          selfNum = p;
        }
      }

      if (selfMessage_) {
        deep_copy_offset(imports, exports, selfReceiveOffset,
                         sendPacketOffsets[selfNum], packetsPerSend[selfNum]);
      }
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myProcID << ": doPosts(4 args, Kokkos, fast) done" << endl;
        *out_ << os.str ();
      }
    }
    else { // data are not blocked by proc, use send buffer
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myProcID << ": doPosts(4 args, Kokkos, slow): posting sends" << endl;
        *out_ << os.str ();
      }

      // FIXME (mfh 05 Mar 2013) This may be broken for Isend.
      typedef typename ExpView::non_const_value_type Packet;
      typedef typename ExpView::array_layout Layout;
      typedef typename ExpView::device_type Device;
      typedef typename ExpView::memory_traits Mem;
      Kokkos::View<Packet*,Layout,Device,Mem> sendArray ("sendArray", maxNumPackets); // send buffer

      TEUCHOS_TEST_FOR_EXCEPTION(
        sendType == Details::DISTRIBUTOR_ISEND,
        std::logic_error,
        "Tpetra::Distributor::doPosts(4 args, Kokkos): "
        "The \"send buffer\" code path may not necessarily work with nonblocking sends.");

      Array<size_t> indicesOffsets (numExportPacketsPerLID.size(), 0);
      size_t ioffset = 0;
      for (int j=0; j<numExportPacketsPerLID.size(); ++j) {
        indicesOffsets[j] = ioffset;
        ioffset += numExportPacketsPerLID[j];
      }

      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + procIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (procsTo_[p] != myProcID) {
          size_t sendArrayOffset = 0;
          size_t j = startsTo_[p];
          size_t numPacketsTo_p = 0;
          for (size_t k = 0; k < lengthsTo_[p]; ++k, ++j) {
            numPacketsTo_p += numExportPacketsPerLID[j];
            deep_copy_offset(sendArray, exports, sendArrayOffset,
                             indicesOffsets[j], numExportPacketsPerLID[j]);
            sendArrayOffset += numExportPacketsPerLID[j];
          }
          if (numPacketsTo_p > 0) {
            ImpView tmpSend =
              subview_offset(sendArray, size_t(0), numPacketsTo_p);

            if (sendType == Details::DISTRIBUTOR_RSEND) {
              readySend<int> (tmpSend,
                              as<int> (tmpSend.size ()),
                              procsTo_[p], tag, *comm_);
            }
            else if (sendType == Details::DISTRIBUTOR_ISEND) {
              exports_view_type tmpSendBuf =
                subview_offset (sendArray, size_t(0), numPacketsTo_p);
              requests_.push_back (isend<int> (tmpSendBuf, procsTo_[p],
                                               tag, *comm_));
            }
            else if (sendType == Details::DISTRIBUTOR_SSEND) {
              ssend<int> (tmpSend,
                          as<int> (tmpSend.size ()),
                          procsTo_[p], tag, *comm_);
            }
            else { // if (sendType == Details::DISTRIBUTOR_SSEND)
              send<int> (tmpSend,
                         as<int> (tmpSend.size ()),
                         procsTo_[p], tag, *comm_);
            }
          }
        }
        else { // "Sending" the message to myself
          selfNum = p;
          selfIndex = startsTo_[p];
        }
      }

      if (selfMessage_) {
        for (size_t k = 0; k < lengthsTo_[selfNum]; ++k) {
          deep_copy_offset(imports, exports, selfReceiveOffset,
                           indicesOffsets[selfIndex],
                           numExportPacketsPerLID[selfIndex]);
          selfReceiveOffset += numExportPacketsPerLID[selfIndex];
          ++selfIndex;
        }
      }
      if (verbose) {
        std::ostringstream os;
        os << "Proc " << myProcID
           << ": doPosts(4 args, Kokkos, slow) done" << endl;
        *out_ << os.str ();
      }
    }
  }

  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
  Distributor::
  doReversePostsAndWaits (const ExpView& exports,
                          size_t numPackets,
                          const ImpView& imports)
  {
    doReversePosts (exports, numPackets, imports);
    doReverseWaits ();
  }

  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
  Distributor::
  doReversePostsAndWaits (const ExpView& exports,
                          const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                          const ImpView& imports,
                          const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(requests_.size() != 0, std::runtime_error,
      "Tpetra::Distributor::doReversePostsAndWaits(4 args): There are "
      << requests_.size() << " outstanding nonblocking messages pending.  It "
      "is incorrect to call this method with posts outstanding.");

    doReversePosts (exports, numExportPacketsPerLID, imports,
                    numImportPacketsPerLID);
    doReverseWaits ();
  }

  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
  Distributor::
  doReversePosts (const ExpView &exports,
                  size_t numPackets,
                  const  ImpView &imports)
  {
    // FIXME (mfh 29 Mar 2012) WHY?
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! indicesTo_.empty (), std::runtime_error,
      "Tpetra::Distributor::doReversePosts(3 args): Can only do "
      "reverse communication when original data are blocked by process.");
    if (reverseDistributor_.is_null ()) {
      createReverseDistributor ();
    }
    reverseDistributor_->doPosts (exports, numPackets, imports);
  }

  template <class ExpView, class ImpView>
  typename std::enable_if<(Kokkos::Impl::is_view<ExpView>::value && Kokkos::Impl::is_view<ImpView>::value)>::type
  Distributor::
  doReversePosts (const ExpView &exports,
                  const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                  const ImpView &imports,
                  const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID)
  {
    // FIXME (mfh 29 Mar 2012) WHY?
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! indicesTo_.empty (), std::runtime_error,
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
  computeSends (const Teuchos::ArrayView<const OrdinalType> & importGIDs,
                const Teuchos::ArrayView<const int> & importProcIDs,
                Teuchos::Array<OrdinalType> & exportGIDs,
                Teuchos::Array<int> & exportProcIDs)
  {
    // NOTE (mfh 19 Apr 2012): There was a note on this code saying:
    // "assumes that size_t >= Ordinal".  The code certainly does
    // assume that sizeof(size_t) >= sizeof(OrdinalType) as well as
    // sizeof(size_t) >= sizeof(int).  This is because it casts the
    // OrdinalType elements of importGIDs (along with their
    // corresponding process IDs, as int) to size_t, and does a
    // doPostsAndWaits<size_t>() to send the packed data.
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using std::endl;
    typedef typename ArrayView<const OrdinalType>::size_type size_type;
    const bool verbose = Tpetra::Details::Behavior::verbose("Distributor");

    Teuchos::OSTab tab (out_);
    const int myRank = comm_->getRank ();
    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myRank << ": computeSends" << endl;
      *out_ << os.str ();
    }

    TEUCHOS_TEST_FOR_EXCEPTION(
      importGIDs.size () != importProcIDs.size (), std::invalid_argument,
      "Tpetra::Distributor::computeSends: On Process " << myRank << ": "
      "importProcIDs.size() = " << importProcIDs.size ()
      << " != importGIDs.size() = " << importGIDs.size () << ".");

    const size_type numImports = importProcIDs.size ();
    Array<size_t> importObjs (2*numImports);
    // Pack pairs (importGIDs[i], my process ID) to send into importObjs.
    for (size_type i = 0; i < numImports; ++i) {
      importObjs[2*i]   = static_cast<size_t> (importGIDs[i]);
      importObjs[2*i+1] = static_cast<size_t> (myRank);
    }
    //
    // Use a temporary Distributor to send the (importGIDs[i], myRank)
    // pairs to importProcIDs[i].
    //
    Distributor tempPlan (comm_, out_);
    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myRank << ": computeSends: tempPlan.createFromSends" << endl;
      *out_ << os.str ();
    }

    // mfh 20 Mar 2014: An extra-cautious cast from unsigned to
    // signed, in order to forestall any possible causes for Bug 6069.
    const size_t numExportsAsSizeT = tempPlan.createFromSends (importProcIDs);
    const size_type numExports = static_cast<size_type> (numExportsAsSizeT);
    TEUCHOS_TEST_FOR_EXCEPTION(
      numExports < 0, std::logic_error, "Tpetra::Distributor::computeSends: "
      "tempPlan.createFromSends() returned numExports = " << numExportsAsSizeT
      << " as a size_t, which overflows to " << numExports << " when cast to "
      << Teuchos::TypeNameTraits<size_type>::name () << ".  "
      "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_type> (tempPlan.getTotalReceiveLength ()) != numExports,
      std::logic_error, "Tpetra::Distributor::computeSends: tempPlan.getTotal"
      "ReceiveLength() = " << tempPlan.getTotalReceiveLength () << " != num"
      "Exports = " << numExports  << ".  Please report this bug to the "
      "Tpetra developers.");

    if (numExports > 0) {
      exportGIDs.resize (numExports);
      exportProcIDs.resize (numExports);
    }

    // exportObjs: Packed receive buffer.  (exportObjs[2*i],
    // exportObjs[2*i+1]) will give the (GID, PID) pair for export i,
    // after tempPlan.doPostsAndWaits(...) finishes below.
    //
    // FIXME (mfh 19 Mar 2014) This only works if OrdinalType fits in
    // size_t.  This issue might come up, for example, on a 32-bit
    // machine using 64-bit global indices.  I will add a check here
    // for that case.
    TEUCHOS_TEST_FOR_EXCEPTION(
      sizeof (size_t) < sizeof (OrdinalType), std::logic_error,
      "Tpetra::Distributor::computeSends: sizeof(size_t) = " << sizeof(size_t)
      << " < sizeof(" << Teuchos::TypeNameTraits<OrdinalType>::name () << ") = "
      << sizeof (OrdinalType) << ".  This violates an assumption of the "
      "method.  It's not hard to work around (just use Array<OrdinalType> as "
      "the export buffer, not Array<size_t>), but we haven't done that yet.  "
      "Please report this bug to the Tpetra developers.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      tempPlan.getTotalReceiveLength () < static_cast<size_t> (numExports),
      std::logic_error,
      "Tpetra::Distributor::computeSends: tempPlan.getTotalReceiveLength() = "
      << tempPlan.getTotalReceiveLength() << " < numExports = " << numExports
      << ".  Please report this bug to the Tpetra developers.");

    Array<size_t> exportObjs (tempPlan.getTotalReceiveLength () * 2);
    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myRank << ": computeSends: tempPlan.doPostsAndWaits" << endl;
      *out_ << os.str ();
    }
    tempPlan.doPostsAndWaits<size_t> (importObjs (), 2, exportObjs ());

    // Unpack received (GID, PID) pairs into exportIDs resp. exportProcIDs.
    for (size_type i = 0; i < numExports; ++i) {
      exportGIDs[i] = static_cast<OrdinalType> (exportObjs[2*i]);
      exportProcIDs[i] = static_cast<int> (exportObjs[2*i+1]);
    }

    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myRank << ": computeSends done" << endl;
      *out_ << os.str ();
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
    const bool verbose = Tpetra::Details::Behavior::verbose("Distributor");

    Teuchos::OSTab tab (out_);
    const int myRank = comm_->getRank();

    if (verbose) {
      *out_ << "Proc " << myRank << ": createFromRecvs" << endl;
    }

#ifdef HAVE_TPETRA_DEBUG
    using Teuchos::outArg;
    using Teuchos::reduceAll;

    // In debug mode, first test locally, then do an all-reduce to
    // make sure that all processes passed.
    const int errProc =
      (remoteGIDs.size () != remoteProcIDs.size ()) ? myRank : -1;
    int maxErrProc = -1;
    reduceAll<int, int> (*comm_, Teuchos::REDUCE_MAX, errProc, outArg (maxErrProc));
    TEUCHOS_TEST_FOR_EXCEPTION(maxErrProc != -1, std::runtime_error,
      Teuchos::typeName (*this) << "::createFromRecvs(): lists of remote IDs "
      "and remote process IDs must have the same size on all participating "
      "processes.  Maximum process ID with error: " << maxErrProc << ".");
#else // NOT HAVE_TPETRA_DEBUG

    // In non-debug mode, just test locally.
    TEUCHOS_TEST_FOR_EXCEPTION(
      remoteGIDs.size () != remoteProcIDs.size (), std::invalid_argument,
      Teuchos::typeName (*this) << "::createFromRecvs<" <<
      Teuchos::TypeNameTraits<OrdinalType>::name () << ">(): On Process " <<
      myRank << ": remoteGIDs.size() = " << remoteGIDs.size () << " != "
      "remoteProcIDs.size() = " << remoteProcIDs.size () << ".");
#endif // HAVE_TPETRA_DEBUG

    computeSends (remoteGIDs, remoteProcIDs, exportGIDs, exportProcIDs);

    const size_t numProcsSendingToMe = createFromSends (exportProcIDs ());

    if (verbose) {
      // NOTE (mfh 20 Mar 2014) If remoteProcIDs could contain
      // duplicates, then its length might not be the right check here,
      // even if we account for selfMessage_.  selfMessage_ is set in
      // createFromSends.
      std::ostringstream os;
      os << "Proc " << myRank << ": {numProcsSendingToMe: "
         << numProcsSendingToMe << ", remoteProcIDs.size(): "
         << remoteProcIDs.size () << ", selfMessage_: "
         << (selfMessage_ ? "true" : "false") << "}" << std::endl;
      *out_ << os.str ();
    }

    if (verbose) {
      *out_ << "Proc " << myRank << ": createFromRecvs done" << endl;
    }

    howInitialized_ = Details::DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_RECVS;
  }


} // namespace Tpetra

#endif // TPETRA_DISTRIBUTOR_HPP
