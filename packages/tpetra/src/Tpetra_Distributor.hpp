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

namespace Tpetra {

  namespace {
    // This is an implementation detail of Distributor.  Please do not
    // rely on these values in your code.  We use this to pick the
    // type of send operation that Distributor uses.
    enum EDistributorSendType {
      DISTRIBUTOR_ISEND, // Use MPI_Isend (Teuchos::isend)
      DISTRIBUTOR_RSEND, // Use MPI_Rsend (Teuchos::readySend)
      DISTRIBUTOR_SEND,  // Use MPI_Send (Teuchos::send)
      DISTRIBUTOR_SSEND  // Use MPI_Ssend (Teuchos::ssend)
    };

    // Convert an EDistributorSendType enum value to a string.
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

  } // namespace (anonymous)

  //! Valid values for Distributor's "Send type" parameter.
  Array<std::string> distributorSendTypes ();

  /// \class Distributor
  /// \brief Sets up and executes a communication plan for a Tpetra DistObject.
  ///
  /// This class encapsulates the general information and services
  /// needed for subclasses of \c DistObject (such as CrsMatrix and
  /// MultiVector) to do data redistribution (Import and Export)
  /// operations.
  class Distributor :
    public Teuchos::Describable,
    public Teuchos::ParameterListAcceptorDefaultBase,
    public Teuchos::VerboseObject<Distributor> {
  public:

    //! @name Constructor/Destructor
    //@{

    /// \brief Construct using the specified communicator and default parameters.
    ///
    /// \param comm [in] Communicator used by the Distributor.
    ///
    /// The constructor doesn't actually set up the distribution
    /// pattern.  You need to call one of the "gather / scatter
    /// 'constructors'" to do that.
    explicit Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm);

    /// \brief Construct using the specified communicator and ParameterList.
    ///
    /// \param comm [in] Communicator used by the Distributor.
    ///
    /// \param plist [in/out] List of parameters controlling how the
    ///   Distributor performs communication.  Must be nonnull.
    ///
    /// The constructor doesn't actually set up the distribution
    /// pattern.  You need to call one of the "gather / scatter
    /// 'constructors'" to do that.
    explicit Distributor (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                          const Teuchos::RCP<Teuchos::ParameterList>& plist);

    //! Copy constructor.
    Distributor(const Distributor &distributor);

    //! Destructor (virtual for memory safety).
    virtual ~Distributor();

    //@}
    //! @name Implementation of ParameterListAcceptorDefaultBase
    //@{

    /// \brief Set Distributor parameters.
    ///
    /// See the documentation of getValidParameters() for the current
    /// list of parameters understood by Distributor.
    void setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist);

    /// \brief List of valid Distributor parameters.
    ///
    /// Current parameters include:
    ///
    /// - "Barrier between receives and sends" (bool): Whether to
    ///   execute a barrier between receives and sends in
    ///   do[Reverse]Posts(). n A barrier is required for correctness
    ///   when the "Send type" parameter is "Rsend".  Otherwise, a
    ///   barrier is correct and may be useful for debugging, but not
    ///   recommended, since it introduces useless synchronization.
    ///
    /// - "Send type" (std::string): When using MPI, the variant of
    ///   MPI_Send to use in do[Reverse]Posts().  Valid values include
    ///   "Isend", "Rsend", "Send", and "Ssend".
    ///
    /// - "VerboseObject" (sublist): Optional sublist for controlling
    ///   behavior of Distributor as a Teuchos::VerboseObject.  This
    ///   is currently useful only for debugging.
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters () const;

    //@}
    //! \name Gather/Scatter Constructors
    //@{

    /// \brief Set up Distributor using list of node IDs to which this node will send.
    ///
    /// Take a list of node IDs and construct a plan for efficiently
    /// scattering to those nodes.  Return the number of nodes which
    /// will send me data.
    ///
    /// \param exportNodeIDs [in] List of nodes that will get the
    ///   exported data.  A node ID greater than or equal to the
    ///   number of nodes will result in a \c std::runtime_error on
    ///   all nodes.  Node IDs less than zero are ignored; their
    ///   placement corresponds to null sends in any future
    ///   exports. That is, if <tt>exportNodeIDs[0] == -1</tt>, then
    ///   the corresponding position in the export array is ignored
    ///   during a call to doPosts() or doPostsAndWaits().  For this
    ///   reason, a negative entry is sufficient to break contiguity.
    ///
    /// \return Number of imports this node will be receiving.
    size_t createFromSends (const ArrayView<const int>& exportNodeIDs);

    /// \brief Set up Distributor using list of node IDs from which to receive.
    ///
    /// Take a list of node IDs and construct a plan for efficiently
    /// scattering to those nodes.  Return the number and list of IDs
    /// being sent by me.
    ///
    /// \c Import invokes this method in order to creating a \c
    /// Distributor from a list of receive neighbors and IDs.  A
    /// common use case for this process is setting up sends and
    /// receives for the remote entries of the source vector in a
    /// distributed sparse matrix-vector multiply.  The Mantevo HPCCG
    /// miniapp shows an annotated and simplified version of this
    /// process for that special case.
    ///
    /// \param remoteIDs [in] List of remote IDs wanted.
    ///
    /// \param remoteNodeIDs [in] List of the nodes that will send the
    ///   remote IDs listed in \remoteIDs. Node IDs less than zero are
    ///   ignored; their placement corresponds to null sends in any
    ///   future exports. A node ID greater than or equal to the
    ///   number of nodes will result in an \c std::runtime_error on
    ///   all nodes.
    ///
    /// \param exportIDs [out] List of IDs that need to be sent from
    ///   this node.
    ///
    /// \param exportNodeIDs [out] List of nodes that will get the
    ///   exported IDs in \c exportIDs.
    ///
    /// The \c exportGIDs and \c exportNodeIDs arrays are allocated by
    /// the Distributor, which is why they are passed in a nonconst
    /// reference to an ArrayRCP.  They may be null on entry.
    template <class Ordinal>
    void createFromRecvs(const ArrayView<const Ordinal> &remoteIDs,
                         const ArrayView<const int> &remoteNodeIDs,
                               ArrayRCP<Ordinal> &exportIDs,
                               ArrayRCP<int> &exportNodeIDs);

    //@}

    //! @name Attribute Accessor Methods
    //@{

    //! The number of nodes from which we will receive data, not include this node ("myself").
    size_t getNumReceives() const;

    //! The number of nodes to which we will send data, not include this node ("myself").
    size_t getNumSends() const;

    //! Indicates whether values are being sent to/recieved from this node.
    /*! If we are sending any elements to ourself, returns true. If we aren't, returns false. */
    bool hasSelfMessage() const;

    //! Maximum number of values that this node is sending to another single node.
    size_t getMaxSendLength() const;

    //! Total number of values that this nodes is receiving from other nodes.
    size_t getTotalReceiveLength() const;

    //! A list of images sending values to this node. (non-persisting view)
    ArrayView<const int> getImagesFrom() const;

    //! A list of images to which this node is sending values. (non-persisting view)
    ArrayView<const int> getImagesTo() const;

    //! Number of values we're receiving from each node. (non-persisting view)
    /*! We will receive <tt>getLengthsFrom[i]</tt> values from node <tt>getImagesFrom[i]</tt>. */
    ArrayView<const size_t> getLengthsFrom() const;

    //! Number of values we're sending to each node. (non-persisting view)
    /*! We will send <tt>getLengthsTo[i]</tt> values to image <tt>getImagesTo[i]</tt>. */
    ArrayView<const size_t> getLengthsTo() const;

    //@}

    //! @name Reverse Communication Methods
    //@{

    //! \brief Returns a Distributor with a reverse plan of this Distributor's plan
    /*! This method creates the reverse Distributor the first time the function
        is called.
    */
    const RCP<Distributor>& getReverse() const;

    //@}

    //! @name Execute Distributor Plan Methods
    //@{

    /// \brief Execute the (forward) communication plan.
    ///
    /// Call this overload when you have the same number of Packets
    /// for each LID to send or receive.
    ///
    /// \param exports [in] Contains the values to be sent by this
    ///   node.  On exit from this method, it's OK to modify the
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
    doPostsAndWaits (const ArrayView<const Packet> &exports,
                     size_t numPackets,
                     const ArrayView<Packet> &imports);

    /// \brief Execute the (forward) communication plan.
    ///
    /// Call this overload when you have possibly different numbers of
    /// Packets for each LID to send or receive.
    ///
    /// \param exports [in] Contains the values to be sent by this
    ///   node.  On exit from this method, it's OK to modify the
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
    doPostsAndWaits (const ArrayView<const Packet> &exports,
                     const ArrayView<size_t> &numExportPacketsPerLID,
                     const ArrayView<Packet> &imports,
                     const ArrayView<size_t> &numImportPacketsPerLID);

    /// \brief Post the data for a forward plan, but do not execute the waits yet.
    ///
    /// Call this overload when you have the same number of Packets
    /// for each LID to send or receive.
    ///
    /// \param exports [in] Contains the values to be sent by this
    ///   node.  This is an ArrayRCP and not an ArrayView so that we
    ///   have the freedom to use nonblocking sends if we wish.  Do
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
    doPosts (const ArrayRCP<const Packet> &exports,
             size_t numPackets,
             const ArrayRCP<Packet> &imports);

    /// \brief Post the data for a forward plan, but do not execute the waits yet.
    ///
    /// Call this overload when you have possibly different numbers of
    /// Packets for each LID to send or receive.
    ///
    /// \param exports [in] Same as in the three-argument version of
    ///   \c doPosts().
    ///
    /// \param numExportPacketsPerLID [in] Same as in the
    ///   four-argument version of \c doPostsAndWaits().
    ///
    /// \param imports [out] Same as in the three-argument version of
    ///   \c doPosts().
    ///
    /// \param numImportPacketsPerLID [in] Same as in the
    ///   four-argument version of \c doPostsAndWaits().
    template <class Packet>
    void
    doPosts (const ArrayRCP<const Packet> &exports,
             const ArrayView<size_t> &numExportPacketsPerLID,
             const ArrayRCP<Packet> &imports,
             const ArrayView<size_t> &numImportPacketsPerLID);

    //! Wait on any outstanding nonblocking message requests to complete.
    void doWaits ();

    /// \brief Execute the reverse communication plan.
    ///
    /// This method takes the same arguments as the three-argument
    /// version of \c doPostsAndWaits().
    template <class Packet>
    void
    doReversePostsAndWaits (const ArrayView<const Packet> &exports,
                            size_t numPackets,
                            const ArrayView<Packet> &imports);

    /// \brief Execute the reverse communication plan.
    ///
    /// This method takes the same arguments as the four-argument
    /// version of \c doPostsAndWaits().
    template <class Packet>
    void
    doReversePostsAndWaits (const ArrayView<const Packet> &exports,
                            const ArrayView<size_t> &numExportPacketsPerLID,
                            const ArrayView<Packet> &imports,
                            const ArrayView<size_t> &numImportPacketsPerLID);

    /// \brief Post the data for a reverse plan, but do not execute the waits yet.
    ///
    /// This method takes the same arguments as the three-argument
    /// version of \c doPosts().
    template <class Packet>
    void
    doReversePosts (const ArrayRCP<const Packet> &exports,
                    size_t numPackets,
                    const ArrayRCP<Packet> &imports);

    /// \brief Post the data for a reverse plan, but do not execute the waits yet.
    ///
    /// This method takes the same arguments as the four-argument
    /// version of \c doPosts().
    template <class Packet>
    void
    doReversePosts (const ArrayRCP<const Packet> &exports,
                    const ArrayView<size_t> &numExportPacketsPerLID,
                    const ArrayRCP<Packet> &imports,
                    const ArrayView<size_t> &numImportPacketsPerLID);

    //! Wait on any outstanding reverse message requests to complete.
    void doReverseWaits ();

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! A simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an \c FancyOStream.
    void describe (Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

    //@}

  private:
    //! The communicator over which to perform distributions.
    RCP<const Comm<int> > comm_;

    //! @name Parameters read in from \c Teuchos::ParameterList
    //@{
    //! The variant of send to use in do[Reverse]Posts().
    EDistributorSendType sendType_;
    //! Whether to do a barrier between receives and sends in do[Reverse]Posts().
    bool barrierBetween_;
    //@}

    /// \brief The number of export process IDs on input to \c createFromSends().
    ///
    /// This may differ from the number of sends.  We always want to
    /// send either zero or one messages to any process.  However, the
    /// user may have specified a process ID twice in \c
    /// createFromSends()'s input array of process IDs (\c
    /// exportNodeIDs).  This is allowed, but may affect whether sends
    /// require a buffer.
    size_t numExports_;

    //! Whether I am supposed to send a message to myself.
    bool selfMessage_;

    /// \brief The number of sends to other nodes.
    ///
    /// This is always less than or equal to the number of nodes.
    /// It does <i>not</i> count self receives.
    ///
    /// This value is computed by the \c createFromSends() method.
    /// That method first includes self receives in the count, but at
    /// the end subtracts one if selfMessage_ is true.
    size_t numSends_;

    /// \brief List of process IDs to which to send.
    ///
    /// This array has length numSends_ + selfMessage_ (that is, it
    /// includes the self message, if there is one).
    Array<int> imagesTo_;

    /// \brief Starting index of the block of Packets to send to each process.
    ///
    /// Given an export buffer that contains all of the data being
    /// sent by this process, the block of Packets to send to process
    /// p will start at position startsTo_[p].
    ///
    /// This array has length numSends_ + selfMessage_ (that is, it
    /// includes the self message, if there is one).
    Array<size_t> startsTo_;

    /// \brief Length (in number of Packets) of my process' send to each process.
    ///
    /// lengthsTo_[p] is the length of my process' send to process p.
    /// This array has length numSends_ + selfMessage_ (that is, it
    /// includes the self message, if there is one).
    Array<size_t> lengthsTo_;

    /// \brief The maximum send length (in number of Packets) to another process.
    ///
    /// maxSendLength_ = max(lengthsTo_[p]) for p != my process ID.
    size_t maxSendLength_;

    /// \brief Offset (by message, not by number of Packets) into exports array.
    ///
    /// This array is used by the three-argument version of \c
    /// doPosts().  In that method, \c indicesTo_[j]*numPackets is the
    /// offset into the \c exports array, where j = startsTo_[p] and p
    /// is an index iterating through the sends in reverse order
    /// (starting with the process ID right before the self message,
    /// if there is a self message, else the largest process ID to
    /// which this process sends).
    ///
    /// This array is only used if export data are not blocked (laid
    /// out) by process ID, i.e., if we need to use a send buffer.
    /// Otherwise, this array has no entries (in fact, Distributor
    /// currently uses this in both overloads of \c doPosts() to test
    /// whether data are laid out by process).
    Array<size_t> indicesTo_;

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

    // imagesFrom_, startsFrom_ and lengthsFrom_ each have size
    //   numReceives_ + selfMessage_
    Array<size_t> lengthsFrom_;
    Array<int> imagesFrom_;
    Array<size_t> startsFrom_;

    /// \brief List that becomes the reverse communicator's indicesTo_.
    ///
    /// Array of length \c totalReceiveLength_.  Allocated and filled
    /// in \c computeReceives().
    Array<size_t> indicesFrom_;

    /// \brief Communication requests associated with nonblocking receives and sends.
    ///
    /// \note To implementers: Distributor uses requests_.size() as
    /// the number of outstanding nonblocking receives and sends.
    /// This means you should always resize to zero after completing
    /// receive and send requests.
    Array<RCP<Teuchos::CommRequest<int> > > requests_;

    /// \brief The reverse distributor.
    ///
    /// This is created on demand in \c getReverse() and cached for
    /// later reuse.  This is why it is declared "mutable".
    mutable RCP<Distributor> reverseDistributor_;

    /// \brief Compute receive info from sends.
    ///
    /// This method computes numReceives_, lengthsFrom_, imagesFrom_,
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
    /// \param importIDs [in] GIDs to receive by my process.
    /// \param importNodeIDs [in] Process IDs from which to receive by
    ///   my process.
    /// \param exportIDs [out] GIDs to send by my process.
    /// \param exportNodeIDs [out] Process IDs to which to send by my
    ///   process.
    template <class Ordinal>
    void computeSends (const ArrayView<const Ordinal> &importIDs,
                       const ArrayView<const int> &importNodeIDs,
                       ArrayRCP<Ordinal> &exportIDs,
                       ArrayRCP<int> &exportNodeIDs);

    //! Create a distributor for the reverse communication pattern.
    void createReverseDistributor() const;

  }; // class Distributor


  template <class Packet>
  void Distributor::
  doPostsAndWaits (const ArrayView<const Packet>& exports,
                   size_t numPackets,
                   const ArrayView<Packet>& imports)
  {
    using Teuchos::as;
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;

    TEUCHOS_TEST_FOR_EXCEPTION(requests_.size() != 0, std::runtime_error,
      Teuchos::typeName(*this) << "::doPostsAndWaits(): There are "
      << requests_.size() << " outstanding nonblocking messages pending.  It is "
      "incorrect to call doPostsAndWaits with posts outstanding.");

    // doPosts() accepts the exports and imports arrays as ArrayRCPs,
    // requiring that the memory location is persisting (as is
    // necessary for nonblocking receives).  However, it need only
    // persist until doWaits() completes, so it is safe for us to use
    // a nonpersisting reference in this case.  The use of a
    // nonpersisting reference is purely a performance optimization.

    typedef typename ArrayRCP<const Packet>::size_type size_type;
    //const Packet* exportsPtr = exports.getRawPtr();
    //ArrayRCP<const Packet> exportsArcp (exportsPtr, as<size_type> (0), exports.size(), false);
    ArrayRCP<const Packet> exportsArcp (exports.getRawPtr(), as<size_type> (0),
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
             arcp<Packet> (imports.getRawPtr(), 0, imports.size(), false));
    doWaits();
  }

  template <class Packet>
  void Distributor::
  doPostsAndWaits (const ArrayView<const Packet>& exports,
                   const ArrayView<size_t> &numExportPacketsPerLID,
                   const ArrayView<Packet> &imports,
                   const ArrayView<size_t> &numImportPacketsPerLID)
  {
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

    TEUCHOS_TEST_FOR_EXCEPTION(requests_.size() != 0, std::runtime_error,
      Teuchos::typeName(*this) << "::doPostsAndWaits(): There are "
      << requests_.size() << " outstanding nonblocking messages pending.  It is "
      "incorrect to call doPostsAndWaits with posts outstanding.");

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
    ArrayRCP<const Packet> exportsArcp (exports.getRawPtr(), as<size_type> (0),
                                        exports.size(), false);
    // mfh 04 Apr 2012: This is the offending code.  This statement
    // would normally be in place of "exportsArcp" in the
    // doPosts() call below.
    //arcp<const Packet> (exports.getRawPtr(), 0, exports.size(), false),
    doPosts (exportsArcp,
             numExportPacketsPerLID,
             arcp<Packet> (imports.getRawPtr(), 0, imports.size(), false),
             numImportPacketsPerLID);
    doWaits();
  }


  template <class Packet>
  void Distributor::
  doPosts (const ArrayRCP<const Packet>& exports,
           size_t numPackets,
           const ArrayRCP<Packet>& imports)
  {
    using Teuchos::as;
    using Teuchos::FancyOStream;
    using Teuchos::includesVerbLevel;
    using Teuchos::ireceive;
    using Teuchos::isend;
    using Teuchos::OSTab;
    using Teuchos::readySend;
    using Teuchos::send;
    using Teuchos::ssend;
    using Teuchos::typeName;
    using std::endl;

    // Run-time configurable parameters that come from the input
    // ParameterList set by setParameterList().
    const EDistributorSendType sendType = sendType_;
    const bool doBarrier = barrierBetween_;

#ifdef HAVE_TEUCHOS_DEBUG
    // Prepare for verbose output, if applicable.
    Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel ();
    RCP<FancyOStream> out = this->getOStream ();
    const bool doPrint = out.get () && (comm_->getRank () == 0) &&
      includesVerbLevel (verbLevel, Teuchos::VERB_EXTREME, true);

    if (doPrint) {
      // Only need one process to print out parameters.
      *out << "Distributor::doPosts" << endl;
    }
    // Add one tab level.  We declare this outside the doPrint scopes
    // so that the tab persists until the end of this method.
    OSTab tab = this->getOSTab();
    if (doPrint) {
      *out << "sendType=" << DistributorSendTypeEnumToString (sendType)
           << ", barrierBetween=" << doBarrier << endl;
    }
#endif // HAVE_TEUCHOS_DEBUG

    TEUCHOS_TEST_FOR_EXCEPTION(sendType == DISTRIBUTOR_RSEND && ! doBarrier,
      std::logic_error, "Ready send implementation requires a barrier between "
      "posting receives and posting ready sends.  This should have been checked "
      "before.  Please report this bug to the Tpetra developers.");

    const int myImageID = comm_->getRank();
    size_t selfReceiveOffset = 0;

#ifdef HAVE_TEUCHOS_DEBUG
    // Each message has the same number of packets.
    const size_t totalNumImportPackets = totalReceiveLength_ * numPackets;
    TEUCHOS_TEST_FOR_EXCEPTION(as<size_t> (imports.size ()) != totalNumImportPackets,
      std::runtime_error, typeName (*this) << "::doPosts(): imports must be "
      "large enough to store the imported data.  imports.size() = "
      << imports.size() << ", but total number of import packets = "
      << totalNumImportPackets << ".");
#endif // HAVE_TEUCHOS_DEBUG

    // Distributor uses requests_.size() as the number of outstanding
    // nonblocking message requests, so we resize to zero to maintain
    // this invariant.
    //
    // numReceives_ does _not_ include the self message, if there is
    // one.  Here, we do actually send a message to ourselves, so we
    // include any self message in the "actual" number of receives to
    // post.
    //
    // QUESTION (mfh 05 Dec 2012): Does Teuchos::Comm's ireceive() do
    // the right thing with self messages?  It should copy from the
    // send buffer into the receive buffer; they need not (and in fact
    // should not, according to MPI semantics) be the same.  This
    // doesn't matter for us, because we do the copying by hand
    // (below) rather than requiring that Teuchos::Comm handle it.
    //
    // NOTE (mfh 19 Mar 2012): Epetra_MpiDistributor::DoPosts()
    // doesn't (re)allocate its array of requests.  That happens in
    // CreateFromSends(), ComputeRecvs_(), DoReversePosts() (on
    // demand), or Resize_().
    const size_t actualNumReceives = numReceives_ + (selfMessage_ ? 1 : 0);
    requests_.resize (0);

    // Post the nonblocking receives.  It's common MPI wisdom to post
    // receives before sends.  In MPI terms, this means favoring
    // adding to the "posted queue" (of receive requests) over adding
    // to the "unexpected queue" (of arrived messages not yet matched
    // with a receive).
    {
      size_t curBufferOffset = 0;
      for (size_t i = 0; i < actualNumReceives; ++i) {
        if (imagesFrom_[i] != myImageID) {
          // If my process is receiving these packet(s) from another
          // process (not a self-receive):
          //
          // 1. Set up the persisting view (recvBuf) of the imports
          //    array, given the offset and size (total number of
          //    packets from process imagesFrom_[i]).
          // 2. Start the Irecv and save the resulting request.
          ArrayRCP<Packet> recvBuf =
            imports.persistingView (curBufferOffset, lengthsFrom_[i]*numPackets);
          requests_.push_back (ireceive<int, Packet> (*comm_, recvBuf, imagesFrom_[i]));
        }
        else { // Receiving from myself
          selfReceiveOffset = curBufferOffset; // Remember the self-recv offset
        }
        curBufferOffset += lengthsFrom_[i]*numPackets;
      }
    }

    if (doBarrier) {
      // If we are using ready sends (MPI_Rsend) below, we need to do
      // a barrier before we post the ready sends.  This is because a
      // ready send requires that its matching receive has already
      // been posted before the send has been posted.  The only way to
      // guarantee that in this case is to use a barrier.
      Teuchos::barrier (*comm_);
    }

    // setup scan through imagesTo_ list starting with higher numbered images
    // (should help balance message traffic)
    size_t numBlocks = numSends_ + selfMessage_;
    size_t imageIndex = 0;
    while ((imageIndex < numBlocks) && (imagesTo_[imageIndex] < myImageID)) {
      ++imageIndex;
    }
    if (imageIndex == numBlocks) {
      imageIndex = 0;
    }

    size_t selfNum = 0;
    size_t selfIndex = 0;

    if (indicesTo_.empty()) {
      // Data are already blocked (laid out) by process, so we don't
      // need a separate send buffer (besides the exports array).
      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + imageIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (imagesTo_[p] != myImageID) {
          ArrayView<const Packet> tmpSend =
            exports.view (startsTo_[p]*numPackets, lengthsTo_[p]*numPackets);
          if (sendType == DISTRIBUTOR_RSEND) {
            readySend<int,Packet> (*comm_, tmpSend, imagesTo_[p]);
          }
          else if (sendType == DISTRIBUTOR_ISEND) {
            ArrayRCP<const Packet> tmpSendBuf =
              exports.persistingView (startsTo_[p] * numPackets,
                                      lengthsTo_[p] * numPackets);
            requests_.push_back (isend<int, Packet> (*comm_, tmpSendBuf, imagesTo_[p]));
          }
          else if (sendType == DISTRIBUTOR_SSEND) {
	    // "ssend" means "synchronous send."
            ssend<int, Packet> (*comm_, tmpSend.size(),
                                tmpSend.getRawPtr(), imagesTo_[p]);

          } else { // if (sendType == DISTRIBUTOR_SEND)
            // We've already validated sendType, so it has to be
            // DISTRIBUTOR_SEND.  If it's not, well, this is a
            // reasonable fallback.
            //
            // FIXME (mfh 23 Mar 2012) Implement a three-argument
            // version of send() that takes an ArrayView instead of a
            // raw array.
            send<int, Packet> (*comm_, tmpSend.size(),
                               tmpSend.getRawPtr(), imagesTo_[p]);
          }
        }
        else { // "Sending" the message to myself
          selfNum = p;
        }
      }

      if (selfMessage_) {
	// This is how we "send a message to ourself": we copy from
	// the export buffer to the import buffer.  That saves
	// Teuchos::Comm implementations the trouble of implementing
	// self messages correctly.  (To do this right, the
	// Teuchos::Comm subclass instance would need internal buffer
	// space for messages, keyed on the message's tag.)
        std::copy (exports.begin()+startsTo_[selfNum]*numPackets,
                   exports.begin()+startsTo_[selfNum]*numPackets+lengthsTo_[selfNum]*numPackets,
                   imports.begin()+selfReceiveOffset);
      }
    }
    else { // data are not blocked by image, use send buffer
      ArrayRCP<Packet> sendArray (maxSendLength_ * numPackets); // send buffer

      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + imageIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (imagesTo_[p] != myImageID) {
          typename ArrayView<const Packet>::iterator srcBegin, srcEnd;
          size_t sendArrayOffset = 0;
          size_t j = startsTo_[p];
          for (size_t k = 0; k < lengthsTo_[p]; ++k, ++j) {
            srcBegin = exports.begin() + indicesTo_[j]*numPackets;
            srcEnd   = srcBegin + numPackets;
            std::copy (srcBegin, srcEnd, sendArray.begin()+sendArrayOffset);
            sendArrayOffset += numPackets;
          }
          ArrayView<const Packet> tmpSend = sendArray.view (0, lengthsTo_[p]*numPackets);

          if (sendType == DISTRIBUTOR_RSEND) {
            readySend<int,Packet> (*comm_, tmpSend, imagesTo_[p]);
          }
          else if (sendType == DISTRIBUTOR_ISEND) {
            ArrayRCP<const Packet> tmpSendBuf =
              sendArray.persistingView (0, lengthsTo_[p] * numPackets);
            requests_.push_back (isend<int, Packet> (*comm_, tmpSendBuf, imagesTo_[p]));
          }
          else if (sendType == DISTRIBUTOR_SSEND) {
            ssend<int,Packet> (*comm_, tmpSend.size(),
                               tmpSend.getRawPtr(), imagesTo_[p]);
          }
          else { // if (sendType == DISTRIBUTOR_SEND)
            // We've already validated sendType, so it has to be
            // DISTRIBUTOR_SEND.  If it's not, well, this is a
            // reasonable fallback.
            send<int,Packet> (*comm_, tmpSend.size(),
                              tmpSend.getRawPtr(), imagesTo_[p]);
          }
        }
        else { // "Sending" the message to myself
          selfNum = p;
          selfIndex = startsTo_[p];
        }
      }

      if (selfMessage_) {
        for (size_t k = 0; k < lengthsTo_[selfNum]; ++k) {
          std::copy (exports.begin()+indicesTo_[selfIndex]*numPackets,
                     exports.begin()+indicesTo_[selfIndex]*numPackets + numPackets,
                     imports.begin() + selfReceiveOffset);
          ++selfIndex;
          selfReceiveOffset += numPackets;
        }
      }
    }
  }

  template <class Packet>
  void Distributor::
  doPosts (const ArrayRCP<const Packet>& exports,
           const ArrayView<size_t>& numExportPacketsPerLID,
           const ArrayRCP<Packet>& imports,
           const ArrayView<size_t>& numImportPacketsPerLID)
  {
    using Teuchos::as;
    using Teuchos::ireceive;
    using Teuchos::isend;
    using Teuchos::readySend;
    using Teuchos::send;
    using Teuchos::ssend;
#ifdef HAVE_TEUCHOS_DEBUG
    using std::endl;
#endif // HAVE_TEUCHOS_DEBUG

    // Run-time configurable parameters that come from the input
    // ParameterList set by setParameterList().
    const EDistributorSendType sendType = sendType_;
    const bool doBarrier = barrierBetween_;

#ifdef HAVE_TEUCHOS_DEBUG
    // Prepare for verbose output, if applicable.
    Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel ();
    RCP<Teuchos::FancyOStream> out = this->getOStream ();
    const bool doPrint = out.get () && (comm_->getRank () == 0) &&
      includesVerbLevel (verbLevel, Teuchos::VERB_EXTREME, true);

    if (doPrint) {
      // Only need one process to print out parameters.
      *out << "Distributor::doPosts" << endl;
    }
    // Add one tab level.  We declare this outside the doPrint scopes
    // so that the tab persists until the end of this method.
    Teuchos::OSTab tab = this->getOSTab();
    if (doPrint) {
      *out << "sendType=" << DistributorSendTypeEnumToString (sendType)
           << ", barrierBetween=" << doBarrier << endl;
    }
#endif // HAVE_TEUCHOS_DEBUG

    TEUCHOS_TEST_FOR_EXCEPTION(sendType == DISTRIBUTOR_RSEND && ! doBarrier,
      std::logic_error, "Ready send implementation requires a barrier between "
      "posting receives and posting ready sends.  This should have been checked "
      "before.  Please report this bug to the Tpetra developers.");

    const int myImageID = comm_->getRank();
    size_t selfReceiveOffset = 0;

#ifdef HAVE_TEUCHOS_DEBUG
    // Different messages may have different numbers of packets.
    size_t totalNumImportPackets = 0;
    for (int ii = 0; ii < numImportPacketsPerLID.size(); ++ii) {
      totalNumImportPackets += numImportPacketsPerLID[ii];
    }
    TEUCHOS_TEST_FOR_EXCEPTION(as<size_t> (imports.size ()) != totalNumImportPackets,
      std::runtime_error, Teuchos::typeName (*this) << "::doPosts(): The "
      "imports array argument must be large enough to store the imported "
      "data.  imports.size() = " << imports.size() << ", but the total number "
      "of packets is " << totalNumImportPackets << ".");
#endif // HAVE_TEUCHOS_DEBUG

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
    const size_t actualNumReceives = numReceives_ + (selfMessage_ ? 1 : 0);
    requests_.resize (0);

    // Post the nonblocking receives.  It's common MPI wisdom to post
    // receives before sends.  In MPI terms, this means favoring
    // adding to the "posted queue" (of receive requests) over adding
    // to the "unexpected queue" (of arrived messages not yet matched
    // with a receive).
    {
      size_t curBufferOffset = 0;
      size_t curLIDoffset = 0;
      for (size_t i = 0; i < actualNumReceives; ++i) {
        size_t totalPacketsFrom_i = 0;
        for (size_t j = 0; j < lengthsFrom_[i]; ++j) {
          totalPacketsFrom_i += numImportPacketsPerLID[curLIDoffset+j];
        }
        curLIDoffset += lengthsFrom_[i];
        if (imagesFrom_[i] != myImageID && totalPacketsFrom_i) {
          // If my process is receiving these packet(s) from another
          // process (not a self-receive), and if there is at least
          // one packet to receive:
          //
          // 1. Set up the persisting view (recvBuf) into the imports
          //    array, given the offset and size (total number of
          //    packets from process imagesFrom_[i]).
          // 2. Start the Irecv and save the resulting request.
          ArrayRCP<Packet> recvBuf =
            imports.persistingView (curBufferOffset, totalPacketsFrom_i);
          requests_.push_back (ireceive<int, Packet> (*comm_, recvBuf, imagesFrom_[i]));
        }
        else { // Receiving these packet(s) from myself
          selfReceiveOffset = curBufferOffset; // Remember the offset
        }
        curBufferOffset += totalPacketsFrom_i;
      }
    }

    if (doBarrier) {
      // Each ready-send below requires that its matching receive has
      // already been posted, so do a barrier to ensure that all the
      // nonblocking receives have posted first.
      Teuchos::barrier (*comm_);
    }

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

    // setup scan through imagesTo_ list starting with higher numbered images
    // (should help balance message traffic)
    size_t numBlocks = numSends_+ selfMessage_;
    size_t imageIndex = 0;
    while ((imageIndex < numBlocks) && (imagesTo_[imageIndex] < myImageID)) {
      ++imageIndex;
    }
    if (imageIndex == numBlocks) {
      imageIndex = 0;
    }

    size_t selfNum = 0;
    size_t selfIndex = 0;

    if (indicesTo_.empty()) {
      // Data are already blocked (laid out) by process, so we don't
      // need a separate send buffer (besides the exports array).
      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + imageIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (imagesTo_[p] != myImageID && packetsPerSend[p] > 0) {
          ArrayView<const Packet> tmpSend =
            exports.view (sendPacketOffsets[p], packetsPerSend[p]);
          if (sendType == DISTRIBUTOR_RSEND) {
            readySend<int,Packet> (*comm_, tmpSend, imagesTo_[p]);
          }
          else if (sendType == DISTRIBUTOR_ISEND) {
            ArrayRCP<const Packet> tmpSendBuf =
              exports.persistingView (sendPacketOffsets[p], packetsPerSend[p]);
            requests_.push_back (isend<int, Packet> (*comm_, tmpSendBuf, imagesTo_[p]));
          }
          else if (sendType == DISTRIBUTOR_SSEND) {
            ssend<int, Packet> (*comm_, tmpSend.size(),
                                tmpSend.getRawPtr(), imagesTo_[p]);
          }
          else { // if (sendType == DISTRIBUTOR_SEND)
            // We've already validated sendType, so it has to be
            // DISTRIBUTOR_SEND.  If it's not, well, this is a
            // reasonable fallback.
            send<int, Packet> (*comm_, tmpSend.size(),
                               tmpSend.getRawPtr(), imagesTo_[p]);
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
    }
    else { // data are not blocked by image, use send buffer
      ArrayRCP<Packet> sendArray (maxNumPackets); // send buffer

      Array<size_t> indicesOffsets (numExportPacketsPerLID.size(), 0);
      size_t ioffset = 0;
      for (int j=0; j<numExportPacketsPerLID.size(); ++j) {
        indicesOffsets[j] = ioffset;
        ioffset += numExportPacketsPerLID[j];
      }

      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + imageIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (imagesTo_[p] != myImageID) {
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
            ArrayView<const Packet> tmpSend = sendArray.view (0, numPacketsTo_p);

            if (sendType == DISTRIBUTOR_RSEND) {
              readySend<int,Packet> (*comm_,tmpSend,imagesTo_[p]);
            }
            else if (sendType == DISTRIBUTOR_ISEND) {
              ArrayRCP<const Packet> tmpSendBuf =
                sendArray.persistingView (0, numPacketsTo_p);
              requests_.push_back (isend<int, Packet> (*comm_, tmpSendBuf,
                                                       imagesTo_[p]));
            }
            else if (sendType == DISTRIBUTOR_SSEND) {
              ssend<int,Packet> (*comm_, tmpSend.size(), tmpSend.getRawPtr(), imagesTo_[p]);
            }
            else { // if (sendType == DISTRIBUTOR_SSEND)
              send<int,Packet> (*comm_, tmpSend.size(), tmpSend.getRawPtr(), imagesTo_[p]);
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
    }
  }

  template <class Packet>
  void Distributor::
  doReversePostsAndWaits (const ArrayView<const Packet>& exports,
                          size_t numPackets,
                          const ArrayView<Packet>& imports)
  {
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
  }

  template <class Packet>
  void Distributor::
  doReversePostsAndWaits (const ArrayView<const Packet>& exports,
                          const ArrayView<size_t> &numExportPacketsPerLID,
                          const ArrayView<Packet> &imports,
                          const ArrayView<size_t> &numImportPacketsPerLID)
  {
    using Teuchos::as;
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;

    TEUCHOS_TEST_FOR_EXCEPTION(requests_.size() != 0, std::runtime_error,
      Teuchos::typeName(*this) << "::doReversePostsAndWaits(): There are "
      << requests_.size() << " outstanding nonblocking messages pending.  It is "
      "incorrect to call doReversePostsAndWaits with posts outstanding.");

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
  }

  template <class Packet>
  void Distributor::
  doReversePosts (const ArrayRCP<const Packet>& exports,
                  size_t numPackets,
                  const ArrayRCP<Packet>& imports)
  {
    // FIXME (mfh 29 Mar 2012) WHY?
    TEUCHOS_TEST_FOR_EXCEPTION(! indicesTo_.empty (), std::runtime_error,
      Teuchos::typeName (*this) << "::doReversePosts(): Can only do reverse "
      "communication when original data are blocked by process.");
    if (reverseDistributor_.is_null ()) {
      createReverseDistributor ();
    }
    reverseDistributor_->doPosts (exports, numPackets, imports);
  }

  template <class Packet>
  void Distributor::
  doReversePosts (const ArrayRCP<const Packet>& exports,
                  const ArrayView<size_t>& numExportPacketsPerLID,
                  const ArrayRCP<Packet>& imports,
                  const ArrayView<size_t>& numImportPacketsPerLID)
  {
    // FIXME (mfh 29 Mar 2012) WHY?
    TEUCHOS_TEST_FOR_EXCEPTION(! indicesTo_.empty (), std::runtime_error,
      Teuchos::typeName (*this) << "::doReversePosts(): Can only do reverse "
      "communication when original data are blocked by process.");
    if (reverseDistributor_.is_null ()) {
      createReverseDistributor ();
    }
    reverseDistributor_->doPosts (exports, numExportPacketsPerLID,
                                  imports, numImportPacketsPerLID);
  }

  template <class OrdinalType>
  void Distributor::
  computeSends (const ArrayView<const OrdinalType> & importIDs,
                const ArrayView<const int> & importNodeIDs,
                ArrayRCP<OrdinalType> & exportIDs,
                ArrayRCP<int> & exportNodeIDs)
  {
    // NOTE (mfh 19 Apr 2012): There was a note on this code saying:
    // "assumes that size_t >= Ordinal".  The code certainly does
    // assume that sizeof(size_t) >= sizeof(OrdinalType) as well as
    // sizeof(size_t) >= sizeof(int).  This is because it casts the
    // OrdinalType elements of importIDs (along with their
    // corresponding process IDs, as int) to size_t, and does a
    // doPostsAndWaits<size_t>() to send the packed data.
    using Teuchos::as;

    const int myRank = comm_->getRank();
    const size_t numImports = importNodeIDs.size();
    TEUCHOS_TEST_FOR_EXCEPTION(as<size_t> (importIDs.size ()) < numImports,
      std::invalid_argument, "Tpetra::Distributor::computeSends: importNodeIDs."
      "size() = " << importNodeIDs.size () << " != importIDs.size() = "
      << importIDs.size () << ".");

    Array<size_t> importObjs (2*numImports);
    // Pack pairs (importIDs[i], my process ID) to send into importObjs.
    for (size_t i = 0; i < numImports; ++i ) {
      importObjs[2*i]   = as<size_t> (importIDs[i]);
      importObjs[2*i+1] = as<size_t> (myRank);
    }
    //
    // Use a temporary Distributor to send the (importIDs[i], myRank)
    // pairs to importNodeIDs[i].
    //
    size_t numExports;
    Distributor tempPlan (comm_);
    numExports = tempPlan.createFromSends (importNodeIDs);
    if (numExports > 0) {
      exportIDs = arcp<OrdinalType> (numExports);
      exportNodeIDs = arcp<int> (numExports);
    }
    Array<size_t> exportObjs (tempPlan.getTotalReceiveLength () * 2);
    tempPlan.doPostsAndWaits<size_t> (importObjs (), 2, exportObjs ());

    // Unpack received (GID, PID) pairs into exportIDs resp. exportNodeIDs.
    for (size_t i = 0; i < numExports; ++i) {
      exportIDs[i]     = as<OrdinalType>(exportObjs[2*i]);
      exportNodeIDs[i] = exportObjs[2*i+1];
    }
  }

  template <class OrdinalType>
  void Distributor::
  createFromRecvs (const ArrayView<const OrdinalType> &remoteIDs,
                   const ArrayView<const int> &remoteImageIDs,
                   ArrayRCP<OrdinalType> &exportGIDs,
                   ArrayRCP<int> &exportNodeIDs)
  {
    using Teuchos::outArg;
    using Teuchos::reduceAll;

    const int myRank = comm_->getRank();
    const int errProc =
      (remoteIDs.size () != remoteImageIDs.size ()) ? myRank : -1;
    int maxErrProc = -1;
    reduceAll (*comm_, Teuchos::REDUCE_MAX, errProc, outArg (maxErrProc));
    TEUCHOS_TEST_FOR_EXCEPTION(maxErrProc != -1, std::runtime_error,
      Teuchos::typeName (*this) << "::createFromRecvs(): lists of remote IDs "
      "and remote process IDs must have the same size on all participating "
      "processes.  Maximum process ID with error: " << maxErrProc << ".");

    computeSends (remoteIDs, remoteImageIDs, exportGIDs, exportNodeIDs);
    (void) createFromSends (exportNodeIDs ());
  }

} // namespace Tpetra

#endif // TPETRA_DISTRIBUTOR_HPP
