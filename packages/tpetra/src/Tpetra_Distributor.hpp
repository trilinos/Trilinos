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


// FINISH: some of the get accessors may not be necessary anymore. clean up.
// FINISH: This class may not be const correct. doPosts() et al. perhaps should be const, with affected members made mutable.

namespace Tpetra {

  /// \class Distributor
  /// \brief Class that sets up gathers and scatters for Tpetra communication.
  ///
  /// This class encapsulates the general information and services
  /// needed for other Tpetra classes to perform gather/scatter
  /// operations on a parallel computer.
  class Distributor : public Teuchos::Describable {
  public:

    //! @name Constructor/Destructor
    //@{ 

    /// \brief Construct the Distributor using the specified communicator.
    ///
    /// This doesn't actually set up the distribution pattern.  You
    /// need to call one of the "gather / scatter 'constructors'" to
    /// do that.
    explicit Distributor(const RCP<const Comm<int> > & comm);

    //! Copy constructor.
    Distributor(const Distributor &distributor);

    //! Destructor.
    ~Distributor();

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

    //! \brief Execute a plan specified by the distributor object.
    /*! 
      \param exports [in]
             Contains the values we're exporting.

      \param numPackets [in]
             Specifies the number of values per export/import.

      \param imports [out]
             On entry, buffer must be large enough to accomodate the data exported to us.
             On exit, contains the values exported to us.
    */
    template <class Packet>
    void doPostsAndWaits(const ArrayView<const Packet> &exports,
                         size_t numPackets,
                         const ArrayView<Packet> &imports);

    //! \brief Execute a plan specified by the distributor object.
    /*! 
      \param exports [in]
             Contains the values we're exporting.

      \param numPackets [in]
             Specifies the number of values per export/import.

      \param imports [out]
             On entry, buffer must be large enough to accomodate the data exported to us.
             On exit, contains the values exported to us.
    */
    template <class Packet>
    void doPostsAndWaits(const ArrayView<const Packet> &exports,
                         const ArrayView<size_t> &numExportPacketsPerLID,
                         const ArrayView<Packet> &imports,
                         const ArrayView<size_t> &numImportPacketsPerLID);

    //! \brief Post the data for a distributor plan, but do not execute the waits yet.
    /*! 
      \param exports [in]
             Contains the values to be sent by this node. 

      \param numPackets [in]
             Specifies the number of scalars per export/import.

      \param imports [out]
             Buffer must be large enough to accomodate the data exported to us. 
             The buffer is not guaranteed to be filled until doWaits() is executed.
    */
    template <class Packet>
    void doPosts(const ArrayView<const Packet> &exports,
                 size_t numPackets,
                 const ArrayRCP<Packet> &imports);

    //! \brief Post the data for a distributor plan, but do not execute the waits yet.
    /*! 
      \param exports [in]
             Contains the values to be sent by this node. 

      \param numPackets [in]
             Specifies the number of scalars per export/import.

      \param imports [out]
             Buffer must be large enough to accomodate the data exported to us. 
             The buffer is not guaranteed to be filled until doWaits() is executed.
    */
    template <class Packet>
    void doPosts(const ArrayView<const Packet> &exports,
                 const ArrayView<size_t> &numExportPacketsPerLID,
                 const ArrayRCP<Packet> &imports,
                 const ArrayView<size_t> &numImportPacketsPerLID);

    //! Wait on any outstanding posts to complete.
    void doWaits();

    //! \brief Execute a reverse plan specified by the distributor object.
    /*! 
      \param exports [in]
             Contains the values to be sent by this node.

      \param numPackets [in]
             Specifies the number of scalars per export/import.

      \param imports [out]
             On entry, buffer must be large enough to accomodate the data exported to us.
             On exit, contains the values exported to us.
    */
    template <class Packet>
    void doReversePostsAndWaits(const ArrayView<const Packet> &exports,
                                size_t numPackets,
                                const ArrayView<Packet> &imports);

    //! \brief Execute a reverse plan specified by the distributor object.
    /*! 
      \param exports [in]
             Contains the values to be sent by this node.

      \param numPackets [in]
             Specifies the number of scalars per export/import.

      \param imports [out]
             On entry, buffer must be large enough to accomodate the data exported to us.
             On exit, contains the values exported to us.
    */
    template <class Packet>
    void doReversePostsAndWaits(const ArrayView<const Packet> &exports,
                                const ArrayView<size_t> &numExportPacketsPerLID,
                                const ArrayView<Packet> &imports,
                                const ArrayView<size_t> &numImportPacketsPerLID);

    //! \brief Post the data for a reverse plan, but do not execute the waits yet.
    /*!
      \param exports [in]
             Contains the values we're exporting.

      \param numPackets [in]
             Specifies the number of scalars per export/import.

      \param imports [out]
             Buffer must be large enough to accomodate the data exported to us. 
             The buffer is not guaranteed to be filled until doWaits() is executed.
    */
    template <class Packet>
    void doReversePosts(const ArrayView<const Packet> &exports,
                        size_t numPackets,
                        const ArrayRCP<Packet> &imports);

    //! \brief Post the data for a reverse plan, but do not execute the waits yet.
    /*!
      \param exports [in]
             Contains the values we're exporting.

      \param numPackets [in]
             Specifies the number of scalars per export/import.

      \param imports [out]
             Buffer must be large enough to accomodate the data exported to us. 
             The buffer is not guaranteed to be filled until doWaits() is executed.
    */
    template <class Packet>
    void doReversePosts(const ArrayView<const Packet> &exports,
                        const ArrayView<size_t> &numExportPacketsPerLID,
                        const ArrayRCP<Packet> &imports,
                        const ArrayView<size_t> &numImportPacketsPerLID);

    //! Wait on any outstanding reverse waits to complete.
    void doReverseWaits();

    //@}

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    std::string description() const;

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

    //@}

  private:
    //! The communicator over which to perform distributions.
    RCP<const Comm<int> > comm_;

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

    // imagesTo_, startsTo_ and lengthsTo_ each have size 
    //   numSends_ + selfMessage_
    Array<int> imagesTo_;

    /// \brief Starting index of the block of values to send to each process.
    ///
    /// Given an export buffer that contains all of the data being
    /// sent by this process, the block of values to send to process i
    /// will start at position startsTo_[i].
    ///
    /// This array has length numSends_ + selfMessage_ (that is, it
    /// includes the self message, if there is one).
    Array<size_t> startsTo_;

    /// \brief Length of my process' send to each process.
    ///
    /// lengthsTo_[i] is the length of my process' send to process i.
    /// This array has length numSends_ + selfMessage_ (that is, it
    /// includes the self message, if there is one).
    Array<size_t> lengthsTo_;

    /// \brief The maximum send length to another node.
    ///
    /// maxSendLength_ = max(lengthsTo_[i]) for i != my process ID.
    size_t maxSendLength_;
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

    // totalReceiveLength_ is the total number of Packet received, used to 
    // allocate the receive buffer
    size_t totalReceiveLength_;
    // imagesFrom_, startsFrom_ and lengthsFrom_ each have size 
    //   numReceives_ + selfMessage_
    Array<size_t> lengthsFrom_;
    Array<int> imagesFrom_;
    Array<size_t> startsFrom_;
    Array<size_t> indicesFrom_;

    //! Communication requests associated with nonblocking receives.
    Array<RCP<Teuchos::CommRequest> > requests_;

    /// \brief The reverse distributor.
    ///
    /// This is created on demand in \c getReverse() and cached for
    /// later reuse.  This is why it is declared "mutable".
    mutable RCP<Distributor> reverseDistributor_;

    //! Compute receive info from sends.
    void computeReceives();

    //! Compute send info from receives.
    template <class Ordinal>
    void computeSends (const ArrayView<const Ordinal> &importIDs,
		       const ArrayView<const int> &importNodeIDs,
		       ArrayRCP<Ordinal> &exportIDs,
		       ArrayRCP<int> &exportNodeIDs);

    //! Create a distributor for the reverse communciation pattern.
    void createReverseDistributor() const;

  }; // class Distributor


  template <class Packet>
  void Distributor::doPostsAndWaits(
      const ArrayView<const Packet>& exports,
      size_t numPackets,
      const ArrayView<Packet>& imports) 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(requests_.size() != 0, std::runtime_error,
        Teuchos::typeName(*this) << "::doPostsAndWaits(): Cannot call with outstanding posts.");
    // doPosts takes imports as an ArrayRCP, requiring that the memory location is persisting
    // however, it need only persist until doWaits is called, so it is safe for us to 
    // use a non-persisting reference in this case
    doPosts(exports, numPackets, arcp<Packet>(imports.getRawPtr(),0,imports.size(),false));
    doWaits();
  }

  template <class Packet>
  void Distributor::doPostsAndWaits(
      const ArrayView<const Packet>& exports,
      const ArrayView<size_t> &numExportPacketsPerLID,
      const ArrayView<Packet> &imports,
      const ArrayView<size_t> &numImportPacketsPerLID)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(requests_.size() != 0, std::runtime_error,
        Teuchos::typeName(*this) << "::doPostsAndWaits(): Cannot call with outstanding posts.");
    // doPosts takes imports as an ArrayRCP, requiring that the memory location is persisting
    // however, it need only persist until doWaits is called, so it is safe for us to 
    // use a non-persisting reference in this case
    doPosts(exports, numExportPacketsPerLID, arcp<Packet>(imports.getRawPtr(),0,imports.size(),false), numImportPacketsPerLID);
    doWaits();
  }


  template <class Packet>
  void Distributor::doPosts(const ArrayView<const Packet>& exports,
                            size_t numPackets,
                            const ArrayRCP<Packet>& imports) {
    // start of actual doPosts function
    const int myImageID = comm_->getRank();
    size_t selfReceiveOffset = 0;

#ifdef HAVE_TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(imports.size()) != totalReceiveLength_ * numPackets, std::runtime_error,
        Teuchos::typeName(*this) << "::doPosts(): imports must be large enough to store the imported data.");
#endif

    // allocate space in requests
    //
    // NOTE (mfh 19 Mar 2012): Epetra_MpiDistributor::DoPosts()
    // doesn't (re)allocate its array of requests.  That happens in
    // CreateFromSends(), ComputeRecvs_(), DoReversePosts() (on
    // demand), or Resize_().
    requests_.resize(0);
    requests_.reserve(numReceives_);

    // start up the Irecv's
    {
      size_t curBufferOffset = 0;
      for (size_t i = 0; i < numReceives_ + (selfMessage_ ? 1 : 0); ++i) {
        if (imagesFrom_[i] != myImageID) { 
          // receiving this one from another image
          // setup reference into imports of the appropriate size and at the appropriate place
          ArrayRCP<Packet> impptr = imports.persistingView(curBufferOffset,lengthsFrom_[i]*numPackets);
          requests_.push_back( Teuchos::ireceive<int,Packet>(*comm_,impptr,imagesFrom_[i]) );
        }
        else {
          // receiving this one from myself 
          // note that offset
          selfReceiveOffset = curBufferOffset;
        }
        curBufferOffset += lengthsFrom_[i]*numPackets;
      }
    }

    // wait for everyone else before posting ready-sends below to ensure that 
    // all non-blocking receives above have been posted
    Teuchos::barrier(*comm_);

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

    if (indicesTo_.empty()) { // data is already blocked by processor
      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + imageIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (imagesTo_[p] != myImageID) {
          // sending it to another image
          ArrayView<const Packet> tmpSend(&exports[startsTo_[p]*numPackets],lengthsTo_[p]*numPackets);
          Teuchos::readySend<int,Packet>(*comm_,tmpSend,imagesTo_[p]);
        }
        else {
          // sending it to ourself
          selfNum = p;
        }
      }

      if (selfMessage_) {
        std::copy(exports.begin()+startsTo_[selfNum]*numPackets, exports.begin()+startsTo_[selfNum]*numPackets+lengthsTo_[selfNum]*numPackets, 
                  imports.begin()+selfReceiveOffset);
      }
    }
    else { // data is not blocked by image, use send buffer
      // allocate sendArray buffer
      Array<Packet> sendArray(maxSendLength_*numPackets); 

      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + imageIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (imagesTo_[p] != myImageID) { 
          // sending it to another image
          typename ArrayView<const Packet>::iterator srcBegin, srcEnd;
          size_t sendArrayOffset = 0;
          size_t j = startsTo_[p];
          for (size_t k = 0; k < lengthsTo_[p]; ++k, ++j) {
            srcBegin = exports.begin() + indicesTo_[j]*numPackets;
            srcEnd   = srcBegin + numPackets;
            std::copy( srcBegin, srcEnd, sendArray.begin()+sendArrayOffset );
            sendArrayOffset += numPackets;
          }
          ArrayView<const Packet> tmpSend = sendArray(0,lengthsTo_[p]*numPackets);
          Teuchos::readySend<int,Packet>(*comm_,tmpSend,imagesTo_[p]);
        }
        else { 
          // sending it to myself
          selfNum = p;
          selfIndex = startsTo_[p];
        }
      }

      if (selfMessage_) {
        for (size_t k = 0; k < lengthsTo_[selfNum]; ++k) {
          std::copy( exports.begin()+indicesTo_[selfIndex]*numPackets,
                     exports.begin()+indicesTo_[selfIndex]*numPackets + numPackets,
                     imports.begin() + selfReceiveOffset );
          ++selfIndex;
          selfReceiveOffset += numPackets;
        }
      }
    }
  }

  template <class Packet>
  void Distributor::doPosts(const ArrayView<const Packet>& exports,
                            const ArrayView<size_t>& numExportPacketsPerLID,
                            const ArrayRCP<Packet>& imports,
                            const ArrayView<size_t>& numImportPacketsPerLID) 
  {
    using Teuchos::as;
    using Teuchos::ireceive;

    const int myImageID = comm_->getRank();
    size_t selfReceiveOffset = 0;

#ifdef HAVE_TEUCHOS_DEBUG
    size_t totalNumPackets = 0;
    for (int ii = 0; ii < numImportPacketsPerLID.size(); ++ii) {
      totalNumPackets += numImportPacketsPerLID[ii];
    }
    TEUCHOS_TEST_FOR_EXCEPTION(as<size_t>(imports.size()) != totalNumPackets, 
      std::runtime_error, Teuchos::typeName(*this) << "::doPosts(): The imports "
      "array argument must be large enough to store the imported data.  imports."
      "size() = " << imports.size() << ", but the total number of packets is " 
      << totalNumPackets << ".");
#endif // HAVE_TEUCHOS_DEBUG

    // allocate space in requests
    //
    // NOTE (mfh 19 Mar 2012): Epetra_MpiDistributor::DoPosts()
    // doesn't (re)allocate its array of requests.  That happens in
    // CreateFromSends(), ComputeRecvs_(), DoReversePosts() (on
    // demand), or Resize_().
    requests_.resize(0);
    requests_.reserve(numReceives_);

    // Post the nonblocking receives (Irecv).
    {
      size_t curBufferOffset = 0;
      size_t curLIDoffset = 0;
      for (size_t i = 0; i < numReceives_ + (selfMessage_ ? 1 : 0); ++i) {
        size_t totalPacketsFrom_i = 0;
        for (size_t j=0; j<lengthsFrom_[i]; ++j) {
          totalPacketsFrom_i += numImportPacketsPerLID[curLIDoffset+j];
        }
        curLIDoffset += lengthsFrom_[i];
        if (imagesFrom_[i] != myImageID && totalPacketsFrom_i) { 
	  // If my process is receiving these packet(s) from another
	  // process (not a self-receive), and if there is at least
	  // one packet to receive: 
	  //
	  // 1. Set up the reference (impptr) into the imports array,
	  //    given the offset and size (total number of packets
	  //    from process imagesFrom_[i]).
	  // 2. Start the Irecv and save the resulting request.
          ArrayRCP<Packet> impptr = 
	    imports.persistingView (curBufferOffset, totalPacketsFrom_i);
          requests_.push_back (ireceive<int, Packet> (*comm_, impptr, imagesFrom_[i]));
        }
        else { // Receiving these packet(s) from myself
          selfReceiveOffset = curBufferOffset; // Remember the offset
        }
        curBufferOffset += totalPacketsFrom_i;
      }
    }

    // NOTE (mfh 19 Mar 2012):
    //
    // The ready-sends below require that each ready-send's matching
    // receive (see above) has already been posted.  We ensure this
    // with a barrier.  (Otherwise, some process that doesn't need to
    // post receives might post its ready-send before the receiving
    // process gets to post its receive.)  If you want to remove the
    // barrier, you'll have to replace the ready-sends below with
    // standard sends or Isends.
    //
    // Epetra_MpiDistributor::DoPosts() uses the same approach
    // (Irecvs, barrier, Rsends).
    Teuchos::barrier(*comm_);

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

    if (indicesTo_.empty()) { // data is already laid out according to processor
      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + imageIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (imagesTo_[p] != myImageID && packetsPerSend[p] > 0) {
          // sending it to another image
          ArrayView<const Packet> tmpSend(&exports[sendPacketOffsets[p]],packetsPerSend[p]);
          Teuchos::readySend<int,Packet>(*comm_,tmpSend,imagesTo_[p]);
        }
        else {
          // sending it to ourself
          selfNum = p;
        }
      }

      if (selfMessage_) {
        std::copy(exports.begin()+sendPacketOffsets[selfNum], exports.begin()+sendPacketOffsets[selfNum]+packetsPerSend[selfNum], 
                  imports.begin()+selfReceiveOffset);
      }
    }
    else { // data is not blocked by image, use send buffer
      // allocate sendArray buffer
      Array<Packet> sendArray(maxNumPackets); 
      Array<size_t> indicesOffsets(numExportPacketsPerLID.size(),0);
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
          // sending it to another image
          typename ArrayView<const Packet>::iterator srcBegin, srcEnd;
          size_t sendArrayOffset = 0;
          size_t j = startsTo_[p];
          size_t numPacketsTo_p = 0;
          for (size_t k = 0; k < lengthsTo_[p]; ++k, ++j) {
            srcBegin = exports.begin() + indicesOffsets[j];
            srcEnd   = srcBegin + numExportPacketsPerLID[j];
            numPacketsTo_p += numExportPacketsPerLID[j];
            std::copy( srcBegin, srcEnd, sendArray.begin()+sendArrayOffset );
            sendArrayOffset += numExportPacketsPerLID[j];
          }
          if (numPacketsTo_p > 0) {
            ArrayView<const Packet> tmpSend = sendArray(0,numPacketsTo_p);
            Teuchos::readySend<int,Packet>(*comm_,tmpSend,imagesTo_[p]);
          }
        }
        else { 
          // sending it to myself
          selfNum = p;
          selfIndex = startsTo_[p];
        }
      }

      if (selfMessage_) {
        for (size_t k = 0; k < lengthsTo_[selfNum]; ++k) {
          std::copy( exports.begin()+indicesOffsets[selfIndex],
                     exports.begin()+indicesOffsets[selfIndex]+numExportPacketsPerLID[selfIndex],
                     imports.begin() + selfReceiveOffset );
          selfReceiveOffset += numExportPacketsPerLID[selfIndex];
          ++selfIndex;
        }
      }
    }
  }


  template <class Packet>
  void Distributor::doReversePostsAndWaits(
      const ArrayView<const Packet>& exports,
      size_t numPackets,
      const ArrayView<Packet>& imports) 
  {
    // doPosts takes imports as an ArrayRCP, requiring that the memory location is persisting
    // however, it need only persist until doWaits is called, so it is safe for us to 
    // use a non-persisting reference in this case
    doReversePosts(exports, numPackets, arcp<Packet>(imports.getRawPtr(),0,imports.size(),false));
    doReverseWaits();
  }

  template <class Packet>
  void Distributor::doReversePostsAndWaits(
       const ArrayView<const Packet>& exports,
       const ArrayView<size_t> &numExportPacketsPerLID,
       const ArrayView<Packet> &imports,
       const ArrayView<size_t> &numImportPacketsPerLID)
  {
    // doPosts takes imports as an ArrayRCP, requiring that the memory location is persisting
    // however, it need only persist until doWaits is called, so it is safe for us to 
    // use a non-persisting reference in this case
    doReversePosts(exports, numExportPacketsPerLID, arcp<Packet>(imports.getRawPtr(),0,imports.size(),false),numImportPacketsPerLID);
    doReverseWaits();
  }


  template <class Packet>
  void Distributor::doReversePosts(
      const ArrayView<const Packet>& exports,
      size_t numPackets,
      const ArrayRCP<Packet>& imports) 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(!indicesTo_.empty(),std::runtime_error,
        Teuchos::typeName(*this) << "::doReversePosts(): Can only do reverse comm when original data is blocked by image.");
    if (reverseDistributor_ == null) {
      createReverseDistributor();
    }
    reverseDistributor_->doPosts(exports,numPackets,imports);
  }

  template <class Packet>
  void Distributor::doReversePosts(
      const ArrayView<const Packet>& exports,
      const ArrayView<size_t>& numExportPacketsPerLID,
      const ArrayRCP<Packet>& imports,
      const ArrayView<size_t>& numImportPacketsPerLID) 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(!indicesTo_.empty(),std::runtime_error,
        Teuchos::typeName(*this) << "::doReversePosts(): Can only do reverse comm when original data is blocked by image.");
    if (reverseDistributor_ == null) {
      createReverseDistributor();
    }
    reverseDistributor_->doPosts(exports,numExportPacketsPerLID,imports,numImportPacketsPerLID);
  }


  // note: assumes that size_t >= Ordinal
  //////////////////////////////////////////////////////////////////////////////////////////
  template <class Ordinal>
  void Distributor::computeSends(
      const ArrayView<const Ordinal> & importIDs,
      const ArrayView<const int> & importNodeIDs,
            ArrayRCP<Ordinal> & exportIDs,
            ArrayRCP<int> & exportNodeIDs)
  {
    int myImageID = comm_->getRank();

    size_t numImports = importNodeIDs.size();
    Array<size_t> importObjs(2*numImports);
    for (size_t i = 0; i < numImports; ++i ) {  
      importObjs[2*i]   = Teuchos::as<size_t>(importIDs[i]);
      importObjs[2*i+1] = Teuchos::as<size_t>(myImageID);
    }

    size_t numExports;
    Distributor tempPlan(comm_);
    numExports = tempPlan.createFromSends(importNodeIDs);
    if (numExports > 0) {
      exportIDs = arcp<Ordinal>(numExports);
      exportNodeIDs = arcp<int>(numExports);
    }

    Array<size_t> exportObjs(tempPlan.getTotalReceiveLength()*2);
    tempPlan.doPostsAndWaits<size_t>(importObjs(),2,exportObjs());

    for (size_t i = 0; i < numExports; ++i) {
      exportIDs[i]     = Teuchos::as<Ordinal>(exportObjs[2*i]);
      exportNodeIDs[i] = exportObjs[2*i+1];
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////
  template <class Ordinal>
  void Distributor::createFromRecvs(
      const ArrayView<const Ordinal> &remoteIDs, 
      const ArrayView<const int> &remoteImageIDs, 
            ArrayRCP<Ordinal> &exportGIDs, 
            ArrayRCP<int> &exportNodeIDs)
  {
    using Teuchos::outArg;
    {
      const int myImageID = comm_->getRank();
      int err_node = (remoteIDs.size() != remoteImageIDs.size()) ? myImageID : -1;
      int gbl_err;
      Teuchos::reduceAll(*comm_,Teuchos::REDUCE_MAX,err_node,outArg(gbl_err));
      TEUCHOS_TEST_FOR_EXCEPTION(gbl_err != -1, std::runtime_error,
          Teuchos::typeName(*this) 
          << "::createFromRecvs(): lists of remote IDs and remote node IDs must have the same size (error on node " 
          << gbl_err << ").");
    }
    computeSends(remoteIDs, remoteImageIDs, exportGIDs, exportNodeIDs);
    (void)createFromSends(exportNodeIDs());
  }

} // namespace Tpetra

#endif // TPETRA_DISTRIBUTOR_HPP
