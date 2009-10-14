// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_DISTRIBUTOR_HPP
#define TPETRA_DISTRIBUTOR_HPP

#include "Tpetra_Util.hpp"
#include <Teuchos_RCP.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>


// FINISH: some of the get accessors may not be necessary anymore. clean up.
// FINISH: This class may not be const correct. doPosts() et al. perhaps should be const, with affected members made mutable.

namespace Tpetra {

  //! \brief The Tpetra gather/scatter setup class.
  /*! The Distributor class is an interface that encapsulates the general
        information and services needed for other Tpetra classes to perform gather/scatter
        operations on a parallel computer.
  */

  class Distributor : public Teuchos::Describable {
  public:

    //! @name Constructor/Destructor
    //@{ 

    //! Construct the Distributor using the specified communicator.
    explicit Distributor(const Teuchos::RCP<const Teuchos::Comm<int> > & comm);

    //! Copy Constructor
    Distributor(const Distributor &distributor);

    //! Destructor.
    ~Distributor();

    //@}


    //! @name Gather/Scatter Constructors
    //@{ 

    //! \brief Create a Distributor object using list of node IDs to send to
    /*! Take a list of node IDs and construct a plan for efficiently scattering to those nodes.
        Return the number of IDs being sent to me.

      \param exportNodeIDs [in]
             List of nodes that will get the exported data. 
             A node ID greater than or equal to the number of nodes will 
             result in a \c std::runtime_error on all nodes.
             Node IDs less than zero
             are ignored; their placement corresponds to null sends in any
             future exports. That is, if <tt>exportNodeIDs[0] == -1</tt>, then 
             the corresponding position in the export array is ignored during a call to
             doPosts() or doPostsAndWaits() is skipped.
             For this reason, a negative entry is sufficient to break contiguity.

      \return Number of imports this node will be receiving.

    */
    size_t createFromSends(const Teuchos::ArrayView<const int> &exportNodeIDs);

    //! \brief Create Distributor object using list of node IDs to receive from
    /*! Take a list of node IDs and construct a plan for efficiently scattering to those nodes.
        Return the number and list of IDs being sent by me.

      \param remoteIDs [in]
             List of remote IDs wanted. 

      \param remoteNodeIDs [in]
             List of nodes that will send the corresponding remote IDs. Node IDs less than zero
             are ignored; their placement corresponds to null sends in any
             future exports. A node ID greater than or equal to the number of nodes will 
             result in a \c std::runtime_error on all nodes.

      \param exportIDs [out]
             List of IDs that need to be sent from this node.

      \param exportNodeIDs [out]
             List of nodes that will get the exported IDs.

      \note \c exportGIDs and \c exportNodeIDs are allocated by the Distributor, but they are reference counted and will be automatically deallocated.
    */
    template <class Ordinal>
    void createFromRecvs(const Teuchos::ArrayView<const Ordinal> &remoteIDs, 
                         const Teuchos::ArrayView<const int> &remoteNodeIDs, 
                               Teuchos::ArrayRCP<Ordinal> &exportIDs, 
                               Teuchos::ArrayRCP<int> &exportNodeIDs);

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
    Teuchos::ArrayView<const int> getImagesFrom() const;

    //! A list of images to which this node is sending values. (non-persisting view)
    Teuchos::ArrayView<const int> getImagesTo() const;

    //! Number of values we're receiving from each node. (non-persisting view)
    /*! We will receive <tt>getLengthsFrom[i]</tt> values from node <tt>getImagesFrom[i]</tt>. */
    Teuchos::ArrayView<const size_t> getLengthsFrom() const;

    //! Number of values we're sending to each node. (non-persisting view)
    /*! We will send <tt>getLengthsTo[i]</tt> values to image <tt>getImagesTo[i]</tt>. */
    Teuchos::ArrayView<const size_t> getLengthsTo() const;

    //@}

    //! @name Reverse Communication Methods
    //@{ 

    //! \brief Returns a Distributor with a reverse plan of this Distributor's plan
    /*! This method creates the reverse Distributor the first time the function
        is called.
    */
    const Teuchos::RCP<Distributor> & getReverse() const;

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
    void doPostsAndWaits(const Teuchos::ArrayView<const Packet> &exports,
                         size_t numPackets,
                         const Teuchos::ArrayView<Packet> &imports);

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
    void doPosts(const Teuchos::ArrayView<const Packet> &exports,
                 size_t numPackets,
                 const Teuchos::ArrayRCP<Packet> &imports);

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
    void doReversePostsAndWaits(const Teuchos::ArrayView<const Packet> &exports,
                                size_t numPackets,
                                const Teuchos::ArrayView<Packet> &imports);

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
    void doReversePosts(const Teuchos::ArrayView<const Packet> &exports,
                        size_t numPackets,
                        const Teuchos::ArrayRCP<Packet> &imports);

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

    // private data members
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;

    size_t numExports_;
    // selfMessage_ is whether I have a send for myself
    bool selfMessage_;
    // numSends_ is number of sends to other nodes; is less than or equal to the number of nodes
    size_t numSends_;
    // imagesTo_, startsTo_ and lengthsTo_ each have size 
    //   numSends_ + selfMessage_
    Teuchos::Array<int> imagesTo_;
    /* Given an export buffer that contains all of the item being sent by this node,
       the block of values for node i will start at position startsTo_[i]  */
    Teuchos::Array<size_t> startsTo_;
    Teuchos::Array<size_t> lengthsTo_;
    // maxSendLength_ is the maximum send to another node: 
    //   max(lengthsTo_[i]) for i != me
    size_t maxSendLength_;
    Teuchos::Array<size_t> indicesTo_;
    // numReceives_ is the number of receives by me from other procs, not
    // counting self receives
    size_t numReceives_;
    // totalReceiveLength_ is the total number of Packet received, used to 
    // allocate the receive buffer
    size_t totalReceiveLength_;
    // imagesFrom_, startsFrom_ and lengthsFrom_ each have size 
    //   numReceives_ + selfMessage_
    Teuchos::Array<size_t> lengthsFrom_;
    Teuchos::Array<int> imagesFrom_;
    Teuchos::Array<size_t> startsFrom_;
    Teuchos::Array<size_t> indicesFrom_;

    // requests associated with non-blocking receives
    Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest> > requests_;

    mutable Teuchos::RCP<Distributor> reverseDistributor_;

    // compute receive info from sends
    void computeReceives();

    // compute send info from receives
    template <class Ordinal>
    void computeSends(const Teuchos::ArrayView<const Ordinal> &importIDs,
                      const Teuchos::ArrayView<const int> &importNodeIDs,
                            Teuchos::ArrayRCP<Ordinal> &exportIDs,
                            Teuchos::ArrayRCP<int> &exportNodeIDs);

    // create a distributor for the reverse communciation pattern (pretty much
    // swap all send and receive info)
    void createReverseDistributor() const;

  }; // class Distributor


  template <class Packet>
  void Distributor::doPostsAndWaits(
      const Teuchos::ArrayView<const Packet>& exports,
      size_t numPackets,
      const Teuchos::ArrayView<Packet>& imports) 
  {
    TEST_FOR_EXCEPTION(requests_.size() != 0, std::runtime_error,
        Teuchos::typeName(*this) << "::doPostsAndWaits(): Cannot call with outstanding posts.");
    // doPosts takes imports as an ArrayRCP, requiring that the memory location is persisting
    // however, it need only persist until doWaits is called, so it is safe for us to 
    // use a non-persisting reference in this case
    doPosts(exports, numPackets, Teuchos::arcp<Packet>(imports.getRawPtr(),0,imports.size(),false));
    doWaits();
  }


  template <class Packet>
  void Distributor::doPosts(const Teuchos::ArrayView<const Packet>& exports,
                            size_t numPackets,
                            const Teuchos::ArrayRCP<Packet>& imports) {
    using Teuchos::ArrayRCP;
    // start of actual doPosts function
    const int myImageID = comm_->getRank();
    size_t selfReceiveOffset = 0;

#ifdef HAVE_TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(Teuchos::as<size_t>(imports.size()) != totalReceiveLength_ * numPackets, std::runtime_error,
        Teuchos::typeName(*this) << "::doPosts(): imports must be large enough to store the imported data.");
#endif

    // allocate space in requests
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
          Teuchos::ArrayView<const Packet> tmpSend(&exports[startsTo_[p]*numPackets],lengthsTo_[p]*numPackets);
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
      Teuchos::Array<Packet> sendArray(maxSendLength_*numPackets); 

      for (size_t i = 0; i < numBlocks; ++i) {
        size_t p = i + imageIndex;
        if (p > (numBlocks - 1)) {
          p -= numBlocks;
        }

        if (imagesTo_[p] != myImageID) { 
          // sending it to another image
          typename Teuchos::ArrayView<const Packet>::iterator srcBegin, srcEnd;
          size_t sendArrayOffset = 0;
          size_t j = startsTo_[p];
          for (size_t k = 0; k < lengthsTo_[p]; ++k, ++j) {
            srcBegin = exports.begin() + indicesTo_[j]*numPackets;
            srcEnd   = srcBegin + numPackets;
            std::copy( srcBegin, srcEnd, sendArray.begin()+sendArrayOffset );
            sendArrayOffset += numPackets;
          }
          Teuchos::ArrayView<const Packet> tmpSend = sendArray(0,lengthsTo_[p]*numPackets);
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
  void Distributor::doReversePostsAndWaits(
      const Teuchos::ArrayView<const Packet>& exports,
      size_t numPackets,
      const Teuchos::ArrayView<Packet>& imports) 
  {
    // doPosts takes imports as an ArrayRCP, requiring that the memory location is persisting
    // however, it need only persist until doWaits is called, so it is safe for us to 
    // use a non-persisting reference in this case
    doReversePosts(exports, numPackets, Teuchos::arcp<Packet>(imports.getRawPtr(),0,imports.size(),false));
    doReverseWaits();
  }


  template <class Packet>
  void Distributor::doReversePosts(
      const Teuchos::ArrayView<const Packet>& exports,
      size_t numPackets,
      const Teuchos::ArrayRCP<Packet>& imports) 
  {
    TEST_FOR_EXCEPTION(!indicesTo_.empty(),std::runtime_error,
        Teuchos::typeName(*this) << "::doReversePosts(): Can only do reverse comm when original data is blocked by image.");
    if (reverseDistributor_ == Teuchos::null) {
      createReverseDistributor();
    }
    reverseDistributor_->doPosts(exports,numPackets,imports);
  }


  // note: assumes that size_t >= Ordinal
  //////////////////////////////////////////////////////////////////////////////////////////
  template <class Ordinal>
  void Distributor::computeSends(
      const Teuchos::ArrayView<const Ordinal> & importIDs,
      const Teuchos::ArrayView<const int> & importNodeIDs,
            Teuchos::ArrayRCP<Ordinal>& exportIDs,
            Teuchos::ArrayRCP<int>& exportNodeIDs)
  {
    int myImageID = comm_->getRank();

    size_t numImports = importNodeIDs.size();
    Teuchos::Array<size_t> importObjs(2*numImports);
    for (size_t i = 0; i < numImports; ++i ) {  
      importObjs[2*i]   = Teuchos::as<size_t>(importIDs[i]);
      importObjs[2*i+1] = Teuchos::as<size_t>(myImageID);
    }

    size_t numExports;
    Distributor tempPlan(comm_);
    numExports = tempPlan.createFromSends(importNodeIDs);
    if (numExports > 0) {
      exportIDs = Teuchos::arcp<Ordinal>(numExports);
      exportNodeIDs = Teuchos::arcp<int>(numExports);
    }

    Teuchos::Array<size_t> exportObjs(tempPlan.getTotalReceiveLength()*2);
    tempPlan.doPostsAndWaits<size_t>(importObjs(),2,exportObjs());

    for (size_t i = 0; i < numExports; ++i) {
      exportIDs[i]     = Teuchos::as<Ordinal>(exportObjs[2*i]);
      exportNodeIDs[i] = exportObjs[2*i+1];
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////
  template <class Ordinal>
  void Distributor::createFromRecvs(
      const Teuchos::ArrayView<const Ordinal> &remoteIDs, 
      const Teuchos::ArrayView<const int> &remoteImageIDs, 
            Teuchos::ArrayRCP<Ordinal> &exportGIDs, 
            Teuchos::ArrayRCP<int> &exportNodeIDs)
  {
    {
      const int myImageID = comm_->getRank();
      int err_node = (remoteIDs.size() != remoteImageIDs.size()) ? myImageID : -1;
      int gbl_err;
      Teuchos::reduceAll(*comm_,Teuchos::REDUCE_MAX,err_node,&gbl_err);
      TEST_FOR_EXCEPTION(gbl_err != -1, std::runtime_error,
          Teuchos::typeName(*this) 
          << "::createFromRecvs(): lists of remote IDs and remote node IDs must have the same size (error on node " 
          << gbl_err << ").");
    }
    computeSends(remoteIDs, remoteImageIDs, exportGIDs, exportNodeIDs);
    (void)createFromSends(exportNodeIDs());
  }

} // namespace Tpetra

#endif // TPETRA_DISTRIBUTOR_HPP
