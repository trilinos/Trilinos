// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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

#include <Teuchos_RCP.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include "Tpetra_Util.hpp"
#include <Teuchos_Object.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

namespace Tpetra {
  
  //! Tpetra::Distributor:  The Tpetra Gather/Scatter Setup Class.
  /*! The Distributor class is an interface that encapsulates the general
        information and services needed for other Tpetra classes to perform gather/scatter
        operations on a parallel computer.
  */

  template<typename Ordinal>
  class Distributor : public Teuchos::Object {
  public:

    //@{ \name Constructor/Destructor

    //! Comm Constuctor (default ctr)
    Distributor(const Teuchos::RCP< Teuchos::Comm<Ordinal> > & comm);

    //! Copy Constructor
    Distributor(const Distributor<Ordinal> & distributor);

    //! Destructor.
    ~Distributor();

    //@}


    //@{ \name Gather/Scatter Constructors

    //! Create Distributor object using list of ImageIDs to send to
    /*! Take a list of ImageIDs and construct a plan for efficiently scattering to these images.
        Return the number of IDs being sent to me.
      \param exportImageIDs In
             List of images that will get the exported IDs.
      \param numRemoteIDs Out
             Number of IDs this image will be receiving.
    */
    void createFromSends(const std::vector<Ordinal> & exportImageIDs,
                         Ordinal & numRemoteIDs);

    //! Create Distributor object using list of Image IDs to receive from
    /*! Take a list of global IDs and construct a plan for efficiently scattering to these images.
        Return the number and list of IDs being sent by me.
      \param remoteGIDs In
             List of IDs that this image wants.
      \param remoteImageIDs In
             List of images that will send the remote IDs.
      \param exportGIDs Out
             List of IDs that need to be sent from this image.
      \param exportImageIDs Out
             List of images that will get the exported IDs.
    */
    void createFromRecvs(const std::vector<Ordinal> & remoteGIDs, 
                         const std::vector<Ordinal> & remoteImageIDs, 
                         std::vector<Ordinal>& exportGIDs, 
                         std::vector<Ordinal>& exportImageIDs);

    //@}

    //@{ \name Attribute Accessor Methods

    //! getTotalReceiveLength
    const Ordinal & getTotalReceiveLength() const;

    //! getNumReceives
    const Ordinal & getNumReceives() const;

    //! getSelfMessage - flag for if we're sending to ourself
    /*! If we are sending any elements to ourself, returns true. If we aren't, returns false. */
    bool getSelfMessage() const;

    //! getNumSends
    const Ordinal & getNumSends() const;

    //! getMaxSendLength - maximum number of elements we're sending to a remote image
    const Ordinal & getMaxSendLength() const;

    //! getImagesFrom - list of images sending elements to us
    const std::vector<Ordinal> & getImagesFrom() const;

    //! getLengthsFrom - number of elements we're receiving from each image
    /*! We will receive lengthsFrom[i] elements from image imagesFrom[i] */
    const std::vector<Ordinal> & getLengthsFrom() const;

    //! getImagesTo - list of images we're sending elements to
    const std::vector<Ordinal> & getImagesTo() const;

    //! getIndicesTo
    /*! (Used only if exportImageIDs was not blocked by image.)
        Gives the order to the export buffer, in order to get
      a version that is sorted by imageID. */
    const std::vector<Ordinal> & getIndicesTo() const;

    //! getStartsTo - list of offsets into export buffer
    /*! Given an export buffer that contains all of the elements we're sending out, 
        image i's block of elements will start at position startsTo[i] */
    const std::vector<Ordinal> & getStartsTo() const;

    //! getLengthsTo - number of elements we're sending to each image
    /*! We will send lengthsTo[i] elements to image imagesTo[i] */
    const std::vector<Ordinal> & getLengthsTo() const;

    //@}

    //@{ \name Reverse Communication Methods

    // getReverse
    //! Returns a Distributor with a reverse plan of this Distributor's plan
    /*! Creates the reverse Distributor if this is the first time this function
        has been called.
    */
    const Distributor<Ordinal> & getReverse() const;

    //@}

    //@{ \name Execute Distributor Plan Methods

    //! doPostsAndWaits
    /*! Execute a plan specified by the distributor object.
      \param exports In
             On entry, contains the values we're exporting.
      \param imports Out
             On exit, contains the values exported to us. (\c imports will be resized
             if necessary, and any existing values will be overwritten.)
    */
    template <typename Packet>
    void doPostsAndWaits(const std::vector<Packet>& exports,
                         const Ordinal numPackets,
                               std::vector<Packet>& imports);

    //! doPosts
    template <typename Packet>
    void doPosts(const std::vector<Packet>& exports,
                 const Ordinal numPackets,
                       std::vector<Packet>& imports);

    //! doWaits
    void doWaits();

    //! doReversePostsAndWaits
    /*! Execute a reverse plan specified by the distributor object.
      \param exports In
             On entry, contains the values we're exporting.
      \param imports Out
             On exit, contains the values exported to us. (imports will be resized
             if necessary, and any existing values will be overwritten.)
    */
    template <typename Packet>
    void doReversePostsAndWaits(const std::vector<Packet>& exports,
                                const Ordinal numPackets,
                                      std::vector<Packet>& imports);

    //! doReversePosts
    template <typename Packet>
    void doReversePosts(const std::vector<Packet>& exports,
                        const Ordinal numPackets,
                              std::vector<Packet>& imports);
    
    //! doReverseWaits
    void doReverseWaits();
    
    //@}

    //@{ \name I/O Methods

    //! Implements Teuchos::Object::print.
    void print(ostream& os) const;

    //@}

  private:

    // private data members
    Teuchos::RCP< Teuchos::Comm<Ordinal> > comm_;

    Ordinal numExports_;
    Ordinal selfMessage_;
    Ordinal numSends_;
    std::vector<Ordinal> imagesTo_;
    std::vector<Ordinal> startsTo_;
    std::vector<Ordinal> lengthsTo_;
    Ordinal maxSendLength_;
    std::vector<Ordinal> indicesTo_;
    Ordinal numReceives_;
    Ordinal totalReceiveLength_;
    std::vector<Ordinal> lengthsFrom_;
    std::vector<Ordinal> imagesFrom_;
    std::vector<Ordinal> indicesFrom_;
    std::vector<Ordinal> startsFrom_;

    mutable Teuchos::RCP< Distributor<Ordinal> > reverseDistributor_;

    void computeReceives();

    void computeSends(const std::vector<Ordinal> & importIDs,
                      const std::vector<Ordinal> & importImageIDs,
                      std::vector<Ordinal>& exportIDs,
                      std::vector<Ordinal>& exportImageIDs);

    void createReverseDistributor() const;

    // requests
    std::vector<Teuchos::RCP<Teuchos::CommRequest> > requests_;

  }; // class Distributor


  template <typename Ordinal>
  Distributor<Ordinal>::Distributor(Teuchos::RCP< Teuchos::Comm<Ordinal> > const& comm) 
    : Teuchos::Object("Tpetra::Distributor")
    , comm_(comm)
    , numExports_(Teuchos::OrdinalTraits<Ordinal>::zero())
    , selfMessage_(false)
    , numSends_(Teuchos::OrdinalTraits<Ordinal>::zero())
    , maxSendLength_(Teuchos::OrdinalTraits<Ordinal>::zero())
    , numReceives_(Teuchos::OrdinalTraits<Ordinal>::zero())
    , totalReceiveLength_(Teuchos::OrdinalTraits<Ordinal>::zero())
  {}

  template <typename Ordinal>
  Distributor<Ordinal>::Distributor(Distributor<Ordinal> const& distributor) 
    : Teuchos::Object(distributor.label())
    , comm_(distributor.comm_)
    , numExports_(distributor.numExports_)
    , selfMessage_(distributor.selfMessage_)
    , numSends_(distributor.numSends_)
    , maxSendLength_(distributor.maxSendLength_)
    , numReceives_(distributor.numReceives_)
    , totalReceiveLength_(distributor.totalReceiveLength_)
    , reverseDistributor_(distributor.reverseDistributor_)
  {}

  template <typename Ordinal>
  Distributor<Ordinal>::~Distributor() 
  {}

  template <typename Ordinal>
  void Distributor<Ordinal>::createFromSends(
      const std::vector<Ordinal> & exportImageIDs,
      Ordinal & numRemoteIDs) 
  {
    const Ordinal zero = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal one  = Teuchos::OrdinalTraits<Ordinal>::one();

    numExports_ = exportImageIDs.size();

    const int myImageID = comm_->getRank();
    const int numImages = comm_->getSize();

    // Check to see if items are grouped by images without gaps
    // If so, indices_to -> 0

    // Setup data structures for quick traversal of arrays
    std::vector<Ordinal> starts(numImages + 1, zero);

    Ordinal numActive = zero;
    bool noSendBuff = true;

    for (int i = 0; i < numExports_; ++i) {
      if (noSendBuff && (i > 0) && (exportImageIDs[i] < exportImageIDs[i-1])) {
        noSendBuff = false;
      }
      if (exportImageIDs[i] >= zero) {
        ++starts[exportImageIDs[i]];
        ++numActive;
      }
    }

    if (starts[myImageID] != zero) {
      selfMessage_ = one;
    }
    else {
      selfMessage_ = zero;
    }

    numSends_ = zero;

    if (noSendBuff) { // grouped by image, no send buffer or indicesTo_ needed
      for (int i=0; i < numImages; ++i) {
        if (starts[i]) ++numSends_;
      }

      imagesTo_.resize(numSends_); // can we change these to reserves?
      startsTo_.resize(numSends_);
      lengthsTo_.resize(numSends_);

      for (Ordinal i = zero, index = zero; i < numSends_; ++i) {
        startsTo_[i] = index;
        Ordinal imageID = exportImageIDs[index];
        imagesTo_[i] = imageID;
        index += starts[imageID];
      }

      if (numSends_ > zero) {
        sortArrays(imagesTo_, startsTo_);
      }

      maxSendLength_ = zero;

      for (int i = 0; i < numSends_; ++i) {
        Ordinal imageID = imagesTo_[i];
        lengthsTo_[i] = starts[imageID];
        if ((imageID != myImageID) && (lengthsTo_[i] > maxSendLength_)) {
          maxSendLength_ = lengthsTo_[i];
        }
      }
    }
    else { // not grouped by image, need send buffer and indicesTo_
      if (starts.front() != zero ) {
        numSends_ = one;
      }

      for (int i = 0; i < numImages; i++) {
        if (starts[i] != zero) {
          ++numSends_;
        }
        starts[i] += starts[i-one];
      }

      for (int i = numImages-1; i != 0; i--) {
        starts[i] = starts[i-one];
      }

      starts.front() = zero;

      indicesTo_.resize(numActive);

      for (Ordinal i = zero; i < numExports_; i++) {
        if (exportImageIDs[i] >= zero) {
          indicesTo_[starts[exportImageIDs[i]]] = i;
          ++starts[exportImageIDs[i]];
        }
      }

      // Reconstuct starts array to index into indicesTo.

      for (int i = numImages-1; i != 0; i--) {
        starts[i] = starts[i-1];
      }
      starts.front() = zero;       
      starts[numImages] = numActive;

      imagesTo_.resize(numSends_); // can we change these to reserves?
      startsTo_.resize(numSends_);
      lengthsTo_.resize(numSends_);

      maxSendLength_ = zero;

      for (int i = 0, j = 0; i < numImages; i++ ) {
        if (starts[i+1] != starts[i]) {
          lengthsTo_[j] = starts[i+1] - starts[i];
          startsTo_[j] = starts[i];
          if ((i != myImageID) && (lengthsTo_[j] > maxSendLength_)) {
            maxSendLength_ = lengthsTo_[j];
          }
          imagesTo_[j] = i;
          j++;
        }
      }
    }

    numSends_ -= selfMessage_;

    // Invert map to see what msgs are received and what length
    computeReceives();

    numRemoteIDs = totalReceiveLength_;
  }

  template <typename Ordinal>
  void Distributor<Ordinal>::createFromRecvs(
      const std::vector<Ordinal> & remoteGIDs, 
      const std::vector<Ordinal> & remoteImageIDs, 
      std::vector<Ordinal> & exportGIDs, 
      std::vector<Ordinal> & exportImageIDs)
  {
    computeSends(remoteGIDs, remoteImageIDs, exportGIDs, exportImageIDs);
    Ordinal testNumRemoteIDs; // dummy
    createFromSends(exportImageIDs, testNumRemoteIDs);
  }

  template <typename Ordinal>
  const Ordinal & Distributor<Ordinal>::getTotalReceiveLength() const 
  { return(totalReceiveLength_); }

  template <typename Ordinal>
  const Ordinal & Distributor<Ordinal>::getNumReceives() const 
  { return(numReceives_); }

  template <typename Ordinal>
  bool Distributor<Ordinal>::getSelfMessage() const 
  { return(selfMessage_); }

  template <typename Ordinal>
  const Ordinal & Distributor<Ordinal>::getNumSends() const 
  { return(numSends_); }

  template <typename Ordinal>
  const Ordinal & Distributor<Ordinal>::getMaxSendLength() const 
  { return(maxSendLength_); }

  template <typename Ordinal>
  const std::vector<Ordinal> & Distributor<Ordinal>::getImagesFrom() const 
  { return(imagesFrom_); }

  template <typename Ordinal>
  const std::vector<Ordinal> & Distributor<Ordinal>::getLengthsFrom() const 
  { return(lengthsFrom_); }

  template <typename Ordinal>
  const std::vector<Ordinal> & Distributor<Ordinal>::getImagesTo() const 
  { return(imagesTo_); }

  template <typename Ordinal>
  const std::vector<Ordinal> & Distributor<Ordinal>::getIndicesTo() const 
  { return(indicesTo_); }

  template <typename Ordinal>
  const std::vector<Ordinal> & Distributor<Ordinal>::getStartsTo() const 
  { return(startsTo_); }

  template <typename Ordinal>
  const std::vector<Ordinal> & Distributor<Ordinal>::getLengthsTo() const 
  { return(lengthsTo_); }

  template <typename Ordinal>
  const Distributor<Ordinal> & Distributor<Ordinal>::getReverse() const
  {
    if (reverseDistributor_ == Teuchos::null) { 
      // need to create reverse distributor
      createReverseDistributor();
    }
    return(*reverseDistributor_);
  }

  template <typename Ordinal>
  void Distributor<Ordinal>::createReverseDistributor() const {
    const Ordinal zero = Teuchos::OrdinalTraits<Ordinal>::zero();

    reverseDistributor_ = Teuchos::rcp(new Distributor<Ordinal>(comm_));

    // compute new totalSendLength
    Ordinal totalSendLength = zero;
    for (Ordinal i = zero; i < (numSends_ + selfMessage_); i++) {
      totalSendLength += lengthsTo_[i];
    }

    // compute new maxReceiveLength
    Ordinal maxReceiveLength = zero;
    const int myImageID = comm_->getRank();
    for (Ordinal i = zero; i < numReceives_; i++) {
      if (imagesFrom_[i] != myImageID) {
        if (lengthsFrom_[i] > maxReceiveLength) {
          maxReceiveLength = lengthsFrom_[i];
        }
      }
    }

    // initialize all of reverseDistributor's data members
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
    // Note: numExports_ was not copied
  }

  template <typename Ordinal>
  template <typename Packet>
  void Distributor<Ordinal>::doPostsAndWaits(
      const std::vector<Packet>& exports,
      const Ordinal numPackets,
            std::vector<Packet>& imports) 
  {
    doPosts(exports, numPackets, imports);
    doWaits();
  }

  template <typename Ordinal>
  template <typename Packet>
  void Distributor<Ordinal>::doPosts(
      const std::vector<Packet>& exports,
      const Ordinal numPackets,
            std::vector<Packet>& imports) 
  {
    using Teuchos::ArrayRCP;

    // start of actual doPosts function
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE  = Teuchos::OrdinalTraits<Ordinal>::one();
    const Ordinal myImageID = comm_->getRank();
    int selfReceiveOffset = 0;

    // FINISH: verify that this resize() doesn't break something
    imports.resize(totalReceiveLength_ * numPackets);

    // allocate space in requests
    requests_.reserve(numReceives_);
    requests_.resize(0);

    // start up the Irecv's
    {
      int curBufferOffset = 0;
      for (int i = 0; i < (numReceives_ + selfMessage_); ++i) {
        if (imagesFrom_[i] != myImageID) { 
          // receiving this one from another image
          // setup reference into imports of the appropriate size and at the appropriate place
          ArrayRCP<Packet> impptr = Teuchos::arcp(&imports[curBufferOffset],0,lengthsFrom_[i]*numPackets,false);
          requests_.push_back( Teuchos::ireceive<Ordinal,Packet>(*comm_,impptr,imagesFrom_[i]) );
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
    Ordinal numBlocks = numSends_+ selfMessage_;
    Ordinal imageIndex = ZERO;
    while ((imageIndex < numBlocks) && (imagesTo_[imageIndex] < myImageID)) {
      imageIndex++;
    }
    if (imageIndex == numBlocks) {
      imageIndex = ZERO;
    }

    Ordinal selfNum = ZERO;
    Ordinal selfIndex = ZERO;

    if (indicesTo_.empty()) { // data is already blocked by processor
      for (Ordinal i = ZERO; i < numBlocks; ++i) {
        Ordinal p = i + imageIndex;
        if (p > (numBlocks - ONE)) {
          p -= numBlocks;
        }

        if (imagesTo_[p] != myImageID) {
          // sending it to another image
          Teuchos::ArrayView<const Packet> tmpSend(&exports[startsTo_[p]*numPackets],lengthsTo_[p]*numPackets);
          Teuchos::readySend<Ordinal,Packet>(*comm_,tmpSend,imagesTo_[p]);
        }
        else {
          // sending it to ourself
          selfNum = p;
        }
      }

      if (selfMessage_ > ZERO) {
        std::copy(exports.begin()+startsTo_[selfNum]*numPackets, exports.begin()+startsTo_[selfNum]*numPackets+lengthsTo_[selfNum]*numPackets, 
                  imports.begin()+selfReceiveOffset);
      }
    }
    else { // data is not blocked by image, use send buffer
      // allocate sendArray buffer
      std::vector<Packet> sendArray(maxSendLength_*numPackets); 

      for (Ordinal i = ZERO; i < numBlocks; i++) {
        Ordinal p = i + imageIndex;
        if (p > (numBlocks - ONE)) {
          p -= numBlocks;
        }

        if (imagesTo_[p] != myImageID) { 
          // sending it to another image
          typename std::vector<Packet>::const_iterator srcBegin, srcEnd;
          int sendArrayOffset = 0;
          int j = startsTo_[p];
          for (Ordinal k = ZERO; k < lengthsTo_[p]; ++k, ++j) {
            srcBegin = exports.begin() + indicesTo_[j]*numPackets;
            srcEnd   = srcBegin + numPackets;
            std::copy( srcBegin, srcEnd, sendArray.begin()+sendArrayOffset );
            sendArrayOffset += numPackets;
          }
          Teuchos::ArrayView<const Packet> tmpSend(&sendArray[0],lengthsTo_[p]*numPackets);
          Teuchos::readySend<Ordinal,Packet>(*comm_,tmpSend,imagesTo_[p]);
        }
        else { 
          // sending it to myself
          selfNum = p;
          selfIndex = startsTo_[p];
        }
      }

      // FINISH: revisit this to make sure it's not stupidly inefficient
      if (selfMessage_ > ZERO) {
        for (Ordinal k = ZERO; k < lengthsTo_[selfNum]; ++k) {
          std::copy( exports.begin()+indicesTo_[selfIndex]*numPackets,
                     exports.begin()+indicesTo_[selfIndex]*numPackets + numPackets,
                     imports.begin() + selfReceiveOffset );
          selfIndex++;
          selfReceiveOffset += numPackets;
        }
      }
    }
  }

  template <typename Ordinal>
  void Distributor<Ordinal>::doWaits() 
  {
    const Ordinal & numReceives = getNumReceives();
    if (numReceives > Teuchos::OrdinalTraits<Ordinal>::zero()) {
      Teuchos::waitAll(*comm_,arrayViewFromVector(requests_));
      // Requests should all be null, clear them
#ifdef TEUCHOS_DEBUG
      for (typename requests_::const_iterator i=requests_.begin(); i != requests_.end(); ++i) {
        TEST_FOR_EXCEPTION(*i != Teuchos::null, std::runtime_error,
            "Tpetra::Distributor<"<<Teuchos::OrdinalTraits<Ordinal>::name()<<">::doWaits(): Requests should be null after call to Teuchos::waitAll().");
      }
#endif
      requests_.clear();
    }
  }

  template <typename Ordinal>
  template <typename Packet>
  void Distributor<Ordinal>::doReversePostsAndWaits(
      const std::vector<Packet>& exports,
      const Ordinal numPackets,
            std::vector<Packet>& imports) 
  {
    doReversePosts(exports, numPackets, imports);
    doReverseWaits();
  }

  template <typename Ordinal>
  template <typename Packet>
  void Distributor<Ordinal>::doReversePosts(
      const std::vector<Packet>& exports,
      const Ordinal numPackets,
            std::vector<Packet>& imports) 
  {
    // FINISH: what does this message mean?
    TEST_FOR_EXCEPTION(getIndicesTo().empty(),std::runtime_error,
        "Tpetra::Distributor<"<<Teuchos::OrdinalTraits<Ordinal>::name()<<">::doReversePosts(): Can only do reverse comm when original data is blocked by image.");
    if (reverseDistributor_ == Teuchos::null) {
      createReverseDistributor();
    }
    reverseDistributor_->doPosts(exports,numPackets,imports);
  }

  template <typename Ordinal>
  void Distributor<Ordinal>::doReverseWaits() 
  {
    // call doWaits() on the reverse Distributor, if it exists
    if (reverseDistributor_ != Teuchos::null) {
      reverseDistributor_->doWaits();
    }
  }

  //! print method inherited from Teuchos::Object
  template <typename Ordinal>
  void Distributor<Ordinal>::print(ostream& os) const 
  {
    int const myImageID = comm_->getRank();
    int const numImages = comm_->getSize();
    for (int i = 0; i < numImages; i++) {
      comm_->barrier();
      if (i == myImageID) {
        os << "[Image " << myImageID << " of " << numImages << "]" << endl;
        os << " numExports: " << numExports_ << endl;
        os << " selfMessage: " << selfMessage_ << endl;
        os << " numSends_: " << numSends_ << endl;
        os << " imagesTo_: " << toString(imagesTo_) << endl;
        os << " startsTo_: " << toString(startsTo_) << endl;
        os << " lengthsTo_: " << toString(lengthsTo_) << endl;
        os << " maxSendLength_: " << maxSendLength_ << endl;
        os << " indicesTo_: " << toString(indicesTo_) << endl;
        os << " numReceives_: " << numReceives_ << endl;
        os << " totalReceiveLength_: " << totalReceiveLength_ << endl;
        os << " lengthsFrom_: " << toString(lengthsFrom_) << endl;
        os << " imagesFrom_: " << toString(imagesFrom_) << endl;
        os << " indicesFrom_: " << toString(indicesFrom_) << endl;
        os << " startsFrom_: " << toString(startsFrom_) << endl;
      }
    }
  }

  template <typename Ordinal>
  void Distributor<Ordinal>::computeReceives()
  {
    int myImageID = comm_->getRank();
    int numImages = comm_->getSize();
    const Ordinal zero = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal  one = Teuchos::OrdinalTraits<Ordinal>::one();
    std::vector<Ordinal> msg_count(numImages, zero);
    std::vector<Ordinal> counts(numImages, 1); // int vector for MPI sumAllAndScatter call

    for (Ordinal i = zero; i < (numSends_ + selfMessage_); i++) {
      msg_count[imagesTo_[i]] = one;
    }
    
    Teuchos::reduceAllAndScatter<Ordinal>(*comm_,Teuchos::REDUCE_SUM,numImages,&msg_count[0],&counts[0],&numReceives_);
    lengthsFrom_.assign(numReceives_, zero);
    imagesFrom_.assign(numReceives_, zero);

    for (Ordinal i = zero; i < (numSends_ + selfMessage_); i++) {
      if (imagesTo_[i] != myImageID ) {
        Teuchos::send(*comm_,lengthsTo_[i],imagesTo_[i]);
      }
      else {
        // set selfMessage_ to end block of recv arrays
        lengthsFrom_[numReceives_-one] = lengthsTo_[i];
        imagesFrom_[numReceives_-one] = myImageID;
      }
    }
    for (Ordinal i = zero; i < (numReceives_ - selfMessage_); i++) {
      // receive 1 Ordinal variable from any sender.
      // store the value in lengthsFrom_[i], and store the sender's ImageID in imagesFrom_[i]
      // imagesFrom_[i] = comm_->receive(&lengthsFrom_[i], 1, -1);
      imagesFrom_[i] = Teuchos::receive(*comm_,-1,&lengthsFrom_[i]);
    }
    comm_->barrier();

    sortArrays(imagesFrom_, lengthsFrom_);

    // Compute indicesFrom_
    totalReceiveLength_ = std::accumulate(lengthsFrom_.begin(), lengthsFrom_.end(), zero);
    indicesFrom_.clear();
    indicesFrom_.reserve(totalReceiveLength_);
    for (Ordinal i = 0; i < totalReceiveLength_; i++) {
      indicesFrom_.push_back(i);
    }

    startsFrom_.reserve(numReceives_);
    for (Ordinal i = zero, j = zero; i < numReceives_; ++i) {
      startsFrom_.push_back(j);
      j += lengthsFrom_[i];
    }

    numReceives_ -= selfMessage_;

    comm_->barrier();
  }

  template <typename Ordinal>
  void Distributor<Ordinal>::computeSends(
      const std::vector<Ordinal> & importIDs,
      const std::vector<Ordinal> & importImageIDs,
      std::vector<Ordinal>& exportIDs,
      std::vector<Ordinal>& exportImageIDs)
  {
    int myImageID = comm_->getRank();
    const Ordinal zero = Teuchos::OrdinalTraits<Ordinal>::zero();

    Distributor<Ordinal> tempPlan(comm_);
    std::vector<Ordinal> imageIDList;
    std::vector<Ordinal> importObjs;

    Ordinal numImports = importImageIDs.size();
    imageIDList = importImageIDs;

    importObjs.resize(2*numImports);

    for (Ordinal i = zero; i < numImports; i++ ) {  
      importObjs[2*i]   = importIDs[i];
      importObjs[2*i+1] = myImageID;
    }

    Ordinal numExports;
    tempPlan.createFromSends(imageIDList, numExports);
    exportIDs.resize(numExports);
    exportImageIDs.resize(numExports);

    std::vector<Ordinal> exportObjs;
    tempPlan.doPostsAndWaits(importObjs,2,exportObjs);

    for (Ordinal i = zero; i < numExports; i++) {
      exportIDs[i]      = exportObjs[2*i];
      exportImageIDs[i] = exportObjs[2*i+1];
    }
  }

} // namespace Tpetra

#endif // TPETRA_DISTRIBUTOR_HPP
