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

  template<typename OrdinalType>
  class Distributor : public Teuchos::Object {
  public:

    //@{ \name Constructor/Destructor

    //! Comm Constuctor (default ctr)
    Distributor(const Teuchos::RCP< Teuchos::Comm<OrdinalType> > & comm);

    //! Copy Constructor
    Distributor(const Distributor<OrdinalType> & distributor);

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
    void createFromSends(const std::vector<OrdinalType> & exportImageIDs,
                         OrdinalType & numRemoteIDs);

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
    void createFromRecvs(const std::vector<OrdinalType> & remoteGIDs, 
                         const std::vector<OrdinalType> & remoteImageIDs, 
                         std::vector<OrdinalType>& exportGIDs, 
                         std::vector<OrdinalType>& exportImageIDs);

    //@}

    //@{ \name Attribute Accessor Methods

    //! getTotalReceiveLength
    const OrdinalType & getTotalReceiveLength() const;

    //! getNumReceives
    const OrdinalType & getNumReceives() const;

    //! getSelfMessage - flag for if we're sending to ourself
    /*! If we are sending any elements to ourself, returns true. If we aren't, returns false. */
    bool getSelfMessage() const;

    //! getNumSends
    const OrdinalType & getNumSends() const;

    //! getMaxSendLength - maximum number of elements we're sending to a remote image
    const OrdinalType & getMaxSendLength() const;

    //! getImagesFrom - list of images sending elements to us
    const std::vector<OrdinalType> & getImagesFrom() const;

    //! getLengthsFrom - number of elements we're receiving from each image
    /*! We will receive lengthsFrom[i] elements from image imagesFrom[i] */
    const std::vector<OrdinalType> & getLengthsFrom() const;

    // getImagesTo - list of images we're sending elements to
    const std::vector<OrdinalType> & getImagesTo() const;

    // getIndicesTo
    /*! (Used only if exportImageIDs was not blocked by image.)
        Gives the order to the export buffer, in order to get
      a version that is sorted by imageID. */
    const std::vector<OrdinalType> & getIndicesTo() const;

    // getStartsTo - list of offsets into export buffer
    /*! Given an export buffer that contains all of the elements we're sending out, 
        image i's block of elements will start at position startsTo[i] */
    const std::vector<OrdinalType> & getStartsTo() const;

    // getLengthsTo - number of elements we're sending to each image
    /*! We will send lengthsTo[i] elements to image imagesTo[i] */
    const std::vector<OrdinalType> & getLengthsTo() const;

    //@}

    //@{ \name Reverse Communication Methods

    // getReverse
    //! Returns a Distributor with a reverse plan of this Distributor's plan
    /*! Creates the reverse Distributor if this is the first time this function
        has been called.
    */
    const Distributor<OrdinalType> & getReverse();

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
    template <typename ScalarType>
    void doPostsAndWaits(std::vector<ScalarType>& exports,
                         OrdinalType packetSize,
                         std::vector<ScalarType>& imports);

    //! doPosts
    template <typename ScalarType>
    void doPosts(std::vector<ScalarType>& exports,
                 OrdinalType packetSize,
                 std::vector<ScalarType>& imports);

    //! doWaits
    void doWaits();

    //! doReversePostsAndWaits
    /*! Execute a reverse plan specified by the distributor object.
      \param exports In
             On entry, contains the values we're exporting.
      \param packetSize In
             On entry, the number of ScalarType variables that make up an element.
      \param imports Out
             On exit, contains the values exported to us. (imports will be resized
             if necessary, and any existing values will be overwritten.)
    */
    template <typename ScalarType>
    void doReversePostsAndWaits(std::vector<ScalarType>& exports,
                                OrdinalType packetSize,
                                std::vector<ScalarType>& imports);

    //! doReversePosts
    template <typename ScalarType>
    void doReversePosts(std::vector<ScalarType>& exports,
                        OrdinalType packetSize,
                        std::vector<ScalarType>& imports);
    
    //! doReverseWaits
    void doReverseWaits();
    
    //@}

    //@{ \name I/O Methods

    //! Implements Teuchos::Object::print.
    void print(ostream& os) const;

    //@}

  private:

    // convenience functions for returning inner data class, both const and nonconst versions.
    Teuchos::Comm<OrdinalType> & getComm();
    const Teuchos::Comm<OrdinalType> & getComm() const;

    // private data members
    Teuchos::RCP< Teuchos::Comm<OrdinalType> > comm_;

    OrdinalType numExports_;
    OrdinalType selfMessage_;
    OrdinalType numSends_;
    std::vector<OrdinalType> imagesTo_;
    std::vector<OrdinalType> startsTo_;
    std::vector<OrdinalType> lengthsTo_;
    OrdinalType maxSendLength_;
    std::vector<OrdinalType> indicesTo_;
    OrdinalType numReceives_;
    OrdinalType totalReceiveLength_;
    std::vector<OrdinalType> lengthsFrom_;
    std::vector<OrdinalType> imagesFrom_;
    std::vector<OrdinalType> indicesFrom_;
    std::vector<OrdinalType> startsFrom_;

    Teuchos::RCP< Distributor<OrdinalType> > reverseDistributor_;

    void computeReceives();

    void computeSends(const std::vector<OrdinalType> & importIDs,
                      const std::vector<OrdinalType> & importImageIDs,
                      std::vector<OrdinalType>& exportIDs,
                      std::vector<OrdinalType>& exportImageIDs);

  }; // class Distributor


  template <typename OrdinalType>
  Distributor<OrdinalType>::Distributor(Teuchos::RCP< Teuchos::Comm<OrdinalType> > const& comm) 
    : Teuchos::Object("Tpetra::Distributor")
    , comm_(comm)
    , numExports_(Teuchos::OrdinalTraits<OrdinalType>::zero())
    , selfMessage_(false)
    , numSends_(Teuchos::OrdinalTraits<OrdinalType>::zero())
    , maxSendLength_(Teuchos::OrdinalTraits<OrdinalType>::zero())
    , numReceives_(Teuchos::OrdinalTraits<OrdinalType>::zero())
    , totalReceiveLength_(Teuchos::OrdinalTraits<OrdinalType>::zero())
  {}

  template <typename OrdinalType>
  Distributor<OrdinalType>::Distributor(Distributor<OrdinalType> const& distributor) 
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

  template <typename OrdinalType>
  Distributor<OrdinalType>::~Distributor() 
  {}

  template <typename OrdinalType>
  void Distributor<OrdinalType>::createFromSends(
      const std::vector<OrdinalType> & exportImageIDs,
      OrdinalType & numRemoteIDs) 
  {
    const OrdinalType zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    const OrdinalType one  = Teuchos::OrdinalTraits<OrdinalType>::one();

    numExports_ = exportImageIDs.size();

    int myImageID = getComm().getRank();
    int numImages = getComm().getSize();

    // Check to see if items are grouped by images without gaps
    // If so, indices_to -> 0

    // Setup data structures for quick traversal of arrays
    std::vector<OrdinalType> starts(numImages + 1, zero);

    OrdinalType numActive = zero;
    bool noSendBuff = true;

    for(OrdinalType i = zero; i < numExports_; i++) {
      if(noSendBuff && (i > zero) && (exportImageIDs[i] < exportImageIDs[i-one]))
        noSendBuff = false;
      if(exportImageIDs[i] >= zero) {
        ++starts[exportImageIDs[i]];
        ++numActive;
      }
    }

    if(starts[myImageID] != zero)
      selfMessage_ = one;
    else
      selfMessage_ = zero;

    numSends_ = zero;

    if(noSendBuff) { // grouped by image, no send buffer or indicesTo_ needed
      for(int i=0; i < numImages; ++i) {
        if(starts[i]) ++numSends_;
      }

      imagesTo_.resize(numSends_); // can we change these to reserves?
      startsTo_.resize(numSends_);
      lengthsTo_.resize(numSends_);

      for(OrdinalType i = zero, index = zero; i < numSends_; ++i) {
        startsTo_[i] = index;
        OrdinalType imageID = exportImageIDs[index];
        imagesTo_[i] = imageID;
        index += starts[imageID];
      }

      if(numSends_ > zero) {
        sortArrays(imagesTo_, startsTo_);
      }

      maxSendLength_ = zero;

      for(OrdinalType i = zero; i < numSends_; ++i) {
        OrdinalType imageID = imagesTo_[i];
        lengthsTo_[i] = starts[imageID];
        if((imageID != myImageID) && (lengthsTo_[i] > maxSendLength_)) {
          maxSendLength_ = lengthsTo_[i];
        }
      }
    }
    else { // not grouped by image, need send buffer and indicesTo_
      if(starts.front() != zero ) {
        numSends_ = one;
      }

      for(int i = 0; i < numImages; i++) {
        if(starts[i] != zero) {
          ++numSends_;
        }
        starts[i] += starts[i-one];
      }

      for(int i = numImages-1; i != 0; i--) {
        starts[i] = starts[i-one];
      }

      starts.front() = zero;

      indicesTo_.resize(numActive);

      for(OrdinalType i = zero; i < numExports_; i++) {
        if(exportImageIDs[i] >= zero) {
          indicesTo_[starts[exportImageIDs[i]]] = i;
          ++starts[exportImageIDs[i]];
        }
      }

      // Reconstuct starts array to index into indicesTo.

      for(int i = numImages-1; i != 0; i--) {
        starts[i] = starts[i-1];
      }
      starts.front() = zero;       
      starts[numImages] = numActive;

      imagesTo_.resize(numSends_); // can we change these to reserves?
      startsTo_.resize(numSends_);
      lengthsTo_.resize(numSends_);

      maxSendLength_ = zero;

      for(int i = 0, j = 0; i < numImages; i++ ) {
        if(starts[i+1] != starts[i]) {
          lengthsTo_[j] = starts[i+1] - starts[i];
          startsTo_[j] = starts[i];
          if((i != myImageID) && (lengthsTo_[j] > maxSendLength_)) {
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

  template <typename OrdinalType>
  void Distributor<OrdinalType>::createFromRecvs(
      const std::vector<OrdinalType> & remoteGIDs, 
      const std::vector<OrdinalType> & remoteImageIDs, 
      std::vector<OrdinalType> & exportGIDs, 
      std::vector<OrdinalType> & exportImageIDs)
  {
    computeSends(remoteGIDs, remoteImageIDs, exportGIDs, exportImageIDs);
    OrdinalType testNumRemoteIDs; // dummy
    createFromSends(exportImageIDs, testNumRemoteIDs);
  }

  template <typename OrdinalType>
  const OrdinalType & Distributor<OrdinalType>::getTotalReceiveLength() const 
  { return(totalReceiveLength_); }

  template <typename OrdinalType>
  const OrdinalType & Distributor<OrdinalType>::getNumReceives() const 
  { return(numReceives_); }

  template <typename OrdinalType>
  bool Distributor<OrdinalType>::getSelfMessage() const 
  { return(selfMessage_); }

  template <typename OrdinalType>
  const OrdinalType & Distributor<OrdinalType>::getNumSends() const 
  { return(numSends_); }

  template <typename OrdinalType>
  const OrdinalType & Distributor<OrdinalType>::getMaxSendLength() const 
  { return(maxSendLength_); }

  template <typename OrdinalType>
  const std::vector<OrdinalType> & Distributor<OrdinalType>::getImagesFrom() const 
  { return(imagesFrom_); }

  template <typename OrdinalType>
  const std::vector<OrdinalType> & Distributor<OrdinalType>::getLengthsFrom() const 
  { return(lengthsFrom_); }

  template <typename OrdinalType>
  const std::vector<OrdinalType> & Distributor<OrdinalType>::getImagesTo() const 
  { return(imagesTo_); }

  template <typename OrdinalType>
  const std::vector<OrdinalType> & Distributor<OrdinalType>::getIndicesTo() const 
  { return(indicesTo_); }

  template <typename OrdinalType>
  const std::vector<OrdinalType> & Distributor<OrdinalType>::getStartsTo() const 
  { return(startsTo_); }

  template <typename OrdinalType>
  const std::vector<OrdinalType> & Distributor<OrdinalType>::getLengthsTo() const 
  { return(lengthsTo_); }

  template <typename OrdinalType>
  const Distributor<OrdinalType> & Distributor<OrdinalType>::getReverse() 
  {
    if (reverseDistributor_ == Teuchos::null) { // need to initialize reverse distributor
      OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();

      reverseDistributor_ = Teuchos::rcp(new Distributor<OrdinalType>(comm_));

      // compute new totalSendLength
      OrdinalType totalSendLength = zero;
      for(OrdinalType i = zero; i < (numSends_ + selfMessage_); i++) {
        totalSendLength += lengthsTo_[i];
      }

      // compute new maxReceiveLength
      OrdinalType maxReceiveLength = zero;
      int const myImageID = getComm().getRank();
      for(OrdinalType i = zero; i < numReceives_; i++)
        if(imagesFrom_[i] != myImageID)
          if(lengthsFrom_[i] > maxReceiveLength)
            maxReceiveLength = lengthsFrom_[i];
      
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

    return(*reverseDistributor_);
  }

  template <typename OrdinalType>
  template <typename ScalarType>
  void Distributor<OrdinalType>::doPostsAndWaits(
      std::vector<ScalarType>& exports,
      OrdinalType packetSize,
      std::vector<ScalarType>& imports) 
  {
    (void)exports;
    (void)packetSize;
    (void)imports;
    // FINISH
  }

  template <typename OrdinalType>
  template <typename ScalarType>
  void Distributor<OrdinalType>::doPosts(
      std::vector<ScalarType>& exports,
      OrdinalType packetSize,
      std::vector<ScalarType>& imports) 
  {
    (void)exports;
    (void)packetSize;
    (void)imports;
    // FINISH
  }

  template <typename OrdinalType>
  void Distributor<OrdinalType>::doWaits() 
  {
    // FINISH
  }

  template <typename OrdinalType>
  template <typename ScalarType>
  void Distributor<OrdinalType>::doReversePostsAndWaits(
      std::vector<ScalarType>& exports,
      OrdinalType packetSize,
      std::vector<ScalarType>& imports) 
  {
    (void)exports;
    (void)packetSize;
    (void)imports;
    // FINISH
  }

  template <typename OrdinalType>
  template <typename ScalarType>
  void Distributor<OrdinalType>::doReversePosts(
      std::vector<ScalarType>& exports,
      OrdinalType packetSize,
      std::vector<ScalarType>& imports) 
  {
    (void)exports;
    (void)packetSize;
    (void)imports;
    // FINISH
  }

  template <typename OrdinalType>
  void Distributor<OrdinalType>::doReverseWaits() 
  {
    // FINISH
  }
  
  //! print method inherited from Teuchos::Object
  template <typename OrdinalType>
  void Distributor<OrdinalType>::print(ostream& os) const 
  {
    int const myImageID = getComm().getRank();
    int const numImages = getComm().getSize();
    for(int i = 0; i < numImages; i++) {
      getComm().barrier();
      if(i == myImageID) {
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

  template <typename OrdinalType>
  Teuchos::Comm<OrdinalType> & Distributor<OrdinalType>::getComm() 
  {return(*comm_);}

  template <typename OrdinalType>
  const Teuchos::Comm<OrdinalType> & Distributor<OrdinalType>::getComm() const 
  {return(*comm_);}

  template <typename OrdinalType>
  void Distributor<OrdinalType>::computeReceives()
  {
    int myImageID = getComm().getRank();
    int numImages = getComm().getSize();
    const OrdinalType zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    const OrdinalType  one = Teuchos::OrdinalTraits<OrdinalType>::one();
    std::vector<OrdinalType> msg_count(numImages, zero);
    std::vector<OrdinalType> counts(numImages, 1); // int vector for MPI sumAllAndScatter call

    for(OrdinalType i = zero; i < (numSends_ + selfMessage_); i++) {
      msg_count[imagesTo_[i]] = one;
    }
    
    Teuchos::reduceAllAndScatter<OrdinalType>(*comm_,Teuchos::REDUCE_SUM,numImages,&msg_count[0],&counts[0],&numReceives_);
    lengthsFrom_.assign(numReceives_, zero);
    imagesFrom_.assign(numReceives_, zero);

    for(OrdinalType i = zero; i < (numSends_ + selfMessage_); i++) {
      if(imagesTo_[i] != myImageID ) {
        Teuchos::send(getComm(),lengthsTo_[i],imagesTo_[i]);
      }
      else {
        // set selfMessage_ to end block of recv arrays
        lengthsFrom_[numReceives_-one] = lengthsTo_[i];
        imagesFrom_[numReceives_-one] = myImageID;
      }
    }
    for(OrdinalType i = zero; i < (numReceives_ - selfMessage_); i++) {
      // receive 1 OrdinalType variable from any sender.
      // store the value in lengthsFrom_[i], and store the sender's ImageID in imagesFrom_[i]
      // imagesFrom_[i] = getComm().receive(&lengthsFrom_[i], 1, -1);
      imagesFrom_[i] = Teuchos::receive(getComm(),-1,&lengthsFrom_[i]);
    }
    getComm().barrier();

    sortArrays(imagesFrom_, lengthsFrom_);

    // Compute indicesFrom_
    totalReceiveLength_ = std::accumulate(lengthsFrom_.begin(), lengthsFrom_.end(), zero);
    indicesFrom_.clear();
    indicesFrom_.reserve(totalReceiveLength_);
    for(OrdinalType i = 0; i < totalReceiveLength_; i++) {
      indicesFrom_.push_back(i);
    }

    startsFrom_.reserve(numReceives_);
    for(OrdinalType i = zero, j = zero; i < numReceives_; ++i) {
      startsFrom_.push_back(j);
      j += lengthsFrom_[i];
    }

    numReceives_ -= selfMessage_;

    getComm().barrier();
  }

  template <typename OrdinalType>
  void Distributor<OrdinalType>::computeSends(
      const std::vector<OrdinalType> & importIDs,
      const std::vector<OrdinalType> & importImageIDs,
      std::vector<OrdinalType>& exportIDs,
      std::vector<OrdinalType>& exportImageIDs)
  {
    int myImageID = getComm().getRank();
    int numImages = getComm().getSize();
    const OrdinalType zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    const OrdinalType one = Teuchos::OrdinalTraits<OrdinalType>::one();
    const OrdinalType two = one + one;

    Distributor<OrdinalType> tempPlan(comm_);
    std::vector<OrdinalType> imageIDList;
    std::vector<OrdinalType> importObjs;

    OrdinalType numImports = importImageIDs.size();
    imageIDList = importImageIDs;

    importObjs.resize(two * numImports);

    for(OrdinalType i = zero; i < numImports; i++ ) {  
      importObjs[two*i] = importIDs[i];
      importObjs[two*i+1] = myImageID;
    }

    OrdinalType numExports;
    tempPlan.createFromSends(imageIDList, numExports);
    exportIDs.resize(numExports);
    exportImageIDs.resize(numExports);

    std::vector<OrdinalType> exportObjs;
    // comm_->doPostsAndWaits(tempPlan, importObjs, two, exportObjs);
    tempPlan.doPostsAndWaits(importObjs,two,exportObjs);

    for(OrdinalType i = zero; i < numExports; i++) {
      exportIDs[i] = exportObjs[two*i];
      exportImageIDs[i] = exportObjs[two*i+one];
    }
  }

} // namespace Tpetra

#endif // TPETRA_DISTRIBUTOR_HPP
