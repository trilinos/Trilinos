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


  Distributor::Distributor(const Teuchos::RCP<const Teuchos::Comm<int> > &comm) 
    : comm_(comm)
    , numExports_(0)
    , selfMessage_(false)
    , numSends_(0)
    , maxSendLength_(0)
    , numReceives_(0)
    , totalReceiveLength_(0)
  {}

  Distributor::Distributor(const Distributor & distributor) 
    : comm_(distributor.comm_)
    , numExports_(distributor.numExports_)
    , selfMessage_(distributor.selfMessage_)
    , numSends_(distributor.numSends_)
    , maxSendLength_(distributor.maxSendLength_)
    , numReceives_(distributor.numReceives_)
    , totalReceiveLength_(distributor.totalReceiveLength_)
    , reverseDistributor_(distributor.reverseDistributor_)
  {}

  Distributor::~Distributor() 
  {
  // we shouldn't have any outstanding requests at this point; verify
    TEST_FOR_EXCEPTION(requests_.size() != 0, std::runtime_error,
        Teuchos::typeName(*this) << "::Distributor~(): Destructor called with outstanding posts.");
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

  Teuchos::ArrayView<const int> Distributor::getImagesFrom() const 
  { return imagesFrom_; }

  Teuchos::ArrayView<const size_t> Distributor::getLengthsFrom() const 
  { return lengthsFrom_; }

  Teuchos::ArrayView<const int> Distributor::getImagesTo() const 
  { return imagesTo_; }

  Teuchos::ArrayView<const size_t> Distributor::getLengthsTo() const 
  { return lengthsTo_; }

  const Teuchos::RCP<Distributor> & 
  Distributor::getReverse() const {
    if (reverseDistributor_ == Teuchos::null) { 
      // need to create reverse distributor
      createReverseDistributor();
    }
    return reverseDistributor_;
  }


  void Distributor::createReverseDistributor() const {

    reverseDistributor_ = Teuchos::rcp(new Distributor(comm_));

    // compute new totalSendLength
    size_t totalSendLength = std::accumulate(lengthsTo_.begin(),lengthsTo_.end(),0);

    // compute new maxReceiveLength
    size_t maxReceiveLength = 0;
    const int myImageID = comm_->getRank();
    for (size_t i=0; i < numReceives_; ++i) {
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
    // Note: technically, I am my reverse distributor's reverse distributor, but 
    //       we will not set this up, as it gives us an opportunity to test 
    //       that reverseDistributor is an inverse operation w.r.t. value semantics of distributors
    // Note: numExports_ was not copied
  }


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


  void Distributor::doWaits() {
    if (getNumReceives() > 0) {
      Teuchos::waitAll(*comm_,requests_());
      // Requests should all be null, clear them
#ifdef HAVE_TEUCHOS_DEBUG
      for (Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest> >::const_iterator i = requests_.begin(); 
           i != requests_.end(); ++i) 
      {
        TEST_FOR_EXCEPTION(*i != Teuchos::null, std::runtime_error,
            Teuchos::typeName(*this) << "::doWaits(): Requests should be null after call to Teuchos::waitAll().");
      }
#endif
      requests_.clear();
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


  void Distributor::doReverseWaits() 
  {
    // call doWaits() on the reverse Distributor, if it exists
    if (reverseDistributor_ != Teuchos::null) {
      reverseDistributor_->doWaits();
    }
  }

  std::string Distributor::description() const
  {
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    return oss.str();
  }

  void Distributor::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
  {
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;
    const int myImageID = comm_->getRank();
    const int numImages = comm_->getSize();
    Teuchos::OSTab tab(out);
    if (vl != VERB_NONE) {
      // VERB_LOW and higher prints description()
      if (myImageID == 0) out << this->description() << std::endl; 
      for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
        if (myImageID == imageCtr) {
          if (vl != VERB_LOW) {
            out << "[Node " << myImageID << " of " << numImages << "]" << endl;
            out << " selfMessage: " << hasSelfMessage() << endl;
            out << " numSends: " << getNumSends() << endl;
            if (vl == VERB_HIGH || vl == VERB_EXTREME) {
              out << " imagesTo: " << toString(imagesTo_) << endl;
              out << " lengthsTo: " << toString(lengthsTo_) << endl;
              out << " maxSendLength: " << getMaxSendLength() << endl;
            }
            if (vl == VERB_EXTREME) {
              out << " startsTo: " << toString(startsTo_) << endl;
              out << " indicesTo: " << toString(indicesTo_) << endl;
            }
            if (vl == VERB_HIGH || vl == VERB_EXTREME) {
              out << " numReceives: " << getNumReceives() << endl;
              out << " totalReceiveLength: " << getTotalReceiveLength() << endl;
              out << " lengthsFrom: " << toString(lengthsFrom_) << endl;
              out << " imagesFrom: " << toString(imagesFrom_) << endl;
            }
          }
        }
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////
  void Distributor::computeReceives()
  {
    int myImageID = comm_->getRank();
    int numImages = comm_->getSize();

    // to_nodes_from_me[i] == number of messages sent by this node to node i
    // the info in numSends_, imagesTo_, lengthsTo_ concerns the contiguous sends
    // therefore, each node will be listed in imagesTo_ at most once
    {
      Teuchos::Array<size_t> to_nodes_from_me(numImages,0);
#     ifdef HAVE_TEUCHOS_DEBUG 
        bool counting_error = false;
#     endif
      for (size_t i=0; i < (numSends_ + (selfMessage_ ? 1 : 0)); ++i) {
#       ifdef HAVE_TEUCHOS_DEBUG
          if (to_nodes_from_me[imagesTo_[i]] != 0) counting_error = true;
#       endif
        to_nodes_from_me[imagesTo_[i]] = 1;
      }
#     ifdef HAVE_TEUCHOS_DEBUG
        SHARED_TEST_FOR_EXCEPTION(counting_error, std::logic_error,
            "Tpetra::Distributor::createFromSends: logic error. Please notify the Tpetra team.",*comm_);
#     endif
      // each proc will get back only one item (hence, counts = ones) from the array of globals sums, 
      // namely that entry corresponding to the node, and detailing how many receives it has.
      // this total includes self sends
      Teuchos::Array<int> counts(numImages, 1);
      Teuchos::reduceAllAndScatter<int,size_t>(*comm_,Teuchos::REDUCE_SUM,numImages,&to_nodes_from_me[0],&counts[0],&numReceives_);
    }

    // assign these to length numReceives, with zero entries
    lengthsFrom_.assign(numReceives_, 0);
    imagesFrom_.assign(numReceives_, 0);

    // FINISH: why do these work? they are blocking sends, and should block until completion, which happens below
    // FINISH: consider switching them to non-blocking
    // NOTE: epetra has both, old (non-blocking) and new (mysterious)

    for (size_t i=0; i < numSends_ + (selfMessage_ ? 1 : 0); ++i) {
      if (imagesTo_[i] != myImageID ) {
        // send a message to imagesTo_[i], telling him that our pattern sends him lengthsTo_[i] blocks of packets
        Teuchos::send(*comm_,lengthsTo_[i],imagesTo_[i]);
      }
      else {
        // set selfMessage_ to end block of recv arrays
        lengthsFrom_[numReceives_-1] = lengthsTo_[i];
        imagesFrom_[numReceives_-1] = myImageID;
      }
    }

    //
    for (size_t i=0; i < numReceives_ - (selfMessage_ ? 1 : 0); ++i) {
      // receive one variable from any sender.
      // store the value in lengthsFrom_[i], and store the sender's ImageID in imagesFrom_[i]
      // imagesFrom_[i] = comm_->receive(&lengthsFrom_[i], 1, -1);
      imagesFrom_[i] = Teuchos::receive(*comm_,-1,&lengthsFrom_[i]);
    }
    comm_->barrier();

    sort2(imagesFrom_.begin(), imagesFrom_.end(), lengthsFrom_.begin());

    // Compute indicesFrom_
    totalReceiveLength_ = std::accumulate(lengthsFrom_.begin(), lengthsFrom_.end(), 0);
    indicesFrom_.clear();
    indicesFrom_.reserve(totalReceiveLength_);
    for (size_t i=0; i < totalReceiveLength_; ++i) {
      indicesFrom_.push_back(i);
    }

    startsFrom_.clear();
    startsFrom_.reserve(numReceives_);
    for (size_t i=0, j = 0; i < numReceives_; ++i) {
      startsFrom_.push_back(j);
      j += lengthsFrom_[i];
    }

    if (selfMessage_) --numReceives_;

    comm_->barrier();
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
  size_t Distributor::createFromSends(const Teuchos::ArrayView<const int> &exportNodeIDs) {
    numExports_ = exportNodeIDs.size();

    const int myImageID = comm_->getRank();
    const int numImages = comm_->getSize();

    // exportNodeIDs tells us the communication pattern for this distributor
    // it dictates the way that the export data will be interpretted in doPosts()
    // we want to perform at most one communication per node; this is for two
    // reasons:
    //   * minimize latency/overhead in the comm routines (nice)
    //   * match the number of receives and sends between nodes (necessary)
    // Teuchos::Comm requires that the data for a send is contiguous in a send
    // buffer.
    // Therefore, if the data in the send buffer for doPosts() is not
    // contiguous, it will need to be copied into a contiguous buffer.
    // 
    // The user has specified this pattern and we can't do anything about it.
    // 
    // However, if they do not provide an efficient pattern, we will warn them 
    // if one of
    //    HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS 
    //    HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS 
    // is on.
    //
    // If the data is contiguous, then we can post the sends in situ.
    // 
    // Determine contiguity. There are a number of ways to do this:
    // * if the export IDs are sorted, then all exports to a particular 
    //   node must be contiguous. This is how Epetra does it. 
    // * if the export ID of the current export already has been listed,
    //   then the previous listing should correspond to the same export.
    //   This tests contiguity, but not sortedness.
    // Both of these tests require O(n), where n is the number of 
    // exports. However, the latter will positively identify a greater
    // portion of contiguous patterns. We will use the latter method.
    // 
    // Check to see if values are grouped by images without gaps
    // If so, indices_to -> 0.

    // Setup data structures for quick traversal of arrays.
    // This contains the number of sends for each image id.
    Teuchos::Array<size_t> starts(numImages + 1, 0);

    // numActive is the number of sends that are not Null
    size_t numActive = 0;
    char needSendBuff = 0;

    int badID = -1;
    for (size_t i = 0; i < numExports_; ++i) {
      int exportID = exportNodeIDs[i];
      if (exportID >= numImages) {
        badID = myImageID;
        break;
      }
      else if (exportID >= 0) {
        // increment starts[exportID]
        ++starts[exportID];
        // if after incrementing it is greater than one, check that the
        // previous export went to this node
        // this is a safe comparison, because starts[exportID] > 1
        // implies that i > 1. 
        // null entries break continuity.
        // e.g.,  [ 0, 0, 0, 1, -99, 1, 2, 2, 2] is not contiguous
        if (needSendBuff==0 && starts[exportID] > 1 && exportID != exportNodeIDs[i-1]) {
          needSendBuff = 1;
        }
        ++numActive;
      }
    }
    {
      int gbl_badID;
      Teuchos::reduceAll(*comm_,Teuchos::REDUCE_MAX,badID,&gbl_badID);
      TEST_FOR_EXCEPTION(gbl_badID >= 0, std::runtime_error,
          Teuchos::typeName(*this) << "::createFromSends(): bad node id listed on node " << gbl_badID << ".");
    }

#   if defined(HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS) || defined(HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS)
    {
      char global_needSendBuff;
      Teuchos::reduceAll(*comm_,Teuchos::REDUCE_MAX,needSendBuff,&global_needSendBuff);
      TPETRA_EFFICIENCY_WARNING(global_needSendBuff,std::runtime_error,
          "::createFromSends(): Grouping export IDs together leads to improved performance.");
    }
#   endif

    if (starts[myImageID] != 0) {
      selfMessage_ = true;
    }
    else {
      selfMessage_ = false;
    }


#ifdef HAVE_TEUCHOS_DEBUG
    bool index_neq_numActive = false;
    bool send_neq_numSends = false;
#endif
    if (!needSendBuff) {
      // grouped by image, no send buffer or indicesTo_ needed
      numSends_ = 0;
      // count total number of sends, i.e., total number of images that we are sending to
      // this includes myself
      for (int i=0; i < numImages; ++i) {
        if (starts[i]) ++numSends_;
      }

      // not only do we not need these, but we must clear it, as empty status of indicesTo 
      // is a flag used later
      indicesTo_.resize(0);
      // size these to numSends_; note, at the moment, numSends_ includes self sends
      // set their values to zeros
      imagesTo_.assign(numSends_,0);
      startsTo_.assign(numSends_,0);
      lengthsTo_.assign(numSends_,0);

      // set startsTo to the offsent for each send (i.e., each image ID)
      // set imagesTo to the image ID for each send
      // in interpretting this code, remember that we are assuming contiguity
      // that is why index skips through the ranks
      {
        size_t index = 0, nodeIndex = 0;
        for (size_t i = 0; i < numSends_; ++i) {
          while (exportNodeIDs[nodeIndex] < 0) {
            ++nodeIndex; // skip all negative node IDs
          }
          startsTo_[i] = nodeIndex;
          int imageID = exportNodeIDs[nodeIndex];
          imagesTo_[i] = imageID;
          index     += starts[imageID];
          nodeIndex += starts[imageID];
        }
#ifdef HAVE_TEUCHOS_DEBUG
        if (index != numActive) {
          index_neq_numActive = true;
        }
#endif
      }
      // sort the startsTo and image IDs together, in ascending order, according
      // to image IDs
      if (numSends_ > 0) {
        sort2(imagesTo_.begin(), imagesTo_.end(), startsTo_.begin());
      }
      // compute the maximum send length
      maxSendLength_ = 0;
      for (size_t i = 0; i < numSends_; ++i) {
        int imageID = imagesTo_[i];
        lengthsTo_[i] = starts[imageID];
        if ((imageID != myImageID) && (lengthsTo_[i] > maxSendLength_)) {
          maxSendLength_ = lengthsTo_[i];
        }
      }
    }
    else {
      // not grouped by image, need send buffer and indicesTo_

      // starts[i] is the number of sends to node i
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
      // starts[i] now contains the number of exports to nodes 0 through i

      for (Teuchos::Array<size_t>::reverse_iterator ip1=starts.rbegin(),
                                                      i=starts.rbegin()+1;
           i != starts.rend(); ++i)
      {
        *ip1 = *i;
        ip1 = i;
      }
      starts[0] = 0;
      // starts[i] now contains the number of exports to nodes 0 through
      // i-1, i.e., all nodes before node i

      indicesTo_.resize(numActive);

      for (size_t i = 0; i < numExports_; ++i) {
        if (exportNodeIDs[i] >= 0) {
          // record the offset to the sendBuffer for this export
          indicesTo_[starts[exportNodeIDs[i]]] = i;
          // now increment the offset for this node
          ++starts[exportNodeIDs[i]];
        }
      }
      // our send buffer will contain the export data for each of the nodes
      // we communicate with, in order by node id
      // sendBuffer = {node_0_data, node_1_data, ..., node_np-1_data}
      // indicesTo now maps each export to the location in our send buffer
      // associated with the export
      // data for export i located at sendBuffer[indicesTo[i]]
      //
      // starts[i] once again contains the number of exports to 
      // nodes 0 through i
      for (int node = numImages-1; node != 0; --node) {
        starts[node] = starts[node-1];
      }
      starts.front() = 0;       
      starts[numImages] = numActive;
      // 
      // starts[node] once again contains the number of exports to 
      // nodes 0 through node-1
      // i.e., the start of my data in the sendBuffer

      // this contains invalid data at nodes we don't care about, that is okay
      imagesTo_.resize(numSends_);
      startsTo_.resize(numSends_);
      lengthsTo_.resize(numSends_);

      // for each group of sends/exports, record the destination node,
      // the length, and the offset for this send into the 
      // send buffer (startsTo_)
      maxSendLength_ = 0;
      size_t snd = 0;
      for (int node = 0; node < numImages; ++node ) {
        if (starts[node+1] != starts[node]) {
          lengthsTo_[snd] = starts[node+1] - starts[node];
          startsTo_[snd] = starts[node];
          // record max length for all off-node sends
          if ((node != myImageID) && (lengthsTo_[snd] > maxSendLength_)) {
            maxSendLength_ = lengthsTo_[snd];
          }
          imagesTo_[snd] = node;
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

    return totalReceiveLength_;
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
