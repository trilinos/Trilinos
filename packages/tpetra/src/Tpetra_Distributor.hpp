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
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>


// FINISH: some of the get accessors may not be necessary anymore. clean up.
// FINISH: This class may not be const correct. doPosts() et al. perhaps should be const, with affected members made mutable.

namespace Tpetra {

  //! Tpetra::Distributor:  The Tpetra Gather/Scatter Setup Class.
  /*! The Distributor class is an interface that encapsulates the general
        information and services needed for other Tpetra classes to perform gather/scatter
        operations on a parallel computer.
  */

  template<typename Ordinal>
  class Distributor : public Teuchos::Object {
  public:

    //! @name Constructor/Destructor
    //@{ 

    //! Comm Constuctor (default ctr)
    Distributor(const Teuchos::RCP<const Teuchos::Comm<Ordinal> > & comm);

    //! Copy Constructor
    Distributor(const Distributor<Ordinal> & distributor);

    //! Destructor.
    ~Distributor();

    //@}


    //! @name Gather/Scatter Constructors
    //@{ 

    //! Create Distributor object using list of ImageIDs to send to
    /*! Take a list of ImageIDs and construct a plan for efficiently scattering to these images.
        Return the number of IDs being sent to me.
      \param exportImageIDs In
             List of images that will get the exported data. Image IDs less than zero
             are ignored; their placement corresponds to null sends in any
             future exports.
      \param numImports Out
             Number of imports this image will be receiving.
    */
    void createFromSends(const Teuchos::ArrayView<const Ordinal> &exportImageIDs,
                         Teuchos_Ordinal &numImports);

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
    void createFromRecvs(const Teuchos::ArrayView<const Ordinal> & remoteGIDs, 
                         const Teuchos::ArrayView<const Ordinal> & remoteImageIDs, 
                         Teuchos::ArrayRCP<Ordinal> &exportGIDs, 
                         Teuchos::ArrayRCP<Ordinal> &exportImageIDs);

    //@}

    //! @name Attribute Accessor Methods
    //@{ 

    //! getNumReceives
    const Ordinal & getNumReceives() const;

    //! getNumSends
    const Ordinal & getNumSends() const;

    //! getSelfMessage - flag for if we're sending to ourself
    /*! If we are sending any elements to ourself, returns true. If we aren't, returns false. */
    bool getSelfMessage() const;

    //! getMaxSendLength - maximum number of elements we're sending to a remote image
    const Ordinal & getMaxSendLength() const;

    //! getTotalReceiveLength
    const Ordinal & getTotalReceiveLength() const;

    //! getImagesFrom - list of images sending elements to us (non-persisting view)
    Teuchos::ArrayView<const Ordinal> getImagesFrom() const;

    //! getImagesTo - list of images we're sending elements to (non-persisting view)
    Teuchos::ArrayView<const Ordinal> getImagesTo() const;

    //! getLengthsFrom - number of elements we're receiving from each image (non-persisting view)
    /*! We will receive lengthsFrom[i] elements from image imagesFrom[i] */
    Teuchos::ArrayView<const Ordinal> getLengthsFrom() const;

    //! getLengthsTo - number of elements we're sending to each image (non-persisting view)
    /*! We will send lengthsTo[i] elements to image imagesTo[i] */
    Teuchos::ArrayView<const Ordinal> getLengthsTo() const;

    //! getStartsTo - list of offsets into export buffer (non-persisting view)
    /*! Given an export buffer that contains all of the elements we're sending out, 
        image i's block of elements will start at position startsTo[i] */
    Teuchos::ArrayView<const Ordinal> getStartsTo() const;

    //! getIndicesTo (non-persisting view)
    /*! (Used only if exportImageIDs was not blocked by image.)
        Gives the order to the export buffer, in order to get
      a version that is sorted by imageID. */
    Teuchos::ArrayView<const Ordinal> getIndicesTo() const;

    //@}

    //! @name Reverse Communication Methods
    //@{ 

    // getReverse
    //! Returns a Distributor with a reverse plan of this Distributor's plan
    /*! Creates the reverse Distributor if this is the first time this function
        has been called.
    */
    Teuchos::RCP<Distributor<Ordinal> > getReverse() const;

    //@}

    //! @name Execute Distributor Plan Methods
    //@{ 

    //! doPostsAndWaits
    /*! Execute a plan specified by the distributor object.
      \param exports In
             Contains the values we're exporting.
      \param numPackets In
             Specifies the number of scalars per export/import.
      \param imports Out
             On entry, buffer must be large enough to accomodate the data exported to us.
             On exit, contains the values exported to us.
    */
    template <typename Packet>
    void doPostsAndWaits(const Teuchos::ArrayView<const Packet> &exports,
                         Ordinal numPackets,
                         const Teuchos::ArrayView<Packet> &imports);

    //! doPosts
    /*! Post the data for a distributor plan, but do not execute the waits yet.
      \param exports In
             Constains the values we're exporting.
      \param numPackets In
             Specifies the number of scalars per export/import.
      \param imports In
             Buffer must be large enough to accomodate the data exported to us. 
             The buffer is not guaranteed to be filled until doWaits() is executed.
    */
    template <typename Packet>
    void doPosts(const Teuchos::ArrayView<const Packet> &exports,
                 Ordinal numPackets,
                 const Teuchos::ArrayRCP<Packet> &imports);

    //! doWaits
    void doWaits();

    //! doReversePostsAndWaits
    /*! Execute a reverse plan specified by the distributor object.
      \param exports In
             Contains the values we're exporting.
      \param numPackets In
             Specifies the number of scalars per export/import.
      \param imports Out
             On entry, buffer must be large enough to accomodate the data exported to us.
             On exit, contains the values exported to us.
    */
    template <typename Packet>
    void doReversePostsAndWaits(const Teuchos::ArrayView<const Packet> &exports,
                                Ordinal numPackets,
                                const Teuchos::ArrayView<Packet> &imports);

    //! doReversePosts
    /*! Post the data for a reverse plan, but do not execute the waits yet.
      \param exports In
             Constains the values we're exporting.
      \param numPackets In
             Specifies the number of scalars per export/import.
      \param imports In
             Buffer must be large enough to accomodate the data exported to us. 
             The buffer is not guaranteed to be filled until doWaits() is executed.
    */
    template <typename Packet>
    void doReversePosts(const Teuchos::ArrayView<const Packet> &exports,
                        Ordinal numPackets,
                        const Teuchos::ArrayRCP<Packet> &imports);

    //! doReverseWaits
    void doReverseWaits();

    //@}

    //! @name I/O Methods
    //@{ 

    //! Implements Teuchos::Object::print.
    void print(std::ostream& os) const;

    //@}

  private:

    // private data members
    Teuchos::RCP<const Teuchos::Comm<Ordinal> > comm_;

    Ordinal numExports_;
    // selfMessage_ is whether I have a send for myself
    bool selfMessage_;
    // numSends_ is number of sends to other nodes
    Ordinal numSends_;
    // imagesTo_, startsTo_ and lengthsTo_ each have size 
    //   numSends_ + selfMessage_
    Teuchos::Array<Ordinal> imagesTo_;
    Teuchos::Array<Ordinal> startsTo_;
    Teuchos::Array<Ordinal> lengthsTo_;
    // maxSendLength_ is the maximum send to another node: 
    //   max(lengthsTo_[i]) for i != me
    Ordinal maxSendLength_;
    Teuchos::Array<Ordinal> indicesTo_;
    // numReceives_ is the number of receives by me from other procs, not
    // counting self receives
    Ordinal numReceives_;
    // totalReceiveLength_ is the total number of Packet received, used to 
    // allocate the receive buffer
    Ordinal totalReceiveLength_;
    // imagesFrom_, startsFrom_ and lengthsFrom_ each have size 
    //   numReceives_ + selfMessage_
    Teuchos::Array<Ordinal> lengthsFrom_;
    Teuchos::Array<Ordinal> imagesFrom_;
    Teuchos::Array<Ordinal> startsFrom_;
    Teuchos::Array<Ordinal> indicesFrom_;

    // requests associated with non-blocking receives
    Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest> > requests_;

    mutable Teuchos::RCP<Distributor<Ordinal> > reverseDistributor_;

    // compute receive info from sends
    void computeReceives();

    // compute send info from receives
    void computeSends(const Teuchos::ArrayView<const Ordinal> &importIDs,
                      const Teuchos::ArrayView<const Ordinal> &importImageIDs,
                            Teuchos::ArrayRCP<Ordinal> &exportIDs,
                            Teuchos::ArrayRCP<Ordinal> &exportImageIDs);

    // create a distributor for the reverse communciation pattern (pretty much
    // swap all send and receive info)
    void createReverseDistributor() const;

  }; // class Distributor


  template <typename Ordinal>
  Distributor<Ordinal>::Distributor(const Teuchos::RCP<const Teuchos::Comm<Ordinal> > &comm) 
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
  {
  // we shouldn't have any outstanding requests at this point; verify
# ifdef TEUCHOS_DEBUG
    for (typename Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest> >::const_iterator i = requests_.begin(); 
         i != requests_.end(); ++i)
    {
      TEST_FOR_EXCEPTION(*i != Teuchos::null, std::runtime_error,
          "Tpetra::Distributor<"<<Teuchos::OrdinalTraits<Ordinal>::name()
          <<">::doWaits(): Requests should be null after call to Teuchos::waitAll().");
    }
#endif
  }

  template <typename Ordinal>
  void Distributor<Ordinal>::createFromSends(
      const Teuchos::ArrayView<const Ordinal> &exportImageIDs,
      Teuchos_Ordinal &numImports) 
  {
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE  = Teuchos::OrdinalTraits<Ordinal>::one();

    numExports_ = exportImageIDs.size();

    const int myImageID = comm_->getRank();
    const int numImages = comm_->getSize();

    // exportImageIDs tells us the communication pattern for this distributor
    // it dictates the way that the export data will be interpretted in doPosts()
    // we want to perform at most one communication per node; this is for two
    // reasons:
    //   * minimize latency/overhead in the comm routines (nice)
    //   * match the number of receives and sends between nodes (necessary)
    // Teuchos::Comm requires that the data for a send is contiguous in a send
    // buffer
    // therefore, if the data in the send buffer for doPosts() is not
    // contiguous, it will need to be copied into a contiguous buffer
    // 
    // the user has specified this pattern and we can't do anything about it,
    // 
    // however, if they do not provide an efficient pattern, we will warn them 
    // if one of
    //    THROW_TPETRA_EFFICIENCY_WARNINGS 
    //    PRINT_TPETRA_EFFICIENCY_WARNINGS 
    // is on.
    //
    // if the data is contiguous, then we can post the sends in situ
    // 
    // determine contiguity. there are a number of ways to do this:
    // * if the export IDs are sorted, then all exports to a particular 
    //   node must contiguous. this is how epetra does it. 
    // * if the export ID of the current export already has been listed,
    //   then the previous listing should correspond to the same export.
    //   this tests contiguity, but not sortedness.
    // both of these tests require O(n), where n is the number of 
    // exports. however, the latter will positively identify a greater
    // portion of contiguous patterns. we will use the latter method.
    // 
    // Check to see if items are grouped by images without gaps
    // If so, indices_to -> 0

    // Setup data structures for quick traversal of arrays
    // this contains the number of sends for each image id
    Teuchos::Array<Ordinal> starts(numImages + 1, ZERO);

    // numActive is the number of sends that are not Null
    Ordinal numActive = ZERO;
    int needSendBuff = 0;

    for (int i = 0; i < numExports_; ++i) {
      Ordinal exportID = exportImageIDs[i];
      if (exportID >= ZERO) {
        // increment starts[exportID]
        ++starts[exportID];
        // if after incrementing it is greater than one, check that the
        // previous export went to this node
        // this is a safe comparison, because starts[exportID] > 1
        // implies that i > 1. 
        // null entries break continuity.
        // e.g.,  [ 0, 0, 0, 1, -99, 1, 2, 2, 2] is not considered contiguous
        if (needSendBuff==0 && starts[exportID]>1 && exportID != exportImageIDs[i-1]) {
          needSendBuff = 1;
        }
        ++numActive;
      }
    }

#   if defined(THROW_TPETRA_EFFICIENCY_WARNINGS) || defined(PRINT_TPETRA_EFFICIENCY_WARNINGS)
    {
      int global_needSendBuff;
      Teuchos::reduceAll(*comm_,Teuchos::REDUCE_SUM,needSendBuff,&global_needSendBuff);
      std::string err;
      err += "Tpetra::Distributor<" + Teuchos::TypeNameTraits<Ordinal>::name() 
        + ">::createFromSends(): Grouping export IDs together leads to improved performance.";
#   if defined(THROW_TPETRA_EFFICIENCY_WARNINGS)
      TEST_FOR_EXCEPTION(global_needSendBuff != 0, std::runtime_error, err);
#   else // print warning
      if (global_needSendBuff) {
        std::cerr << err << std::endl;
      }
    }
#   endif
#   endif

    if (starts[myImageID] != ZERO) {
      selfMessage_ = true;
    }
    else {
      selfMessage_ = false;
    }


#ifdef TEUCHOS_DEBUG
    bool index_neq_numActive = false;
    bool send_neq_numSends = false;
#endif
    if (!needSendBuff) {
      // grouped by image, no send buffer or indicesTo_ needed
      numSends_ = ZERO;
      for (int i=0; i < numImages; ++i) {
        if (starts[i]) ++numSends_;
      }

      // not only do we not need these, but empty indicesTo is also a flag for
      // later
      indicesTo_.resize(0);
      // size these to numSends_; note, at the moment, numSends_ includes sends
      // to myself; set their values to zeros
      imagesTo_.assign(numSends_,ZERO);
      startsTo_.assign(numSends_,ZERO);
      lengthsTo_.assign(numSends_,ZERO);

      // set startsTo to the offsent for each send (i.e., each image ID)
      // set imagesTo to the image ID for each send
      // in interpretting this code, remember that we are assuming contiguity
      // that is why index skips through the ranks
      {
        Ordinal index = ZERO;
        for (Ordinal i = ZERO; i < numSends_; ++i) {
          startsTo_[i] = index;
          Ordinal imageID = exportImageIDs[index];
          imagesTo_[i] = imageID;
          index += starts[imageID];
        }
#ifdef TEUCHOS_DEBUG
        if (index != numActive) {
          index_neq_numActive = true;
        }
#endif
      }

      // sort the startsTo and image IDs together, in ascending order, according
      // to image IDs
      if (numSends_ > ZERO) {
        sortArrays(imagesTo_(), startsTo_());
      }

      // compute the maximum send length
      maxSendLength_ = ZERO;
      for (int i = 0; i < numSends_; ++i) {
        Ordinal imageID = imagesTo_[i];
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
      if (starts[0] == ZERO ) {
        numSends_ = ZERO;
      }
      else {
        numSends_ = ONE;
      }
      for (typename Teuchos::Array<Ordinal>::iterator i=starts.begin()+1,
                                                 im1=starts.begin();
           i != starts.end(); ++i) 
      {
        if (*i != ZERO) ++numSends_;
        *i += *im1;
        im1 = i;
      }
      // starts[i] now contains the number of exports to nodes 0 through i

      for (typename Teuchos::Array<Ordinal>::reverse_iterator ip1=starts.rbegin(),
                                                             i=starts.rbegin()+1;
           i != starts.rend(); ++i)
      {
        *ip1 = *i;
        ip1 = i;
      }
      starts[0] = ZERO;
      // starts[i] now contains the number of exports to nodes 0 through
      // i-1, i.e., all nodes before node i

      indicesTo_.resize(numActive);

      for (Ordinal i = ZERO; i < numExports_; ++i) {
        if (exportImageIDs[i] >= ZERO) {
          // record the offset to the sendBuffer for this export
          indicesTo_[starts[exportImageIDs[i]]] = i;
          // now increment the offset for this node
          ++starts[exportImageIDs[i]];
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
      starts.front() = ZERO;       
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
      maxSendLength_ = ZERO;
      int snd = 0;
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
#ifdef TEUCHOS_DEBUG
      if (snd != numSends_) {
        send_neq_numSends = true;
      }
#endif
    }
#ifdef TEUCHOS_DEBUG
        SHARED_TEST_FOR_EXCEPTION(index_neq_numActive, std::logic_error,
            "Tpetra::Distributor::createFromSends: logic error. Please notify the Tpetra team.",*comm_);
        SHARED_TEST_FOR_EXCEPTION(send_neq_numSends, std::logic_error,
            "Tpetra::Distributor::createFromSends: logic error. Please notify the Tpetra team.",*comm_);
#endif

    if (selfMessage_) --numSends_;

    // Invert map to see what msgs are received and what length
    computeReceives();

    numImports = totalReceiveLength_;
  }

  template <typename Ordinal>
  void Distributor<Ordinal>::createFromRecvs(
      const Teuchos::ArrayView<const Ordinal> &remoteGIDs, 
      const Teuchos::ArrayView<const Ordinal> &remoteImageIDs, 
            Teuchos::ArrayRCP<Ordinal> &exportGIDs, 
            Teuchos::ArrayRCP<Ordinal> &exportImageIDs)
  {
    computeSends(remoteGIDs, remoteImageIDs, exportGIDs, exportImageIDs);
    Teuchos_Ordinal testNumRemoteIDs; // dummy
    createFromSends(exportImageIDs(), testNumRemoteIDs);
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
  Teuchos::ArrayView<const Ordinal> Distributor<Ordinal>::getImagesFrom() const 
  { return imagesFrom_; }

  template <typename Ordinal>
  Teuchos::ArrayView<const Ordinal> Distributor<Ordinal>::getLengthsFrom() const 
  { return lengthsFrom_; }

  template <typename Ordinal>
  Teuchos::ArrayView<const Ordinal> Distributor<Ordinal>::getImagesTo() const 
  { return imagesTo_; }

  template <typename Ordinal>
  Teuchos::ArrayView<const Ordinal> Distributor<Ordinal>::getIndicesTo() const 
  { return indicesTo_; }

  template <typename Ordinal>
  Teuchos::ArrayView<const Ordinal> Distributor<Ordinal>::getStartsTo() const 
  { return startsTo_; }

  template <typename Ordinal>
  Teuchos::ArrayView<const Ordinal> Distributor<Ordinal>::getLengthsTo() const 
  { return lengthsTo_; }

  template <typename Ordinal>
  Teuchos::RCP<Distributor<Ordinal> > Distributor<Ordinal>::getReverse() const
  {
    if (reverseDistributor_ == Teuchos::null) { 
      // need to create reverse distributor
      createReverseDistributor();
    }
    return reverseDistributor_;
  }

  template <typename Ordinal>
  void Distributor<Ordinal>::createReverseDistributor() const {
    const Ordinal zero = Teuchos::OrdinalTraits<Ordinal>::zero();

    reverseDistributor_ = Teuchos::rcp(new Distributor<Ordinal>(comm_));

    // compute new totalSendLength
    Ordinal totalSendLength = zero;
    for (Ordinal i = zero; i < (numSends_ + (selfMessage_ ? 1 : 0)); ++i) {
      totalSendLength += lengthsTo_[i];
    }

    // compute new maxReceiveLength
    Ordinal maxReceiveLength = zero;
    const int myImageID = comm_->getRank();
    for (Ordinal i = zero; i < numReceives_; ++i) {
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

  template <typename Ordinal>
  template <typename Packet>
  void Distributor<Ordinal>::doPostsAndWaits(
      const Teuchos::ArrayView<const Packet>& exports,
      Ordinal numPackets,
      const Teuchos::ArrayView<Packet>& imports) 
  {
    // doPosts takes imports as an ArrayRCP, requiring that the memory location is persisting
    // however, it need only persist until doWaits is called, so it is safe for us to 
    // use a non-persisting reference in this case
    doPosts(exports, numPackets, Teuchos::arcp<Packet>(imports.getRawPtr(),0,imports.size(),false));
    doWaits();
  }

  template <typename Ordinal>
  template <typename Packet>
  void Distributor<Ordinal>::doPosts(
      const Teuchos::ArrayView<const Packet>& exports,
      Ordinal numPackets,
      const Teuchos::ArrayRCP<Packet>& imports) 
  {
    using Teuchos::ArrayRCP;

    // start of actual doPosts function
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE  = Teuchos::OrdinalTraits<Ordinal>::one();
    const Ordinal myImageID = comm_->getRank();
    int selfReceiveOffset = 0;

#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(imports.size() != totalReceiveLength_ * numPackets, std::runtime_error,
        "Tpetra::Distributor<"<<Teuchos::OrdinalTraits<Ordinal>::name()
        <<">::doPosts(): imports must be large enough to store the imported data.");
#endif

    // allocate space in requests
    requests_.resize(0);
    requests_.reserve(numReceives_);

    // start up the Irecv's
    {
      int curBufferOffset = 0;
      for (int i = 0; i < (numReceives_ + (selfMessage_ ? 1 : 0)); ++i) {
        if (imagesFrom_[i] != myImageID) { 
          // receiving this one from another image
          // setup reference into imports of the appropriate size and at the appropriate place
          ArrayRCP<Packet> impptr = imports.persistingView(curBufferOffset,lengthsFrom_[i]*numPackets);
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
      ++imageIndex;
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

      if (selfMessage_) {
        std::copy(exports.begin()+startsTo_[selfNum]*numPackets, exports.begin()+startsTo_[selfNum]*numPackets+lengthsTo_[selfNum]*numPackets, 
                  imports.begin()+selfReceiveOffset);
      }
    }
    else { // data is not blocked by image, use send buffer
      // allocate sendArray buffer
      Teuchos::Array<Packet> sendArray(maxSendLength_*numPackets); 

      for (Ordinal i = ZERO; i < numBlocks; ++i) {
        Ordinal p = i + imageIndex;
        if (p > (numBlocks - ONE)) {
          p -= numBlocks;
        }

        if (imagesTo_[p] != myImageID) { 
          // sending it to another image
          typename Teuchos::ArrayView<const Packet>::iterator srcBegin, srcEnd;
          int sendArrayOffset = 0;
          int j = startsTo_[p];
          for (Ordinal k = ZERO; k < lengthsTo_[p]; ++k, ++j) {
            srcBegin = exports.begin() + indicesTo_[j]*numPackets;
            srcEnd   = srcBegin + numPackets;
            std::copy( srcBegin, srcEnd, sendArray.begin()+sendArrayOffset );
            sendArrayOffset += numPackets;
          }
          Teuchos::ArrayView<const Packet> tmpSend = sendArray(0,lengthsTo_[p]*numPackets);
          Teuchos::readySend<Ordinal,Packet>(*comm_,tmpSend,imagesTo_[p]);
        }
        else { 
          // sending it to myself
          selfNum = p;
          selfIndex = startsTo_[p];
        }
      }

      if (selfMessage_) {
        for (Ordinal k = ZERO; k < lengthsTo_[selfNum]; ++k) {
          std::copy( exports.begin()+indicesTo_[selfIndex]*numPackets,
                     exports.begin()+indicesTo_[selfIndex]*numPackets + numPackets,
                     imports.begin() + selfReceiveOffset );
          ++selfIndex;
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
      Teuchos::waitAll(*comm_,requests_());
      // Requests should all be null, clear them
#ifdef TEUCHOS_DEBUG
      for (typename Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest> >::const_iterator i = requests_.begin(); 
           i != requests_.end(); ++i) 
      {
        TEST_FOR_EXCEPTION(*i != Teuchos::null, std::runtime_error,
            "Tpetra::Distributor<"<<Teuchos::OrdinalTraits<Ordinal>::name()
            <<">::doWaits(): Requests should be null after call to Teuchos::waitAll().");
      }
#endif
      requests_.empty();
    }
  }

  template <typename Ordinal>
  template <typename Packet>
  void Distributor<Ordinal>::doReversePostsAndWaits(
      const Teuchos::ArrayView<const Packet>& exports,
      Ordinal numPackets,
      const Teuchos::ArrayView<Packet>& imports) 
  {
    // doPosts takes imports as an ArrayRCP, requiring that the memory location is persisting
    // however, it need only persist until doWaits is called, so it is safe for us to 
    // use a non-persisting reference in this case
    doReversePosts(exports, numPackets, Teuchos::arcp<Packet>(imports.getRawPtr(),0,imports.size(),false));
    doReverseWaits();
  }

  template <typename Ordinal>
  template <typename Packet>
  void Distributor<Ordinal>::doReversePosts(
      const Teuchos::ArrayView<const Packet>& exports,
      Ordinal numPackets,
      const Teuchos::ArrayRCP<Packet>& imports) 
  {
    TEST_FOR_EXCEPTION(!indicesTo_.empty(),std::runtime_error,
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
  void Distributor<Ordinal>::print(std::ostream& os) const 
  {
    using std::endl;
    int const myImageID = comm_->getRank();
    int const numImages = comm_->getSize();
    for (int i = 0; i < numImages; ++i) {
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
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal  ONE = Teuchos::OrdinalTraits<Ordinal>::one();

    // to_nodes_from_me[i] == number of messages sent by this node to node i
    // the info in numSends_, imagesTo_, lengthsTo_ concerns the contiguous sends
    // therefore, each node will be listed in imagesTo_ at most once
    {
      std::vector<Ordinal> to_nodes_from_me(numImages, ZERO);
#     ifdef TEUCHOS_DEBUG 
        bool counting_error = false;
#     endif
      for (int i = ZERO; i < (numSends_ + (selfMessage_ ? 1 : 0)); ++i) {
#       ifdef TEUCHOS_DEBUG
          if (to_nodes_from_me[imagesTo_[i]] != 0) counting_error = true;
#       endif
        to_nodes_from_me[imagesTo_[i]] = ONE;
      }
#     ifdef TEUCHOS_DEBUG
        SHARED_TEST_FOR_EXCEPTION(counting_error, std::logic_error,
            "Tpetra::Distributor::createFromSends: logic error. Please notify the Tpetra team.",*comm_);
#     endif
      // each proc will get back only one item (hence, counts = ones) from the array of globals sums, 
      // namely that entry corresponding to the node, and detailing how many receives it has.
      // this total includes self sends
      std::vector<Ordinal> counts(numImages, 1);
      Teuchos::reduceAllAndScatter<Ordinal>(*comm_,Teuchos::REDUCE_SUM,numImages,&to_nodes_from_me[0],&counts[0],&numReceives_);
    }

    // assign these to length numReceives, with zero entries
    lengthsFrom_.assign(numReceives_, ZERO);
    imagesFrom_.assign(numReceives_, ZERO);

    // FINISH: why do these work? they are blocking sends, and should block until completion, which happens below
    // FINISH: consider switching them to non-blocking
    // NOTE: epetra has both, old (non-blocking) and new (mysterious)

    for (Ordinal i = ZERO; i < (numSends_ + (selfMessage_ ? 1 : 0)); ++i) {
      if (imagesTo_[i] != myImageID ) {
        // send a message to imagesTo_[i], telling him that our pattern sends him lengthsTo_[i] blocks of packets
        Teuchos::send(*comm_,lengthsTo_[i],imagesTo_[i]);
      }
      else {
        // set selfMessage_ to end block of recv arrays
        lengthsFrom_[numReceives_-ONE] = lengthsTo_[i];
        imagesFrom_[numReceives_-ONE] = myImageID;
      }
    }

    //
    for (Ordinal i = ZERO; i < (numReceives_ - (selfMessage_ ? 1 : 0)); ++i) {
      // receive one Ordinal variable from any sender.
      // store the value in lengthsFrom_[i], and store the sender's ImageID in imagesFrom_[i]
      // imagesFrom_[i] = comm_->receive(&lengthsFrom_[i], 1, -1);
      imagesFrom_[i] = Teuchos::receive(*comm_,-1,&lengthsFrom_[i]);
    }
    comm_->barrier();

    sortArrays(imagesFrom_(), lengthsFrom_());

    // Compute indicesFrom_
    totalReceiveLength_ = std::accumulate(lengthsFrom_.begin(), lengthsFrom_.end(), ZERO);
    indicesFrom_.clear();
    indicesFrom_.reserve(totalReceiveLength_);
    for (Ordinal i = 0; i < totalReceiveLength_; ++i) {
      indicesFrom_.push_back(i);
    }

    startsFrom_.clear();
    startsFrom_.reserve(numReceives_);
    for (Ordinal i = ZERO, j = ZERO; i < numReceives_; ++i) {
      startsFrom_.push_back(j);
      j += lengthsFrom_[i];
    }

    if (selfMessage_) --numReceives_;

    comm_->barrier();
  }

  template <typename Ordinal>
  void Distributor<Ordinal>::computeSends(
      const Teuchos::ArrayView<const Ordinal> & importIDs,
      const Teuchos::ArrayView<const Ordinal> & importImageIDs,
            Teuchos::ArrayRCP<Ordinal>& exportIDs,
            Teuchos::ArrayRCP<Ordinal>& exportImageIDs)
  {
    int myImageID = comm_->getRank();
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();

    Ordinal numImports = importImageIDs.size();
    Teuchos::Array<Ordinal> importObjs(2*numImports);
    for (Ordinal i = ZERO; i < numImports; ++i ) {  
      importObjs[2*i]   = importIDs[i];
      importObjs[2*i+1] = myImageID;
    }

    Teuchos_Ordinal numExports;
    Distributor<Ordinal> tempPlan(comm_);
    tempPlan.createFromSends(importImageIDs, numExports);
    exportIDs = Teuchos::arcp<Ordinal>(numExports);
    exportImageIDs = Teuchos::arcp<Ordinal>(numExports);

    Teuchos::Array<Ordinal> exportObjs(tempPlan.getTotalReceiveLength()*2);
    tempPlan.doPostsAndWaits(importObjs().getConst(),2,exportObjs());

    for (Teuchos_Ordinal i = 0; i < numExports; ++i) {
      exportIDs[i]      = exportObjs[2*i];
      exportImageIDs[i] = exportObjs[2*i+1];
    }
  }

} // namespace Tpetra

#endif // TPETRA_DISTRIBUTOR_HPP
