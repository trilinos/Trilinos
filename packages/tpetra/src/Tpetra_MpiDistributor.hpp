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

#ifndef TPETRA_MPIDISTRIBUTOR_HPP
#define TPETRA_MPIDISTRIBUTOR_HPP

#include <mpi.h>
#include <Teuchos_RefCountPtr.hpp>
#include "Tpetra_Object.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_MpiTraits.hpp"
#include "Tpetra_Util.hpp"

namespace Tpetra {

  // forward declaration of MpiData, needed to prevent circular inclusions
  // actual #include statement at the end of this file
  class MpiData;
  
  //! Tpetra::MpiDistributor:  The Tpetra MPI implementation of the Tpetra::Distributor Gather/Scatter Setup Class.
  /*! The MpiDistributor class is an MPI implement of Tpetra::Distributor.
		  An MpiDistributor object is actually produced by calling a method in the Tpetra::MpiPlatform class.
  */
  
	template<typename OrdinalType>
	class MpiDistributor : public Object, public virtual Distributor<OrdinalType> {
  public:
    
    //@{ \name Constructor/Destructor
    
    //! Platform Constuctor (default ctr)
		MpiDistributor(Teuchos::RefCountPtr<MpiData> const& data) 
      : Object("Tpetra::MpiDistributor")
      , MpiData_(data)
      , numExports_(Teuchos::OrdinalTraits<OrdinalType>::zero())
      , selfMessage_(Teuchos::OrdinalTraits<OrdinalType>::zero())
      , numSends_(Teuchos::OrdinalTraits<OrdinalType>::zero())
      , maxSendLength_(Teuchos::OrdinalTraits<OrdinalType>::zero())
      , numReceives_(Teuchos::OrdinalTraits<OrdinalType>::zero())
      , request_(0)
      , status_(0)
      , totalReceiveLength_(Teuchos::OrdinalTraits<OrdinalType>::zero())
      , tag_(26000)
    {};
    
    //! Copy Constructor
    MpiDistributor(Distributor<OrdinalType> const& distributor) 
      : Object(distributor.label())
      , MpiData_(distributor.MpiData_)
      , numExports_(distributor.numExports_)
      , selfMessage_(distributor.selfMessage_)
      , numSends_(distributor.numSends_)
      , maxSendLength_(distributor.maxSendLength_)
      , numReceives_(distributor.numReceives_)
      , request_(distributor.request_) // Is this correct?
      , status_(distributor.status_) // Is this correct?
      , totalReceiveLength_(distributor.totalReceiveLength_)
      , tag_(distributor.tag_) // Is this correct?
    {};

    //! Clone constructor
    Distributor<OrdinalType>* clone() {
      MpiDistributor<OrdinalType>* distributor = new MpiDistributor<OrdinalType>(*this); 
      return(distributor); 
    };
    
    //! Destructor.
    ~MpiDistributor() {
      delete[] request_; // C++ guarantees that deleting a null pointer will have no effect.
      delete[] status_;  // So we don't need to check to see that they're != 0.
    };

    //@}
    

    //@{ \name Gather/Scatter Constructors

    //! Create Distributor object using list of ImageIDs to send to
    /*! Take a list of ImageIDs and construct a plan for efficiently scattering to these images.
        Return the number of IDs being sent to me.
        \param numExportIDs In
               Number of IDs that need to be sent from this image.
        \param exportImageIDs In
               List of images that will get the exported IDs (should be of length numExportIDs).
        \param deterministic In
               No Op.
        \param numRemoteIDs Out
               Number of IDs this image will be receiving.
    */
    void createFromSends(OrdinalType const& numExportIDs, 
                         std::vector<OrdinalType> const& exportImageIDs,
                         bool const& deterministic,
                         OrdinalType& numRemoteIDs) 
		{
      OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
      OrdinalType const one = Teuchos::OrdinalTraits<OrdinalType>::one();

      numExports_ = numExportIDs;
      
      OrdinalType myImageID = data().getMyImageID();
      OrdinalType numImages = data().getNumImages();
      
      // Check to see if items are grouped by images without gaps
      // If so, indices_to -> 0
      
      // Setup data structures for quick traversal of arrays
      std::vector<OrdinalType> starts(numImages + one, zero);
      
      OrdinalType numActive = zero;
      bool noSendBuff = true;
      
      for(OrdinalType i = zero; i < numExportIDs; i++) {
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
        for(OrdinalType i = zero; i < numImages; ++i)
          if(starts[i]) 
            ++numSends_;
        
        imagesTo_.resize(numSends_); // can we change these to reserves?
        startsTo_.resize(numSends_);
        lengthsTo_.resize(numSends_);
        
        for(OrdinalType i = zero, index = zero; i < numSends_; ++i) {
          startsTo_[i] = index;
          OrdinalType imageID = exportImageIDs[index];
          imagesTo_[i] = imageID;
          index += starts[imageID];
        }
        
        if(numSends_ > zero)
          sortArrays(imagesTo_, startsTo_);
        
        maxSendLength_ = zero;
        
        for(OrdinalType i = zero; i < numSends_; ++i) {
          OrdinalType imageID = imagesTo_[i];
          lengthsTo_[i] = starts[imageID];
          if((imageID != myImageID) && (lengthsTo_[i] > maxSendLength_))
            maxSendLength_ = lengthsTo_[i];
        }
      }
      else { // not grouped by image, need send buffer and indicesTo_
        if(starts.front() != zero ) 
          numSends_ = one;
          
        for(OrdinalType i = one; i < numImages; i++) {
          if(starts[i] != zero) 
            ++numSends_;
          starts[i] += starts[i-one];
        }
        
        for(OrdinalType i = (numImages - one); i != zero; i--)
          starts[i] = starts[i-one];
        
        starts.front() = zero;
          
        indicesTo_.resize(numActive);
          
        for(OrdinalType i = zero; i < numExportIDs; i++)
          if(exportImageIDs[i] >= zero) {
            indicesTo_[starts[exportImageIDs[i]]] = i;
            ++starts[exportImageIDs[i]];
          }
        
        // Reconstuct starts array to index into indicesTo.
        
        for(OrdinalType i = (numImages - one); i != zero; i--)
          starts[i] = starts[i-1];
        starts.front() = zero;       
        starts[numImages] = numActive;
        
        imagesTo_.resize(numSends_); // can we change these to reserves?
        startsTo_.resize(numSends_);
        lengthsTo_.resize(numSends_);
        
        maxSendLength_ = zero;
        
        for(OrdinalType i = zero, j = 0; i < numImages; i++ )
          if(starts[i+1] != starts[i]) {
            lengthsTo_[j] = starts[i+1] - starts[i];
            startsTo_[j] = starts[i];
            if((i != myImageID) && (lengthsTo_[j] > maxSendLength_))
              maxSendLength_ = lengthsTo_[j];
            imagesTo_[j] = i;
            j++;
          }
      }
      
      numSends_ -= selfMessage_;
      
      // Invert map to see what msgs are received and what length
      computeReceives(myImageID, numImages);
      
      if(numReceives_ > zero) {
        request_ = new MPI_Request[numReceives_];
        status_ = new MPI_Status[numReceives_];
      }
      
      numRemoteIDs = totalReceiveLength_;
    };
    
    //! Create Distributor object using list of Image IDs to receive from
    /*! Take a list of global IDs and construct a plan for efficiently scattering to these images.
        Return the number and list of IDs being sent by me.
        \param numRemoteIDs In
               Number of IDs this image will be receiving.
        \param remoteGIDs In
               List of IDs that this image wants.
        \param remoteImageIDs In
               List of images that will send the remote IDs.
        \param deterministic In
               No Op.
        \param numExportIDs Out
               Number of IDs that need to be sent from this image.
        \param exportGIDs Out
               List of IDs that need to be sent from this image.
        \param exportImageIDs Out
               List of images that will get the exported IDs.
  */
    void createFromRecvs(OrdinalType const& numRemoteIDs, 
                         std::vector<OrdinalType> const& remoteGIDs, 
                         std::vector<OrdinalType> const& remoteImageIDs, 
                         bool const& deterministic, 
                         OrdinalType& numExportIDs, 
                         std::vector<OrdinalType>& exportGIDs, 
                         std::vector<OrdinalType>& exportImageIDs)
		{
      computeSends(numRemoteIDs, remoteGIDs, remoteImageIDs, numExportIDs, exportGIDs, exportImageIDs, data().getMyImageID());
      
      OrdinalType testNumRemoteIDs; // dummy-ish variable
      createFromSends(numExportIDs, exportImageIDs, deterministic, testNumRemoteIDs);
    };
    
    //@}
    

    //@{ \name Execute Gather/Scatter Operations (Constant size objects)

    //! Execute plan on buffer of export objects in a single step
    void doPostsAndWaits(char* export_objs, OrdinalType const& obj_size, OrdinalType& len_import_objs, char* import_objs) {
      doPosts(export_objs, obj_size, len_import_objs, import_objs);
      doWaits();
    };
    
    //! Post buffer of export objects (can do other local work before executing Waits)
    void doPosts(char* export_objs, OrdinalType const& obj_size, OrdinalType& len_import_objs, char* import_objs) {
      OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
      OrdinalType const one = Teuchos::OrdinalTraits<OrdinalType>::one();
      OrdinalType myImageID = data().getMyImageID();
      OrdinalType selfReceiveAddress = zero;
        
      if(len_import_objs < (totalReceiveLength_ * obj_size)) {
        if(len_import_objs) 
          delete [] import_objs;
        len_import_objs = totalReceiveLength_ * obj_size;
        import_objs = new char[len_import_objs];
        for(OrdinalType i = zero; i < len_import_objs; i++) 
          import_objs[i] = 0;
      }
      
      for(OrdinalType i = zero, j = zero, k = zero; i < (numReceives_ + selfMessage_); i++) {
        if(imagesFrom_[i] != myImageID) {
          MPI_Irecv(&(import_objs[j]), (lengthsFrom_[i] * obj_size), MPI_CHAR, imagesFrom_[i],
                    tag_, data().getMpiComm(), &(request_[k]));
          k++;
        }
        else
          selfReceiveAddress = j;
        
        j += (lengthsFrom_[i] * obj_size);
      }
      
      MPI_Barrier(data().getMpiComm());
      
      //setup scan through procs_to list starting w/ higher numbered procs 
      //Should help balance msg traffic
      OrdinalType numBlocks = numSends_ + selfMessage_; 
      OrdinalType imageIndex = zero; 
      while((imageIndex < numBlocks) && (imagesTo_[imageIndex] < myImageID))
        imageIndex++;                    
      if(imageIndex == numBlocks) 
        imageIndex = zero;
      
      OrdinalType selfNum, selfIndex;
      
      if(indicesTo_.empty()) { //data already blocked by processor
        for(OrdinalType i = zero; i < numBlocks; i++) {
          OrdinalType p = i + imageIndex;
          if(p > (numBlocks - one)) 
            p -= numBlocks;
              
          if(imagesTo_[p] != myImageID)
            MPI_Rsend(&export_objs[startsTo_[p]*obj_size], lengthsTo_[p]*obj_size, MPI_CHAR,
                      imagesTo_[p], tag_, data().getMpiComm());
          else
            selfNum = p;
        }
        
        if(selfMessage_ > zero)
          memcpy(&import_objs[selfReceiveAddress], &export_objs[startsTo_[selfNum]*obj_size],
                 lengthsTo_[selfNum]*obj_size);
      }
      else { //data not blocked by proc, use send buffer
        if(sendArraySize_ < (maxSendLength_ * obj_size)) {
          if(sendArraySize_ > zero) 
            delete[] sendArray_;
                sendArraySize_ = maxSendLength_ * obj_size;
                sendArray_ = new char[sendArraySize_];
        }
        
        for(OrdinalType i = zero, j = zero; i < numBlocks; i++ ) {
          OrdinalType p = i + imageIndex;
          if(p > (numBlocks - one)) 
            p -= numBlocks;

          if(imagesTo_[p] != myImageID) {
            OrdinalType offset = zero;
            j = startsTo_[p];
            for(OrdinalType k = zero; k < lengthsTo_[p]; k++) {
              memcpy(&(sendArray_[offset]), &(export_objs[indicesTo_[j]*obj_size]), obj_size);
              j++;
              offset += obj_size;
            }
            MPI_Rsend(sendArray_, lengthsTo_[p] * obj_size, MPI_CHAR, imagesTo_[p], tag_, data().getMpiComm());
          }
          else {
            selfNum = p;
            selfIndex = startsTo_[p];
          }
        }
        
        if(selfMessage_ > zero)
          for(OrdinalType k = zero; k < lengthsTo_[selfNum]; k++) {
            memcpy(&(import_objs[selfReceiveAddress]), &(export_objs[indicesTo_[selfIndex]*obj_size]), obj_size);
            selfIndex++;
            selfReceiveAddress += obj_size;
          }
      }
    };
    
    //! Wait on a set of posts
    void doWaits() {
      if(numReceives_ > Teuchos::OrdinalTraits<OrdinalType>::zero())
        MPI_Waitall(numReceives_, request_, status_);
    };

    //! Execute reverse of plan on buffer of export objects in a single step
    void doReversePostsAndWaits(char* export_objs, OrdinalType const& obj_size, OrdinalType& len_import_objs, char* import_objs) {};
    
    //! Do reverse post of buffer of export objects (can do other local work before executing Waits)
    void doReversePosts(char* export_objs, OrdinalType const& obj_size, OrdinalType& len_import_objs, char* import_objs) {};
    
    //! Wait on a reverse set of posts
    void doReverseWaits() {};
    
    //@}
    

    //@{ \name I/O Methods

    //! print method inherited from Object
    void print(ostream& os) const {};

    //! printInfo method inherited from Distributor
    void printInfo(ostream& os) const {os << *this;};

    //@}
    
  private:
    
    // convenience functions for returning inner data class, both const and nonconst versions.
    MpiData& data() {return(*MpiData_);};
    MpiData const& data() const {return(*MpiData_);};

    // private data members
    Teuchos::RefCountPtr<MpiData> MpiData_;

    OrdinalType numExports_;
    OrdinalType selfMessage_;
    OrdinalType numSends_;
    std::vector<OrdinalType> imagesTo_;
    std::vector<OrdinalType> startsTo_;
    std::vector<OrdinalType> lengthsTo_;
    OrdinalType maxSendLength_;
    std::vector<OrdinalType> indicesTo_;
    OrdinalType numReceives_;
    MPI_Request* request_;
    MPI_Status* status_;
    OrdinalType totalReceiveLength_;
    std::vector<OrdinalType> lengthsFrom_;
    std::vector<OrdinalType> imagesFrom_;
    std::vector<OrdinalType> indicesFrom_;
    std::vector<OrdinalType> startsFrom_;
    int tag_;
    OrdinalType sendArraySize_;
    char* sendArray_;

    void computeReceives(OrdinalType myImageID, OrdinalType numImages) {
      OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
      OrdinalType const one = Teuchos::OrdinalTraits<OrdinalType>::one();
      std::vector<OrdinalType> msg_count(numImages, zero);
      std::vector<OrdinalType> counts(numImages, one);
      
      MPI_Status status; // MPI
      
      for(OrdinalType i = zero; i < (numSends_ + selfMessage_); i++)
        msg_count[imagesTo_[i]] = one;
      
      MPI_Reduce_scatter(&msg_count.front(), &numReceives_, &counts.front(), MpiTraits<OrdinalType>::datatype(), MpiTraits<OrdinalType>::sumOp(), data().getMpiComm()); // MPI

      lengthsFrom_.assign(numReceives_, zero);
      imagesFrom_.assign(numReceives_, zero);
      
      for(OrdinalType i = zero; i < (numSends_ + selfMessage_); i++)
        if(imagesTo_[i] != myImageID )
          MPI_Send(&lengthsTo_[i], 1, MpiTraits<OrdinalType>::datatype(), imagesTo_[i], tag_, data().getMpiComm()); // MPI
        else {
          // set selfMessage_ to end block of recv arrays
          lengthsFrom_[numReceives_-one] = lengthsTo_[i];
          imagesFrom_[numReceives_-one] = myImageID;
        }
      
      for(OrdinalType i = zero; i < (numReceives_ - selfMessage_); i++) {
        MPI_Recv(&lengthsFrom_[i], 1, MpiTraits<OrdinalType>::datatype(), MPI_ANY_SOURCE, tag_, data().getMpiComm(), &status); // MPI
        imagesFrom_[i] = status.MPI_SOURCE; // MPI
      }
      
      MPI_Barrier(data().getMpiComm()); // MPI
      
      sortArrays(imagesFrom_, lengthsFrom_);
      
      // Compute indicesFrom_
      totalReceiveLength_ = std::accumulate(lengthsFrom_.begin(), lengthsFrom_.end(), zero);
      indicesFrom_.clear();
      indicesFrom_.reserve(totalReceiveLength_);
      for(OrdinalType i = 0; i < totalReceiveLength_; i++)
        indicesFrom_.push_back(i);

      
      startsFrom_.reserve(numReceives_);
      for(OrdinalType i = zero, j = zero; i < numReceives_; ++i) {
        startsFrom_.push_back(j);
        j += lengthsFrom_[i];
      }
      
      numReceives_ -= selfMessage_;
      
      MPI_Barrier(data().getMpiComm()); // MPI
    };

    void computeSends(OrdinalType const numImports, 
                      std::vector<OrdinalType> const& importIDs,
                      std::vector<OrdinalType> const& importImageIDs,
                      OrdinalType& numExports,
                      std::vector<OrdinalType>& exportIDs,
                      std::vector<OrdinalType>& exportImageIDs,
                      OrdinalType myImageID) 
    {
      OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
      OrdinalType const one = Teuchos::OrdinalTraits<OrdinalType>::one();
      OrdinalType const two = one + one;

      MpiDistributor<OrdinalType> tempPlan(MpiData_);
      std::vector<OrdinalType> imageIDList;
      std::vector<OrdinalType> importObjs;

      imageIDList = importImageIDs; // assumes importImageIDs is of length numImports
      if(importImageIDs.size() != numImports) // hence the error-checking
        throw reportError("Internal error in MpiDistributor::computeSends", -99);
      
      importObjs.resize(two * numImports);
      
      for(OrdinalType i = zero; i < numImports; i++ ) {  
        importObjs[two*i] = importIDs[i];
        importObjs[two*i+1] = myImageID;
      }
      
      tempPlan.createFromSends(numImports, imageIDList, true, numExports);
      exportIDs.resize(numExports);
      exportImageIDs.resize(numExports);
      
      char* c_export_objs = 0;
      OrdinalType len_c_export_objs = zero;
      tempPlan.doPostsAndWaits(reinterpret_cast<char*>(&importObjs.front()), 
                               (two * sizeof(OrdinalType)),
                               len_c_export_objs,
                               c_export_objs);
      OrdinalType* exportObjs = reinterpret_cast<OrdinalType*>(c_export_objs);

      for(OrdinalType i = zero; i < numExports; i++) {
        exportIDs[i] = exportObjs[two*i];
        exportImageIDs[i] = exportObjs[two*i+one];
      }

      delete[] c_export_objs;
      
    };
    
  }; // class MpiDistributor
  
} // namespace Tpetra

#include "Tpetra_MpiData.hpp"

#endif // TPETRA_MPIDISTRIBUTOR_HPP
