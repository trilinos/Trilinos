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

#include <Teuchos_RefCountPtr.hpp>
#include "Tpetra_Object.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Util.hpp"

namespace Tpetra {
  
	//! Tpetra::MpiDistributor:  The Tpetra MPI implementation of the Tpetra::Distributor Gather/Scatter Setup Class.
	/*! The MpiDistributor class is an MPI implement of Tpetra::Distributor.
	    An MpiDistributor object is actually produced by calling a method in the Tpetra::MpiPlatform class.
	*/
  
	template<typename OrdinalType>
	class MpiDistributor : public Object, public virtual Distributor<OrdinalType> {
	public:
    
		//@{ \name Constructor/Destructor
    
		//! Platform Constuctor (default ctr)
		MpiDistributor(Teuchos::RefCountPtr< Comm<OrdinalType, OrdinalType> > const& comm) 
			: Object("Tpetra::MpiDistributor")
			, Comm_(comm)
			, numExports_(Teuchos::OrdinalTraits<OrdinalType>::zero())
			, selfMessage_(Teuchos::OrdinalTraits<OrdinalType>::zero())
			, numSends_(Teuchos::OrdinalTraits<OrdinalType>::zero())
			, maxSendLength_(Teuchos::OrdinalTraits<OrdinalType>::zero())
			, numReceives_(Teuchos::OrdinalTraits<OrdinalType>::zero())
			, totalReceiveLength_(Teuchos::OrdinalTraits<OrdinalType>::zero())
		{};
    
		//! Copy Constructor
		MpiDistributor(Distributor<OrdinalType> const& distributor) 
			: Object(distributor.label())
			, Comm_(distributor.Comm_)
			, numExports_(distributor.numExports_)
			, selfMessage_(distributor.selfMessage_)
			, numSends_(distributor.numSends_)
			, maxSendLength_(distributor.maxSendLength_)
			, numReceives_(distributor.numReceives_)
			, totalReceiveLength_(distributor.totalReceiveLength_)
		{};

		//! Clone constructor
		Distributor<OrdinalType>* clone() {
			MpiDistributor<OrdinalType>* distributor = new MpiDistributor<OrdinalType>(*this); 
			return(distributor); 
		};
    
		//! Destructor.
		~MpiDistributor() {};

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
      
			OrdinalType myImageID = comm().getMyImageID();
			OrdinalType numImages = comm().getNumImages();
      
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
			computeSends(numRemoteIDs, remoteGIDs, remoteImageIDs, 
						 numExportIDs, exportGIDs, exportImageIDs, comm().getMyImageID());
			OrdinalType testNumRemoteIDs; // dummy-ish variable
			createFromSends(numExportIDs, exportImageIDs, deterministic, testNumRemoteIDs);
		};
    
		//@}

		//@{ \name Attribute Accessor Methods

		OrdinalType const& getTotalReceiveLength() const {
			return(totalReceiveLength_);
		};
		OrdinalType const& getNumReceives() const {
			return(numReceives_);
		};
		OrdinalType const& getSelfMessage() const {
			return(selfMessage_);
		};
		OrdinalType const& getNumSends() const {
			return(numSends_);
		};
		OrdinalType const& getMaxSendLength() const {
			return(maxSendLength_);
		};
		std::vector<OrdinalType> const& getImagesFrom() const {
			return(imagesFrom_);
		};
		std::vector<OrdinalType> const& getLengthsFrom() const {
			return(lengthsFrom_);
		};
		std::vector<OrdinalType> const& getImagesTo() const {
			return(imagesTo_);
		};
		std::vector<OrdinalType> const& getIndicesTo() const {
			return(indicesTo_);
		};
		std::vector<OrdinalType> const& getStartsTo() const {
			return(startsTo_);
		};
		std::vector<OrdinalType> const& getLengthsTo() const {
			return(lengthsTo_);
		};

		//@}
    

		//@{ \name I/O Methods

		//! print method inherited from Object
		void print(ostream& os) const {
			int const myImageID = comm().getMyImageID();
			int const numImages = comm().getNumImages();
			for(int i = 0; i < numImages; i++) {
				comm().barrier();
				if(i == myImageID) {
					cout << "[Image " << myImageID << " of " << numImages << "]" << endl;
					cout << " numExports: " << numExports_ << endl;
					cout << " selfMessage: " << selfMessage_ << endl;
					cout << " numSends_: " << numSends_ << endl;
					cout << " imagesTo_: " << toString(imagesTo_) << endl;
					cout << " startsTo_: " << toString(startsTo_) << endl;
					cout << " lengthsTo_: " << toString(lengthsTo_) << endl;
					cout << " maxSendLength_: " << maxSendLength_ << endl;
					cout << " indicesTo_: " << toString(indicesTo_) << endl;
					cout << " numReceives_: " << numReceives_ << endl;
					cout << " totalReceiveLength_: " << totalReceiveLength_ << endl;
					cout << " lengthsFrom_: " << toString(lengthsFrom_) << endl;
					cout << " imagesFrom_: " << toString(imagesFrom_) << endl;
					cout << " indicesFrom_: " << toString(indicesFrom_) << endl;
					cout << " startsFrom_: " << toString(startsFrom_) << endl;
				}
			}
		};

		//! printInfo method inherited from Distributor
		void printInfo(ostream& os) const {os << *this;};

		//@}
    
	private:

		// convenience functions for returning inner data class, both const and nonconst versions.
		Comm<OrdinalType, OrdinalType>& comm() {return(*Comm_);};
		Comm<OrdinalType, OrdinalType> const& comm() const {return(*Comm_);};

		// private data members
		Teuchos::RefCountPtr< Comm<OrdinalType, OrdinalType> > Comm_;

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

		void computeReceives(OrdinalType myImageID, OrdinalType numImages) {
			OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
			OrdinalType const one = Teuchos::OrdinalTraits<OrdinalType>::one();
			std::vector<OrdinalType> msg_count(numImages, zero);
			std::vector<int> counts(numImages, 1); // int vector for MPI sumAllAndScatter call
      
			for(OrdinalType i = zero; i < (numSends_ + selfMessage_); i++)
				msg_count[imagesTo_[i]] = one;
			
			comm().sumAllAndScatter(&msg_count.front(), &numReceives_, &counts.front());
			lengthsFrom_.assign(numReceives_, zero);
			imagesFrom_.assign(numReceives_, zero);
      
			for(OrdinalType i = zero; i < (numSends_ + selfMessage_); i++)
				if(imagesTo_[i] != myImageID )
					comm().send(&lengthsTo_[i], 1, imagesTo_[i]);
				else {
					// set selfMessage_ to end block of recv arrays
					lengthsFrom_[numReceives_-one] = lengthsTo_[i];
					imagesFrom_[numReceives_-one] = myImageID;
				}
			for(OrdinalType i = zero; i < (numReceives_ - selfMessage_); i++) {
				// receive 1 OrdinalType variable from any sender.
				// store the value in lengthsFrom_[i], and store the sender's ImageID in imagesFrom_[i]
				imagesFrom_[i] = comm().receive(&lengthsFrom_[i], 1, -1);
			}
			comm().barrier();
      
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
      
			comm().barrier();
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

			MpiDistributor<OrdinalType> tempPlan(Comm_);
			std::vector<OrdinalType> imageIDList;
			std::vector<OrdinalType> importObjs;

			imageIDList = importImageIDs; // assumes importImageIDs is of length numImports
			if(static_cast<OrdinalType>(importImageIDs.size()) != numImports) // hence the error-checking
				throw reportError("Internal error in MpiDistributor::computeSends", -99);
      
			importObjs.resize(two * numImports);
      
			for(OrdinalType i = zero; i < numImports; i++ ) {  
				importObjs[two*i] = importIDs[i];
				importObjs[two*i+1] = myImageID;
			}
      
			tempPlan.createFromSends(numImports, imageIDList, true, numExports);
			exportIDs.resize(numExports);
			exportImageIDs.resize(numExports);

			std::vector<OrdinalType> exportObjs;
			Comm_->doPostsAndWaits(tempPlan, importObjs, two, exportObjs);

			for(OrdinalType i = zero; i < numExports; i++) {
				exportIDs[i] = exportObjs[two*i];
				exportImageIDs[i] = exportObjs[two*i+one];
			}
		};
    
	}; // class MpiDistributor
  
} // namespace Tpetra

#endif // TPETRA_MPIDISTRIBUTOR_HPP
