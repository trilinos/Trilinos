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

#ifndef TPETRA_MPICOMM_HPP
#define TPETRA_MPICOMM_HPP

#include <mpi.h>
#include <Teuchos_RefCountPtr.hpp>
#include "Tpetra_Object.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_MpiTraits.hpp"
#include "Tpetra_PacketTraits.hpp"
#include "Tpetra_MpiData.hpp"
#include "Tpetra_Distributor.hpp"

namespace Tpetra {

	//! Tpetra::MpiComm:  The Tpetra MPI Communication Class.
	/*! The MpiComm class is an implementation of Tpetra::Comm, providing the general
	    information and services needed for other Tpetra classes to run on 
		a parallel computer using MPI.

		MpiComm error codes (positive for non-fatal, negative for fatal):
		<ol>
		<li> +1  Pointer to buffer is null.
		<li> +2  Specified source ImageID does not exist.
		<li> +3  Specified destination ImageID does not exist.
		<li> -99 Internal MpiComm error.  Contact developer.
		</ol>
	*/

	template<typename PacketType, typename OrdinalType>
	class MpiComm : public Object, public virtual Comm<PacketType, OrdinalType> {
	public:
    
		//@{ \name Constructor/Destructor Methods

		//! default constuctor
		/*! This is used when the user wants to create an MpiComm instance directly.
		  \param Comm In
		         MPI_Comm communicator we will use.
		*/
		MpiComm(MPI_Comm Comm) 
			: Object("Tpetra::MpiComm") 
			, MpiData_()
			, tag_(-1)
		{
			MpiData_ = Teuchos::rcp(new MpiData(Comm));
			tag_ = data().getMpiTag();
		}

		//! platform constructor
		/*! This is used by MpiPlatform to create an MpiComm instance. It should
		  not be called directly by the user.
		  \param MpiData In
		         MpiData inner data class passed in by MpiPlatform.
		*/
		MpiComm(Teuchos::RefCountPtr<MpiData> const& mpidata)
			: Object("Tpetra::MpiComm")
			, MpiData_(mpidata)
			   //, tag_(-1)
			, tag_(data().getMpiTag())
		{
			//tag_ = data().getMpiTag();
		}

		//! copy constructor
		MpiComm(MpiComm<PacketType, OrdinalType> const& comm) 
			: Object(comm.label())
			, MpiData_(comm.MpiData_)
			   //, tag_(-1)
			, tag_(data().getMpiTag())
		{
			tag_ = data().getMpiTag();
		}

		~MpiComm() {};

		//@}
    
		//@{ \name Barrier Methods

		//! MpiComm barrier function.
		/*! Causes each image in the communicator to wait until all images have arrived. */
		void barrier() const {
			MPI_Barrier(getMpiComm());
		}

		//@}
    
		//@{ \name Broadcast Methods

		//! MpiComm Broadcast function.
		/*!Takes list of input values from the root image and sends to all other images.
		  \param myVals InOut
		         On entry, the root image contains the list of values.  On exit,
				 all images will have the same list of values.  Note that values must be
				 allocated on all images before the broadcast.
		  \param count In
		         On entry, contains the length of myVals.
		  \param root In
		         On entry, contains the image from which all images will receive a copy of myVals.
		*/
		void broadcast(PacketType* myVals, OrdinalType const count, int const root) const {
			MPI_Bcast(myVals, MpiTraits<PacketType>::count(count), MpiTraits<PacketType>::datatype(), root, getMpiComm());
		}

		//@}
    
		//@{ \name Gather Methods

		//! MpiComm All Gather function.
		/*! Takes list of values from all images in the communicator and creates an ordered contiguous list of 
		    those values on each image.
		  \param myVals In
		         On entry, contains the list of values, to be sent to all images.
		  \param allVals Out
		         On exit, contains the list of values from all images. Must by of size numImages*count.
		  \param count In
		         On entry, contains the length of myVals.
		*/
		void gatherAll(PacketType* myVals, PacketType* allVals, OrdinalType const count) const {
			MPI_Allgather(myVals, MpiTraits<PacketType>::count(count), 
						  MpiTraits<PacketType>::datatype(), 
						  allVals, MpiTraits<PacketType>::count(count), 
						  MpiTraits<PacketType>::datatype(), getMpiComm());
		}

		//@}
    
		//@{ \name Sum Methods

		//! MpiComm Global Sum function.
		/*! Takes list of input values from all images in the communicator, computes the sum and returns the
		    sum to all images.
		  \param partialSums In
		         On entry, contains the list of values, usually partial sums computed locally,
				 to be summed across all images.
		  \param globalSums Out
		         On exit, contains the list of values summed across all images.
		  \param count In
		         On entry, contains the length of partialSums.
		*/
		void sumAll(PacketType* partialSums, PacketType* globalSums, OrdinalType const count) const {
			MPI_Allreduce(partialSums, globalSums, MpiTraits<PacketType>::count(count), 
						  MpiTraits<PacketType>::datatype(), MpiTraits<PacketType>::sumOp(), getMpiComm());
		}

		//! Scattered Global Sum function
		/*! Take a list of input values from each image, compute the global sums, and distribute 
		    the list of sums across all images.
		  \param sendVals In
		         On entry, contains the list of values to sum from this image.
		  \param recvVals Out
		         On exit, contains the list of sums distributed to this image.
		  \param recvCounts In
		         On entry, contains a list of sizes. On exit, recvVals on image i will contain
				 recvCounts[i] entries.
		*/
		void sumAllAndScatter(PacketType* sendVals, PacketType* recvVals, OrdinalType* recvCounts) const {
			int const numImages = getNumImages();
			std::vector<int> mpiRecvCounts(numImages, 0);
			for(int i = 0; i < numImages; i++)
				mpiRecvCounts[i] = MpiTraits<PacketType>::count(recvCounts[i]);

			MPI_Reduce_scatter(sendVals, recvVals, &mpiRecvCounts.front(), MpiTraits<PacketType>::datatype(),
							   MpiTraits<PacketType>::sumOp(), getMpiComm());
		}

		//@}
    
		//@{ \name Max/Min Methods

		//! MpiComm Global Max function.
		/*! Takes list of input values from all images in the communicator, computes the max and returns the 
		    max to all images.
		  \param partialMaxs In
		         On entry, contains the list of values to compute the max for.
		  \param globalMaxs Out
		         On exit, contains the list of max values across all images.
		  \param count In
		         On entry, contains the length of partialMaxs.
		*/
		void maxAll(PacketType* partialMaxs, PacketType* globalMaxs, OrdinalType const count) const {
			MPI_Allreduce(partialMaxs, globalMaxs, MpiTraits<PacketType>::count(count), 
						  MpiTraits<PacketType>::datatype(), MpiTraits<PacketType>::maxOp(), getMpiComm());
		}

		//! MpiComm Global Min function.
		/*! Takes list of input values from all images in the communicator, computes the min and returns the 
		    min to all images.
		  \param partialMins In
		         On entry, contains the list of values to compute the min for.
		  \param globalMins Out
		         On exit, contains the list of min values across all images.
		  \param count In
		         On entry, contains the length of partialMins.
		*/
		void minAll(PacketType* partialMins, PacketType* globalMins, OrdinalType const count) const {
			MPI_Allreduce(partialMins, globalMins, MpiTraits<PacketType>::count(count), 
						  MpiTraits<PacketType>::datatype(), MpiTraits<PacketType>::minOp(), getMpiComm());
		}

		//@}
    
		//@{ \name Parallel Prefix Methods

		//! MpiComm Scan Sum function.
		/*! Takes list of input values from all images in the communicator, computes the scan sum and returns it
		    to all images such that image i receives the sum of values from image 0 up to and including image i.
		  \param myVals In
		         On entry, contains the list of values to be summed across all images.
		  \param scanSums Out
		         On exit, contains the list of values summed across images 0 through i.
		  \param count In
		         On entry, contains the length of myVals.
		*/
		void scanSum(PacketType* myVals, PacketType* scanSums, OrdinalType const count) const {
			MPI_Scan(myVals, scanSums, MpiTraits<PacketType>::count(count), 
					 MpiTraits<PacketType>::datatype(), MpiTraits<PacketType>::sumOp(), getMpiComm());
		}

		//@}

		//@{ \name Point-to-Point Methods

		//! MpiComm Blocking Send function
		/*! Sends a list of values to another image. Control will not exit this function until
		    the data have been acknowledged by the receiving image.
		  \param myVals In
		         On entry, contains the list of values to be sent.
		  \param count In
		         On entry, contains the length of myVals.
		  \param destImageID In
		         On entry, contains the ImageID of the image to send the values to.
		*/
		void send(PacketType* myVals, OrdinalType const count, int destImageID) const {
			// Throw an exception if myVals is null.
			if(myVals == 0)
				throw reportError("myVals is null.", 1);
			// Throw an exception if destImageID is not a valid ID.
			if(destImageID < 0 || destImageID >= data().getNumImages())
				throw reportError("Invalid destImageID specified. Should be in the range [0," +
								  toString(data().getNumImages()) + ").", 3);
			
			int err = MPI_Send(myVals, MpiTraits<PacketType>::count(count), MpiTraits<PacketType>::datatype(),
							   destImageID, tag_, data().getMpiComm());
			if(err != 0)
				cerr << "MpiComm error on image " << data().getMyImageID() << ", code = " << err << endl;
			
		}

		//! MpiComm Blocking Receive function
		/*! Receive a list of values from another image. Control will not exit this function until
		    data has been received from the source image specified.
		  \param myVals Out
			     On exit, contains the list of values received.
		  \param count In
		         On entry, contains the length (capacity) of myVals.
		  \param sourceImageID In
		         On entry, contains the ImageID of the image to receive the values from.
				 (A sourceImageID of -1 means receive values from any image.)
		  \return If successful, the ImageID of the sender (>= 0).
		          If not successful, aminteger return code (< 0).
		*/
		int receive(PacketType* myVals, OrdinalType const count, int sourceImageID) const {
			// Throw an exception if myVals is null.
			if(myVals == 0)
				throw reportError("myVals is null.", 1);
			// Throw an exception if sourceImageID is not a valid ID.
			if(sourceImageID < -1 || sourceImageID >= data().getNumImages())
				throw reportError("Invalid sourceImageID specified. Should be in the range [-1," +
								  toString(data().getNumImages()) + ").", 2);
			
			// Because MPI_ANY_SOURCE is an MPI constant, we have the caller pass 
			// the generic value -1 instead. Here we convert that to the value MPI expects.
			if(sourceImageID == -1)
				sourceImageID = MPI_ANY_SOURCE;

			MPI_Status status; // A dummy MPI_Status object, needed for the MPI_Recv call.

			MPI_Recv(myVals, MpiTraits<PacketType>::count(count), MpiTraits<PacketType>::datatype(),
							   sourceImageID, tag_, data().getMpiComm(), &status);

			if(status.MPI_ERROR != 0) {
				cerr << "MpiComm error on image " << data().getMyImageID() << ", code = " << status.MPI_ERROR << endl;
				return(-1);
			}

			return(status.MPI_SOURCE); // if still here, successful. return sender's ImageID.
		}

		//@}

		//@{ \name Execute Distributor Plan Methods

		//! doPostsAndWaits
		/*! Execute a plan specified by the distributor object passed in. 
		  \param distributor In
			     Contains the specifications of the plan we're executing.
		  \param exports In
		         On entry, contains the values we're exporting.
		  \param packetSize In
		         On entry, the number of PacketType variables that make up an element.
		  \param imports Out
		         On exit, contains the values exported to us. (imports will be resized
				 if necessary, and any existing values will be overwritten.)
		*/
		void doPostsAndWaits(Distributor<OrdinalType> const& distributor,
							 std::vector<PacketType>& exports,
							 OrdinalType packetSize,
							 std::vector<PacketType>& imports) {
			doPosts(distributor, exports, packetSize, imports);
			doWaits(distributor);
		}

		//! doPosts
		void doPosts(Distributor<OrdinalType> const& distributor,
					 std::vector<PacketType>& exports,
					 OrdinalType elementSize,
					 std::vector<PacketType>& imports) {
			// setup references to Distributor's data
			OrdinalType const& totalReceiveLength = distributor.getTotalReceiveLength();
			OrdinalType const& numReceives = distributor.getNumReceives();
			OrdinalType const& selfMessage = distributor.getSelfMessage();
			OrdinalType const& numSends = distributor.getNumSends();
			OrdinalType const& maxSendLength = distributor.getMaxSendLength();
			std::vector<OrdinalType> const& imagesFrom = distributor.getImagesFrom();
			std::vector<OrdinalType> const& lengthsFrom = distributor.getLengthsFrom();
			std::vector<OrdinalType> const& imagesTo = distributor.getImagesTo();
			std::vector<OrdinalType> const& indicesTo = distributor.getIndicesTo();
			std::vector<OrdinalType> const& startsTo = distributor.getStartsTo();
			std::vector<OrdinalType> const& lengthsTo = distributor.getLengthsTo();

			// start of actual doPosts function
			OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
			OrdinalType const one = Teuchos::OrdinalTraits<OrdinalType>::one();
			OrdinalType const myImageID = getMyImageID();
			OrdinalType selfReceiveAddress = zero;

			// allocate space in imports if needed
			if(static_cast<OrdinalType>(imports.size()) < (totalReceiveLength * elementSize))
				imports.resize(totalReceiveLength * elementSize);

			// allocate space in requests array if needed
			if(static_cast<OrdinalType>(request_.size()) < numReceives)
				request_.resize(numReceives);

			// start up the Irecv's
			for(OrdinalType i = zero, j = zero, k = zero; i < (numReceives + selfMessage); i++) {
				if(imagesFrom[i] != myImageID) { // receiving this one from another image
					irecv(imports, j, (lengthsFrom[i] * elementSize), imagesFrom[i], request_[k]);
					k++;
				}
				else // receiving this one from ourself
					selfReceiveAddress = j;
				j += (lengthsFrom[i] * elementSize);
			}
			barrier();

			// setup scan through imagesTo list starting with higher numbered images
			// (should help balance message traffic)
			OrdinalType numBlocks = numSends+ selfMessage;
			OrdinalType imageIndex = zero;
			while((imageIndex < numBlocks) && (imagesTo[imageIndex] < myImageID))
				imageIndex++;
			if(imageIndex == numBlocks)
				imageIndex = zero;

			OrdinalType selfNum = zero;
			OrdinalType selfIndex = zero;

			if(indicesTo.empty()) { // data is already blocked by processor
				for(OrdinalType i = zero; i < numBlocks; i++) {
					OrdinalType p = i + imageIndex;
					if(p > (numBlocks - one))
						p -= numBlocks;

					if(imagesTo[p] != myImageID) // sending it to another image
						rsend(exports, (startsTo[p] * elementSize), (lengthsTo[p] * elementSize),
							  imagesTo[p]);
					else // sending it to ourself
						selfNum = p;
				}
				
				if(selfMessage > zero)
					memcopy(exports, (startsTo[selfNum] * elementSize), (lengthsTo[selfNum] * elementSize),
							imports, selfReceiveAddress);
			}
			else { // data is not blocked by image, use send buffer
				std::vector<PacketType> sendArray(maxSendLength * elementSize); // allocate sendArray buffer

				for(OrdinalType i = zero, j = zero; i < numBlocks; i++) {
					OrdinalType p = i + imageIndex;
					if(p > (numBlocks - one)) 
						p -= numBlocks;

					if(imagesTo[p] != myImageID) { // sending it to another image
						OrdinalType offset = zero;
						j = startsTo[p];
						for(OrdinalType k = zero; k < lengthsTo[p]; k++) {
							memcopy(exports, (indicesTo[j] * elementSize), elementSize, sendArray, offset);
							j++;
							offset += elementSize;
						}
						rsend(sendArray, zero, (lengthsTo[p] * elementSize), imagesTo[p]);
					}
					else { // sending it to ourself
						selfNum = p;
						selfIndex = startsTo[p];
					}
				}

				if(selfMessage > zero)
					for(OrdinalType k = zero; k < lengthsTo[selfNum]; k++) {
						memcopy(exports, (indicesTo[selfIndex] * elementSize), elementSize, imports, selfReceiveAddress);
						selfIndex++;
						selfReceiveAddress += elementSize;
					}
			
			}
		}

		//! doWaits
		void doWaits(Distributor<OrdinalType> const& distributor) {
			OrdinalType const& numReceives = distributor.getNumReceives();
			std::vector<MPI_Status> status(numReceives);
			if(numReceives > Teuchos::OrdinalTraits<OrdinalType>::zero())
				MPI_Waitall(numReceives, &request_.front(), &status.front());
		}

		//! doReversePostsAndWaits
		/*! Execute a reverse plan specified by the distributor object passed in. 
		  \param distributor In
			     Contains the specifications of the plan we're reverse-executing.
		  \param exports In
		         On entry, contains the values we're exporting.
		  \param packetSize In
		         On entry, the number of PacketType variables that make up an element.
		  \param imports Out
		         On exit, contains the values exported to us. (imports will be resized
				 if necessary, and any existing values will be overwritten.)
		*/
		void doReversePostsAndWaits(Distributor<OrdinalType> const& distributor,
									std::vector<PacketType>& exports,
									OrdinalType packetSize,
									std::vector<PacketType>& imports) 
		{
			doReversePosts(distributor, exports, packetSize, imports);
			doReverseWaits(distributor);
		}

		//! doReversePosts
		void doReversePosts(Distributor<OrdinalType> const& distributor,
							std::vector<PacketType>& exports,
							OrdinalType packetSize,
							std::vector<PacketType>& imports)
		{
			if(!distributor.getIndicesTo().empty())
				throw reportError("Can only do reverse comm when original data is blocked by image", -99);
			
			// Distributor.getReverse() is logically const but not bitwise const.
			Distributor<OrdinalType>& nc_distrib = const_cast<Distributor<OrdinalType>&>(distributor);
			doPosts(nc_distrib.getReverse(), exports, packetSize, imports);
		}

		//! doReverseWaits
		void doReverseWaits(Distributor<OrdinalType> const& distributor)
		{
			// Distributor.getReverse() is logically const but not bitwise const.
			Distributor<OrdinalType>& nc_distrib = const_cast<Distributor<OrdinalType>&>(distributor);
			doWaits(nc_distrib.getReverse());
		}

		//@}

		//@{ \name Image Info Methods

		//! getMyImageID
		/*! returns the rank of the calling image in the MPI communicator we are using. 
		    (Obtained by calling MPI_Comm_rank.)
		 */
		int getMyImageID() const {return(data().getMyImageID());}

		//! getNumImages - returns the MPI size
		/*! returns the size of the MPI communicator we are using. 
		    (Obtained by calling MPI_Comm_size.)
		 */
		int getNumImages() const {return(data().getNumImages());}

		//@}
    
		//@{ \name I/O Methods

		//! Print methods
		void print(ostream& os) const {}
		void printInfo(ostream& os) const {os << *this;}

		//@}
    
		//@{ \name MPI-specific methods, not inherited from Tpetra::Comm

		//! Access method to the MPI Communicator we're using.
		MPI_Comm getMpiComm() const {
			return(data().getMpiComm());
		}
    
		//@}
    
	private:
    
		// convenience functions for returning inner data class, both const and nonconst versions.
		MpiData& data() {return(*MpiData_);}
		MpiData const& data() const {return(*MpiData_);}

		// private data members
		Teuchos::RefCountPtr<MpiData> MpiData_;
		int tag_;
		std::vector<MPI_Request> request_;

		// templated MPI_Rsend functions
		void rsend(PacketType& myVal, int destinationImageID) const {
			// this one's for a single value
			MPI_Rsend(&myVal, MpiTraits<PacketType>::count(1), MpiTraits<PacketType>::datatype(), 
					  destinationImageID, tag_, data().getMpiComm());
		}

		void rsend(std::vector<PacketType>& myVals, OrdinalType startIndex, OrdinalType count, int destinationImageID) const {
			// this one's for a vector or partial vector
			MPI_Rsend(&myVals[startIndex], MpiTraits<PacketType>::count(count), MpiTraits<PacketType>::datatype(), 
					  destinationImageID, tag_, data().getMpiComm());
		}

		// templated MPI_Irecv functions
		void irecv(PacketType& myVal, int sourceImageID, MPI_Request& request) const {
			MPI_Irecv(&myVal, MpiTraits<PacketType>::count(1), MpiTraits<PacketType>::datatype(), 
					  sourceImageID, tag_, data().getMpiComm(), &request);
		}

		void irecv(std::vector<PacketType>& myVals, OrdinalType startIndex, OrdinalType count, 
				   int sourceImageID, MPI_Request& request) const {
			MPI_Irecv(&myVals[startIndex], MpiTraits<PacketType>::count(count), MpiTraits<PacketType>::datatype(), 
					  sourceImageID, tag_, data().getMpiComm(), &request);
		}

		// replacement for previously used memcpy
		// same functionality as std::copy, but takes array offsets instead of iterators
		// source = STL vector to copy from. dest = STL vector to copy into.
		void memcopy(std::vector<PacketType>& source, OrdinalType sourceStartPos, OrdinalType const length,
					 std::vector<PacketType>& dest, OrdinalType destStartPos) {
			OrdinalType sPos = sourceStartPos;
			OrdinalType dPos = destStartPos;
			for(OrdinalType i = Teuchos::OrdinalTraits<OrdinalType>::zero(); i < length; i++) {
				dest[dPos] = source[sPos];
				sPos++;
				dPos++;
			}
		}
		
    
	}; // MpiComm class
  
} // namespace Tpetra

#endif // TPETRA_MPICOMM_HPP
