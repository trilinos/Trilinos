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
			, tag_(26050)
		{
			MpiData_ = Teuchos::rcp(new MpiData(Comm));
			//tag_ = data().getMpiTag();
		};

		//! platform constructor
		/*! This is used by MpiPlatform to create an MpiComm instance. It should
		  not be called directly by the user.
		  \param MpiData In
		         MpiData inner data class passed in by MpiPlatform.
		*/
		MpiComm(Teuchos::RefCountPtr<MpiData> const& mpidata)
			: Object("Tpetra::MpiComm")
			, MpiData_(mpidata)
			, tag_(26050)
		{
			//tag_ = data().getMpiTag();
		};

		//! copy constructor
		MpiComm(MpiComm<PacketType, OrdinalType> const& comm) 
			: Object(comm.label())
			, MpiData_(comm.MpiData_)
			, tag_(26050)
		{
			//tag_ = data().getMpiTag();
		};

		~MpiComm() {};

		//@}
    
		//@{ \name Barrier Methods

		//! MpiComm barrier function.
		/*! Causes each image in the communicator to wait until all images have arrived. */
		void barrier() const {
			MPI_Barrier(getMpiComm());
		};

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
		};

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
		};

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
		};

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
		};

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
		};

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
		};

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
		void doPostsAndWaits(Distributor<OrdinalType>& distributor,
							 std::vector<PacketType>& exports,
							 OrdinalType packetSize,
							 std::vector<PacketType>& imports) {
			doPosts(distributor, exports, packetSize, imports);
			doWaits(distributor);
		}

		//! doPosts
		void doPosts(Distributor<OrdinalType>& distributor,
					 std::vector<PacketType>& exports,
					 OrdinalType packetSize,
					 std::vector<PacketType>& imports) {

		}

		//! doWaits
		void doWaits(Distributor<OrdinalType>& distributor) {

		}

		//@}

		//@{ \name Image Info Methods

		//! getMyImageID
		/*! returns the rank of the calling image in the MPI communicator we are using. 
		    (Obtained by calling MPI_Comm_rank.)
		 */
		int getMyImageID() const {return(data().getMyImageID());};

		//! getNumImages - returns the MPI size
		/*! returns the size of the MPI communicator we are using. 
		    (Obtained by calling MPI_Comm_size.)
		 */
		int getNumImages() const {return(data().getNumImages());};

		//@}
    
		//@{ \name I/O Methods

		//! Print methods
		void print(ostream& os) const {};
		void printInfo(ostream& os) const {os << *this;};

		//@}
    
		//@{ \name MPI-specific methods, not inherited from Tpetra::Comm

		//! Access method to the MPI Communicator we're using.
		MPI_Comm getMpiComm() const {
			return(data().getMpiComm());
		};
    
		//@}
    
	private:
    
		// convenience functions for returning inner data class, both const and nonconst versions.
		MpiData& data() {return(*MpiData_);};
		MpiData const& data() const {return(*MpiData_);};

		// private data members
		Teuchos::RefCountPtr<MpiData> MpiData_;
		int tag_;

		// templated MPI_Rsend functions
		int rsend(PacketType& myVal, int destinationImageID) const {
			MPI_Rsend(&myVal, MpiTraits<PacketType>::count(1), MpiTraits<PacketType>::datatype(), 
					  destinationImageID, tag_, data().getMpiComm());
		}

		int rsend(std::vector<PacketType>& myVals, int destinationImageID) const {
			MPI_Rsend(&myVals.front(), MpiTraits<PacketType>::count(myVals.size()), MpiTraits<PacketType>::datatype(), 
					  destinationImageID, tag_, data().getMpiComm());
		}

		// templated MPI_Irecv functions
		int irecv(PacketType& myVal, int sourceImageID) const {
			MPI_Request request;
			MPI_Irecv(&myVal, MpiTraits<PacketType>::count(1), MpiTraits<PacketType>::datatype(), 
					  sourceImageID, tag_, data().getMpiComm(), &request);
		}

		int irecv(std::vector<PacketType>& myVals, int sourceImageID) const {
			MPI_Request request;
			MPI_Irecv(&myVals.front(), MpiTraits<PacketType>::count(myVals.size()), MpiTraits<PacketType>::datatype(), 
					  sourceImageID, tag_, data().getMpiComm(), &request);
		}
    
	}; // MpiComm class
  
} // namespace Tpetra

#endif // TPETRA_MPICOMM_HPP
