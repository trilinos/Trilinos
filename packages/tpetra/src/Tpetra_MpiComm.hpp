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

namespace Tpetra {

// forward declaration of MpiData, needed to prevent circular inclusions
// actual #include statement at the end of this file
class MpiData;

	template<typename PacketType, typename OrdinalType>
	class MpiComm : public Object, public virtual Comm<PacketType, OrdinalType> {
	public:
    
    //@{ \name Constructor/Destructor Methods
    //! default constuctor
    /*! This is used when the user wants to create an MpiComm instance directly.
      \param In
      Comm - MPI_Comm communicator we will use.
     */
		MpiComm(MPI_Comm Comm) 
      : Object("Tpetra::MpiComm") 
      , MpiData_()
    {
      MpiData_ = Teuchos::rcp(new MpiData(Comm));
    };

    //! platform constructor
    /*! This is used by MpiPlatform to create an MpiComm instance. It should
      not be called directly by the user.
      \param In
      MpiData - MpiData inner data class passed in by MpiPlatform.
     */
    MpiComm(Teuchos::RefCountPtr<MpiData> const& data)
      : Object("Tpetra::MpiComm")
      , MpiData_(data)
    {};

    //! copy constructor
    MpiComm(MpiComm<PacketType, OrdinalType> const& comm) 
      : Object(comm.label())
      , MpiData_(comm.MpiData_)
    {};

		~MpiComm() {};
    //@}
    
    //@{ \name Barrier Methods
    //! MpiComm barrier function.
    /*! Causes each image in the communicator to wait until all images have arrived.
     */
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
      MPI_Allgather(myVals, MpiTraits<PacketType>::count(count), MpiTraits<PacketType>::datatype(), 
                    allVals, MpiTraits<PacketType>::count(count), MpiTraits<PacketType>::datatype(), getMpiComm());
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
      MPI_Allreduce(partialSums, globalSums, MpiTraits<PacketType>::count(count), MpiTraits<PacketType>::datatype(), MpiTraits<PacketType>::sumOp(), getMpiComm());
    };
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
      MPI_Allreduce(partialMaxs, globalMaxs, MpiTraits<PacketType>::count(count), MpiTraits<PacketType>::datatype(), MpiTraits<PacketType>::maxOp(), getMpiComm());
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
      MPI_Allreduce(partialMins, globalMins, MpiTraits<PacketType>::count(count), MpiTraits<PacketType>::datatype(), MpiTraits<PacketType>::minOp(), getMpiComm());
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
      MPI_Scan(myVals, scanSums, MpiTraits<PacketType>::count(count), MpiTraits<PacketType>::datatype(), MpiTraits<PacketType>::sumOp(), getMpiComm());
    };
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
    
	}; // MpiComm class
  
} // namespace Tpetra

#endif // TPETRA_MPICOMM_HPP
