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

#ifndef _TPETRA_MPICOMM_HPP_
#define _TPETRA_MPICOMM_HPP_

#include "Tpetra_Object.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_PacketTraits.hpp"
#include <mpi.h>

namespace Tpetra {

	template<typename PacketType, typename OrdinalType>
	class MpiComm : public Object, public virtual Comm<PacketType, OrdinalType> {
	public:
    
    //@{ \name Constructor/Destructor Methods
		MpiComm(MPI_Comm Comm) 
      : Object("Tpetra::Comm[MPI]") 
      , MpiComm_(Comm)
    {
        MPI_Comm_size(Comm, &size_);
        MPI_Comm_rank(Comm, &rank_);
    };

    MpiComm(MpiComm<PacketType, OrdinalType> const& comm) 
      : Object(comm.label())
      , MpiComm_(comm.MpiComm_)
      , size_(comm.size_)
      , rank_(comm.rank_)
    {};

		~MpiComm() {};
    //@}
    
    //@{ \name Image Info Methods
    int getMyImageID() const {return(rank_);};
    int getNumImages() const {return(size_);};
    //@}
    
    //@{ \name Barrier Methods
    void barrier() const {
      MPI_Barrier(MpiComm_);
    };
    //@}
    
    //@{ \name Broadcast Methods
    void broadcast(PacketType* myVals, OrdinalType const count, int const root) const {
      MPI_Bcast(myVals, count, PacketTraits<PacketType>::mpiDataType(), root, MpiComm_);
    };
    //@}
    
    //@{ \name Gather Methods
    void gatherAll(PacketType* myVals, PacketType* allVals, OrdinalType const count) const {
      MPI_Allgather(myVals, count, PacketTraits<PacketType>::mpiDataType(), allVals, count, PacketTraits<PacketType>::mpiDataType(), MpiComm_);
    };
    //@}
    
    //@{ \name Sum Methods
    void sumAll(PacketType* partialSums, PacketType* globalSums, OrdinalType const count) const {
      MPI_Allreduce(partialSums, globalSums, count, PacketTraits<PacketType>::mpiDataType(), MPI_SUM, MpiComm_);
    };
    //@}
    
    //@{ \name Max/Min Methods
    void maxAll(PacketType* partialMaxs, PacketType* globalMaxs, OrdinalType const count) const {
      MPI_Allreduce(partialMaxs, globalMaxs, count, PacketTraits<PacketType>::mpiDataType(), MPI_MAX, MpiComm_);
    };
    void minAll(PacketType* partialMins, PacketType* globalMins, OrdinalType const count) const {
      MPI_Allreduce(partialMins, globalMins, count, PacketTraits<PacketType>::mpiDataType(), MPI_MIN, MpiComm_);
    };
    //@}
    
    //@{ \name Parallel Prefix Methods
    void scanSum(PacketType* myVals, PacketType* scanSums, OrdinalType const count) const {
      MPI_Scan(myVals, scanSums, count, PacketTraits<PacketType>::mpiDataType(), MPI_SUM, MpiComm_);
    };
    //@}
    
    //@{ \name I/O Methods
    //! Print methods
    void print(ostream& os) const {os << "Image " << getMyImageID() << " of " << getNumImages() << " total images." << endl;};
    void printInfo(ostream& os) const {print(os);};
    //@}
    
    //@{ \name MPI-specific methods, not inherited from Tpetra::Comm
    //! Access method to the MPI Communicator we're using.
    MPI_Comm getMpiComm() const {
      return(MpiComm_);
    };
    
    //@}
    
private:
    MPI_Comm MpiComm_; // The MPI Communicator passed in at construction (actually a new one cpy ctr'd from it).
    int rank_;
    int size_;
    
	}; // MpiComm class
  
} // namespace Tpetra

#endif // _TPETRA_MPICOMM_HPP_
