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

namespace Tpetra {

	template<typename PacketType, typename OrdinalType>
	class MpiComm : public Object, public virtual Comm<PacketType, OrdinalType> {
	public:
    
    //@{ \name Constructor/Destructor Methods
		MpiComm() : Object("Tpetra::Comm[MPI]") {};

    MpiComm(MpiComm<PacketType, OrdinalType> const& comm) {};

		~MpiComm() {};
    //@}
    
    //@{ \name Image Info Methods
    int getMyImageID() const {return(-1);};
    int getNumImages() const {return(-1);};
    //@}
    
    //@{ \name Barrier Methods
    void barrier() const {};
    //@}
    
    //@{ \name Broadcast Methods
    void broadcast(PacketType* myVals, OrdinalType const count, int const root) const {};
    //@}
    
    //@{ \name Gather Methods
    void gatherAll(PacketType* myVals, PacketType* allVals, OrdinalType const count) const {};
    //@}
    
    //@{ \name Sum Methods
    void sumAll(PacketType* partialSums, PacketType* globalSums, OrdinalType const count) const {};
    //@}
    
    //@{ \name Max/Min Methods
    void maxAll(PacketType* partialMaxs, PacketType* globalMaxs, OrdinalType const count) const {};
    void minAll(PacketType* partialMins, PacketType* globalMins, OrdinalType const count) const {};
    //@}
    
    //@{ \name Parallel Prefix Methods
    void scanSum(PacketType* myVals, PacketType* scanSums, OrdinalType const count) const {};
    //@}
    
    //@{ \name I/O Methods
    //! Print methods
    void print(ostream& os) const {os << "MpiComm print function." << endl;};
    void printInfo(ostream& os) const {print(os);};
    //@}
    
	}; // MpiComm class
  
} // namespace Tpetra

#endif // _TPETRA_MPICOMM_HPP_
