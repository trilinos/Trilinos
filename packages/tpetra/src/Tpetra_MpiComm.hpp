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

#include "Tpetra_Comm.hpp"
#include "Tpetra_Object.hpp"

namespace Tpetra {

	//! Tpetra::MpiComm:  The Tpetra MPI Communication Class.
	/*! The MpiComm class is an implementation of Tpetra::Comm, providing the general
		information and services needed for other Tpetra classes to run on a parallel 
		computer using MPI.
	*/
	
	template<typename PacketType, typename OrdinalType>
	class MpiComm : public Object, public virtual Comm<PacketType, OrdinalType> {
	public:
		//@{ \name Constructor/Destructor Methods
		//! Constructor
		/*! Builds an instance of an MPI communicator.  
		 */
		MpiComm() : Object("Tpetra::Comm[MPI]") {};
		
		//! Destructor.
		/*! Completely deletes a MpiComm object.  
			\warning Note:  All objects that depend on a MpiComm instance 
			should be destroyed prior to calling this function.
		*/
		~MpiComm() {};
		//@}
		
		//@{ \name I/O Methods
		//! print - implements Tpetra::Object virtual print method.
		void print(ostream& os) const { os << label() << endl;};
		
		//! printInfo - implements Tpetra::Platform virtual printInfo method.
		void printInfo(ostream& os) const {print(os);};
		//@}
		
	}; // class MpiComm
	
} // namespace Tpetra

#endif // _TPETRA_MPICOMM_HPP_
