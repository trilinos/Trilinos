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

#ifndef _TPETRA_MPIPLATFORM_HPP_
#define _TPETRA_MPIPLATFORM_HPP_

#include "Tpetra_Object.hpp"
#include "Tpetra_Platform.hpp"
#include "Tpetra_MpiComm.hpp"

namespace Tpetra {

	//! Tpetra::MpiPlatform: Mpi Implementation of the Platform class.

	template<typename OrdinalType, typename ScalarType>
	class MpiPlatform : public Object, public virtual Platform<OrdinalType, ScalarType> {
	public:

		//@{ \name Constructor/Destructor Methods

		//! Constructor
		MpiPlatform() : Object("Tpetra::Platform[MPI]") {};
		//! Destructor
		~MpiPlatform() {};
		
		//@}
		
		//@{ \name Class Creation and Accessor Methods

		//! Comm Instances
		Comm<ScalarType, OrdinalType>* createScalarComm() const {
			// static_cast casts MpiComm* to Comm*
			Comm<ScalarType, OrdinalType>* comm = static_cast<Comm<ScalarType, OrdinalType>*>(new MpiComm<ScalarType, OrdinalType>());
			return(comm);
		};
		Comm<OrdinalType, OrdinalType>* createOrdinalComm() const {
			// static_cast casts MpiComm* to Comm*
			Comm<OrdinalType, OrdinalType>* comm = static_cast<Comm<OrdinalType, OrdinalType>*>(new MpiComm<OrdinalType, OrdinalType>());
			return(comm);
		};

		//@}

		//@{ \name I/O Methods
		//! print - implements Tpetra::Object virtual print method.
		void print(ostream& os) const { os << label() << endl;};

		//! printInfo - implements Tpetra::Platform virtual printInfo method.
		void printInfo(ostream& os) const {print(os);};
		//@}

	}; // MpiPlatform class

} // namespace Tpetra

#endif // _TPETRA_MPIPLATFORM_HPP_
