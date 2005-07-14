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

#ifndef TPETRA_OMNIPLATFORM_HPP
#define TPETRA_OMNIPLATFORM_HPP

#include "Tpetra_ConfigDefs.hpp"
#include <Teuchos_RefCountPtr.hpp>
#include "Tpetra_Object.hpp"
#include "Tpetra_SerialComm.hpp"
#ifdef TPETRA_MPI
#include "Tpetra_MpiComm.hpp"
#include "Tpetra_MpiData.hpp"
#endif

namespace Tpetra {

	//! Tpetra::OmniPlatform - Experimental new implementation of the Tpetra Platform class.

	class OmniPlatform : public Object {
	public:

		enum CommType { GENERIC, SERIAL, MPI };
	
		//@{ \name Constructor/Destructor Methods
		
		//! Constructor (serial)
		OmniPlatform() 
			: Object("Tpetra::OmniPlatform(Serial)")
		{}

#ifdef TPETRA_MPI
		//! Constructor (MPI)
		OmniPlatform(MPI_Comm Comm)
			: Object("Tpetra::OmniPlatform(MPI)")
			, MpiData_()
		{
			MpiData_ = Teuchos::rcp(new MpiData(Comm));
		}
#endif

		//! Copy constructor
		OmniPlatform(OmniPlatform const& rhs)
			: Object(rhs.label())
#ifdef TPETRA_MPI
			, MpiData_(rhs.MpiData_)
#endif
		{}

		//! Destructor
		~OmniPlatform() {}

		//@}
	
		//@{ \name Class Creation and Accessor Methods

		//! Comm instances
		/*! Creates a Tpetra::Comm object that will communicate over whatever
		    system was passed to this Platform at construction (Serial, MPI, etc.)
			If you would like a specific type of Comm object, specify a CommType
			argument (i.e. "SERIAL", "MPI", etc.)

			The new Comm instance will be pointed to by the RefCountPtr passed in.
		*/
		template <typename PacketType, typename OrdinalType>
		void createComm(Teuchos::RefCountPtr< Comm<PacketType, OrdinalType> >& comm, CommType ct = GENERIC) const {
			switch(ct) {

			case SERIAL:
				createSerialComm(comm);
				break;

			case MPI:
				createMpiComm(comm);
				break;

			case GENERIC:
#ifdef TPETRA_MPI
				createMpiComm(comm);
#else
				createSerialComm(comm);
#endif
				break;

			default:
				throw reportError("Unknown CommType", -99);
				break;
			}
		}

		//@}
	
		//@{ \name I/O Methods

		//! print - implements Tpetra::Object virtual print method.
		//void print(ostream& os) const {
		//	os << label();
		//}

		//@}

	private:
		// private member functions - used by createComm

		template <typename PacketType, typename OrdinalType>
		void createSerialComm(Teuchos::RefCountPtr< Comm<PacketType, OrdinalType> >& comm) const {
			// serial is always enabled - no need for an #ifdef
			comm = Teuchos::rcp(new SerialComm<PacketType, OrdinalType>());
		}

		template <typename PacketType, typename OrdinalType>
		void createMpiComm(Teuchos::RefCountPtr< Comm<PacketType, OrdinalType> >& comm) const {
#ifdef TPETRA_MPI
			comm = Teuchos::rcp(new MpiComm<PacketType, OrdinalType>(MpiData_));
#else
			throw reportError("Mpi is not enabled.", -1);
#endif
		}

		// private data members

#ifdef TPETRA_MPI
		Teuchos::RefCountPtr<MpiData> MpiData_;
#endif
	
	}; // OmniPlatform class
	
} // namespace Tpetra

#endif // TPETRA_OMNIPLATFORM_HPP
