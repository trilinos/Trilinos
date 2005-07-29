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
#include "Tpetra_OmniPlatformData.hpp"
#include "Tpetra_SerialComm.hpp"
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiComm.hpp"
#endif
#ifdef TPETRA_THREADED_MPI
// -- includes for a threaded MPI would go here --
#endif

namespace Tpetra {

	//! Tpetra::OmniPlatform - Experimental new implementation of the Tpetra Platform class.
	/*! Platform is a sort of Factory that generates Tpetra::Comm instances. It is meant
	    to insulate Tpetra classes from the specifics of how communication is being done.

		For example, having been passed a Platform instance at construction, a Tpetra class
		can call Platform.createComm, and receive a Tpetra::Comm object that it can use to do
		global communications with. This may actually be a Tpetra::SerialComm, a Tpetra::MpiComm,
		or some other implementation of the Tpetra::Comm interface.
		
		If they desire, a caller can also create a Comm instance of a specific type. For instance,
		they may want a SerialComm object. This can be done by calling Platform::createComm with an
		additional argument specifying the type of Comm they want. This is an enumerable type (CommType),
		with options "SERIAL", "MPI", "THREADED_MPI".

		To determine if a specific type of Comm is available, the user can call Platform::isEnabled.
		Given a CommType as an argument, it will return a boolean specifying if that CommType is enabled
		in this Platform.

		It is also possible to find out what the "generic" CommType is, by calling the 
		Platform::getGeneric method.
	*/

	class OmniPlatform : public Object {
	public:

		enum CommType { GENERIC, SERIAL, MPI, THREADED_MPI };
	
		//@{ \name Constructor/Destructor Methods
		
		//! Constructor (serial)
		OmniPlatform() 
			: Object("Tpetra::OmniPlatform(Serial)")
			, OmniPlatformData_()
		{
			OmniPlatformData_ = Teuchos::rcp(new OmniPlatformData());
		}

#ifdef TPETRA_MPI
		//! Constructor (MPI)
		OmniPlatform(MPI_Comm Comm)
			: Object("Tpetra::OmniPlatform(MPI)")
		    , OmniPlatformData_()
		{
			OmniPlatformData_ = Teuchos::rcp(new OmniPlatformData(Comm));
		}
#endif

#ifdef TPETRA_THREADED_MPI
		// -- A constructor for a threaded MPI would go here --
#endif

		//! Copy constructor
		OmniPlatform(OmniPlatform const& rhs)
			: Object(rhs.label())
			, OmniPlatformData_(rhs.OmniPlatformData_)
		{}

		//! Destructor
		~OmniPlatform() {}

		//@}

		//@{ \name Attribute Accessor Methods

		//! isEnabled - Returns true if the specified CommType is enabled in the Platform.
		bool isEnabled(CommType const& ct) const {
			bool enabled = false;

			switch (ct) {
			case GENERIC:
				enabled = isEnabled(getGeneric());
				break;		
			case SERIAL:
				enabled = data().serialEnabled_;
				break;
			case MPI:
				enabled = data().mpiEnabled_;
				break;
			case THREADED_MPI:
				enabled = data().threadedMpiEnabled_;
				break;
			default:
				throw reportError("Unknown CommType", -99);
				break;
			}

			return(enabled);
		}

		//! getGeneric - Returns the CommType of the generic communications library we're using.
		/*! This is how we decide which type of Comm to create if more than one is enabled. */
		CommType getGeneric() const {
			if(isEnabled(THREADED_MPI))
				return(THREADED_MPI);
			else if(isEnabled(MPI))
				return(MPI);
			else
				return(SERIAL);
		}

		//@}
	
		//@{ \name Class Creation and Accessor Methods

		//! Comm instances
		/*! Creates a Tpetra::Comm object that will communicate over whatever
		    system was passed to this Platform at construction (Serial, MPI, etc.)
			If you would like a specific type of Comm object, specify a CommType
			argument (i.e. "SERIAL", "MPI", etc.)

			The new Comm instance will be pointed to by the RefCountPtr passed in.
		*/
		template <typename OrdinalType, typename ScalarType>
		void createComm(Teuchos::RefCountPtr< Comm<OrdinalType, ScalarType> >& comm, CommType ct = GENERIC) const {
			if(!isEnabled(ct))
				throw reportError("That CommType is not enabled", -1);

			switch(ct) {

			case SERIAL:
				createSerialComm(comm);
				break;

			case MPI:
				createMpiComm(comm);
				break;

			case THREADED_MPI:
				createThreadedMpiComm(comm);
				break;

			case GENERIC:
				createComm(comm, getGeneric());
				break;

			default:
				throw reportError("Unknown CommType", -99);
				break;
			}
		}

		//@}
	
		//@{ \name I/O Methods

		//! print - implements Tpetra::Object virtual print method.
		void print(ostream& os) const {
			// ..
		}

		//@}

	private:
		// private member functions - used by createComm

		template <typename OrdinalType, typename ScalarType>
		void createSerialComm(Teuchos::RefCountPtr< Comm<OrdinalType, ScalarType> >& comm) const {
			// serial is always enabled - no need for an #ifdef
			comm = Teuchos::rcp(new SerialComm<OrdinalType, ScalarType>());
		}

		template <typename OrdinalType, typename ScalarType>
		void createMpiComm(Teuchos::RefCountPtr< Comm<OrdinalType, ScalarType> >& comm) const {
#ifdef TPETRA_MPI
			comm = Teuchos::rcp(new MpiComm<OrdinalType, ScalarType>(data().MpiComm_));
#else
			throw reportError("MPI is not enabled.", -1);
#endif
		}

		template <typename OrdinalType, typename ScalarType>
		void createThreadedMpiComm(Teuchos::RefCountPtr< Comm<OrdinalType, ScalarType> >& comm) const {
#ifdef TPETRA_THREADED_MPI
			// -- creation of a Tpetra::ThreadedMpiComm would go here --
#else
			throw reportError("Threaded MPI is not enabled.", -1);
#endif
		}

		// convenience functions for returning inner data class, both const and nonconst versions.
		OmniPlatformData& data() {return(*OmniPlatformData_);};
		OmniPlatformData const& data() const {return(*OmniPlatformData_);};

		// private data members

		Teuchos::RefCountPtr<OmniPlatformData> OmniPlatformData_;
	
	}; // OmniPlatform class
	
} // namespace Tpetra

#endif // TPETRA_OMNIPLATFORM_HPP
