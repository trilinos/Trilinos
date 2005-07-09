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
	
		//@{ \name Constructor/Destructor Methods
		
		//! Constructor (serial)
		OmniPlatform() 
			: Object("Tpetra::OmniPlatform(Serial)")
			, myImageID_(0)
			, numImages_(1)
		{}

#ifdef TPETRA_MPI
		//! Constructor (MPI)
		OmniPlatform(MPI_Comm Comm)
			: Object("Tpetra::OmniPlatform(MPI)")
			, myImageID_(-99)
			, numImages_(-99)
			, MpiData_()
		{
			MpiData_ = Teuchos::rcp(new MpiData(Comm));
			MPI_Comm_size(MPI_COMM_WORLD, &numImages_);
			MPI_Comm_rank(MPI_COMM_WORLD, &myImageID_);
		}
#endif

		//! Copy constructor
		OmniPlatform(OmniPlatform const& rhs)
			: Object(rhs.label())
			, myImageID_(rhs.getMyImageID())
			, numImages_(rhs.getNumImages())
#ifdef TPETRA_MPI
			, MpiData_(rhs.MpiData_)
#endif
		{}

		//! Destructor
		~OmniPlatform() {}

		//@}

		//@{ \name Image Info Methods

		//! getMyImageID - returns my rank on this machine
		/*! ImageIDs are always in the range [0, numImages), and are returned as an int.
		 */
		int getMyImageID() const { 
			return(myImageID_); 
		}

		//! getNumImages - returns the number of images on this machine
		/*! The number of images on this machine is returned as an int, and should always be greater than zero.
		 */
		int getNumImages() const { 
			return(numImages_); 
		}

		//@}
	
		//@{ \name Class Creation and Accessor Methods

		//! Comm instances
		/*! Returns a Tpetra::Comm object that will communicate over whatever
		    system was passed to this Platform at construction (Serial, MPI, etc.)
		*/
		template <typename PacketType, typename OrdinalType>
		Teuchos::RefCountPtr< Comm<PacketType, OrdinalType> > createComm() const {
#ifdef TPETRA_MPI
			// We are running in MPI, create a Tpetra::MpiComm.
			Teuchos::RefCountPtr< MpiComm<PacketType, OrdinalType> > comm;
			comm = Teuchos::rcp(new MpiComm<PacketType, OrdinalType>(MpiData_));
			return(comm);
#else
			// We are running in serial, create a Tpetra::SerialComm.
			Teuchos::RefCountPtr< SerialComm<PacketType, OrdinalType> > comm;
			comm = Teuchos::rcp(new SerialComm<PacketType, OrdinalType>());
			return(comm);
#endif
		}
		
		//! SerialComm instances
		/*! This will always return a Tpetra::Comm object that will communicate serially,
		    regardless of what communication system we are using. (In other words, even
			if createComm() returns a Tpetra::MpiComm object, createSerialComm() 
			will still return a Tpetra::SerialComm object.
		*/
		template <typename PacketType, typename OrdinalType>
		Teuchos::RefCountPtr< Comm<PacketType, OrdinalType> > createSerialComm() const {
			Teuchos::RefCountPtr< SerialComm<PacketType, OrdinalType> > comm;
			comm = Teuchos::rcp(new SerialComm<PacketType, OrdinalType>());
			return(comm);
		}

		//@}
	
		//@{ \name I/O Methods

		//! print - implements Tpetra::Object virtual print method.
		void print(ostream& os) const {
			os << label() << ",Image " << getMyImageID() << " of " << getNumImages();
		}

		//@}

	private:
		int myImageID_;
		int numImages_;
#ifdef TPETRA_MPI
		Teuchos::RefCountPtr<MpiData> MpiData_;
#endif
	
	}; // OmniPlatform class
	
} // namespace Tpetra

#endif // TPETRA_OMNIPLATFORM_HPP
