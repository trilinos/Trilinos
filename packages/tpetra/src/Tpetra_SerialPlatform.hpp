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

#ifndef TPETRA_SERIALPLATFORM_HPP
#define TPETRA_SERIALPLATFORM_HPP

#include <Teuchos_RefCountPtr.hpp>
#include "Tpetra_Object.hpp"
#include "Tpetra_Platform.hpp"
#include "Tpetra_SerialComm.hpp"
#include "Tpetra_BasicDirectory.hpp"
#include "Tpetra_SerialDistributor.hpp"

namespace Tpetra {

// forward definition
template<typename OrdinalType> class ElementSpace;

	//! Tpetra::SerialPlatform: Serial Implementation of the Platform class.

 template<typename OrdinalType, typename ScalarType>
	class SerialPlatform : public Object, public virtual Platform<OrdinalType, ScalarType> {
	public:

		//@{ \name Constructor/Destructor Methods
		//! Constructor
		SerialPlatform() : Object("Tpetra::SerialPlatform") {};
		//! Copy constructor
		SerialPlatform(SerialPlatform<OrdinalType, ScalarType> const& platform) : Object(platform.label()) {};
		//! Destructor
		~SerialPlatform() {};
		//! Clone constructor
		Teuchos::RefCountPtr< Platform<OrdinalType, ScalarType> > clone() const {
      Teuchos::RefCountPtr< Platform<OrdinalType, ScalarType> > platform;
      platform = Teuchos::rcp(new SerialPlatform<OrdinalType, ScalarType>(*this));
			return(platform);
		};
		//@}

    //@{ \name Image Info Methods

    //! getMyImageID - In serial mode, always returns 0.
    int getMyImageID() const {return(0);};

    //! getNumImages - In serial mode, always returns 1.
    int getNumImages() const {return(1);};

    //@}

		//@{ \name Class Creation and Accessor Methods

		//! Comm Instances
		Teuchos::RefCountPtr< Comm<ScalarType, OrdinalType> > createScalarComm() const {
			Teuchos::RefCountPtr< SerialComm<ScalarType, OrdinalType> > comm;
      comm = Teuchos::rcp(new SerialComm<ScalarType, OrdinalType>());
			return(comm);
		};
		Teuchos::RefCountPtr< Comm<OrdinalType, OrdinalType> > createOrdinalComm() const {
			Teuchos::RefCountPtr< SerialComm<OrdinalType, OrdinalType> > comm;
      comm = Teuchos::rcp(new SerialComm<OrdinalType, OrdinalType>());
			return(comm);
		};

		//! Distributor Instance
		Teuchos::RefCountPtr< Distributor<OrdinalType> > createDistributor() const {
			Teuchos::RefCountPtr< SerialDistributor<OrdinalType> > distributor;
      distributor = Teuchos::rcp(new SerialDistributor<OrdinalType>());
			return(distributor);
		};

		//! Directory Instance
		Teuchos::RefCountPtr< Directory<OrdinalType> > createDirectory(ElementSpace<OrdinalType> const& elementSpace) const {
      Teuchos::RefCountPtr< BasicDirectory<OrdinalType> > directory;
      directory = Teuchos::rcp(new BasicDirectory<OrdinalType>(elementSpace)); 
			return(directory);
		};

		//@}

		//@{ \name I/O Methods
		//! print - implements Tpetra::Object virtual print method.
		void print(ostream& os) const {};

		//! printInfo - implements Tpetra::Platform virtual printInfo method.
		void printInfo(ostream& os) const {os << *this;};
		//@}

	}; // SerialPlatform class

} // namespace Tpetra

#endif // TPETRA_SERIALPLATFORM_HPP
