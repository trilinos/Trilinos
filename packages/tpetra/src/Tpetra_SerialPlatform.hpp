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

#include "Tpetra_Object.hpp"
#include "Tpetra_Platform.hpp"
#include "Tpetra_SerialComm.hpp"
#include "Tpetra_SerialDirectory.hpp"
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
		SerialPlatform() : Object("Tpetra::Platform[Serial]") {};
		//! Copy constructor
		SerialPlatform(SerialPlatform<OrdinalType, ScalarType> const& platform) : Object(platform.label()) {};
		//! Destructor
		~SerialPlatform() {};
		//! Clone constructor
		Platform<OrdinalType, ScalarType>* clone() const {
			Platform<OrdinalType, ScalarType>* platform = static_cast<Platform<OrdinalType, ScalarType>*>
				(new SerialPlatform<OrdinalType, ScalarType>(*this));
			return(platform);
		};
		//@}

		//@{ \name Class Creation and Accessor Methods

		//! Comm Instances
		Comm<ScalarType, OrdinalType>* createScalarComm() const {
			// static_cast casts SerialComm* to Comm*
			Comm<ScalarType, OrdinalType>* comm = static_cast<Comm<ScalarType, OrdinalType>*>(new SerialComm<ScalarType, OrdinalType>());
			return(comm);
		};
		Comm<OrdinalType, OrdinalType>* createOrdinalComm() const {
			// static_cast casts SerialComm* to Comm*
			Comm<OrdinalType, OrdinalType>* comm = static_cast<Comm<OrdinalType, OrdinalType>*>(new SerialComm<OrdinalType, OrdinalType>());
			return(comm);
		};

		//! Distributor Instances
		Distributor<ScalarType, OrdinalType>* createScalarDistributor() const {
			// static_cast casts SerialDistributor* to Distributor*
			Distributor<ScalarType, OrdinalType>* distributor = 
				static_cast<Distributor<ScalarType, OrdinalType>*>(new SerialDistributor<ScalarType, OrdinalType>());
			return(distributor);
		};
		Distributor<OrdinalType, OrdinalType>* createOrdinalDistributor() const {
			// static_cast casts SerialDistributor* to Distributor*
			Distributor<OrdinalType, OrdinalType>* distributor = 
				static_cast<Distributor<OrdinalType, OrdinalType>*>(new SerialDistributor<OrdinalType, OrdinalType>());
			return(distributor);
		};

		//! Directory Instance
		Directory<OrdinalType>* createDirectory(ElementSpace<OrdinalType> const& elementSpace) const {
		  // static_cast casts SerialDirectory* to Directory*
		  Directory<OrdinalType>* dir = static_cast<Directory<OrdinalType>*>(new SerialDirectory<OrdinalType>(elementSpace)); 
			return(dir);
		};

		//@}

		//@{ \name I/O Methods
		//! print - implements Tpetra::Object virtual print method.
		void print(ostream& os) const { os << label() << endl;};

		//! printInfo - implements Tpetra::Platform virtual printInfo method.
		void printInfo(ostream& os) const {print(os);};
		//@}

	}; // SerialPlatform class

} // namespace Tpetra

#endif // TPETRA_SERIALPLATFORM_HPP
