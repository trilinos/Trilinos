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

#ifndef TPETRA_PLATFORM_HPP
#define TPETRA_PLATFORM_HPP

#include "Tpetra_ConfigDefs.hpp"
#include <Teuchos_RefCountPtr.hpp>
#include "Tpetra_Comm.hpp"
#include "Tpetra_Directory.hpp"
#include "Tpetra_Distributor.hpp"

namespace Tpetra {

template<typename OrdinalType> class ElementSpace;
// Comm, Directory, and Distributor are not forward declared because they are used as return types. 

//! Tpetra::Platform: The Tpetra Platform Abstract Base Class
/*! Platform is an abstract base class. It should never be called directly.
		Rather, an implementation of Platform, such as SerialPlatform, should be used instead.
		Platform is used to generate Comm, Distributor, and Directory instances. It also manages 
		platform-specific information, such as how inter-image communication is implemented.
		An implementation of Platform, such as SerialPlatform, will create corresponding classes,
		such as SerialComm and SerialDistributor. These will then be cast to their base class,
		and passed back to other Tpetra modules. As a result, other Tpetra modules don't need to know
		anything about the platform they're running on, or any implementation-specific details.

    NOTE: Methods that return a new object (such as clone, createComm, etc.) return them 
    encapsulated in a Teuchos RefCountPtr object. This is done whenever the new object is
    allocated on the heap.
*/

template<typename OrdinalType, typename ScalarType>
class Platform {
public:
	
	//@{ \name Constructor/Destructor Methods
	//! Destructor
	virtual ~Platform() {};
	//! Clone method
  /*! Returns a copy of this Platform instance. It is allocated on the heap and
      encapsulated in a Teuchos RefCountPtr.
  */
	virtual Teuchos::RefCountPtr< Platform<OrdinalType, ScalarType> > clone() const = 0;
	//@}

	//@{ \name Image Info Methods
	//! getMyImageID - returns my rank on this machine
  /*! ImageIDs are always in the range [0, numImages), and are returned as an int,
   */
	virtual int getMyImageID() const = 0;
	//! getNumImages - returns the number of images on this machine
  /*! The number of images on this machine is returned as an int, and should always be greater than zero.
   */
	virtual int getNumImages() const = 0;
	//@}
	
	//@{ \name Class Creation and Accessor Methods
	//! Comm Instances
	virtual Teuchos::RefCountPtr< Comm<ScalarType, OrdinalType> > createScalarComm() const = 0;
	virtual Teuchos::RefCountPtr< Comm<OrdinalType, OrdinalType> > createOrdinalComm() const = 0;
	//! Distributor Instances
	virtual Teuchos::RefCountPtr< Distributor<ScalarType, OrdinalType> > createScalarDistributor() const = 0;
	virtual Teuchos::RefCountPtr< Distributor<OrdinalType, OrdinalType> > createOrdinalDistributor() const = 0;
	//! Directory Instance
	virtual Teuchos::RefCountPtr< Directory<OrdinalType> > createDirectory(ElementSpace<OrdinalType> const& elementSpace) const = 0;
	//@}
	
	//@{ \name I/O Methods
	//! printInfo
	virtual void printInfo(ostream& os) const = 0;
	//@}
	
}; // Platform class
	
} // namespace Tpetra

#endif // TPETRA_PLATFORM_HPP
