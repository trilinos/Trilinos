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
#include "Tpetra_MpiDirectory.hpp"
#include "Tpetra_MpiDistributor.hpp"
#include <mpi.h>

namespace Tpetra {

	//! Tpetra::MpiPlatform: MPI Implementation of the Platform class.

	template<typename OrdinalType, typename ScalarType>
	class MpiPlatform : public Object, public virtual Platform<OrdinalType, ScalarType> {
	public:

    //@{ \name Constructor/Destructor Methods

		//! Constructor
		MpiPlatform() : Object("Tpetra::Platform[MPI]") {};

    //! Copy Constructor
    MpiPlatform(MpiPlatform<OrdinalType, ScalarType> const& platform) : Object(platform.label()) {};

    //! Destructor
		~MpiPlatform() {};

    //! Clone Constructor - implements Tpetra::Platform virtual clone method.
		Platform<OrdinalType, ScalarType>* clone() const {
      MpiPlatform<OrdinalType, ScalarType>* platform = new MpiPlatform<OrdinalType, ScalarType>(*this);
      return(platform);
		};

    //@}

    //@{ \name Class Creation and Accessor Methods

		//! Comm Instances
    Comm<ScalarType, OrdinalType>* createScalarComm() const {
			MpiComm<ScalarType, OrdinalType>* comm = new MpiComm<ScalarType, OrdinalType>(MPI_COMM_WORLD);
      return(comm);
		};
		Comm<OrdinalType, OrdinalType>* createOrdinalComm() const {
			MpiComm<OrdinalType, OrdinalType>* comm = new MpiComm<OrdinalType, OrdinalType>(MPI_COMM_WORLD);
      return(comm);
		};

    //! Distributor Instances
		Distributor<ScalarType, OrdinalType>* createScalarDistributor() const {
			MpiDistributor<ScalarType, OrdinalType>* distributor = new MpiDistributor<ScalarType, OrdinalType>();
      return(distributor);
		};
		Distributor<OrdinalType, OrdinalType>* createOrdinalDistributor() const {
			MpiDistributor<OrdinalType, OrdinalType>* distributor = new MpiDistributor<OrdinalType, OrdinalType>();
      return(distributor);
		};

		//! Directory Instance
		Directory<OrdinalType>* createDirectory(ElementSpace<OrdinalType> const& elementSpace) const {
		  MpiDirectory<OrdinalType>* dir = new MpiDirectory<OrdinalType>(elementSpace); 
      return(dir);
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
