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

#ifndef TPETRA_MPIPLATFORM_HPP
#define TPETRA_MPIPLATFORM_HPP

#include <mpi.h>
#include <Teuchos_RefCountPtr.hpp>
#include "Tpetra_Object.hpp"
#include "Tpetra_Platform.hpp"
#include "Tpetra_MpiComm.hpp"
#include "Tpetra_BasicDirectory.hpp"
#include "Tpetra_MpiDistributor.hpp"

namespace Tpetra {

// forward declaration of MpiData, needed to prevent circular inclusions
// actual #include statement at the end of this file
class MpiData;

	//! Tpetra::MpiPlatform: MPI Implementation of the Platform class.

	template<typename OrdinalType, typename ScalarType>
	class MpiPlatform : public Object, public virtual Platform<OrdinalType, ScalarType> {
	public:

    //@{ \name Constructor/Destructor Methods

		//! Constructor
		MpiPlatform(MPI_Comm Comm) 
      : Object("Tpetra::MpiPlatform")
      , MpiData_()
    {
      MpiData_ = Teuchos::rcp(new MpiData(Comm));
    };

    //! Copy Constructor
    MpiPlatform(MpiPlatform<OrdinalType, ScalarType> const& platform) 
      : Object(platform.label())
      , MpiData_(platform.MpiData_)
    {};

    //! Destructor
		~MpiPlatform() {};

    //! Clone Constructor - implements Tpetra::Platform virtual clone method.
		Teuchos::RefCountPtr< Platform<OrdinalType, ScalarType> > clone() const {
      Teuchos::RefCountPtr< MpiPlatform<OrdinalType, ScalarType> > platform;
      platform = Teuchos::rcp(new MpiPlatform<OrdinalType, ScalarType>(*this));
      return(platform);
		};

    //@}

    //@{ \name Image Info Methods

    //! getMyImageID
    /*! returns the rank of the calling image in the MPI communicator we are using. (Obtained by calling MPI_Comm_rank.)
     */
    int getMyImageID() const {return(data().getMyImageID());};

    //! getNumImages - returns the MPI size
    /*! returns the size of the MPI communicator we are using. (Obtained by calling MPI_Comm_size.)
     */
    int getNumImages() const {return(data().getNumImages());};

    //@}

    //@{ \name Class Creation and Accessor Methods

		//! Comm Instances
    Teuchos::RefCountPtr< Comm<ScalarType, OrdinalType> > createScalarComm() const {
			Teuchos::RefCountPtr< MpiComm<ScalarType, OrdinalType> > comm;
      comm = Teuchos::rcp(new MpiComm<ScalarType, OrdinalType>(MpiData_));
      return(comm);
		};
		Teuchos::RefCountPtr< Comm<OrdinalType, OrdinalType> > createOrdinalComm() const {
			Teuchos::RefCountPtr< MpiComm<OrdinalType, OrdinalType> > comm;
      comm = Teuchos::rcp(new MpiComm<OrdinalType, OrdinalType>(MpiData_));
      return(comm);
		};
    
		//! Distributor Instance
		Teuchos::RefCountPtr< Distributor<OrdinalType> > createDistributor() const {
		  Teuchos::RefCountPtr< MpiDistributor<OrdinalType> > distributor;
      distributor = Teuchos::rcp(new MpiDistributor<OrdinalType>(MpiData_)); 
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
    
    //@{ \name MPI-specific methods, not inherited from Tpetra::Platform

    //! Access method to the MPI Communicator we're using.
    MPI_Comm getMpiComm() const {
      return(data().MpiComm_);
    };
    
    //@}
    
  private:
    Teuchos::RefCountPtr<MpiData> MpiData_;
    
    // convenience functions for returning inner data class, both const and nonconst versions.
    MpiData& data() {return(*MpiData_);};
    MpiData const& data() const {return(*MpiData_);};
    
	}; // MpiPlatform class
  
} // namespace Tpetra

#include "Tpetra_MpiData.hpp"

#endif // TPETRA_MPIPLATFORM_HPP
