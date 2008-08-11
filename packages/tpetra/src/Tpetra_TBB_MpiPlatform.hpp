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

#ifndef TPETRA_TBB_MPIPLATFORM_HPP
#define TPETRA_TBB_MPIPLATFORM_HPP

#include <mpi.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Object.hpp>
#include "Tpetra_Platform.hpp"

#ifdef HAVE_TPETRA_TBB

#include "Tpetra_TBB_MpiComm.hpp"
#include "Tpetra_MpiData.hpp"
#include "Tpetra_ElementSpace.hpp"

namespace Tpetra {

  //! Tpetra::TBB_MpiPlatform: MPI Implementation of the Platform class, that can
  //  produce instances of the TBB-capable TBB_MpiComm object.

  template<typename OrdinalType, typename ScalarType>
  class TBB_MpiPlatform : public Teuchos::Object, public virtual Platform<OrdinalType, ScalarType> {
  public:

    //@{ \name Constructor/Destructor Methods

    //! Constructor
    TBB_MpiPlatform(MPI_Comm Comm, int num_threads = 0) 
      : Teuchos::Object("Tpetra::MpiPlatform")
      , MpiData_()
      , num_threads_(num_threads)
    {
      MpiData_ = Teuchos::rcp(new MpiData(Comm));
    };

    //! Copy Constructor
    TBB_MpiPlatform(TBB_MpiPlatform<OrdinalType, ScalarType> const& platform) 
      : Teuchos::Object(platform.label())
      , MpiData_(platform.MpiData_)
      , num_threads_(platform.num_threads_)
    {};

    //! Destructor
    ~TBB_MpiPlatform() {};

    //! Clone Constructor - implements Tpetra::Platform virtual clone method.
    Teuchos::RCP< Platform<OrdinalType, ScalarType> > clone() const {
      Teuchos::RCP< TBB_MpiPlatform<OrdinalType, ScalarType> > platform;
      platform = Teuchos::rcp(new TBB_MpiPlatform<OrdinalType, ScalarType>(*this));
      return(platform);
    };

    //@}

    //@{ \name Class Creation and Accessor Methods

    //! Comm Instances
    Teuchos::RCP< Comm<OrdinalType, ScalarType> > createScalarComm() const {
      Teuchos::RCP< TBB_MpiComm<OrdinalType, ScalarType> > comm;
      comm = Teuchos::rcp(new TBB_MpiComm<OrdinalType, ScalarType>(MpiData_, num_threads_));
      return(comm);
    };
    Teuchos::RCP< Comm<OrdinalType, OrdinalType> > createOrdinalComm() const {
      Teuchos::RCP< TBB_MpiComm<OrdinalType, OrdinalType> > comm;
      comm = Teuchos::rcp(new TBB_MpiComm<OrdinalType, OrdinalType>(MpiData_, num_threads_));
      return(comm);
    };

    //@}

    //@{ \name I/O Methods

    //! print - implements Teuchos::Object virtual print method.
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
    Teuchos::RCP<MpiData> MpiData_;
    
    // convenience functions for returning inner data class, both const and nonconst versions.
    MpiData& data() {return(*MpiData_);};
    MpiData const& data() const {return(*MpiData_);};
    
    int num_threads_;
  }; // TBB_MpiPlatform class
  
} // namespace Tpetra

#endif //HAVE_TPETRA_TBB
#endif // TPETRA_TBB_MPIPLATFORM_HPP
