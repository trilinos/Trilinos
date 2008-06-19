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
#include <Teuchos_RCP.hpp>
#include <Teuchos_Object.hpp>
#include "Tpetra_Platform.hpp"
#include "Tpetra_MpiComm.hpp"
#include "Tpetra_ElementSpace.hpp"

namespace Tpetra {

  //! Tpetra::MpiPlatform: MPI Implementation of the Platform class.

  template<typename OrdinalType, typename ScalarType>
  class MpiPlatform : public Teuchos::Object, public virtual Platform<OrdinalType, ScalarType> {
  public:

    //@{ \name Constructor/Destructor Methods

    //! Constructor
    MpiPlatform(MPI_Comm Comm) 
      : Teuchos::Object("Tpetra::MpiPlatform")
      , MpiComm_(Comm)
    {
    }

    //! Copy Constructor
    MpiPlatform(MpiPlatform<OrdinalType, ScalarType> const& platform) 
      : Teuchos::Object(platform.label())
      , MpiComm_(platform.MpiComm_)
    {}

    //! Destructor
    ~MpiPlatform() {}

    //! Clone Constructor - implements Tpetra::Platform virtual clone method.
    Teuchos::RCP< Platform<OrdinalType, ScalarType> > clone() const {
      Teuchos::RCP< MpiPlatform<OrdinalType, ScalarType> > platform;
      platform = Teuchos::rcp(new MpiPlatform<OrdinalType, ScalarType>(*this));
      return(platform);
    }

    //@}

    //@{ \name Class Creation and Accessor Methods

    //! Comm Instances
    Teuchos::RCP< Comm<OrdinalType, ScalarType> > createScalarComm() const {
      Teuchos::RCP< MpiComm<OrdinalType, ScalarType> > comm;
      comm = Teuchos::rcp(new MpiComm<OrdinalType, ScalarType>(MpiComm_));
      return(comm);
    };
    Teuchos::RCP< Comm<OrdinalType, OrdinalType> > createOrdinalComm() const {
      Teuchos::RCP< MpiComm<OrdinalType, OrdinalType> > comm;
      comm = Teuchos::rcp(new MpiComm<OrdinalType, OrdinalType>(MpiComm_));
      return(comm);
    };

    //@}

    //@{ \name I/O Methods

    //! print - implements Teuchos::Object virtual print method.
    void print(ostream& os) const {};

    //! printInfo - implements Tpetra::Platform virtual printInfo method.
    void printInfo(ostream& os) const {os << *this;};

    //@}

  private:
    MPI_Comm MpiComm_;

  }; // MpiPlatform class

} // namespace Tpetra

#endif // TPETRA_MPIPLATFORM_HPP

