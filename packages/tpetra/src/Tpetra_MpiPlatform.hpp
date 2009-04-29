// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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

#include <Teuchos_DefaultMpiComm.hpp>
#include "Tpetra_Platform.hpp"

namespace Tpetra {

	//! \brief A implementation of the Platform class for MPI-based platforms.
	template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal>
  class MpiPlatform : public virtual Platform<Scalar, LocalOrdinal, GlobalOrdinal> 
  {
  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructor
    MpiPlatform(const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm);

    //! Constructor with object label
    MpiPlatform(const std::string &label);

    //! Destructor
    ~MpiPlatform();

    //! Clone constructor - implements Tpetra::Platform clone() method.
    Teuchos::RCP< Platform<Scalar,LocalOrdinal,GlobalOrdinal> > clone() const;

    //@}

    //! @name Class Creation and Accessor Methods
    //@{ 

    //! Comm Instance
    Teuchos::RCP< Teuchos::Comm<int> > getComm() const;

    //@}

    private:
    Teuchos::RCP<Teuchos::MpiComm<int> > comm_;
    MpiPlatform(const MpiPlatform<Scalar,LocalOrdinal,GlobalOrdinal> &platform);

  };

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  MpiPlatform<Scalar,LocalOrdinal,GlobalOrdinal>::MpiPlatform(const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm)
  {
    comm_ = Teuchos::createMpiComm<int>(rawMpiComm);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  MpiPlatform<Scalar,LocalOrdinal,GlobalOrdinal>::MpiPlatform(const std::string &label) 
  : Teuchos::LabeledObject(label) 
  {
    comm_ = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));
  } 

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  MpiPlatform<Scalar,LocalOrdinal,GlobalOrdinal>::~MpiPlatform() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  MpiPlatform<Scalar,LocalOrdinal,GlobalOrdinal>::MpiPlatform(const MpiPlatform<Scalar,LocalOrdinal,GlobalOrdinal> &platform)
  {
    comm_ = platform.comm_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP< Platform<Scalar,LocalOrdinal,GlobalOrdinal> > 
  MpiPlatform<Scalar,LocalOrdinal,GlobalOrdinal>::clone() const 
  {
    return Teuchos::rcp(new MpiPlatform<Scalar,LocalOrdinal,GlobalOrdinal>(comm_->getRawMpiComm()));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP< Teuchos::Comm<int> > 
  MpiPlatform<Scalar,LocalOrdinal,GlobalOrdinal>::getComm() const 
  {
    return comm_;
  }

} // namespace Tpetra

#endif // TPETRA_MPIPLATFORM_HPP

