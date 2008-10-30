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

#include <Teuchos_RCP.hpp>
#include <Teuchos_Object.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include "Tpetra_Platform.hpp"

// FINISH: this should have a getComm() method; Comm should then have a clone() method
//         because of MPI tag issues, we want to always use the same Comm, unless the user 
//         specifies a different Comm

namespace Tpetra {

  //! Tpetra::MpiPlatform: MPI Implementation of the Platform class.
  template<typename Ordinal>
  class MpiPlatform : public virtual Platform<Ordinal> {
  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructor
    MpiPlatform(const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm);

    //! Copy Constructor
    MpiPlatform(const MpiPlatform<Ordinal> & platform);

    //! Destructor
    ~MpiPlatform();

    //! Clone Constructor - implements Tpetra::Platform virtual clone method.
    Teuchos::RCP< Platform<Ordinal> > clone() const;

    //@}

    //! @name Class Creation and Accessor Methods
    //@{ 

    //! Comm Instance
    Teuchos::RCP< Teuchos::Comm<Ordinal> > createComm() const;

    //@}

  private:
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > MpiComm_;

  }; // MpiPlatform class

  template <typename Ordinal>
  MpiPlatform<Ordinal>::MpiPlatform(const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm)
  : Platform<Ordinal>("Tpetra::MpiPlatform<"+Teuchos::TypeNameTraits<Ordinal>::name()+">") 
  , MpiComm_(rawMpiComm)
  {}

  template <typename Ordinal>
  MpiPlatform<Ordinal>::MpiPlatform(const MpiPlatform<Ordinal> & platform) 
  : Platform<Ordinal>("Tpetra::MpiPlatform<"+Teuchos::TypeNameTraits<Ordinal>::name()+">") 
  , MpiComm_(platform.MpiComm_)
  {}

  template <typename Ordinal>
  MpiPlatform<Ordinal>::~MpiPlatform() 
  {}

  template <typename Ordinal>
  Teuchos::RCP< Platform<Ordinal> > 
  MpiPlatform<Ordinal>::clone() const 
  {
    Teuchos::RCP< MpiPlatform<Ordinal> > platform;
    platform = Teuchos::rcp(new MpiPlatform<Ordinal>(*this));
    return platform;
  }

  template <typename Ordinal>
  Teuchos::RCP< Teuchos::Comm<Ordinal> > 
  MpiPlatform<Ordinal>::createComm() const 
  {
    return Teuchos::createMpiComm<Ordinal>(MpiComm_);
  }

} // namespace Tpetra

#endif // TPETRA_MPIPLATFORM_HPP

