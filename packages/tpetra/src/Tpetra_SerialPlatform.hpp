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

#include <Teuchos_RCP.hpp>
#include <Teuchos_Object.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
#include "Tpetra_Platform.hpp"

namespace Tpetra {

	//! Tpetra::SerialPlatform: Serial Implementation of the Platform class.
	template<typename Ordinal>
	class SerialPlatform : public virtual Platform<Ordinal> {
	public:

    //! @name Constructor/Destructor Methods
		//@{ 

		//! Constructor
		SerialPlatform();

		//! Copy constructor
		SerialPlatform(const SerialPlatform<Ordinal> & platform);

		//! Destructor
		~SerialPlatform();

		//! Clone constructor
		Teuchos::RCP< Platform<Ordinal> > clone() const;

		//@}

    //! @name Class Creation and Accessor Methods
		//@{ 

		//! Comm Instance
		Teuchos::RCP< Teuchos::Comm<Ordinal> > createComm() const;

		//@}

	}; // SerialPlatform class

  template <typename Ordinal>
  SerialPlatform<Ordinal>::SerialPlatform() 
    : Platform<Ordinal>("Tpetra::SerialPlatform<"+Teuchos::TypeNameTraits<Ordinal>::name()+">") 
  {}

  template <typename Ordinal>
  SerialPlatform<Ordinal>::SerialPlatform(const SerialPlatform<Ordinal> & /*platform*/) 
    : Platform<Ordinal>("Tpetra::SerialPlatform<"+Teuchos::TypeNameTraits<Ordinal>::name()+">") 
  {}

  template <typename Ordinal>
  SerialPlatform<Ordinal>::~SerialPlatform() 
  {}

  template <typename Ordinal>
  Teuchos::RCP< Platform<Ordinal> > 
  SerialPlatform<Ordinal>::clone() const 
  {
    Teuchos::RCP< Platform<Ordinal> > platform;
    platform = Teuchos::rcp(new SerialPlatform<Ordinal>(*this));
    return platform;
  }

  template <typename Ordinal>
  Teuchos::RCP< Teuchos::Comm<Ordinal> > 
  SerialPlatform<Ordinal>::createComm() const 
  {
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
  }

} // namespace Tpetra

#endif // TPETRA_SERIALPLATFORM_HPP
