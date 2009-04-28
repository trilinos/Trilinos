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

#ifndef TPETRA_SERIALPLATFORM_HPP
#define TPETRA_SERIALPLATFORM_HPP

#include <Teuchos_DefaultSerialComm.hpp>
#include "Tpetra_Platform.hpp"

namespace Tpetra {

	//! Tpetra::SerialPlatform: Serial Implementation of the Platform class.
	template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal>
	class SerialPlatform : public virtual Platform<Scalar, LocalOrdinal, GlobalOrdinal> 
  {
  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructor
    SerialPlatform();

    //! Constructor with object label
    SerialPlatform(const std::string &label);

    //! Destructor
    ~SerialPlatform();

    //! Clone constructor - implements Tpetra::Platform clone() method.
    Teuchos::RCP< Platform<Scalar,LocalOrdinal,GlobalOrdinal> > clone() const;

    //@}

    //! @name Class Creation and Accessor Methods
    //@{ 

    //! Comm Instance
    Teuchos::RCP< Teuchos::Comm<int> > getComm() const;

    //@}

  };

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  SerialPlatform<Scalar,LocalOrdinal,GlobalOrdinal>::SerialPlatform() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  SerialPlatform<Scalar,LocalOrdinal,GlobalOrdinal>::SerialPlatform(const std::string &label) 
  : Teuchos::LabeledObject(label) {} 

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  SerialPlatform<Scalar,LocalOrdinal,GlobalOrdinal>::~SerialPlatform() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP< Platform<Scalar,LocalOrdinal,GlobalOrdinal> > 
  SerialPlatform<Scalar,LocalOrdinal,GlobalOrdinal>::clone() const 
  {
    return Teuchos::rcp(new SerialPlatform<Scalar,LocalOrdinal,GlobalOrdinal>());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP< Teuchos::Comm<int> > 
  SerialPlatform<Scalar,LocalOrdinal,GlobalOrdinal>::getComm() const 
  {
    return Teuchos::rcp(new Teuchos::SerialComm<int>() );
  }

} // namespace Tpetra

#endif // TPETRA_SERIALPLATFORM_HPP
