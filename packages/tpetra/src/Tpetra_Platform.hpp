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

#ifndef TPETRA_PLATFORM_HPP
#define TPETRA_PLATFORM_HPP

#include "Tpetra_ConfigDefs.hpp"
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Describable.hpp>

namespace Tpetra {

  //! \brief The Tpetra platform abstract base class.
  /*!
     This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
     The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
     type, if omitted, defaults to the \c LocalOrdinal type.
   */
  template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal>
  class Platform : public Teuchos::Describable {
  public:
  
    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructor
    Platform();

    //! Constructor with object label
    Platform(const std::string &label);

    //! Destructor
    virtual ~Platform() {};

    //! Clone method
    /*! Returns a copy of this Platform instance. It is dynamically allocated and 
        encapsulated in a Teuchos RCP.
    */
    virtual Teuchos::RCP< Platform<Scalar, LocalOrdinal, GlobalOrdinal> > clone() const = 0;

    //@}
  
    //! @name Class Creation and Accessor Methods
    //@{ 

    //! Create a Comm instance for global communication between nodes.
    virtual Teuchos::RCP< Teuchos::Comm<int> > getComm() const = 0;

    //@}

  }; // Platform class

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Platform<Scalar, LocalOrdinal, GlobalOrdinal>::Platform() {}
  
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Platform<Scalar, LocalOrdinal, GlobalOrdinal>::Platform(const std::string &str)
  : Teuchos::LabeledObject(str) {}
  
} // namespace Tpetra

#endif // TPETRA_PLATFORM_HPP

