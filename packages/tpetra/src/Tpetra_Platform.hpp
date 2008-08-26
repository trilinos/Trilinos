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
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Object.hpp>

namespace Tpetra {

  //! Tpetra::Platform: The Tpetra Platform Abstract Base Class
  /*! Platform is an abstract base class. 
      Platform is used to generate Comm instances. It also manages platform-specific information, 
      such as how inter-image communication is implemented.
      An implementation of Platform (e.g., SerialPlatform) will create corresponding classes,
      (e.g., SerialComm). These will then be cast to their base class, and passed back to 
      other Tpetra modules. As a result, other Tpetra modules don't need to know anything about 
      the platform they're running on, or any implementation-specific details.
  */

  template<typename Ordinal>
  class Platform : public Teuchos::Object {
  public:
  
    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructor
    Platform(const std::string &str);

    //! Destructor
    virtual ~Platform() {};

    //! Clone method
    /*! Returns a copy of this Platform instance. It is allocated on the heap and
        encapsulated in a Teuchos RCP.
    */
    virtual Teuchos::RCP< Platform<Ordinal> > clone() const = 0;

    //@}
  
    //! @name Class Creation and Accessor Methods
    //@{ 

    //! Create a Comm instance for global communication between nodes.
    virtual Teuchos::RCP< Teuchos::Comm<Ordinal> > createComm() const = 0;

    //@}

  }; // Platform class

  template<typename Ordinal>
  Platform<Ordinal>::Platform(const std::string &str)
  : Teuchos::Object(str.c_str()) {}
  
} // namespace Tpetra

#endif // TPETRA_PLATFORM_HPP

