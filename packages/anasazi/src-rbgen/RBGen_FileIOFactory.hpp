// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef RBGEN_FILEIO_FACTORY_HPP
#define RBGEN_FILEIO_FACTORY_HPP

#include "Teuchos_RCP.hpp"

// Forward declarations for Teuchos.
namespace Teuchos {
  class ParameterList;
}

namespace RBGen {

  // Forward declarations for FileIOHandler.
  template< class DataSetType >
  class FileIOHandler;

  //! Abstract factory for instantiating FileIOHandler concrete classes.
  /*!
   */
  template< class DataSetType > 
  class FileIOFactory {

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    FileIOFactory() {};

    //! Destructor.
    virtual ~FileIOFactory() {};
    //@}

    //! @name Factory methods
    //@{

    virtual Teuchos::RCP< FileIOHandler< DataSetType > > create( const Teuchos::ParameterList& params ) = 0;

    //@}

  };

} // end of RBGen namespace

#endif // RBGEN_FILEIO_FACTORY_HPP
