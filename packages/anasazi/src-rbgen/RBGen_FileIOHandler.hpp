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

#ifndef RBGEN_FILE_IO_HANDLER_HPP
#define RBGEN_FILE_IO_HANDLER_HPP

#include "RBGen_ConfigDefs.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace RBGen {

  //! Abstract base class for reading datasets from files.
  /*!
   */
  template<class DataSetType>  
  class FileIOHandler {  
  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    FileIOHandler() {};

    //! Destructor.
    virtual ~FileIOHandler() {};
    //@}

    //! @name Initialization/Reset Methods
    //@{

    //! Initialize file reader using 
    virtual void Initialize( const Teuchos::RCP<Teuchos::ParameterList>& params ) = 0;
    
    void Reset() {};
    //@}
    
    //! @name File Reading Methods
    //@{
    
    //! Method for reading multiple files and putting them into an data set.
    virtual Teuchos::RCP< DataSetType > Read( const std::vector<std::string>& filenames ) = 0;

    //@}

    //! @name Writing Methods
    //@{

    //! Method for writing one data set into a file.
    virtual void Write( const Teuchos::RCP< const DataSetType >& MV, const std::string& filename ) = 0;

    //@}

    //! @name Handler Status Methods
    //@{

    //! Return initialized status of the handler
    virtual bool isInitialized() const = 0;

    //@}
  };

} // end of RBGen namespace

#endif // RBGEN_FILE_IO_HANDLER_HPP
