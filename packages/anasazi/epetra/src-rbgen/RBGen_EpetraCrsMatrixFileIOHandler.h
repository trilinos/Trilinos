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

#ifndef EPETRA_CRSMATRIX_FILE_IO_HANDLER_H
#define EPETRA_CRSMATRIX_FILE_IO_HANDLER_H


#include "RBGen_FileIOHandler.hpp"

// Forward declaration of Epetra_MultiVector class
class Epetra_Operator;

namespace RBGen {
  
  //! FileIOHandler for reading EpetraCrsMatrix data from a file using EpetraExt.
  class EpetraCrsMatrixFileIOHandler : public virtual FileIOHandler< Epetra_Operator > {
    
  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    EpetraCrsMatrixFileIOHandler();

    //! Destructor.
    virtual ~EpetraCrsMatrixFileIOHandler() {};

    //@}

    //! @name Initialization/Reset Methods
    //@{

    void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params );

    void Reset() { isInit = false; };

    //@}

    //! @name File Reading Methods
    //@{

    //! Method for reading a file and constructing an Epetra_CrsMatrix.
    Teuchos::RCP<Epetra_Operator> Read( const std::vector<std::string>& filenames );

    //@}

    //! @name Writing Methods
    //@{

    //! Method for writing one Epetra_CrsMatrix into a file using the same type as was.
    void Write( const Teuchos::RCP<const Epetra_Operator>& MTX, const std::string& filename );

    //@}
    //! @name Handler Status Methods
    //@{

    //! Return initialized status of the handler
    bool isInitialized() const { return isInit; };

    //@}

  private:

    // Whether or not we know the file format.
    bool isInit;

    // File input / output paths
    std::string in_path, out_path;

    // ParameterList that this file handler was initialized with.
    Teuchos::RCP< Teuchos::ParameterList > params_;

  };
  
} // namespace RBGen

#endif // EPETRA_CRSMATRIX_FILE_IO_HANDLER_H

