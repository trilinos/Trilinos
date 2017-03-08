// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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

