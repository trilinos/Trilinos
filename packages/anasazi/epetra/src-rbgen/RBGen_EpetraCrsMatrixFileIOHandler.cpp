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


#include "RBGen_EpetraCrsMatrixFileIOHandler.h"
#include "RBGen_ConfigDefs.h"

#include "Epetra_BLAS.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"

#include "EpetraExt_readEpetraLinearSystem.h"
#include "EpetraExt_RowMatrixOut.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "AnasaziGlobalComm.hpp"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_Assert.hpp"

namespace RBGen {

  EpetraCrsMatrixFileIOHandler::EpetraCrsMatrixFileIOHandler()
    : isInit(false)
  {
  }

  void EpetraCrsMatrixFileIOHandler::Initialize( const Teuchos::RCP<Teuchos::ParameterList>& params )
  {
    // Get the "File I/O" sublist.
    Teuchos::ParameterList& fileio_params = params->sublist( "File IO" );

    // Get the input path.
    in_path = "";
    if ( fileio_params.isParameter( "Data Input Path" ) ) {
      in_path = Teuchos::getParameter<std::string>( fileio_params, "Data Input Path" );
    }

    // Get the output path.
    out_path = "";
    if ( fileio_params.isParameter( "Data Output Path" ) ) {
      out_path = Teuchos::getParameter<std::string>( fileio_params, "Data Output Path" );
    }

    // This file i/o handler is now initialized.
    isInit = true;
  }

  Teuchos::RCP<Epetra_Operator> EpetraCrsMatrixFileIOHandler::Read( const std::vector<std::string>& filenames )
  {

    Teuchos::RCP<Epetra_CrsMatrix> newMTX;

    if (isInit) {

#ifdef EPETRA_MPI
      Epetra_MpiComm comm( Anasazi::get_global_comm() );
#else
      Epetra_SerialComm comm;
#endif
      // This reader is only equipped to read in one matrix file.
      TEUCHOS_TEST_FOR_EXCEPTION(filenames.size() > 1, std::invalid_argument, "File I/O handler cannot read more than one file!");

      // Open the data file
      std::string temp_filename = in_path + filenames[0];

      // Create a null pointer to the Epetra_Map
      Teuchos::RCP<Epetra_Map> Map;

      // Read in the matrix from file
      EpetraExt::readEpetraLinearSystem( temp_filename, comm, &newMTX, &Map );

    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "File I/O handler is not initialized!");
    }
    // Return.
    return newMTX;
  }

  void EpetraCrsMatrixFileIOHandler::Write( const Teuchos::RCP<const Epetra_Operator>& MTX, const std::string& filename )
  {
    if (isInit) {

      std::string temp_filename = out_path + filename;
      EpetraExt::RowMatrixToMatrixMarketFile( temp_filename.c_str(), *(Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(MTX)) );

    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "File I/O handler is not initialized!");
    }
  }

} // namespace RBGen


