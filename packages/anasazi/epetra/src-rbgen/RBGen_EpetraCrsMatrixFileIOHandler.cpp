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
      Epetra_MpiComm comm( MPI_COMM_WORLD );
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


