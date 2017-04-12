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

#include "RBGen_EpetraMVFileIOFactory.h"
#include "Teuchos_Assert.hpp"

namespace RBGen {

  EpetraMVFileIOFactory::EpetraMVFileIOFactory() 
  { 
     // Insert the acceptable input file types for this factory
     file_formats.push_back("Burkardt");
#ifdef HAVE_ANASAZI_NETCDF
     file_formats.push_back("NetCDF");
#endif
#ifdef HAVE_ANASAZI_EPETRAEXT
     file_formats.push_back("Matrix Market");
#endif
  }

  Teuchos::RCP< FileIOHandler< Epetra_MultiVector > >
  EpetraMVFileIOFactory::create( const Teuchos::ParameterList& params )
  {
    // See if the "File I/O" sublist exists
    TEUCHOS_TEST_FOR_EXCEPTION(!params.isSublist( "File IO" ), std::invalid_argument, "File IO sublist does not exist!");

    // Get the "File I/O" sublist.
    const Teuchos::ParameterList& fileio_params = params.sublist( "File IO" );

    // Get the file format type
    std::string file_format = Teuchos::getParameter<std::string>( const_cast<Teuchos::ParameterList&>(fileio_params),
                                                                  "Type" );

    Teuchos::RCP< FileIOHandler< Epetra_MultiVector > > RBFileIO;

    // File input format based on Burkardt's input files
    if ( file_format == "Burkardt" ) {
      RBFileIO = Teuchos::rcp( new BurkardtFileIOHandler() );
    } else
    // File input format for NetCDF files
#ifdef HAVE_ANASAZI_NETCDF
    if ( file_format == "NetCDF" ) {
      RBFileIO = Teuchos::rcp( new NetCDFFileIOHandler() );
    } else
#endif
    // File input format for Matrix Market files
#ifdef HAVE_ANASAZI_EPETRAEXT
    if ( file_format == "Matrix Market" ) {
      RBFileIO = Teuchos::rcp( new MatrixMarketFileIOHandler() );
    } else 
#endif
    {
    // Throw an exception because the format type is not recognized by this factory
    }
    //
    // Return the method created
    //
    return RBFileIO;
  }
  
} // end of RBGen namespace

