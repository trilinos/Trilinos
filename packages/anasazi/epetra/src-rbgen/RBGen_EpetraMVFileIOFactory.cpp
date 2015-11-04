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

