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


#ifndef NETCDF_FILE_IO_HANDLER_H
#define NETCDF_FILE_IO_HANDLER_H


#include "RBGen_FileIOHandler.hpp"
#include "RBGen_ConfigDefs.h"

#include <netcdf.h>

// Forward declaration of Epetra_MultiVector class
class Epetra_MultiVector;

namespace RBGen {
 
  //! FileIOHandler for reading an Epetra_MultiVector from a NetCDF file.
  class NetCDFFileIOHandler : public virtual FileIOHandler< Epetra_MultiVector > 
  {  
  public:
    //! @name Constructor/Destructor.
    //@{
    
    //! Default constructor.
    NetCDFFileIOHandler();
    
    //! Destructor.
    virtual ~NetCDFFileIOHandler();
    
    //@}
    
    //! @name Initialization/Reset Methods
    //@{
    
    void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params );
    
    void Reset();
    
    //@}
    
    //! @name File Reading Methods
    //@{
    
    //! Method for reading multiple files and putting them into an Epetra_MultiVector.
    Teuchos::RCP< Epetra_MultiVector > Read( const std::vector<std::string>& filenames );

    //@}

    //! @name Writing Methods
    //@{

    //! Method for writing one Epetra_MultiVector into a file.
    void Write( const Teuchos::RCP< const Epetra_MultiVector >& MV, const std::string& filename );

    //@}

    //! @name Handler Status Methods
    //@{

    //! Return initialized status of the handler
    bool isInitialized() const { return isInit; };

    //@}

  private:

    bool isInit;
    int num_nodes, num_nod_var, len_string;
    char **var_name;
    string in_path, out_path;
    Teuchos::RCP< Teuchos::ParameterList > params_;

    // Method for handling error from NetCDF.
    void handle_error( int status ) {
      if (status != NC_NOERR) {
	fprintf(stderr,"%s\n", nc_strerror(status));
	exit(-1);
      }
    }
  };
  
} // namespace RBGen

#endif // NETCDF_FILE_IO_HANDLER_H

