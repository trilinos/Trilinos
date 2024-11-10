// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
