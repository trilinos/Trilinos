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

