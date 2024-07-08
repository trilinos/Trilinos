// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_EPETRAMV_FILEIO_FACTORY_HPP
#define RBGEN_EPETRAMV_FILEIO_FACTORY_HPP

#include "RBGen_FileIOFactory.hpp"
#include "RBGen_BurkardtFileIOHandler.h"
#include "RBGen_ConfigDefs.h"

#ifdef HAVE_ANASAZI_EPETRAEXT
#include "RBGen_MatrixMarketFileIOHandler.h"
#endif

#ifdef HAVE_ANASAZI_NETCDF
#include "RBGen_NetCDFFileIOHandler.h"
#endif

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace RBGen {

  //! Abstract factory for creating a concrete FileIOFactory for reading an Epetra_MultiVector.
  class EpetraMVFileIOFactory : public virtual FileIOFactory<Epetra_MultiVector> {
 
  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    EpetraMVFileIOFactory();

    //! Destructor.
    virtual ~EpetraMVFileIOFactory() {};
    //@}

    //! @name Factory methods
    //@{

    Teuchos::RCP< FileIOHandler< Epetra_MultiVector > > create( const Teuchos::ParameterList& params );

    //@}

  private:

    // Available file formats
    std::vector<std::string> file_formats;

  };

} // end of RBGen namespace

#endif
