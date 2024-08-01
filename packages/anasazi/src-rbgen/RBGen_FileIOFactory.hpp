// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_FILEIO_FACTORY_HPP
#define RBGEN_FILEIO_FACTORY_HPP

#include "Teuchos_RCP.hpp"

// Forward declarations for Teuchos.
namespace Teuchos {
  class ParameterList;
}

namespace RBGen {

  // Forward declarations for FileIOHandler.
  template< class DataSetType >
  class FileIOHandler;

  //! Abstract factory for instantiating FileIOHandler concrete classes.
  /*!
   */
  template< class DataSetType > 
  class FileIOFactory {

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    FileIOFactory() {};

    //! Destructor.
    virtual ~FileIOFactory() {};
    //@}

    //! @name Factory methods
    //@{

    virtual Teuchos::RCP< FileIOHandler< DataSetType > > create( const Teuchos::ParameterList& params ) = 0;

    //@}

  };

} // end of RBGen namespace

#endif // RBGEN_FILEIO_FACTORY_HPP
