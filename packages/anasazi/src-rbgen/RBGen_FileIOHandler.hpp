// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_FILE_IO_HANDLER_HPP
#define RBGEN_FILE_IO_HANDLER_HPP

#include "RBGen_ConfigDefs.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace RBGen {

  //! Abstract base class for reading datasets from files.
  /*!
   */
  template<class DataSetType>  
  class FileIOHandler {  
  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    FileIOHandler() {};

    //! Destructor.
    virtual ~FileIOHandler() {};
    //@}

    //! @name Initialization/Reset Methods
    //@{

    //! Initialize file reader using 
    virtual void Initialize( const Teuchos::RCP<Teuchos::ParameterList>& params ) = 0;
    
    void Reset() {};
    //@}
    
    //! @name File Reading Methods
    //@{
    
    //! Method for reading multiple files and putting them into an data set.
    virtual Teuchos::RCP< DataSetType > Read( const std::vector<std::string>& filenames ) = 0;

    //@}

    //! @name Writing Methods
    //@{

    //! Method for writing one data set into a file.
    virtual void Write( const Teuchos::RCP< const DataSetType >& MV, const std::string& filename ) = 0;

    //@}

    //! @name Handler Status Methods
    //@{

    //! Return initialized status of the handler
    virtual bool isInitialized() const = 0;

    //@}
  };

} // end of RBGen namespace

#endif // RBGEN_FILE_IO_HANDLER_HPP
