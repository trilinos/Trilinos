// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_PREPROCESSOR_HPP
#define RBGEN_PREPROCESSOR_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace RBGen {
 
  template< class DataSetType > 
  class FileIOHandler;

  //! Abstract base class for encapsulating dataset preprocessing.
  template< class DataSetType >
  class Preprocessor {
    
  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    Preprocessor() {};

    //! Destructor.
    virtual ~Preprocessor() {};
    //@}

    //! @name Initialization/Reset Methods
    //@{

    //! Initialize preprocessor
    virtual void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params, 
			     const Teuchos::RCP<FileIOHandler<DataSetType> >& fileio ) = 0;
    
    //! Reset preprocessor
    virtual void Reset() = 0;
    //@}

    //! @name Preprocess Methods
    //@{

    //! Preprocess the snapshot set passed in
    virtual void Preprocess( Teuchos::RCP< DataSetType >& ss ) = 0;
    //@}

    //! @name Status Methods
    //@{

    //! Return initialized status of the preprocessor
    virtual bool isInitialized() const = 0;

    //@}
  };
  
} // end of RBGen namespace

#endif // RBGEN_PREPROCESSOR_HPP
