// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_NO_PREPROCESSOR_HPP
#define RBGEN_NO_PREPROCESSOR_HPP

#include "RBGen_Preprocessor.hpp"

namespace RBGen {

  //! Null preprocessor concrete class.
  template< class DataSetType > 
  class NoPreprocessor : public Preprocessor< DataSetType > {

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    NoPreprocessor() {};

    //! Destructor.
    virtual ~NoPreprocessor() {};
    //@}

    //! @name Initialization/Reset Methods
    //@{

    //! Initialize preprocessor
    void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params, 
                     const Teuchos::RCP< FileIOHandler<DataSetType> >& fileio ) {};

    //! Reset preprocessor
    void Reset() {};
    //@}

    //! @name Preprocess Methods
    //@{

    //! Preprocess the snapshot set passed in
    void Preprocess( Teuchos::RCP<DataSetType>& ss ) {};
    //@}

    //! @name Status Methods
    //@{

    //! Return initialized status of the preprocessor
    bool isInitialized() const { return true; };
    //@}
  };

} // end of RBGen namespace

#endif // RBGEN_NO_PREPROCESSOR_HPP
