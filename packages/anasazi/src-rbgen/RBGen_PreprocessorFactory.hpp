// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_PREPROCESSOR_FACTORY_HPP
#define RBGEN_PREPROCESSOR_FACTORY_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace RBGen {

  template< class DataSetType >
  class Preprocessor;

  //! Abstract factory for instantiating Preprocessor concrete classes.
  template< class DataSetType >
  class PreprocessorFactory {
    
  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    PreprocessorFactory() {};

    //! Destructor.
    virtual ~PreprocessorFactory() {};
    //@}

    //! @name Factory methods
    //@{

    virtual Teuchos::RCP<Preprocessor<DataSetType> > create( const Teuchos::ParameterList& params ) = 0;
    
    //@}

  };
  
} // end of RBGen namespace

#endif // RBGEN_PREPROCESSOR_FACTORY_HPP
