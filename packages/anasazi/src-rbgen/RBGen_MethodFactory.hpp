// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_METHOD_FACTORY_HPP
#define RBGEN_METHOD_FACTORY_HPP

#include "Teuchos_ParameterList.hpp"
#include "RBGen_Method.hpp"

namespace RBGen {

  //! Abstract factory for creating basis generation methods.
  template< class DataSetType, class OperatorType > 
  class MethodFactory {

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    MethodFactory() {};

    //! Destructor.
    virtual ~MethodFactory() {};
    //@}

    //! @name Factory methods
    //@{

    virtual Teuchos::RCP<Method< DataSetType, OperatorType > > create( const Teuchos::ParameterList& params ) = 0;

    //@}

  };

} // end of RBGen namespace

#endif // RBGEN_METHOD_FACTORY_HPP
