// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef RBGEN_EPETRAMV_METHOD_FACTORY_H
#define RBGEN_EPETRAMV_METHOD_FACTORY_H

#include "Teuchos_ParameterList.hpp"
#include "RBGen_MethodFactory.hpp"
#include "RBGen_ConfigDefs.h"

// Forward declaration of Epetra_Multivector
class Epetra_MultiVector;
class Epetra_Operator;

namespace RBGen {
 
  //! Specialization of MethodFactory for Epetra_MultiVector datasets.
  class EpetraMVMethodFactory : public virtual MethodFactory<Epetra_MultiVector,Epetra_Operator> {
    
  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    EpetraMVMethodFactory() {};

    //! Destructor.
    virtual ~EpetraMVMethodFactory() {};
    //@}

    //! @name Factory methods
    //@{

    Teuchos::RCP<Method< Epetra_MultiVector,Epetra_Operator > > create( const Teuchos::ParameterList& params );
    
    //@}

  };
  
} // end of RBGen namespace

#endif // RBGEN_EPETRAMV_METHOD_FACTORY_H
