// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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
