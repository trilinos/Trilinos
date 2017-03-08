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
