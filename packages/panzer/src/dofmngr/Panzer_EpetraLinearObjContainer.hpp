// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_EpetraLinearObjContainer_hpp__
#define __Panzer_EpetraLinearObjContainer_hpp__

#include "Panzer_config.hpp"

#include <map>

// Epetra includes
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Panzer_LinearObjFactory.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {

class EpetraLinearObjContainer : public LinearObjContainer {
public:
   typedef LinearObjContainer::Members Members;

   virtual void initialize() 
   {
      if(x!=Teuchos::null) x->PutScalar(0.0);
      if(dxdt!=Teuchos::null) dxdt->PutScalar(0.0);
      if(f!=Teuchos::null) f->PutScalar(0.0);
      if(A!=Teuchos::null) A->PutScalar(0.0);
   }

   //! Wipe out stored data.
   void clear()
   {
      x = Teuchos::null;
      dxdt = Teuchos::null;
      f = Teuchos::null;
      A = Teuchos::null;
   }
   
   Teuchos::RCP<Epetra_Vector> x, dxdt, f;
   Teuchos::RCP<Epetra_CrsMatrix> A;
};

}

#endif
