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

   virtual void initialize() 
   {
      if(x!=Teuchos::null) x->PutScalar(0.0);
      if(dxdt!=Teuchos::null) dxdt->PutScalar(0.0);
      if(f!=Teuchos::null) f->PutScalar(0.0);
      if(A!=Teuchos::null) A->PutScalar(0.0);
   }
   
   Teuchos::RCP<Epetra_Vector> x, dxdt, f;
   Teuchos::RCP<Epetra_CrsMatrix> A;
};

}

#endif
