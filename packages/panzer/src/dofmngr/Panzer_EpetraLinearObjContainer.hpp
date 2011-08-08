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
   typedef enum { X=0x1, DxDt=0x2, F=0x4, Mat=0x8} Members;

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
