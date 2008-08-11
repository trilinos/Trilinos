
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterList.hpp"

int simpleStratimikosSolve( 
   Epetra_CrsMatrix                          const& epetra_A,  // non-persisting, non-changeable
   Epetra_MultiVector                        const& epetra_B,  // non-persisting, non-changeable
   Epetra_MultiVector                             * epetra_X,  // non-persisting, changeable
   Teuchos::ParameterList                         * paramList  // non-persisting, changeable
   ) ;
