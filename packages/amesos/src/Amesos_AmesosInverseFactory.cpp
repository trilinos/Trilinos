#include "Amesos_ConfigDefs.h"
#include "Amesos_InverseFactory.h"
#include "Amesos_AmesosInverseFactory.h"
#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Amesos_Preconditioner.h"
#include "Teuchos_ParameterList.hpp"

Epetra_Operator* Amesos_AmesosInverseFactory::
Create(char* Type, Epetra_RowMatrix* OverlappingMatrix,
       Epetra_RowMatrix* Matrix, Teuchos::ParameterList& List)
{

  Amesos_Preconditioner* APrec;

  // matrix is already local
  APrec = new Amesos_Preconditioner(Type, OverlappingMatrix, List, false);

  return(APrec);

}
