#ifndef AMESOS_INVERSEFACTORY_H
#define AMESOS_INVERSEFACTORY_H

class Epetra_RowMatrix;
class Epetra_Operator;
#include "Teuchos_ParameterList.hpp"

class Amesos_InverseFactory {

public:
  virtual Epetra_Operator* Create(char* Type,
				  Epetra_RowMatrix* OverlappingMatrix,
				  Epetra_RowMatrix* Matrix,
				  Teuchos::ParameterList& List) = 0;

};

#endif // AMESOS_INVERSEFACTORY_H
