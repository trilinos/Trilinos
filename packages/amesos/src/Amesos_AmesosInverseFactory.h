#ifndef AMESOS_AMESOSINVERSEFACTORY_H
#define AMESOS_AMESOSINVERSEFACTORY_H
#include "Teuchos_ParameterList.hpp"
#include "Amesos_InverseFactory.h"

class Epetra_Operator;
class Epetra_RowMatrix;

class Amesos_AmesosInverseFactory : public Amesos_InverseFactory {

public:
  virtual Epetra_Operator* Create(char* Type,
				  Epetra_RowMatrix* OverlappingMatrix,
				  Epetra_RowMatrix* Matrix,
				  Teuchos::ParameterList& List);

};

#endif // AMESOS_AMESOSINVERSEFACTORY_H
