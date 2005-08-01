#include <iostream>
#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

class PyOperator : public Epetra_Operator
{
public:
  PyOperator(const Epetra_Comm& Comm) :
    Comm_(Comm)
  {}

  ~PyOperator() {}

  int SetUseTranspose(bool UseTranspose) = 0;

  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

  virtual double NormInf() const = 0;

  virtual const char * Label() const
  {
    return("PySerialOperator");
  }

  virtual bool UseTranspose() const = 0;

  virtual bool HasNormInf() const = 0;

  virtual const Epetra_Comm & Comm() const
  {
    return(Comm_);
  }

  virtual const Epetra_Map & OperatorDomainMap() const = 0;

  virtual const Epetra_Map & OperatorRangeMap() const = 0;

  virtual const Epetra_Map & Map() const = 0;

private:

  const Epetra_Comm& Comm_;
};
