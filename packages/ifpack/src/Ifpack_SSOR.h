#ifndef IFPACK_SSOR_H
#define IFPACK_SSOR_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointPreconditioner.h"
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Map;
class Epetra_Comm;
class Epetra_Vector;
class Epetra_RowMatrix;

//! Ifpack_SSOR: a class to define point Gauss-Seidel preconditioners of Epetra_RowMatrix's.

class Ifpack_SSOR : public Ifpack_PointPreconditioner {

public:

  //@{ \name Constructor/Destructor
  //! Constructs a Gauss-Seidel preconditioner object for the input Epetra_RowMatrix.
  Ifpack_SSOR(const Epetra_RowMatrix* Matrix) :
    Ifpack_PointPreconditioner(Matrix),
    FirstTime_(true)
  {}

  //! Destructor.  
  ~Ifpack_SSOR()
  {}

  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  virtual int SetLabel();

private:

  virtual int ApplySSOR(Epetra_MultiVector& X) const;
  mutable bool FirstTime_;
};
#endif
