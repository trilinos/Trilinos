#ifndef IFPACK_SOR_H
#define IFPACK_SOR_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointPreconditioner.h"
#include <vector>
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Map;
class Epetra_Comm;
class Epetra_Vector;
class Epetra_RowMatrix;

//! Ifpack_SOR: a class to define point SOR preconditioners of Epetra_RowMatrix's.

class Ifpack_SOR : public Ifpack_PointPreconditioner {

public:

  //@{ \name Constructor/Destructor
  //! Constructs a SOR preconditioner object for the input Epetra_RowMatrix.
  Ifpack_SOR(const Epetra_RowMatrix* Matrix) :
    Ifpack_PointPreconditioner(Matrix),
    FirstTime_(true)
  {}

  //! Destructor.  
  ~Ifpack_SOR()
  {}

  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  virtual int SetLabel();

private:

  mutable bool FirstTime_;
};

#endif 
