#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointPreconditioner.h"
#include "Ifpack_Jacobi.h"
#include "Ifpack_Utils.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

//==============================================================================
int Ifpack_Jacobi::SetLabel()
{

  sprintf(Label_,"Ifpack Jacobi, sweeps = %d, damping = %e",
	  NumSweeps(), DampingFactor());

  return(0);
}

//==============================================================================
int Ifpack_Jacobi::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  if (IsComputed() == false)
    IFPACK_CHK_ERR(-4);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-3);

  // first time this method is called, compute the inverse of 
  // each diagonal element
  if (FirstTime_) {
    for (int i = 0 ; i < Matrix()->NumMyRows() ; ++i) {
      double diag = (*Diagonal_)[i];
      if (diag != 0.0)
	(*Diagonal_)[i] = 1.0 / diag;
      else
	(*Diagonal_)[i] = 0.0;
    }
    FirstTime_ = false;
  }

  // single sweep, easy to do. No output,
  // as I wanne be efficient and clean

  if (NumSweeps() == 1) {
    
    IFPACK_CHK_ERR(Y.Multiply(DampingFactor(),Y,*Diagonal_,0.0));
    return(0);
  }

  Epetra_MultiVector Xtmp(X);

  // general case (solver)
  // need an additional vector for AztecOO preconditioning
  // (as X and Y both point to the same memory space)

  Epetra_MultiVector AX(X);
  Y.PutScalar(0.0);
  
  if (PrintLevel() > 5)
    Ifpack_PrintResidual(Label(),*Matrix(),Y,Xtmp);

  for (int j = 0; j < NumSweeps() ; j++) {

    // compute the residual
    IFPACK_CHK_ERR(Apply(Y,AX));
    AX.Update(1.0,Xtmp,-1.0);

    // apply the inverse of the diagonal
    AX.Multiply(1.0, AX, *Diagonal_, 0.0);

    // update the solution at step `j'
    Y.Update(DampingFactor(), AX, 1.0);

    if (PrintLevel() > 5)
      Ifpack_PrintResidual(j + 1,*Matrix(),Y,Xtmp);

  }

  return(0);

}

