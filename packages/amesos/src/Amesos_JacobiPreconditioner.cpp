#include "Amesos_ConfigDefs.h"
#include "Amesos_JacobiPreconditioner.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

//==============================================================================
Amesos_JacobiPreconditioner::
Amesos_JacobiPreconditioner(const Epetra_RowMatrix* Matrix,
			    Teuchos::ParameterList& List)
{
  IsComputed_ = false;
  Matrix_ = Matrix;
  UseTranspose_ = false;
  NumMyRows_ = Matrix_->NumMyRows();

  // fix parameters

  NumApplications_ = List.get("sweeps",1);
  Omega_ = List.get("omega", 1.0);
  DebugSmoother_ = List.get("debug", false);

  sprintf(Label_,"Amesos_Jacobi prec, sweeps = %d, omega = %e",
	  NumApplications_, Omega_);

  // sanity checks

  if (Matrix_->NumGlobalRows() != Matrix_->NumGlobalCols())
    AMESOS_CHK_ERRV(-3);

  if (NumApplications_ <= 0)
    AMESOS_CHK_ERRV(-3);

  // extract diagonal, store it in vector
  
  InvDiagonal_ = new Epetra_Vector(Matrix_->RowMatrixRowMap());

  if (InvDiagonal_ == 0)
    AMESOS_CHK_ERRV(-11);

  Matrix_->ExtractDiagonalCopy(*InvDiagonal_);

  // check that no zero diagonal element are found
  // and store the inverse of each element in InvDiagonal_

  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    if ((*InvDiagonal_)[i] != 0.0)
      (*InvDiagonal_)[i] = 1.0 / (*InvDiagonal_)[i];
    else
      (*InvDiagonal_)[i] = 1.0;
  }

  // allocate room for an auxiliary vector

  AX_ = new Epetra_MultiVector(Matrix_->RowMatrixRowMap(),1);

  if (AX_ == 0)
    AMESOS_CHK_ERRV(-11);


  IsComputed_ = true;

}

//==============================================================================
int Amesos_JacobiPreconditioner::
ApplyInverse(const Epetra_MultiVector& RHS, Epetra_MultiVector& LHS) const
{

  // sanity checks

  if (IsComputed() == false)
        AMESOS_CHK_ERR(-4);

  if (!AX_->Map().SameAs(RHS.Map())) 
    AMESOS_CHK_ERR(-3);

  if (RHS.NumVectors() != LHS.NumVectors())
    AMESOS_CHK_ERR(-3);

  // may need to reallocate AX_

  if (AX_->NumVectors() != RHS.NumVectors()) {

    delete AX_;
    AX_ = new Epetra_MultiVector(RHS.Map(), RHS.NumVectors());

    if (AX_ == 0)
      AMESOS_CHK_ERR(-11);
  }

  Epetra_MultiVector* RHS2;
  vector<double> Norm2;

  if (DebugSmoother_) {
    RHS2 = new Epetra_MultiVector(RHS);
    Norm2.resize(RHS.NumVectors());

    LHS.Norm2(&Norm2[0]);
    if (Comm().MyPID() == 0) {
      cout << "Jacobi preconditioner for the solution of Ax = b" << endl;
      cout << "||x||_2 = " << Norm2[0] << endl;
    }

    RHS.Norm2(&Norm2[0]);
    if (Comm().MyPID() == 0)
      cout << "||b||_2 = " << Norm2[0] << endl;
  }

  Epetra_MultiVector RHStmp(RHS);
  
  // Jacobi loop

  for (int j = 0; j < NumApplications_ ; j++) {

    AMESOS_CHK_ERR(Apply(LHS,*AX_));

    AX_->Update(1.0,RHStmp,-1.0);

    AX_->Multiply(1.0, *AX_, *InvDiagonal_, 0.0);

#ifdef FIXME
    if (smooth_ptr->omega == ML_ONE_STEP_CG) {
      /* Compute damping parameter that corresonds to one step of CG. */
      r_z_dot = 0.;
      for (i = 0; i < n; i++) r_z_dot += res[i]*res[i]*diagonal[i];
      r_z_dot = ML_gsum_double(r_z_dot, smooth_ptr->my_level->comm);
      ML_Operator_Apply(Amat, n, res, n, res2);
      p_ap_dot = ML_gdot(n, res, res2, smooth_ptr->my_level->comm);
      if (p_ap_dot != 0.0) omega = r_z_dot/p_ap_dot;
      else omega = 1.;
    }
#endif

    LHS.Update(Omega_, *AX_, 1.0);

    if (DebugSmoother_) {

      AMESOS_CHK_ERR(Apply(LHS,*RHS2));
      RHS2->Update(1.0, RHS, -1.0);

      RHS2->Norm2(&Norm2[0]);

      cout << "sweep " << j << ":  ||Ax - b||_2 = " 
	   << Norm2[0] << endl;

    }

  }

  if (DebugSmoother_)
    delete RHS2;

  return(0);

}

//==============================================================================
Amesos_JacobiPreconditioner::~Amesos_JacobiPreconditioner()
{
  if (InvDiagonal_)
    delete InvDiagonal_;

  if (AX_)
    delete AX_;

}

//==============================================================================
const Epetra_Comm & Amesos_JacobiPreconditioner::Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
const Epetra_Map & Amesos_JacobiPreconditioner::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map & Amesos_JacobiPreconditioner::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}

//==============================================================================
int Amesos_JacobiPreconditioner::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  AMESOS_CHK_ERR(Matrix_->Multiply(UseTranspose(),X,Y));
}
