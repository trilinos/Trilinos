#include "Amesos_ConfigDefs.h"
#include "Amesos_SSORPreconditioner.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_BlockMap.h"
#include <vector>

//==============================================================================
Amesos_SSORPreconditioner::
Amesos_SSORPreconditioner(const Epetra_RowMatrix* Matrix,
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

  sprintf(Label_,"Amesos_SSOR prec, sweeps = %d, omega = %e",
	  NumApplications_, Omega_);

  // sanity checks

  if (Matrix_->NumGlobalRows() != Matrix_->NumGlobalCols())
    AMESOS_CHK_ERRV(-3);

  if (NumApplications_ <= 0)
    AMESOS_CHK_ERRV(-3);

  // extract diagonal, store it in vector

  Diagonal_ = new Epetra_Vector(Matrix_->RowMatrixRowMap());

  if (Diagonal_ == 0)
    AMESOS_CHK_ERRV(-11);

  Matrix_->ExtractDiagonalCopy(*Diagonal_);

  // allocate room for an auxiliary vector

  AX_ = new Epetra_MultiVector(Matrix_->RowMatrixRowMap(),1);

  if (AX_ == 0)
    AMESOS_CHK_ERRV(-11);


  IsComputed_ = true;

}

//==============================================================================
int Amesos_SSORPreconditioner::
ApplyInverse(const Epetra_MultiVector& RHS, Epetra_MultiVector& LHS) const
{

  // sanity checks

  if (Matrix_->Filled() == false)
    AMESOS_CHK_ERR(-4);

  if (IsComputed() == false)
    AMESOS_CHK_ERR(-4);

  if (!AX_->Map().SameAs(RHS.Map())) 
    AMESOS_CHK_ERR(-3);

  if (RHS.NumVectors() != LHS.NumVectors())
    AMESOS_CHK_ERR(-3);

  // may need to reallocate AX_

  int NumVectors = RHS.NumVectors();

  if (AX_->NumVectors() != NumVectors) {

    delete AX_;
    AX_ = new Epetra_MultiVector(RHS.Map(), NumVectors);

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
      cout << "SSOR preconditioner for the solution of Ax = b" << endl;
      cout << "||x||_2 = " << Norm2[0] << endl;
    }

    RHS.Norm2(&Norm2[0]);
    if (Comm().MyPID() == 0)
      cout << "||b||_2 = " << Norm2[0] << endl;
  }

  Epetra_MultiVector RHStmp(RHS);
  
  // SSOR loop

  // the case NumApplications_ == 1 can be coded more
  // efficiently 

  if (NumApplications_ == 1) {

    ApplySSOR(LHS);

    return(0);
  }
  
  // now the more general case

  for (int j = 0; j < NumApplications_ ; j++) {

    // compute the residual

    AMESOS_CHK_ERR(Apply(LHS,*AX_));

    AX_->Update(1.0,RHStmp,-1.0);

    // apply the lower triangular part of A

    ApplySSOR(*AX_);

    // update the residual
 
    LHS.Update(Omega_, *AX_, 1.0);

    if (DebugSmoother_) {

      AMESOS_CHK_ERR(Apply(LHS,*RHS2));
      RHS2->Update(1.0, RHStmp, -1.0);

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
int Amesos_SSORPreconditioner::
ApplySSOR(Epetra_MultiVector& X) const
{

  int NumVectors = X.NumVectors();
  int Length = Matrix_->MaxNumEntries();
  vector<int> Indices;
  vector<double> Values;
  Indices.resize(Length);
  Values.resize(Length);

  // Apply (D - \omega E)^{-1}

  for (int i = 0 ; i < NumMyRows_ ; ++i) {

    int NumEntries;
    int col;
    double diag;
    double dtemp = 0.0;

    AMESOS_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
					     &Values[0], &Indices[0]));

    for (int m = 0 ; m < NumVectors ; ++m) {

      for (int k = 0 ; k < NumEntries ; ++k) {

	col = Indices[k];
	if (col > i) 
	  continue;

	if (col == i) 
	  diag = Values[k];
	else
	  dtemp += Values[k] * X[m][col];
      }

      X[m][i] = (X[m][i] - Omega_ * dtemp) / diag;

    }
  }

  // Apply D

  X.Multiply(1.0,X,*Diagonal_,0.0);

  // Apply (D - \omega F)^{-1}

  for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

    int NumEntries;
    int col;
    double diag;
    double dtemp = 0.0;

    AMESOS_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
					     &Values[0], &Indices[0]));

    for (int m = 0 ; m < NumVectors ; ++m) {

      for (int k = 0 ; k < NumEntries ; ++k) {

	col = Indices[k];
	
	if (col >= NumMyRows_)
	  continue;

	if (col < i) 
	  continue;

	if (col == i) 
	  diag = Values[k];
	else
	  dtemp += Values[k] * X[m][col];
      }

      X[m][i] = (X[m][i] - Omega_ * dtemp) / diag;

    }
  }
}

//==============================================================================
Amesos_SSORPreconditioner::~Amesos_SSORPreconditioner()
{
  if (AX_)
    delete AX_;

  if (Diagonal_)
    delete Diagonal_;
}

//==============================================================================
const Epetra_Comm & Amesos_SSORPreconditioner::Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
const Epetra_Map & Amesos_SSORPreconditioner::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map & Amesos_SSORPreconditioner::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}

//==============================================================================
int Amesos_SSORPreconditioner::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  AMESOS_CHK_ERR(Matrix_->Multiply(UseTranspose(),X,Y));
}
