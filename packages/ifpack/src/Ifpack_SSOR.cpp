#include "Ifpack_ConfigDefs.h"
#include "Ifpack_SSOR.h"
#include "Ifpack_Utils.h"
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
int Ifpack_SSOR::SetLabel()
{

  sprintf(Label_,"Ifpack SSOR, sweeps = %d, damping = %e",
	  NumSweeps(), DampingFactor());

  return(0);
}

//==============================================================================
int Ifpack_SSOR::
ApplyInverse(const Epetra_MultiVector& RHS, Epetra_MultiVector& LHS) const
{

  // sanity checks

  if (Matrix()->Filled() == false)
    IFPACK_CHK_ERR(-4);

  if (IsComputed() == false)
    IFPACK_CHK_ERR(-4);

  if (RHS.NumVectors() != LHS.NumVectors())
    IFPACK_CHK_ERR(-3);

  int NumVectors = RHS.NumVectors();

  // first time this method is called, replace zero diagonal
  // elements with 1.0
  if (FirstTime_) {
    for (int i = 0 ; i < Matrix()->NumMyRows() ; ++i) {
      double diag = (*Diagonal_)[i];
      if (diag == 0.0)
        (*Diagonal_)[i] = 1.0;
    }
    FirstTime_ = false;
  }
  
  if (PrintLevel() > 5)
    Ifpack_PrintResidual(Label(),*Matrix(),RHS,LHS);

  // the case NumSweeps_ == 1 can be coded more
  // efficiently 
  if (NumSweeps() == 1) {

    IFPACK_CHK_ERR(ApplySSOR(LHS));

    if (PrintLevel() > 5)
      Ifpack_PrintResidual(1,*Matrix(),RHS,LHS);

    return(0);
  }
  
  // now the more general case
  Epetra_MultiVector RHStmp(RHS);
  Epetra_MultiVector AX(RHS);

  for (int j = 0; j < NumSweeps() ; j++) {

    // compute the residual
    IFPACK_CHK_ERR(Apply(LHS,AX));
    AX.Update(1.0,RHStmp,-1.0);

    // apply the lower triangular part of A
    IFPACK_CHK_ERR(ApplySSOR(AX));

    // update the residual
    LHS.Update(DampingFactor(), AX, 1.0);

    if (PrintLevel() > 5)
      Ifpack_PrintResidual(j + 1,*Matrix(),RHS,LHS);
  }

  return(0);

}

//==============================================================================
int Ifpack_SSOR::
ApplySSOR(Epetra_MultiVector& X) const
{

  int NumVectors = X.NumVectors();
  int Length = Matrix()->MaxNumEntries();
  vector<int> Indices;
  vector<double> Values;
  Indices.resize(Length);
  Values.resize(Length);

  // Apply (D - \omega E)^{-1}

  for (int i = 0 ; i < Matrix()->NumMyRows() ; ++i) {

    int NumEntries;
    int col;
    double diag;
    double dtemp = 0.0;

    IFPACK_CHK_ERR(Matrix()->ExtractMyRowCopy(i, Length,NumEntries,
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

      X[m][i] = (X[m][i] - DampingFactor() * dtemp) / diag;

    }
  }

  // Apply D

  X.Multiply(1.0,X,*Diagonal_,0.0);

  // Apply (D - \omega F)^{-1}

  for (int i = Matrix()->NumMyRows()  - 1 ; i > -1 ; --i) {

    int NumEntries;
    int col;
    double diag;
    double dtemp = 0.0;

    IFPACK_CHK_ERR(Matrix()->ExtractMyRowCopy(i, Length,NumEntries,
					     &Values[0], &Indices[0]));

    for (int m = 0 ; m < NumVectors ; ++m) {

      for (int k = 0 ; k < NumEntries ; ++k) {

	col = Indices[k];
	
	if (col >= Matrix()->NumMyRows())
	  continue;

	if (col <= i) 
	  continue;

	dtemp += Values[k] * X[m][col];
      }

      X[m][i] = (X[m][i] - DampingFactor() * dtemp) / ((*Diagonal_)[i]);

    }
  }
}

