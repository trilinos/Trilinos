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
void Ifpack_SSOR::SetLabel()
{
  sprintf(Label_,"Ifpack SSOR, sweeps = %d, damping = %e",
	  NumSweeps(), DampingFactor());
}

//==============================================================================
int Ifpack_SSOR::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  // sanity checks

  if (Matrix().Filled() == false)
    IFPACK_CHK_ERR(-4);

  if (IsComputed() == false)
    IFPACK_CHK_ERR(-4);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-3);

  int NumVectors = X.NumVectors();
  int NumMyRows = Matrix().NumMyRows();

  // first time this method is called, replace zero diagonal
  // elements with 1.0
  if (FirstTime_) {
    for (int i = 0 ; i < NumMyRows ; ++i) {
      double diag = (*Diagonal_)[i];
      if (diag == 0.0)
        (*Diagonal_)[i] = 1.0;
    }
    FirstTime_ = false;
  }

  // ---------------- //
  // single sweep can //
  // ---------------- //
 
  if ((NumSweeps() == 1) && ZeroStartingSolution_
      && (PrintFrequency() != 0)) {
    IFPACK_CHK_ERR(ApplySSOR(Y));
    return(0);
  }
  
  // --------------------- //
  // general case (solver) //
  // --------------------- //
  
  Epetra_MultiVector Xtmp(X);
  Epetra_MultiVector AX(X);

  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  if (PrintFrequency())
    Ifpack_PrintResidual(Label(),Matrix(),Y,Xtmp);
  
  for (int j = 0; j < NumSweeps() ; j++) {

    // compute the residual
    IFPACK_CHK_ERR(Apply(Y,AX));
    AX.Update(1.0,Xtmp,-1.0);

    // apply the lower triangular part of A
    IFPACK_CHK_ERR(ApplySSOR(AX));

    // update the residual
    Y.Update(DampingFactor(), AX, 1.0);

    if (PrintFrequency() && (j != 0) && (j % PrintFrequency() == 0))
      Ifpack_PrintResidual(j,Matrix(),Y,Xtmp);
    
  }

  if (PrintFrequency())
    Ifpack_PrintResidual(NumSweeps(),Matrix(),Y,Xtmp);

  return(0);

}

//==============================================================================
int Ifpack_SSOR::
ApplySSOR(Epetra_MultiVector& X) const
{

  int NumVectors = X.NumVectors();
  int Length = Matrix().MaxNumEntries();
  vector<int> Indices;
  vector<double> Values;
  Indices.resize(Length);
  Values.resize(Length);
  int NumMyRows = Matrix().NumMyRows();

  // Apply (D - \omega E)^{-1}

  for (int i = 0 ; i < NumMyRows ; ++i) {

    int NumEntries;
    int col;
    double diag;
    double dtemp = 0.0;

    IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(i, Length,NumEntries,
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

  for (int i = NumMyRows  - 1 ; i > -1 ; --i) {

    int NumEntries;
    int col;
    double diag;
    double dtemp = 0.0;

    IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(i, Length,NumEntries,
					     &Values[0], &Indices[0]));

    for (int m = 0 ; m < NumVectors ; ++m) {

      for (int k = 0 ; k < NumEntries ; ++k) {

	col = Indices[k];
	
	if (col >= NumMyRows)
	  continue;

	if (col <= i) 
	  continue;

	dtemp += Values[k] * X[m][col];
      }

      X[m][i] = (X[m][i] - DampingFactor() * dtemp) / ((*Diagonal_)[i]);

    }
  }
  return(0);
}

