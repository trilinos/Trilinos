#include "Ifpack_ConfigDefs.h"
#include "Ifpack_GaussSeidel.h"
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
void Ifpack_GaussSeidel::SetLabel()
{

  sprintf(Label_,"Ifpack Gauss-Seidel, sweeps = %d, damping = %e",
	  NumSweeps(), DampingFactor());
}

//==============================================================================
int Ifpack_GaussSeidel::
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

  // first time this method is called, compute the inverse of
  // each diagonal element
  if (FirstTime_) {
    for (int i = 0 ; i < Matrix().NumMyRows() ; ++i) {
      double diag = (*Diagonal_)[i];
      if (diag == 0.0)
	diag = 1.0;
      (*Diagonal_)[i] = 1.0/diag;
    }
    FirstTime_ = false;
  }
  
  int Length = Matrix().MaxNumEntries();
  vector<int> Indices;
  vector<double> Values;
  Indices.resize(Length);
  Values.resize(Length);

  int NumMyRows = Matrix().NumMyRows();

  // ============ //
  // single sweep //
  // ------------ //

  if ((NumSweeps() == 1)&& ZeroStartingSolution_
      && (PrintFrequency() != 0)) {

    for (int i = 0 ; i < NumMyRows ; ++i) {

      int NumEntries;
      int col;
      IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(i, Length,NumEntries,
						&Values[0], &Indices[0]));

      for (int m = 0 ; m < NumVectors ; ++m) {

	double dtemp = 0.0;
	for (int k = 0 ; k < NumEntries ; ++k) {

	  col = Indices[k];
	  if (col >= i || (col >= NumMyRows)) continue;

	  dtemp += Values[k] * Y[m][col];
	}

	Y[m][i] = DampingFactor() * ((*Diagonal_)[i])
	  * (Y[m][i] - dtemp);
      }
    }

    return(0);
  }

  // --------------------- //
  // general case (solver) //
  // --------------------- //

  // need an additional vector for AztecOO preconditioning
  // (as X and Y both point to the same memory space)
  Epetra_MultiVector Xtmp(X);

  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  if (PrintFrequency())
    Ifpack_PrintResidual(Label(),Matrix(),Y,Xtmp);

  for (int j = 0; j < NumSweeps() ; j++) {

    for (int i = 0 ; i < NumMyRows ; ++i) {

      int NumEntries;
      int col;
      IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(i, Length,NumEntries,
						&Values[0], &Indices[0]));

      for (int m = 0 ; m < NumVectors ; ++m) {

	double dtemp = 0.0;
	for (int k = 0 ; k < NumEntries ; ++k) {

	  col = Indices[k];

	  if ((col == i) || (col >= NumMyRows)) continue;

	  dtemp += Values[k] * Y[m][col];
	}

	Y[m][i] = DampingFactor() * ((*Diagonal_)[i])
	  * (Xtmp[m][i] - dtemp);
      }
    }

    if (PrintFrequency() && (j != 0) && (j % PrintFrequency() == 0))
      Ifpack_PrintResidual(j,Matrix(),Y,Xtmp);
  }

  if (PrintFrequency())
    Ifpack_PrintResidual(NumSweeps(),Matrix(),Y,Xtmp);

  return(0);

}
