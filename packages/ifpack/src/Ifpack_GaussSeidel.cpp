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
int Ifpack_GaussSeidel::SetLabel()
{

  sprintf(Label_,"Ifpack Gauss-Seidel, sweeps = %d, damping = %e",
	  NumSweeps(), DampingFactor());

  return(0);
}

//==============================================================================
int Ifpack_GaussSeidel::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  // sanity checks

  if (Matrix()->Filled() == false)
    IFPACK_CHK_ERR(-4);

  if (IsComputed() == false)
    IFPACK_CHK_ERR(-4);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-3);

  int NumVectors = X.NumVectors();

  // first time this method is called, compute the inverse of
  // each diagonal element
  if (FirstTime_) {
    for (int i = 0 ; i < Matrix()->NumMyRows() ; ++i) {
      double diag = (*Diagonal_)[i];
      if (diag == 0.0)
	diag = 1.0;
      (*Diagonal_)[i] = 1.0/diag;
    }
    FirstTime_ = false;
  }
  
  int Length = Matrix()->MaxNumEntries();
  vector<int> Indices;
  vector<double> Values;
  Indices.resize(Length);
  Values.resize(Length);

  // NumSweeps() == 1 can he handled more efficiently
  // than the general case, as I can overwrite the Y vector
  // directly
  // no print out for this, as I want to save one vector

  if (NumSweeps() == 1) {

    for (int i = 0 ; i < Matrix()->NumMyRows() ; ++i) {

      int NumEntries;
      int col;
      IFPACK_CHK_ERR(Matrix()->ExtractMyRowCopy(i, Length,NumEntries,
						&Values[0], &Indices[0]));

      for (int m = 0 ; m < NumVectors ; ++m) {

	double dtemp = 0.0;
	for (int k = 0 ; k < NumEntries ; ++k) {

	  col = Indices[k];
	  if (col >= i) continue;

	  dtemp += Values[k] * Y[m][col];
	}

	Y[m][i] = DampingFactor() * ((*Diagonal_)[i])
	  * (Y[m][i] - dtemp);
      }
    }

    return(0);
  }

  // now the general case

  Epetra_MultiVector Xtmp(X);
  Y.PutScalar(0.0);

  if (PrintLevel() > 5)
    Ifpack_PrintResidual(Label(),*Matrix(),Y,Xtmp);

  for (int j = 0; j < NumSweeps() ; j++) {

    for (int i = 0 ; i < Matrix()->NumMyRows() ; ++i) {

      int NumEntries;
      int col;
      IFPACK_CHK_ERR(Matrix()->ExtractMyRowCopy(i, Length,NumEntries,
						&Values[0], &Indices[0]));

      for (int m = 0 ; m < NumVectors ; ++m) {

	double dtemp = 0.0;
	for (int k = 0 ; k < NumEntries ; ++k) {

	  col = Indices[k];

	  if (col == i) continue;

	  dtemp += Values[k] * Y[m][col];
	}

	Y[m][i] = DampingFactor() * ((*Diagonal_)[i])
	  * (Xtmp[m][i] - dtemp);

      }
    }

    if (PrintLevel() > 5)
      Ifpack_PrintResidual(j + 1,*Matrix(),Y,Xtmp);
  }

  return(0);

}
