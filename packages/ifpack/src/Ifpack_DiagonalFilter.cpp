#include "Ifpack_ConfigDefs.h"
#include "Ifpack_DiagonalFilter.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

//==============================================================================
int Ifpack_DiagonalFilter::
ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, 
		 double* Values, int* Indices) const
{

  IFPACK_CHK_ERR(A_.ExtractMyRowCopy(MyRow, Length, NumEntries,
				     Values,Indices));

  // modify the diagonal value
  for (int i = 0 ; i < NumEntries ; ++i) {
    if (Indices[i] == MyRow) {
      Values[i] = Values[i] * RelativeThreshold_ +
        AbsoluteThreshold_ * EPETRA_SGN(Values[i]);
    }
  }

  return(0);
}

//==============================================================================
int Ifpack_DiagonalFilter::
Multiply(bool TransA, const Epetra_MultiVector& X, 
	 Epetra_MultiVector& Y) const
{

  int NumVectors = X.NumVectors();
  if (NumVectors != Y.NumVectors())
    IFPACK_CHK_ERR(-2);

  Y.PutScalar(0.0);

  vector<int> Indices(MaxNumEntries());
  vector<double> Values(MaxNumEntries());

  for (int i = 0 ; i < NumMyRows() ; ++i) {

    int Nnz;
    ExtractMyRowCopy(i,MaxNumEntries(),Nnz,
		     &Values[0], &Indices[0]);
    if (!TransA) {
      // no transpose first
      for (int j = 0 ; j < NumVectors ; ++j) {
	for (int k = 0 ; k < Nnz ; ++k) {
	  Y[j][i] += Values[k] * X[j][Indices[k]];
	}
      }
    }
    else {
      // transpose here
      for (int j = 0 ; j < NumVectors ; ++j) {
	for (int k = 0 ; k < Nnz ; ++k) {
	  Y[j][Indices[k]] += Values[k] * X[j][i];
	}
      }
    }
  }

  return(0);
}
