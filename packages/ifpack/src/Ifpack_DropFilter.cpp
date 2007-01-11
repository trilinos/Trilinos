#include "Ifpack_ConfigDefs.h"
#include "Ifpack_DropFilter.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

//==============================================================================
Ifpack_DropFilter::Ifpack_DropFilter(const Teuchos::RefCountPtr<Epetra_RowMatrix>& Matrix,
				     double DropTol) :
  A_(Matrix),
  DropTol_(DropTol),
  MaxNumEntries_(0),
  MaxNumEntriesA_(0),
  NumNonzeros_(0)
{
  // use this filter only on serial matrices
  if (A_->Comm().NumProc() != 1) {
    cerr << "Ifpack_DropFilter can be used with Comm().NumProc() == 1" << endl;
    cerr << "only. This class is a tool for Ifpack_AdditiveSchwarz," << endl;
    cerr << "and it is not meant to be used otherwise." << endl;
    exit(EXIT_FAILURE);
  }
  
  if ((A_->NumMyRows() != A_->NumGlobalRows()) ||
      (A_->NumMyRows() != A_->NumMyCols()))
    IFPACK_CHK_ERRV(-2);

  NumRows_ = A_->NumMyRows();
  MaxNumEntriesA_ = A_->MaxNumEntries();

  NumEntries_.resize(NumRows_);
  Indices_.resize(MaxNumEntriesA_);
  Values_.resize(MaxNumEntriesA_);

  vector<int>    Ind(MaxNumEntriesA_);
  vector<double> Val(MaxNumEntriesA_);

  for (int i = 0 ; i < NumRows_ ; ++i) {
    NumEntries_[i] = MaxNumEntriesA_;
    int Nnz;
    IFPACK_CHK_ERRV(ExtractMyRowCopy(i,MaxNumEntriesA_,Nnz,
				     &Val[0], &Ind[0]));
 
    NumEntries_[i] = Nnz;
    NumNonzeros_ += Nnz;
    if (Nnz > MaxNumEntries_)
      MaxNumEntries_ = Nnz;
  }

}

//==============================================================================
int Ifpack_DropFilter::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
		 double *Values, int * Indices) const
{
  if (Length < NumEntries_[MyRow])
    IFPACK_CHK_ERR(-1);

  int Nnz;

  IFPACK_CHK_ERR(A_->ExtractMyRowCopy(MyRow,MaxNumEntriesA_,Nnz,
				     &Values_[0],&Indices_[0]));

  // loop over all nonzero elements of row MyRow,
  // and drop elements below specified threshold.
  // Also, add value to zero diagonal
  int count = 0;
  for (int i = 0 ; i < Nnz ; ++i) {

    // if element is above specified tol, add to the
    // user's defined arrays. Check that we are not
    // exceeding allocated space. Do not drop any diagonal entry.
    if ((Indices_[i] == MyRow) || (IFPACK_ABS(Values_[i]) >= DropTol_)) {
      if (count == Length)
	IFPACK_CHK_ERR(-1);
      Values[count] = Values_[i];
      Indices[count] = Indices_[i];
      count++;
    }
  }

  NumEntries = count;
  return(0);
}

//==============================================================================
int Ifpack_DropFilter::
ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  IFPACK_CHK_ERR(A_->ExtractDiagonalCopy(Diagonal));
  return(0);
}

//==============================================================================
int Ifpack_DropFilter::
Multiply(bool TransA, const Epetra_MultiVector& X, 
	 Epetra_MultiVector& Y) const
{
  // NOTE: I suppose that the matrix has been localized,
  // hence all maps are trivial.
  int NumVectors = X.NumVectors();
  if (NumVectors != Y.NumVectors())
    IFPACK_CHK_ERR(-1);

  Y.PutScalar(0.0);

  vector<int> Indices(MaxNumEntries_);
  vector<double> Values(MaxNumEntries_);

  for (int i = 0 ; i < NumRows_ ; ++i) {

    int Nnz;
    ExtractMyRowCopy(i,MaxNumEntries_,Nnz,
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

//==============================================================================
int Ifpack_DropFilter::
Solve(bool Upper, bool Trans, bool UnitDiagonal, 
      const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(-99);
}

//==============================================================================
int Ifpack_DropFilter::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_RETURN(Multiply(UseTranspose(),X,Y));
}

//==============================================================================
int Ifpack_DropFilter::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(-99);
}

//==============================================================================
int Ifpack_DropFilter::InvRowSums(Epetra_Vector& x) const
{
  IFPACK_CHK_ERR(-1);
}

//==============================================================================
int Ifpack_DropFilter::InvColSums(Epetra_Vector& x) const
{
  IFPACK_CHK_ERR(-1);
}
