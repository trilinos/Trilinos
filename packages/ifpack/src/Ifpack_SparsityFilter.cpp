#include "Ifpack_ConfigDefs.h"
#include "Ifpack_SparsityFilter.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include <algorithm>

//==============================================================================
Ifpack_SparsityFilter::Ifpack_SparsityFilter(Epetra_RowMatrix* Matrix,
					     int MaxRowEntries, 
					     int MaxBandwidth) :
  A_(*Matrix),
  MaxBandwidth_(MaxBandwidth),
  MaxRowEntries_(MaxRowEntries),
  MaxNumEntries_(0),
  UpperTriangular_(true),
  LowerTriangular_(true),
  InvRowSum_(0),
  InvColSum_(0),
  NumMyNonzeros_(0),
  NumGlobalNonzeros_(0),
  NormInf_(0.0),
  NormOne_(0.0)
{
  MaxNumEntries_ = A_.MaxNumEntries();
  Indices_.resize(MaxNumEntries_);
  Values_.resize(MaxNumEntries_);

  // allocates space for row and colum sum
  InvRowSum_ = new Epetra_Vector(A_.RowMatrixRowMap());
  InvColSum_ = new Epetra_Vector(A_.RowMatrixColMap());
  assert (InvRowSum_ != 0);
  assert (InvColSum_ != 0);
  InvRowSum_->PutScalar(0.0);
  InvColSum_->PutScalar(0.0);

  // computes the number of nonzeros, the maximum nonzeros 
  // per row, the inverse of row
  // and column sum, and checks whether the dropped matrix
  // is lower or upper triangular.

  vector<int>    Indices(A_.MaxNumEntries());
  vector<double> Values(A_.MaxNumEntries());

  for (int i = 0 ; i < A_.NumMyRows() ; ++i) {
    int Nnz;
    IFPACK_CHK_ERRV(ExtractMyRowCopy(i,MaxNumEntries_,Nnz,
				     &Values[0], &Indices[0]));
    NumMyNonzeros_ += Nnz;
    if (Nnz > MaxNumEntries_)
      MaxNumEntries_ = Nnz;
    // check whether the dropped matrix is lower/upper triangular
    for (int j = 0 ; j < Nnz ; ++j) {
      if (Indices[j] > i)
	LowerTriangular_ = false;
      if (Indices[j] < i)
	UpperTriangular_ = false;
      // add contribution to the row and col sum
      (*InvRowSum_)[i] += Values[0];
      (*InvColSum_)[Indices[j]] += Values[0];
    }

  }
  // sum global values
  A_.Comm().SumAll(&NumMyNonzeros_,&NumGlobalNonzeros_,1);

  IFPACK_CHK_ERRV(InvRowSum_->MaxValue(&NormInf_));
  IFPACK_CHK_ERRV(InvColSum_->MaxValue(&NormOne_));
  // compute the inverse of row and col sums
  (InvRowSum_->Reciprocal(*InvRowSum_));
  (InvColSum_->Reciprocal(*InvColSum_));
}

//==============================================================================
Ifpack_SparsityFilter::~Ifpack_SparsityFilter()
{
  if (InvRowSum_)
    delete InvRowSum_;
  if (InvColSum_)
    delete InvColSum_;
}

//==============================================================================
int Ifpack_SparsityFilter::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
		 double *Values, int * Indices) const
{
  int NumMyRows = A_.NumMyRows();
  int Nnz;

  IFPACK_CHK_ERR(A_.ExtractMyRowCopy(MyRow,A_.MaxNumEntries(),Nnz,
				     &Values_[0],&Indices_[0]));

  // if MaxBandwidth is not trivial, drop all off-bandwidth
  // elements
  if (MaxBandwidth_ != -1) {
    for (int i = 0 ; i < Nnz ; ++i) {
      if (IFPACK_ABS(Indices_[i] - MyRow) >= MaxBandwidth_)
	Indices_[i] = NumMyRows + 1; // it will be ignored
    }
  }

  double Threshold = 0.0;
  // if MaxEntries is not trivial, keep only largest-magnitude
  // elements, only if I have more than MaxRowEntries elements.
  if ((MaxRowEntries_ != -1) && (MaxRowEntries_ < Nnz)) {
    vector<double> Values2(Nnz);
    for (int i = 0 ; i < Nnz ; ++i) {
      if (Indices_[i] >= NumMyRows)
	Values2[i] = 0;
      else
	Values2[i] = IFPACK_ABS(Values_[i]);
    }
    sort(Values2.rbegin(),Values2.rend());
    Threshold = Values2[MaxRowEntries_ - 1];
  }

  // loop over all nonzero elements of row MyRow,
  // and drop elements below specified threshold.
  // Also, add value to zero diagonal
  int count = 0;
  bool FoundDiag = false;
  for (int i = 0 ; i < Nnz ; ++i) {
    // drop off-process connections
    if (Indices_[i] >= NumMyRows)
      continue;

    if (IFPACK_ABS(Values_[i]) < Threshold)
      continue;

    Values[count] = Values_[i];
    Indices[count] = Indices_[i];
    count++;
  }

  NumEntries = count;
  return(0);
}

//==============================================================================
int Ifpack_SparsityFilter::
ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  IFPACK_CHK_ERR(A_.ExtractDiagonalCopy(Diagonal));
}

//==============================================================================
int Ifpack_SparsityFilter::
Multiply(bool TransA, const Epetra_MultiVector& X, 
	 Epetra_MultiVector& Y) const
{
  // FIXME: I do not work with funky Range and Domain maps

  int NumVectors = X.NumVectors();
  if (NumVectors != Y.NumVectors())
    IFPACK_CHK_ERR(-1);

  Y.PutScalar(0.0);

  vector<int> Indices(MaxNumEntries_);
  vector<double> Values(MaxNumEntries_);

  // NOTE: at this point the off-process elements are ignored.
  // TO BE FIXED???
  for (int i = 0 ; i < A_.NumMyRows() ; ++i) {

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
int Ifpack_SparsityFilter::
Solve(bool Upper, bool Trans, bool UnitDiagonal, 
      const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (!Filled()) {
    IFPACK_CHK_ERR(-1); // Matrix must be filled.
  }
  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-1);
  if ((Upper) && (!UpperTriangular()))
    IFPACK_CHK_ERR(-2);
  if ((!Upper) && (!LowerTriangular()))
    IFPACK_CHK_ERR(-3);

  int NumVectors = X.NumVectors();
  vector<int> Indices(MaxNumEntries_);
  vector<double> Values(MaxNumEntries_);

  // extract a copy of the diagonal
  Epetra_Vector Diagonal(RowMatrixRowMap());
  IFPACK_CHK_ERR(ExtractDiagonalCopy(Diagonal));
  // check that no diagonal is zero
  // FIXME: do be done

  Y.PutScalar(0.0);

  if (Upper) {
    // solve the upper triangular linear system
    for (int i = NumMyRows() - 1 ; i >= 0 ; --i) {
      int Nnz;
      ExtractMyRowCopy(i,MaxNumEntries_,Nnz,&Values[0],&Indices[0]);
      // backsolve for all vectors
      // Diagonal value for Y is now zero, so I can
      // skip an annoying if at this stage
      for (int j = 0 ; j < NumVectors ; ++j) {
	for (int k = 0 ; k < Nnz ; ++k) {
	  Y[j][i] -= Values[k] * X[j][Indices[j]];
	}
	Y[j][i] /= Diagonal[i];
      }
    }	
  }
  else {
    // solve the lower triangular linear system
    for (int i = 0 ; i < NumMyRows() ; ++i) {
      int Nnz;
      ExtractMyRowCopy(i,MaxNumEntries_,Nnz,&Values_[0],&Indices[0]);
      // backsolve for all vectors
      // Diagonal value for Y is now zero, so I can
      // skip an annoying if at this stage
      for (int j = 0 ; j < NumVectors ; ++j) {
	for (int k = 0 ; k < Nnz ; ++k) {
	  Y[j][i] -= Values[k] * X[j][Indices[j]];
	}
	Y[j][i] /= Diagonal[i];
      }
    }	
  }

  return(0);
}

//==============================================================================
int Ifpack_SparsityFilter::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(Multiply(false,X,Y));
  return(0);
}

//==============================================================================
int Ifpack_SparsityFilter::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return(-1); // NOT IMPLEMENTED AT THIS STAGE
}
