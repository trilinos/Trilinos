#include "Ifpack_ConfigDefs.h"

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Map.h"
#include "Ifpack_LocalRowMatrix.h"

//==============================================================================
Ifpack_LocalRowMatrix::Ifpack_LocalRowMatrix(const Epetra_RowMatrix* Matrix) :
  Matrix_(Matrix)
{
  Label_ = new char[80];
  sprintf(Label_,"Ifpack_LocalRowMatrix");

#ifdef HAVE_MPI
  SerialComm_ = new Epetra_MpiComm(MPI_COMM_SELF);
#else
  SerialComm_ = new Epetra_SerialComm;
#endif

  NumRows_ = Matrix->NumMyRows();
  // FIXME: compute the local ones
  NumCols_ = NumRows_;
  Map_ = new Epetra_Map(NumRows_,0,*SerialComm_);

  NumEntries_.resize(NumRows_);

  // FIXME: not all combinations of RowMatrixRowMap() and RangeMap() are
  // allowed... (Same for columns)
  // loop over all local rows to compute the nonzeros in each row

  MaxNumEntries_ = Matrix_->MaxNumEntries();
  NumDiagonals_ = 0;
  NumNonzeros_ = 0;

  Indices_ = new int[MaxNumEntries_];
  Values_ = new double[MaxNumEntries_];

  for (int i = 0 ; i < NumRows_ ; ++i) {
    
    NumEntries_[i] = 0;
    int Nnz;
    Matrix_->ExtractMyRowCopy(i,MaxNumEntries_,Nnz,Values_,Indices_);

    for (int j = 0 ; j < Nnz ; ++j) {
      if (Indices_[j] == i)
	++NumDiagonals_;
      if (Indices_[j] < NumRows_ ) {
	++NumNonzeros_;
	NumEntries_[i]++;
      }
    }
  }

  // compute MaxNumEntries_ for the local matrix
  MaxNumLocalEntries_ = 0;
  for (int i = 0 ; i < NumRows_ ; ++i) {
    if (NumEntries_[i] > MaxNumLocalEntries_)
      MaxNumLocalEntries_ = NumEntries_[i];
  }

}

//==============================================================================
Ifpack_LocalRowMatrix::~Ifpack_LocalRowMatrix()
{
  if (Label_)
    delete Label_;
}

//==============================================================================
int Ifpack_LocalRowMatrix::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
		 double *Values, int * Indices) const
{
  if ((MyRow < 0) || (MyRow >= NumRows_)) {
    IFPACK_CHK_ERR(-1); // range not valid
  }

  if (Length < NumEntries_[MyRow])
    return(-1);

  // always extract using the object Values_ and Indices_.
  // This is because I need more space than that given by
  // the user (for the external nodes)
  int Nnz;
  int ierr = Matrix_->ExtractMyRowCopy(MyRow,MaxNumEntries_,Nnz,
				       Values_,Indices_);

  // this should never happen
  if (ierr < 0) {
    IFPACK_CHK_ERR(ierr);
  }

  // populate the user's vectors
  NumEntries = 0;
  for (int j = 0 ; j < Nnz ; ++j) {
    if (Indices_[j] < NumRows_ ) {
      Indices[NumEntries] = Indices_[j];
      Values[NumEntries] = Values_[j];
      ++NumEntries;
    }
  }
    
  assert (NumEntries == NumEntries_[MyRow]);

  return(0);

}

//==============================================================================
int Ifpack_LocalRowMatrix::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  if (!Map_->SameAs(Diagonal.Map())) {
    IFPACK_CHK_ERR(-1); 
  }

  for (int i = 0 ; i < NumRows_ ; ++i) {
    
    int Nnz;
    Matrix_->ExtractMyRowCopy(i,MaxNumEntries_,Nnz,
			      Values_,Indices_);

    assert (Nnz = NumEntries_[i]);

    for (int j = 0 ; j < Nnz ; ++j) {
      if (Indices_[j] == i )
	Diagonal[i] = Values_[j];
    }
  }

  return(0);

}

//==============================================================================
int Ifpack_LocalRowMatrix::Apply(const Epetra_MultiVector& X,
	  Epetra_MultiVector& Y) const 
{

  if (!Map_->SameAs(X.Map())) {
    IFPACK_CHK_ERR(-1); 
  }

  if (!Map_->SameAs(Y.Map())) {
    IFPACK_CHK_ERR(-1); 
  }

  if (X.NumVectors() != Y.NumVectors()) {
    IFPACK_CHK_ERR(-2);
  }

  Y.PutScalar(0.0);
  int NumVectors = Y.NumVectors();

  for (int i = 0 ; i < NumRows_ ; ++i) {
    
    int Nnz;
    Matrix_->ExtractMyRowCopy(i,MaxNumEntries_,Nnz,Values_,Indices_);

    for (int j = 0 ; j < Nnz ; ++j) {
      if (Indices_[j] < NumRows_ ) {
	for (int k = 0 ; k < NumVectors ; ++k)
	  Y[k][i] = Y[k][i] + Values_[j] * X[k][Indices_[j]];
      }
    }
  }

  return(0);
}

//==============================================================================
int Ifpack_LocalRowMatrix::ApplyInverse(const Epetra_MultiVector& X,
		 Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(-1); // not implemented
}

