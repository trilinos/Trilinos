
#include "Epetra_CrsKundertSparse.h"

extern "C" {
#include "spmatrix.h"
  int spFactorAndSolve(char *eMatrix, double *RHS); // Sparse has no prototype for this function

}

#include "Epetra_LinearProblem.h"
#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"


Epetra_CrsKundertSparse::Epetra_CrsKundertSparse( Epetra_LinearProblem * Problem,
			      const double RelThreshold,
			      const double AbsThreshold,
			      const int DiagPivoting)

  : RelThreshold_(RelThreshold),
    AbsThreshold_(AbsThreshold),
    DiagPivoting_(DiagPivoting),
    FirstSolve_(true)
{

  A_ = dynamic_cast<Epetra_CrsMatrix *> (Problem->GetOperator());
  X_ = Problem->GetLHS();
  B_ = Problem->GetRHS();
  if (A_->Comm().NumProc()!=1) 
    A_->ReportError("Can only use Spice Sparse solver single PE", -1);
  NumMyRows_ = A_->NumMyRows();
  NumMyCols_ = A_->NumMyCols();
  NumGlobalRows_ = A_->NumGlobalRows();
  NumGlobalCols_ = A_->NumGlobalCols();


  if (NumGlobalRows_ != NumGlobalCols_) A_->ReportError("Matrix must be square", -2);

  // Create a Sparse matrix
  int err = 0;
  Matrix_ = (char *) spCreate (NumGlobalRows_, 0, &err);
  if (err!=0) A_->ReportError("Error occurred in Spice Sparse spCreate", err);

  int NumEntries;
  int * Indices;
  double * Values;
    
  int NumMyNonzeros = A_->NumMyNonzeros();
  
  int curValue = 0;
  addr_list_ = new double *[NumMyNonzeros];
  for (int i=0; i<NumMyRows_; i++) {
    // int curGRID = A_->RowMap().GID(i); // Needed for parallel (later)
    // View of current row
    int ierr = A_->ExtractMyRowView(i, NumEntries, Values, Indices); 
    if (ierr!=0) A_->ReportError("Error occurred in ExtractMyRowView", ierr);
    for (int j=0; j<NumEntries; j++) {
      // int columnIndex = A_->ImportMap().GID(Indices[j]);// parallel (later)
      int columnIndex = Indices[j];
      // Register this entry into Sparse matrix
      double * p = (double *) spGetElement(Matrix_, i+1, columnIndex+1);
      *p = Values[j];
      addr_list_[curValue++] = p;
    }
  }

}

Epetra_CrsKundertSparse::~Epetra_CrsKundertSparse() {
  deleteArrays();
}

void Epetra_CrsKundertSparse::deleteArrays() {

  if (Matrix_!=0) {
    spDestroy(Matrix_);
    Matrix_ = 0;
  }
  if (addr_list_!=0) {
    delete [] addr_list_;
    addr_list_ = 0;
  }

}

int Epetra_CrsKundertSparse::Solve() {

  // If not first call to solver, we need to copy values to solver matrix.
  // NOTE: We are proceeding through the matrix in the same order as it was
  //       constructed.  As a result, we do not need to access index information.
  if (!FirstSolve_) {
    spClear (Matrix_); // Clear previous factorization and matrxi values
    int curValue = 0;
    int NumEntries;
    double * Values;
    for (int i=0; i<NumMyRows_; i++) {
      // int curGRID = A_->RowMap().GID(i); // Needed for parallel (later)
      // View of current row
      EPETRA_CHK_ERR(A_->ExtractMyRowView(i, NumEntries, Values)); 
      for (int j=0; j<NumEntries; j++)
	*(addr_list_[curValue++]) = Values[j];
    }
  }

 /* Create right-hand side matrix B. */
  double * rhs;
  double * solution;
  int LDA_x, LDA_b;
  B_->ExtractView( &rhs, &LDA_b );
  X_->ExtractView( &solution, &LDA_x );
  rhs--; solution--; // adjust for 1-based indexing

  if (B_->NumVectors()>1) EPETRA_CHK_ERR(1); // Can only handle one RHS at this time

  if (FirstSolve_) {
   spOrderAndFactor (Matrix_, rhs, RelThreshold_, AbsThreshold_, DiagPivoting_);
   spSolve (Matrix_, rhs, solution, NULL, NULL);
  }
  else {
    *X_ = *B_; // Copy B to X
    spFactorAndSolve (Matrix_, solution);
  }

  FirstSolve_ = false;
  return 0;
}

