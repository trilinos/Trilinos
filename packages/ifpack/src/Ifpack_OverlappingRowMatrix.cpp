#include "Ifpack_ConfigDefs.h"
#include "Ifpack_OverlappingRowMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"

// ====================================================================== 
Ifpack_OverlappingRowMatrix::
Ifpack_OverlappingRowMatrix(const Teuchos::RefCountPtr<const Epetra_RowMatrix>& Matrix,
                            int OverlapLevel)  :
  Matrix_(Matrix),
  OverlapLevel_(OverlapLevel)
{
  // should not be here if no overlap
  if (OverlapLevel == 0)
    IFPACK_CHK_ERRV(-1);

  // nothing to do as well with one process
  if (Comm().NumProc() == 1)
    IFPACK_CHK_ERRV(-1);
  
  NumMyRowsA_ = A().NumMyRows();

  // construct the external matrix
  vector<int> ExtElements; 

  Epetra_Map* TmpMap = 0;
  Epetra_CrsMatrix* TmpMatrix = 0; 
  Epetra_Import* TmpImporter = 0; 

  // importing rows corresponding to elements that are 
  // in ColMap, but not in RowMap 
  const Epetra_Map *RowMap; 
  const Epetra_Map *ColMap; 

  for (int overlap = 0 ; overlap < OverlapLevel ; ++overlap) {
    if (TmpMatrix) {
      RowMap = &(TmpMatrix->RowMatrixRowMap()); 
      ColMap = &(TmpMatrix->RowMatrixColMap()); 
    }
    else {
      RowMap = &(A().RowMatrixRowMap()); 
      ColMap = &(A().RowMatrixColMap()); 
    }

    int size = ColMap->NumMyElements() - RowMap->NumMyElements(); 
    vector<int> list(size); 

    int count = 0; 

    // define the set of rows that are in ColMap but not in RowMap 
    for (int i = 0 ; i < ColMap->NumMyElements() ; ++i) { 
      int GID = ColMap->GID(i); 
      if (A().RowMatrixRowMap().LID(GID) == -1) { 
        vector<int>::iterator pos 
          = find(ExtElements.begin(),ExtElements.end(),GID); 
        if (pos == ExtElements.end()) { 
          ExtElements.push_back(GID);
          list[count] = GID; 
          ++count; 
        } 
      } 
    } 

    if (TmpMap)
      delete TmpMap;
    TmpMap = new Epetra_Map(-1,count, &list[0],0,Comm()); 

    if (TmpMatrix)
      delete TmpMatrix;
    TmpMatrix = new Epetra_CrsMatrix(Copy,*TmpMap,0); 

    if (TmpImporter)
      delete TmpImporter;
    TmpImporter = new Epetra_Import(*TmpMap,A().RowMatrixRowMap()); 

    TmpMatrix->Import(A(),*TmpImporter,Insert); 
    TmpMatrix->FillComplete(A().OperatorDomainMap(),*TmpMap); 

  }

  // clean up memory
  delete TmpMap;
  delete TmpMatrix;
  delete TmpImporter;

  // build the map containing all the nodes (original
  // matrix + extended matrix)
  vector<int> list(NumMyRowsA_ + ExtElements.size());
  for (int i = 0 ; i < NumMyRowsA_ ; ++i)
    list[i] = A().RowMatrixRowMap().GID(i);
  for (int i = 0 ; i < (int)ExtElements.size() ; ++i)
    list[i + NumMyRowsA_] = ExtElements[i];

  Map_ = Teuchos::rcp( new Epetra_Map(-1, NumMyRowsA_ + ExtElements.size(),
				      &list[0], 0, Comm()) );
  // now build the map corresponding to all the external nodes
  // (with respect to A().RowMatrixRowMap().
  ExtMap_ = Teuchos::rcp( new Epetra_Map(-1,ExtElements.size(),
					 &ExtElements[0],0,A().Comm()) );
  ExtMatrix_ = Teuchos::rcp( new Epetra_CrsMatrix(Copy,*ExtMap_,*Map_,0) ); 

  ExtImporter_ = Teuchos::rcp( new Epetra_Import(*ExtMap_,A().RowMatrixRowMap()) ); 
  ExtMatrix_->Import(A(),*ExtImporter_,Insert); 
  ExtMatrix_->FillComplete(A().OperatorDomainMap(),*Map_);

  Importer_ = Teuchos::rcp( new Epetra_Import(*Map_,A().RowMatrixRowMap()) );

  // fix indices for overlapping matrix
  NumMyRowsB_ = B().NumMyRows();
  NumMyRows_ = NumMyRowsA_ + NumMyRowsB_;
  NumMyCols_ = NumMyRows_;
  
  NumMyDiagonals_ = A().NumMyDiagonals() + B().NumMyDiagonals();
  
  NumMyNonzeros_ = A().NumMyNonzeros() + B().NumMyNonzeros();
  Comm().SumAll(&NumMyNonzeros_,&NumGlobalNonzeros_,1);
  MaxNumEntries_ = A().MaxNumEntries();
  
  if (MaxNumEntries_ < B().MaxNumEntries())
    MaxNumEntries_ = B().MaxNumEntries();

}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
NumMyRowEntries(int MyRow, int & NumEntries) const
{
  if (MyRow < NumMyRowsA_)
    return(A().NumMyRowEntries(MyRow,NumEntries));
  else
    return(B().NumMyRowEntries(MyRow - NumMyRowsA_, NumEntries));
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, 
                 int * Indices) const
{
  int ierr;
  if (MyRow < NumMyRowsA_)
    ierr = A().ExtractMyRowCopy(MyRow,Length,NumEntries,Values,Indices);
  else
    ierr = B().ExtractMyRowCopy(MyRow - NumMyRowsA_,Length,NumEntries,
                                Values,Indices);

  IFPACK_RETURN(ierr);
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  IFPACK_CHK_ERR(-1);
}


// ======================================================================
int Ifpack_OverlappingRowMatrix::
Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();
  vector<int> Ind(MaxNumEntries_);
  vector<double> Val(MaxNumEntries_);

  Y.PutScalar(0.0);

  // matvec with A (local rows)
  for (int i = 0 ; i < NumMyRowsA_ ; ++i) {
    for (int k = 0 ; k < NumVectors ; ++k) {
      int Nnz;
      IFPACK_CHK_ERR(A().ExtractMyRowCopy(i,MaxNumEntries_,Nnz, 
                                          &Val[0], &Ind[0]));
      for (int j = 0 ; j < Nnz ; ++j) {
        Y[k][i] += Val[j] * X[k][Ind[j]];
      }
    }
  }

  // matvec with B (overlapping rows)
  for (int i = 0 ; i < NumMyRowsB_ ; ++i) {
    for (int k = 0 ; k < NumVectors ; ++k) {
      int Nnz;
      IFPACK_CHK_ERR(B().ExtractMyRowCopy(i,MaxNumEntries_,Nnz, 
                                          &Val[0], &Ind[0]));
      for (int j = 0 ; j < Nnz ; ++j) {
        Y[k][i + NumMyRowsA_] += Val[j] * X[k][Ind[j]];
      }
    }
  }
  return(0);
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(Multiply(UseTranspose(),X,Y));
  return(0);
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(-1);
}

// ======================================================================
Epetra_RowMatrix& Ifpack_OverlappingRowMatrix::B() const
{
  return(*ExtMatrix_);
}

// ======================================================================
const Epetra_BlockMap& Ifpack_OverlappingRowMatrix::Map() const
{
  return(*Map_);
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ImportMultiVector(const Epetra_MultiVector& X, Epetra_MultiVector& OvX,
                  Epetra_CombineMode CM)
{
  OvX.Import(X,*Importer_,CM);
  return(0);
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ExportMultiVector(const Epetra_MultiVector& OvX, Epetra_MultiVector& X,
                  Epetra_CombineMode CM)
{
  X.Export(OvX,*Importer_,CM);
  return(0);
}

