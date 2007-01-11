#include "Ifpack_ConfigDefs.h"
#include "Ifpack_ReorderFilter.h"
#include "Ifpack_Reordering.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

//==============================================================================
Ifpack_ReorderFilter::Ifpack_ReorderFilter(const Teuchos::RefCountPtr<Epetra_RowMatrix>& Matrix,
                                           const Teuchos::RefCountPtr<Ifpack_Reordering>& Reordering) :
  A_(Matrix),
  Reordering_(Reordering),
  NumMyRows_(Matrix->NumMyRows()),
  MaxNumEntries_(Matrix->MaxNumEntries())
{
}

//==============================================================================
Ifpack_ReorderFilter::Ifpack_ReorderFilter(const Ifpack_ReorderFilter& RHS) :
  A_(Matrix()),
  Reordering_(Reordering()),
  NumMyRows_(RHS.NumMyRows()),
  MaxNumEntries_(RHS.MaxNumEntries())
{
  strcpy(Label_,RHS.Label());
}

//==============================================================================
Ifpack_ReorderFilter& 
Ifpack_ReorderFilter::operator=(const Ifpack_ReorderFilter& RHS)
{
  if (this == &RHS)
    return (*this);

  A_ = RHS.Matrix();

  Reordering_ = RHS.Reordering();
  MaxNumEntries_ = RHS.MaxNumEntries();
  NumMyRows_ = RHS.NumMyRows();

  strcpy(Label_,RHS.Label());
  return(*this);
}

//==============================================================================
int Ifpack_ReorderFilter::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
                 double *Values, int * Indices) const
{
  int MyReorderdRow = Reordering_->InvReorder(MyRow);

  IFPACK_CHK_ERR(Matrix()->ExtractMyRowCopy(MyReorderdRow,MaxNumEntries_,
                                           NumEntries, Values,Indices));

  // suppose all elements are local. Note that now
  // Indices can have indices in non-increasing order.
  for (int i = 0 ; i < NumEntries ; ++i) {
    Indices[i] = Reordering_->Reorder(Indices[i]);
  }

  return(0);
}

//==============================================================================
int Ifpack_ReorderFilter::
ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  Epetra_Vector DiagonalTilde(Diagonal.Map());
  IFPACK_CHK_ERR(Matrix()->ExtractDiagonalCopy(DiagonalTilde));
  IFPACK_CHK_ERR((Reordering_->P(DiagonalTilde,Diagonal)));
  return(0);
}

//==============================================================================
int Ifpack_ReorderFilter::
Multiply(bool TransA, const Epetra_MultiVector& X, 
         Epetra_MultiVector& Y) const
{
  // need two additional vectors
  Epetra_MultiVector Xtilde(X.Map(),X.NumVectors());
  Epetra_MultiVector Ytilde(Y.Map(),Y.NumVectors());
  // bring X back to original ordering
  Reordering_->Pinv(X,Xtilde);
  // apply original matrix
  IFPACK_CHK_ERR(Matrix()->Multiply(TransA,Xtilde,Ytilde));
  // now reorder result
  Reordering_->P(Ytilde,Y);


  return(0);
}

//==============================================================================
int Ifpack_ReorderFilter::
Solve(bool Upper, bool Trans, bool UnitDiagonal, 
      const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(-98);
}

//==============================================================================
int Ifpack_ReorderFilter::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_RETURN(Multiply(UseTranspose(),X,Y));
}
