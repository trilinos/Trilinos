#include "Ifpack_ConfigDefs.h"
#include "Ifpack_ReorderFilter.h"
#include "Ifpack_Reordering.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include <algorithm>

//==============================================================================
Ifpack_ReorderFilter::Ifpack_ReorderFilter(Epetra_RowMatrix* Matrix,
					   Ifpack_Reordering* Reordering,
					   const bool LightWeigth) :
  A_(Matrix),
  Reordering_(Reordering),
  MaxNumEntries_(Matrix->MaxNumEntries()),
  NumMyRows_(Matrix->NumMyRows())
{
  Indices_.resize(MaxNumEntries_);
  Values_.resize(MaxNumEntries_);
}

//==============================================================================
Ifpack_ReorderFilter::Ifpack_ReorderFilter(const Ifpack_ReorderFilter& RHS) :
  A_(&RHS.Matrix()),
  Reordering_(&RHS.Reordering()),
  MaxNumEntries_(RHS.MaxNumEntries()),
  NumMyRows_(RHS.NumMyRows())
{
  Indices_.resize(MaxNumEntries_);
  Values_.resize(MaxNumEntries_);
  strcpy(Label_,RHS.Label());
}

//==============================================================================
Ifpack_ReorderFilter::~Ifpack_ReorderFilter()
{
}

//==============================================================================
Ifpack_ReorderFilter& 
Ifpack_ReorderFilter::operator=(const Ifpack_ReorderFilter& RHS)
{
  if (this == &RHS)
    return (*this);

  A_ = &RHS.Matrix();

  Reordering_ = &RHS.Reordering();
  MaxNumEntries_ = RHS.MaxNumEntries();
  NumMyRows_ = RHS.NumMyRows();

  Indices_.resize(MaxNumEntries_);
  Values_.resize(MaxNumEntries_);
  strcpy(Label_,RHS.Label());
}

//==============================================================================
int Ifpack_ReorderFilter::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
		 double *Values, int * Indices) const
{
  int MyReorderdRow = Reordering().InvReorder(MyRow);

  // FIXME: delete Indices_ Values_
  IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(MyReorderdRow,MaxNumEntries_,
					   NumEntries, Values,Indices));

  for (int i = 0 ; i < NumEntries ; ++i) {
    // ignore off-process elements, leave them as they are
    if (Indices[i] >= NumMyRows_)
      continue;
    Indices[i] = Reordering().Reorder(Indices[i]);
  }

  return(0);
}

//==============================================================================
int Ifpack_ReorderFilter::
ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  Epetra_Vector DiagonalTilde(Diagonal);
  IFPACK_CHK_ERR(Matrix().ExtractDiagonalCopy(DiagonalTilde));
  Reordering().P(DiagonalTilde,Diagonal);
}

//==============================================================================
int Ifpack_ReorderFilter::
Multiply(bool TransA, const Epetra_MultiVector& X, 
	 Epetra_MultiVector& Y) const
{
  Epetra_MultiVector Xtilde(X);
  Epetra_MultiVector Ytilde(Y);
  // bring X back to original ordering
  Reordering().Pinv(X,Xtilde);
  // apply original matrix
  IFPACK_CHK_ERR(Matrix().Multiply(TransA,Xtilde,Ytilde));
  // now reorder result
  Reordering().P(Ytilde,Y);


  return(0);
}

//==============================================================================
int Ifpack_ReorderFilter::
Solve(bool Upper, bool Trans, bool UnitDiagonal, 
      const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(-1);
}

//==============================================================================
int Ifpack_ReorderFilter::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(Multiply(false,X,Y));
  return(0);
}

//==============================================================================
int Ifpack_ReorderFilter::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return(-1); // NOT IMPLEMENTED AT THIS STAGE
}
