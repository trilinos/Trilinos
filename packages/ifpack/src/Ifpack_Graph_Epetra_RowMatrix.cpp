#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Graph.h"
#include "Ifpack_Graph_Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"

//==============================================================================
Ifpack_Graph_Epetra_RowMatrix::Ifpack_Graph_Epetra_RowMatrix(const Epetra_RowMatrix* RowMatrix) :
RowMatrix_(RowMatrix)
{
  NumMyRows_ = RowMatrix_->NumMyRows();
  NumMyCols_ = RowMatrix_->NumMyCols();
  NumGlobalRows_ = RowMatrix_->NumGlobalRows();
  NumGlobalCols_ = RowMatrix_->NumGlobalCols();
  MaxNumIndices_ = RowMatrix_->MaxNumEntries();

  Values_.resize(MaxNumIndices_);
}

//==============================================================================
Ifpack_Graph_Epetra_RowMatrix::~Ifpack_Graph_Epetra_RowMatrix()
{
}

//==============================================================================
const Epetra_Comm& Ifpack_Graph_Epetra_RowMatrix::Comm() const
{
  return(RowMatrix_->Comm());
}

//==============================================================================
bool Ifpack_Graph_Epetra_RowMatrix::Filled() const
{
  return(RowMatrix_->Filled());
}
 
//==============================================================================
int Ifpack_Graph_Epetra_RowMatrix::GRID(int LRID) const
{
  return(RowMatrix_->RowMatrixRowMap().GID(LRID));
}

//==============================================================================
int Ifpack_Graph_Epetra_RowMatrix::GCID(int LCID) const
{
  return(RowMatrix_->RowMatrixColMap().GID(LCID));
}

//==============================================================================
int Ifpack_Graph_Epetra_RowMatrix::LRID(int GRID) const
{
  return(RowMatrix_->RowMatrixRowMap().LID(GRID));
}

//==============================================================================
int Ifpack_Graph_Epetra_RowMatrix::LCID(int GCID) const
{
  return(RowMatrix_->RowMatrixColMap().LID(GCID));
}

//==============================================================================
int Ifpack_Graph_Epetra_RowMatrix::
ExtractMyRowCopy(int MyRow, int LenOfIndices, 
		 int &NumIndices, int *Indices) const
{
  return(RowMatrix_->ExtractMyRowCopy(MyRow, LenOfIndices,
				      NumIndices, &Values_[0],
				      Indices));
}

//==============================================================================
int Ifpack_Graph_Epetra_RowMatrix::NumMyNonzeros() const
{
  return(RowMatrix_->NumMyNonzeros());
}


