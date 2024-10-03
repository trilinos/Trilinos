// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Graph.hpp"
#include "Ifpack2_Graph_Tpetra_RowMatrix.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"

//==============================================================================
Ifpack2_Graph_Tpetra_RowMatrix::Ifpack2_Graph_Tpetra_RowMatrix(const Teuchos::RCP<const Tpetra_RowMatrix>& RowMatrix) :
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
const Tpetra_Comm& Ifpack2_Graph_Tpetra_RowMatrix::Comm() const
{
  return(RowMatrix_->Comm());
}

//==============================================================================
bool Ifpack2_Graph_Tpetra_RowMatrix::Filled() const
{
  return(RowMatrix_->Filled());
}
 
//==============================================================================
int Ifpack2_Graph_Tpetra_RowMatrix::GRID(int LRID_in) const
{
  return(RowMatrix_->RowMatrixRowMap().GID(LRID_in));
}

//==============================================================================
int Ifpack2_Graph_Tpetra_RowMatrix::GCID(int LCID_in) const
{
  return(RowMatrix_->RowMatrixColMap().GID(LCID_in));
}

//==============================================================================
int Ifpack2_Graph_Tpetra_RowMatrix::LRID(int GRID_in) const
{
  return(RowMatrix_->RowMatrixRowMap().LID(GRID_in));
}

//==============================================================================
int Ifpack2_Graph_Tpetra_RowMatrix::LCID(int GCID_in) const
{
  return(RowMatrix_->RowMatrixColMap().LID(GCID_in));
}

//==============================================================================
int Ifpack2_Graph_Tpetra_RowMatrix::
ExtractMyRowCopy(int MyRow, int LenOfIndices, 
		 int &NumIndices, int *Indices) const
{
  return(RowMatrix_->ExtractMyRowCopy(MyRow, LenOfIndices,
				      NumIndices, &Values_[0],
				      Indices));
}

//==============================================================================
int Ifpack2_Graph_Tpetra_RowMatrix::NumMyNonzeros() const
{
  return(RowMatrix_->NumMyNonzeros());
}

// ======================================================================
ostream& Ifpack2_Graph_Tpetra_RowMatrix::Print(std::ostream& os) const
{

  if (Comm().MyPID())
    return(os);

  os << "================================================================================" << endl;
  os << "Ifpack2_Graph_Tpetra_RowMatrix" << endl;
  os << "Number of local rows  = " << NumMyRows_ << endl;
  os << "Number of global rows = " << NumGlobalRows_ << endl;
  os << "================================================================================" << endl;

  return(os);

}

