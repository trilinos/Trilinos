/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Graph.h"
#include "Ifpack_Graph_Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"

//==============================================================================
Ifpack_Graph_Epetra_RowMatrix::Ifpack_Graph_Epetra_RowMatrix(const Teuchos::RefCountPtr<const Epetra_RowMatrix>& RowMatrix) :
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

// ======================================================================
ostream& Ifpack_Graph_Epetra_RowMatrix::Print(std::ostream& os) const
{

  if (Comm().MyPID())
    return(os);

  os << "================================================================================" << endl;
  os << "Ifpack_Graph_Epetra_RowMatrix" << endl;
  os << "Number of local rows  = " << NumMyRows_ << endl;
  os << "Number of global rows = " << NumGlobalRows_ << endl;
  os << "================================================================================" << endl;

  return(os);

}

