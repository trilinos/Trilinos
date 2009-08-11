/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_Graph.hpp"
#include "Tifpack_Graph_Tpetra_RowMatrix.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"

//==============================================================================
Tifpack_Graph_Tpetra_RowMatrix::Tifpack_Graph_Tpetra_RowMatrix(const Teuchos::RefCountPtr<const Tpetra_RowMatrix>& RowMatrix) :
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
const Tpetra_Comm& Tifpack_Graph_Tpetra_RowMatrix::Comm() const
{
  return(RowMatrix_->Comm());
}

//==============================================================================
bool Tifpack_Graph_Tpetra_RowMatrix::Filled() const
{
  return(RowMatrix_->Filled());
}
 
//==============================================================================
int Tifpack_Graph_Tpetra_RowMatrix::GRID(int LRID_in) const
{
  return(RowMatrix_->RowMatrixRowMap().GID(LRID_in));
}

//==============================================================================
int Tifpack_Graph_Tpetra_RowMatrix::GCID(int LCID_in) const
{
  return(RowMatrix_->RowMatrixColMap().GID(LCID_in));
}

//==============================================================================
int Tifpack_Graph_Tpetra_RowMatrix::LRID(int GRID_in) const
{
  return(RowMatrix_->RowMatrixRowMap().LID(GRID_in));
}

//==============================================================================
int Tifpack_Graph_Tpetra_RowMatrix::LCID(int GCID_in) const
{
  return(RowMatrix_->RowMatrixColMap().LID(GCID_in));
}

//==============================================================================
int Tifpack_Graph_Tpetra_RowMatrix::
ExtractMyRowCopy(int MyRow, int LenOfIndices, 
		 int &NumIndices, int *Indices) const
{
  return(RowMatrix_->ExtractMyRowCopy(MyRow, LenOfIndices,
				      NumIndices, &Values_[0],
				      Indices));
}

//==============================================================================
int Tifpack_Graph_Tpetra_RowMatrix::NumMyNonzeros() const
{
  return(RowMatrix_->NumMyNonzeros());
}

// ======================================================================
ostream& Tifpack_Graph_Tpetra_RowMatrix::Print(std::ostream& os) const
{

  if (Comm().MyPID())
    return(os);

  os << "================================================================================" << endl;
  os << "Tifpack_Graph_Tpetra_RowMatrix" << endl;
  os << "Number of local rows  = " << NumMyRows_ << endl;
  os << "Number of global rows = " << NumGlobalRows_ << endl;
  os << "================================================================================" << endl;

  return(os);

}

