/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
  NumGlobalRows_ = RowMatrix_->NumGlobalRows64();
  NumGlobalCols_ = RowMatrix_->NumGlobalCols64();
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
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Ifpack_Graph_Epetra_RowMatrix::GRID(int LRID_in) const
{
  return(RowMatrix_->RowMatrixRowMap().GID(LRID_in));
}

//==============================================================================
int Ifpack_Graph_Epetra_RowMatrix::GCID(int LCID_in) const
{
  return(RowMatrix_->RowMatrixColMap().GID(LCID_in));
}
#endif

long long Ifpack_Graph_Epetra_RowMatrix::GRID64(int LRID_in) const
{
  return(RowMatrix_->RowMatrixRowMap().GID64(LRID_in));
}

//==============================================================================
long long Ifpack_Graph_Epetra_RowMatrix::GCID64(int LCID_in) const
{
  return(RowMatrix_->RowMatrixColMap().GID64(LCID_in));
}

//==============================================================================
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Ifpack_Graph_Epetra_RowMatrix::LRID(int GRID_in) const
{
  return(RowMatrix_->RowMatrixRowMap().LID(GRID_in));
}

//==============================================================================
int Ifpack_Graph_Epetra_RowMatrix::LCID(int GCID_in) const
{
  return(RowMatrix_->RowMatrixColMap().LID(GCID_in));
}
#endif
//==============================================================================
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Ifpack_Graph_Epetra_RowMatrix::LRID(long long GRID_in) const
{
  return(RowMatrix_->RowMatrixRowMap().LID(GRID_in));
}

//==============================================================================
int Ifpack_Graph_Epetra_RowMatrix::LCID(long long GCID_in) const
{
  return(RowMatrix_->RowMatrixColMap().LID(GCID_in));
}
#endif
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
std::ostream& Ifpack_Graph_Epetra_RowMatrix::Print(std::ostream& os) const
{
  using std::endl;

  if (Comm().MyPID())
    return(os);

  os << "================================================================================" << endl;
  os << "Ifpack_Graph_Epetra_RowMatrix" << endl;
  os << "Number of local rows  = " << NumMyRows_ << endl;
  os << "Number of global rows = " << NumGlobalRows_ << endl;
  os << "================================================================================" << endl;

  return(os);

}

