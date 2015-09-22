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
#include "Ifpack_Graph_Epetra_CrsGraph.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"

//==============================================================================
Ifpack_Graph_Epetra_CrsGraph::
Ifpack_Graph_Epetra_CrsGraph(const Teuchos::RefCountPtr<const Epetra_CrsGraph>& CrsGraph) :
CrsGraph_(CrsGraph)
{
  NumMyRows_ = CrsGraph_->NumMyRows();
  NumMyCols_ = CrsGraph_->NumMyCols();
  NumGlobalRows_ = CrsGraph_->NumGlobalRows64();
  NumGlobalCols_ = CrsGraph_->NumGlobalCols64();
  MaxNumIndices_ = CrsGraph_->MaxNumIndices();
}

//==============================================================================
const Epetra_Comm& Ifpack_Graph_Epetra_CrsGraph::Comm() const
{
  return(CrsGraph_->Comm());
}

//==============================================================================
bool Ifpack_Graph_Epetra_CrsGraph::Filled() const
{
  return(CrsGraph_->Filled());
}

//==============================================================================
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Ifpack_Graph_Epetra_CrsGraph::GRID(int LRID_in) const
{
  return(CrsGraph_->GRID(LRID_in));
}

//==============================================================================
int Ifpack_Graph_Epetra_CrsGraph::GCID(int LCID_in) const
{
  return(CrsGraph_->GCID(LCID_in));
}
#endif

long long Ifpack_Graph_Epetra_CrsGraph::GRID64(int LRID_in) const
{
  return(CrsGraph_->GRID64(LRID_in));
}

//==============================================================================
long long Ifpack_Graph_Epetra_CrsGraph::GCID64(int LCID_in) const
{
  return(CrsGraph_->GCID64(LCID_in));
}

//==============================================================================
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Ifpack_Graph_Epetra_CrsGraph::LRID(int GRID_in) const
{
  return(CrsGraph_->LRID(GRID_in));
}

//==============================================================================
int Ifpack_Graph_Epetra_CrsGraph::LCID(int GCID_in) const
{
  return(CrsGraph_->LCID(GCID_in));
}
#endif

//==============================================================================
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Ifpack_Graph_Epetra_CrsGraph::LRID(long long GRID_in) const
{
  return(CrsGraph_->LRID(GRID_in));
}

//==============================================================================
int Ifpack_Graph_Epetra_CrsGraph::LCID(long long GCID_in) const
{
  return(CrsGraph_->LCID(GCID_in));
}
#endif
//==============================================================================
int Ifpack_Graph_Epetra_CrsGraph::
ExtractMyRowCopy(int MyRow, int LenOfIndices,
                 int &NumIndices, int *Indices) const
{
  return(CrsGraph_->ExtractMyRowCopy(MyRow, LenOfIndices,
                                         NumIndices, Indices));
}

//==============================================================================
int Ifpack_Graph_Epetra_CrsGraph::NumMyNonzeros() const
{
  return(CrsGraph_->NumMyEntries());
}

// ======================================================================
std::ostream& Ifpack_Graph_Epetra_CrsGraph::Print(std::ostream& os) const
{
  using std::endl;

  if (Comm().MyPID())
    return(os);

  os << "================================================================================" << endl;
  os << "Ifpack_Graph_Epetra_CrsGraph" << endl;
  os << "Number of local rows  = " << NumMyRows_ << endl;
  os << "Number of global rows = " << NumGlobalRows_ << endl;
  os << "================================================================================" << endl;

  return(os);

}
