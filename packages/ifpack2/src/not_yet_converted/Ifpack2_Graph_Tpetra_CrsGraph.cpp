/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Graph.hpp"
#include "Ifpack2_Graph_Tpetra_CrsGraph.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_CrsGraph.hpp"

//==============================================================================
Ifpack2_Graph_Tpetra_CrsGraph::
Ifpack2_Graph_Tpetra_CrsGraph(const Teuchos::RCP<const Tpetra_CrsGraph>& CrsGraph) :
CrsGraph_(CrsGraph)
{
  NumMyRows_ = CrsGraph_->NumMyRows();
  NumMyCols_ = CrsGraph_->NumMyCols();
  NumGlobalRows_ = CrsGraph_->NumGlobalRows();
  NumGlobalCols_ = CrsGraph_->NumGlobalCols();
  MaxNumIndices_ = CrsGraph_->MaxNumIndices();
}

//==============================================================================
const Tpetra_Comm& Ifpack2_Graph_Tpetra_CrsGraph::Comm() const
{
  return(CrsGraph_->Comm());
}

//==============================================================================
bool Ifpack2_Graph_Tpetra_CrsGraph::Filled() const
{
  return(CrsGraph_->Filled());
}
 
//==============================================================================
int Ifpack2_Graph_Tpetra_CrsGraph::GRID(int LRID_in) const
{
  return(CrsGraph_->GRID(LRID_in));
}

//==============================================================================
int Ifpack2_Graph_Tpetra_CrsGraph::GCID(int LCID_in) const
{
  return(CrsGraph_->GCID(LCID_in));
}

//==============================================================================
int Ifpack2_Graph_Tpetra_CrsGraph::LRID(int GRID_in) const
{
  return(CrsGraph_->LRID(GRID_in));
}

//==============================================================================
int Ifpack2_Graph_Tpetra_CrsGraph::LCID(int GCID_in) const
{
  return(CrsGraph_->LCID(GCID_in));
}

//==============================================================================
int Ifpack2_Graph_Tpetra_CrsGraph::
ExtractMyRowCopy(int MyRow, int LenOfIndices, 
		 int &NumIndices, int *Indices) const
{
  return(CrsGraph_->ExtractMyRowCopy(MyRow, LenOfIndices,
					 NumIndices, Indices));
}

//==============================================================================
int Ifpack2_Graph_Tpetra_CrsGraph::NumMyNonzeros() const
{
  return(CrsGraph_->NumMyEntries());
}

// ======================================================================
ostream& Ifpack2_Graph_Tpetra_CrsGraph::Print(std::ostream& os) const
{

  if (Comm().MyPID())
    return(os);

  os << "================================================================================" << endl;
  os << "Ifpack2_Graph_Tpetra_CrsGraph" << endl;
  os << "Number of local rows  = " << NumMyRows_ << endl;
  os << "Number of global rows = " << NumGlobalRows_ << endl;
  os << "================================================================================" << endl;

  return(os);

}
