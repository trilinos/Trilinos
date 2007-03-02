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
  NumGlobalRows_ = CrsGraph_->NumGlobalRows();
  NumGlobalCols_ = CrsGraph_->NumGlobalCols();
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
int Ifpack_Graph_Epetra_CrsGraph::GRID(int LRID) const
{
  return(CrsGraph_->GRID(LRID));
}

//==============================================================================
int Ifpack_Graph_Epetra_CrsGraph::GCID(int LCID) const
{
  return(CrsGraph_->GCID(LCID));
}

//==============================================================================
int Ifpack_Graph_Epetra_CrsGraph::LRID(int GRID) const
{
  return(CrsGraph_->LRID(GRID));
}

//==============================================================================
int Ifpack_Graph_Epetra_CrsGraph::LCID(int GCID) const
{
  return(CrsGraph_->LCID(GCID));
}

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
ostream& Ifpack_Graph_Epetra_CrsGraph::Print(std::ostream& os) const
{

  if (Comm().MyPID())
    return(os);

  os << "================================================================================" << endl;
  os << "Ifpack_Graph_Epetra_CrsGraph" << endl;
  os << "Number of local rows  = " << NumMyRows_ << endl;
  os << "Number of global rows = " << NumGlobalRows_ << endl;
  os << "================================================================================" << endl;

  return(os);

}
