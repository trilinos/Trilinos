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
#include "Tifpack_Graph_Tpetra_CrsGraph.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_CrsGraph.hpp"

//==============================================================================
Tifpack_Graph_Tpetra_CrsGraph::
Tifpack_Graph_Tpetra_CrsGraph(const Teuchos::RCP<const Tpetra_CrsGraph>& CrsGraph) :
CrsGraph_(CrsGraph)
{
  NumMyRows_ = CrsGraph_->NumMyRows();
  NumMyCols_ = CrsGraph_->NumMyCols();
  NumGlobalRows_ = CrsGraph_->NumGlobalRows();
  NumGlobalCols_ = CrsGraph_->NumGlobalCols();
  MaxNumIndices_ = CrsGraph_->MaxNumIndices();
}

//==============================================================================
const Tpetra_Comm& Tifpack_Graph_Tpetra_CrsGraph::Comm() const
{
  return(CrsGraph_->Comm());
}

//==============================================================================
bool Tifpack_Graph_Tpetra_CrsGraph::Filled() const
{
  return(CrsGraph_->Filled());
}
 
//==============================================================================
int Tifpack_Graph_Tpetra_CrsGraph::GRID(int LRID_in) const
{
  return(CrsGraph_->GRID(LRID_in));
}

//==============================================================================
int Tifpack_Graph_Tpetra_CrsGraph::GCID(int LCID_in) const
{
  return(CrsGraph_->GCID(LCID_in));
}

//==============================================================================
int Tifpack_Graph_Tpetra_CrsGraph::LRID(int GRID_in) const
{
  return(CrsGraph_->LRID(GRID_in));
}

//==============================================================================
int Tifpack_Graph_Tpetra_CrsGraph::LCID(int GCID_in) const
{
  return(CrsGraph_->LCID(GCID_in));
}

//==============================================================================
int Tifpack_Graph_Tpetra_CrsGraph::
ExtractMyRowCopy(int MyRow, int LenOfIndices, 
		 int &NumIndices, int *Indices) const
{
  return(CrsGraph_->ExtractMyRowCopy(MyRow, LenOfIndices,
					 NumIndices, Indices));
}

//==============================================================================
int Tifpack_Graph_Tpetra_CrsGraph::NumMyNonzeros() const
{
  return(CrsGraph_->NumMyEntries());
}

// ======================================================================
ostream& Tifpack_Graph_Tpetra_CrsGraph::Print(std::ostream& os) const
{

  if (Comm().MyPID())
    return(os);

  os << "================================================================================" << endl;
  os << "Tifpack_Graph_Tpetra_CrsGraph" << endl;
  os << "Number of local rows  = " << NumMyRows_ << endl;
  os << "Number of global rows = " << NumGlobalRows_ << endl;
  os << "================================================================================" << endl;

  return(os);

}
