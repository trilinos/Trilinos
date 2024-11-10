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
