#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Graph.h"
#include "Ifpack_Graph_Epetra_CrsGraph.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"

//==============================================================================
Ifpack_Graph_Epetra_CrsGraph::
Ifpack_Graph_Epetra_CrsGraph(const Epetra_CrsGraph* CrsGraph) :
CrsGraph_(CrsGraph)
{
  NumMyRows_ = CrsGraph_->NumMyRows();
  NumMyCols_ = CrsGraph_->NumMyCols();
  NumGlobalRows_ = CrsGraph_->NumGlobalRows();
  NumGlobalCols_ = CrsGraph_->NumGlobalCols();
  MaxNumIndices_ = CrsGraph_->MaxNumIndices();
}

//==============================================================================
Ifpack_Graph_Epetra_CrsGraph::~Ifpack_Graph_Epetra_CrsGraph()
{
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

