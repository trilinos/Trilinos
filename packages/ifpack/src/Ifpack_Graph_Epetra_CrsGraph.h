#ifndef IFPACK_EPETRA_CRSGRAPH_H
#define IFPACK_EPETRA_CRSGRAPH_H
#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Graph.h"
class Epetra_Comm;

class Epetra_CrsGraph;

class Ifpack_Graph_Epetra_CrsGraph : public Ifpack_Graph {

public:
    
  Ifpack_Graph_Epetra_CrsGraph(const Epetra_CrsGraph* CrsGraph);

  ~Ifpack_Graph_Epetra_CrsGraph();

  int NumMyRows() const
  {
    return(NumMyRows_);
  }

  int NumMyCols() const
  {
    return(NumMyCols_);
  }

  int NumGlobalRows() const
  {
    return(NumGlobalRows_);
  }

  int NumGlobalCols() const
  {
    return(NumGlobalCols_);
  }

  int MaxMyNumEntries() const
  {
    return(MaxNumIndices_);
  }

  int NumMyNonzeros() const;

  bool Filled() const;

  int GRID(int) const;

  int GCID(int) const;
  
  int LRID(int) const;

  int LCID(int) const;

  int ExtractMyRowCopy(int GlobalRow, int LenOfIndices, 
		       int &NumIndices, int *Indices) const;

  const Epetra_Comm& Comm() const;
  
  
private:

  int NumMyRows_;
  int NumMyCols_;
  int NumGlobalRows_;
  int NumGlobalCols_;
  int MaxNumIndices_;
  const Epetra_CrsGraph* CrsGraph_;
};

#endif
