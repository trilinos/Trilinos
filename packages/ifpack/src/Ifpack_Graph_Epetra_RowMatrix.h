#ifndef IFPACK_GRAPH_EPETRA_ROWMATRIX_H
#define IFPACK_GRAPH_EPETRA_ROWMATRIX_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Graph.h"
#include <vector>
class Epetra_Comm;
class Epetra_RowMatrix;

class Ifpack_Graph_Epetra_RowMatrix : public Ifpack_Graph {

public:
    
  Ifpack_Graph_Epetra_RowMatrix(const Epetra_RowMatrix* RowMatrix);

  ~Ifpack_Graph_Epetra_RowMatrix();

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
  const Epetra_RowMatrix* RowMatrix_;
  mutable std::vector<double> Values_;
};

#endif
