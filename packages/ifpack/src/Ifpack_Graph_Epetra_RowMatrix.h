#ifndef IFPACK_GRAPH_EPETRA_ROWMATRIX_H
#define IFPACK_GRAPH_EPETRA_ROWMATRIX_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Graph.h"
#include <vector>
class Epetra_Comm;
class Epetra_RowMatrix;

class Ifpack_Graph_Epetra_RowMatrix : public Ifpack_Graph {

public:
    
  //! Constructor.
  Ifpack_Graph_Epetra_RowMatrix(const Epetra_RowMatrix* RowMatrix);

  //! Destructor.
  ~Ifpack_Graph_Epetra_RowMatrix();

  //! Returns the number of local rows.  
  int NumMyRows() const
  {
    return(NumMyRows_);
  }

  //! Returns the number of local columns.
  int NumMyCols() const
  {
    return(NumMyCols_);
  }

  //! Returns the number of global rows.
  int NumGlobalRows() const
  {
    return(NumGlobalRows_);
  }

  //! Returns the number of global columns.
  int NumGlobalCols() const
  {
    return(NumGlobalCols_);
  }

  //! Returns the maximun number of entries for row.
  int MaxMyNumEntries() const
  {
    return(MaxNumIndices_);
  }

  //! Returns the number of local nonzero entries.
  int NumMyNonzeros() const;

  //! Returns \c true is graph is filled.
  bool Filled() const;

  //! Returns the global row ID of input local row.
  int GRID(int) const;

  //! Returns the global column ID of input local column.
  int GCID(int) const;
  
  //! Returns the local row ID of input global row.
  int LRID(int) const;

  //! Returns the local column ID of input global column.
  int LCID(int) const;

  //! Extracts a copy of input local row.
  int ExtractMyRowCopy(int GlobalRow, int LenOfIndices, 
		       int &NumIndices, int *Indices) const;

  //! Returns the communicator object of the graph.
  const Epetra_Comm& Comm() const;  
  
private:

  //! Number of local rows.
  int NumMyRows_;
  //! Number of local columns.
  int NumMyCols_;
  //! Number of global rows.
  int NumGlobalRows_;
  //! Number of global columns.
  int NumGlobalCols_;
  //! Maximum number of indices per row.
  int MaxNumIndices_;
  //! Pointer to the wrapped Epetra_CrsGraph.
  const Epetra_RowMatrix* RowMatrix_;
  //! Vectors that can be used in calls to ExtractMyRowView of the Row matrix.
  mutable std::vector<double> Values_;
};

#endif
