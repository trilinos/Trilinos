#ifndef IFPACK_EPETRA_CRSGRAPH_H
#define IFPACK_EPETRA_CRSGRAPH_H
#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Graph.h"
class Epetra_Comm;
class Epetra_CrsGraph;

//! Ifpack_Graph_Epetra_CrsGraph: a class to define Ifpack_Graph as a light-weight conversion of Epetra_CrsGraph.

/*! 
Class Ifpack_Graph_Epetra_CrsGraph enables the construction of an
Ifpack_Graph based on the input Epetra_CrsGraph. Note that data are
not copied to \e this object; instead, wrappers are furnished.

\date Set-04.
*/
class Ifpack_Graph_Epetra_CrsGraph : public Ifpack_Graph {

public:
    
  //! Constructor.
  Ifpack_Graph_Epetra_CrsGraph(const Epetra_CrsGraph* CrsGraph);

  //! Destructor.
  ~Ifpack_Graph_Epetra_CrsGraph();

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
  const Epetra_CrsGraph* CrsGraph_;
};

#endif
