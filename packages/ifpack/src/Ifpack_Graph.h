#ifndef IFPACK_GRAPH_H
#define IFPACK_GRAPH_H
class Epetra_Comm;

//! Ifpack_Graph: a pure virtual class that defines graphs for IFPACK.
/*!
Class Ifpack_Graph defines the abstract interface to use graphs in
IFPACK. This class contains all the functions that are required by
IFPACK classes.

\date Sep-04.

*/

class Ifpack_Graph {

public:
    
  //! Destructor.
  virtual ~Ifpack_Graph() {};

  //! Returns the number of local rows.
  virtual int NumMyRows() const = 0;

  //! Returns the number of local columns.
  virtual int NumMyCols() const = 0;

  //! Returns the number of global rows.
  virtual int NumGlobalRows() const = 0;

  //! Returns the number of global columns.
  virtual int NumGlobalCols() const = 0;

  //! Returns the maximun number of entries for row.
  virtual int MaxMyNumEntries() const = 0;

  //! Returns the number of local nonzero entries.
  virtual int NumMyNonzeros() const = 0;

  //! Returns \c true is graph is filled.
  virtual bool Filled() const = 0;

  //! Returns the global row ID of input local row.
  virtual int GRID(int) const = 0;

  //! Returns the global column ID of input local column.
  virtual int GCID(int) const = 0;
  
  //! Returns the local row ID of input global row.
  virtual int LRID(int) const = 0;

  //! Returns the local column ID of input global column.
  virtual int LCID(int) const = 0;

  //! Extracts a copy of input local row.
  virtual int ExtractMyRowCopy(int MyRow, int LenOfIndices, 
			       int &NumIndices, int *Indices) const = 0;

  //! Returns the communicator object of the graph.
  virtual const Epetra_Comm& Comm() const = 0;
};

#endif // iFPACK_GRAPH_H
