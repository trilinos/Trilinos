#ifndef IFPACK_GRAPH_H
#define IFPACK_GRAPH_H
class Epetra_Comm;

class Ifpack_Graph {

public:
    
  virtual ~Ifpack_Graph() {};

  virtual int NumMyRows() const = 0;

  virtual int NumMyCols() const = 0;

  virtual int NumGlobalRows() const = 0;

  virtual int NumGlobalCols() const = 0;

  virtual int MaxMyNumEntries() const = 0;

  virtual int NumMyNonzeros() const = 0;

  virtual bool Filled() const = 0;

  virtual int GRID(int) const = 0;

  virtual int GCID(int) const = 0;
  
  virtual int LRID(int) const = 0;

  virtual int LCID(int) const = 0;

  virtual int ExtractMyRowCopy(int GlobalRow, int LenOfIndices, 
			       int &NumIndices, int *Indices) const = 0;

  virtual const Epetra_Comm& Comm() const = 0;
};

#endif // iFPACK_GRAPH_H
