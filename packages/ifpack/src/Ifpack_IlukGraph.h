#ifndef _IFPACK_ILUK_GRAPH_H_
#define _IFPACK_ILUK_GRAPH_H_
#include "Epetra_Object.h" 
#include "Epetra_CrsGraph.h"
//! Ifpack_IlukGraph: A class for constructing level filled graphs for use with ILU(k) class preconditioners.

/*! The Ifpack_IlukGraph class enable the construction matrix graphs using level-fill algorithms.  The only function
    required for construction is an ExtractRowView capability, i.e., the matrix that is passed in to the constructor
    must implement the Ifpack_CrsGraph interface defined in Ifpack_CrsMatrix.h 


<b>Constructing Ifpack_IlukGraph objects</b>

Constructing Ifpack_IlukGraph objects is usually a two step process of passing in a Ifpack_CrsGraph object and 
an integer indicating the desired level of fill and then calling the ConstructFilledGraph function to complete the
process.  This allows warning error codes to be returned to the calling routine.

It is worth noting that an Ifpack_IlukGraph object has two
Epetra_CrsGraph objects containing L and U, the graphs for the lower and upper triangular parts of the ILU(k) graph.
Thus, it is possible to manually insert and delete graph entries in L and U via the Epetra_CrsGraph
InsertIndices and RemoveIndices functions.  However, in this case FillComplete must be
called before the graph is used for subsequent operations.

*/    


class Ifpack_IlukGraph {
      
  // Give ostream << function some access to private and protected data/functions.

  friend ostream& operator << (ostream& os, const Ifpack_IlukGraph& A);

 public:

  //! Ifpack_IlukGraph constuctor.
  /*! Creates a Ifpack_IlukGraph object using the input graph and specified level of fill.  
    
    \param In
           Graph - An existing Ifpack_CrsGraph.  This object must implement the Ifpack_CrsGraph functions
	   that provide graph dimension and pattern information.
    \param In
           LevelFill - The level of fill to compute via ILU(k) algorithm.
    \param In
           LevelOverlap - The level of between subdomains.

	   \warning Actual construction occurs in ConstructFilledGraph.  This allows error codes to 
	            be passed back to the user.
  */
  Ifpack_IlukGraph(const Epetra_CrsGraph & Graph, int LevelFill, int LevelOverlap);
  
  //! Copy constructor.
  Ifpack_IlukGraph(const Ifpack_IlukGraph & Graph);

  //! Ifpack_IlukGraph Destructor
  virtual ~Ifpack_IlukGraph();
  

  //! Does the actual construction of the graph.
  /* 
    \return Integer error code, set to 0 if successful.
    
  */
  virtual int ConstructFilledGraph();

  //! Does the actual construction of the overlap matrix graph.
  /* 
    \return Integer error code, set to 0 if successful.
    
  */
  virtual int ConstructOverlapGraph();
  
  //! Returns the level of fill used to construct this graph.
  virtual int LevelFill() const {return(LevelFill_);};
  
  //! Returns the level of overlap used to construct this graph.
  virtual int LevelOverlap() const {return(LevelOverlap_);};
    
  //! Returns the number of global matrix rows.
  int NumGlobalRows() const {return(NumGlobalRows_);};
  
  //! Returns the number of global matrix columns.
  int NumGlobalCols() const {return(NumGlobalCols_);};
  
  //! Returns the number of nonzero entries in the global graph.
  int NumGlobalNonzeros() const {return(NumGlobalNonzeros_);};
  
  //! Returns the number of diagonal entries found in the global input graph.
  virtual int NumGlobalDiagonals() const {return(NumGlobalDiagonals_);};
  
  //! Returns the number of local matrix rows.
  int NumMyRows() const {return(NumMyRows_);};
  
  //! Returns the number of local matrix columns.
  int NumMyCols() const {return(NumMyCols_);};
  
  //! Returns the number of nonzero entries in the local graph.
  int NumMyNonzeros() const {return(NumMyNonzeros_);};
  
  //! Returns the number of diagonal entries found in the local input graph.
  virtual int NumMyDiagonals() const {return(NumMyDiagonals_);};
  
  //! Returns the index base for row and column indices for this graph.
  int IndexBase() const {return(IndexBase_);};
  
  //! Returns the graph of lower triangle of the ILU(k) graph as a Epetra_CrsGraph.
  virtual Epetra_CrsGraph & L_Graph() {return(*L_Graph_);};
  
  //! Returns the graph of upper triangle of the ILU(k) graph as a Epetra_CrsGraph.
  virtual Epetra_CrsGraph & U_Graph() {return(*U_Graph_);};
  
  //! Returns the graph of lower triangle of the ILU(k) graph as a Epetra_CrsGraph.
  virtual Epetra_CrsGraph & L_Graph() const {return(*L_Graph_);};
  
  //! Returns the graph of upper triangle of the ILU(k) graph as a Epetra_CrsGraph.
  virtual Epetra_CrsGraph & U_Graph() const {return(*U_Graph_);};
  
  //! Returns the importer used to create the overlapped graph.
  virtual Epetra_Import * OverlapImporter() const  {return(OverlapImporter_);};
  
  //! Returns the the overlapped graph.
  virtual Epetra_CrsGraph * OverlapGraph() const  {return(OverlapGraph_);};

    //! Returns the Epetra_BlockMap object associated with the domain of this matrix operator.
    virtual const Epetra_BlockMap & DomainMap() const {return(DomainMap_);};

    //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
    virtual const Epetra_BlockMap & RangeMap() const{return(RangeMap_);};

    //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
    virtual const Epetra_Comm & Comm() const{return(Comm_);};
  
 private:
    

  const Epetra_CrsGraph & Graph_;
  const Epetra_BlockMap & DomainMap_;
  const Epetra_BlockMap & RangeMap_;
  const Epetra_Comm & Comm_;
  Epetra_CrsGraph * OverlapGraph_;
  Epetra_Import * OverlapImporter_;
  Epetra_BlockMap * OverlapRowMap_;
  int LevelFill_;
  int LevelOverlap_;
  Epetra_CrsGraph * L_Graph_;
  Epetra_CrsGraph * U_Graph_;
  int IndexBase_;
  int NumGlobalRows_;
  int NumGlobalCols_;
  int NumGlobalDiagonals_;
  int NumGlobalNonzeros_;
  int NumMyRows_;
  int NumMyCols_;
  int NumMyDiagonals_;
  int NumMyNonzeros_;
  

 };

//! << operator will work for Ifpack_IlukGraph.
ostream& operator << (ostream& os, const Ifpack_IlukGraph& A);

#endif /* _IFPACK_ILUK_GRAPH_H_ */
