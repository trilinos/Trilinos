#ifndef _PETRA_CRS_GRAPH_H_
#define _PETRA_CRS_GRAPH_H_
//! Petra_CRS_Graph: A class for constructing and using sparse compressed row graphs.

/*! The Petra_CRS_Graph enable the piecewise construction and use of sparse matrix graphs (the integer structure without
    values) where entries are intended for row access.  

    Petra_CRS_Graph is a base class for all row-based matrix classes.

<b>Constructing Petra_CRS_Graph objects</b>

Constructing Petra_CRS_Graph objects is a multi-step process.  The basic steps are as follows:
<ol>
  <li> Create Petra_CRS_Graph instance, including storage,  via constructor.
  <li> Enter row and column entry information via calls to the InsertIndices function.
  <li> Complete construction via TransformToLocal call.
  <li> (Optional) Optimize the graph storage via a call to StorageOptimize.
</ol>

Note that, even after a matrix is constructed, it is possible to add or remove entries.  However, 
TransformToLocal must be
called again before the graph is used for subsequent operations.



*/    

#include "Petra_Petra.h" 
#include "Petra_Util.h"
#include "Petra_Time.h" 
#include "Petra_Map.h" 
#include "Petra_BlockMap.h" 
#include "Petra_Import.h" 
#include "Petra_Export.h" 


class Petra_CRS_Graph {
      
  // Give ostream << function some access to private and protected data/functions.
  
  friend ostream& operator << (ostream& os, const Petra_CRS_Graph& A);
  
 public:

  //! Petra_CRS_Graph constuctor with variable number of indices per row.
  /*! Creates a Petra_CRS_Graph object and allocates storage.  
    
    \param In
           CV - A Petra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - A Petra_BlockMap (or Petra_Map or Petra_LocalMap).
    \param In
           NumIndicesPerRow - An integer array of length NumMyRows
	   such that NumIndicesPerRow[i] indicates the (approximate) number of entries in the ith row.
  */
  Petra_CRS_Graph(Petra_DataAccess CV, const Petra_BlockMap& RowMap, int *NumIndicesPerRow);
  
  //! Petra_CRS_Graph constuctor with fixed number of indices per row.
  /*! Creates a Petra_CRS_Graph object and allocates storage.  
    
    \param In
           CV - A Petra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - A Petra_BlockMap (or Petra_Map or Petra_LocalMap).
    \param In
           NumIndicesPerRow - An integer that indicates the (approximate) number of entries in the each row.
	   Note that it is possible to use 0 for this value and let fill occur during the insertion phase.
	   
  */
  Petra_CRS_Graph(Petra_DataAccess CV, const Petra_BlockMap& RowMap, int NumIndicesPerRow);
  
  //! Petra_CRS_Graph constuctor with variable number of indices per row.
  /*! Creates a Petra_CRS_Graph object and allocates storage.  
    
    \param In
           CV - A Petra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - A Petra_BlockMap (or Petra_Map or Petra_LocalMap).
    \param In 
           ColMap - A Petra_BlockMap (or Petra_Map or Petra_LocalMap).
    \param In
           NumIndicesPerRow - An integer array of length NumMyRows
	   such that NumIndicesPerRow[i] indicates the (approximate) number of entries in the ith row.
  */
  Petra_CRS_Graph(Petra_DataAccess CV, const Petra_BlockMap& RowMap, const Petra_BlockMap& ColMap, int *NumIndicesPerRow);
  
  //! Petra_CRS_Graph constuctor with fixed number of indices per row.
  /*! Creates a Petra_CRS_Graph object and allocates storage.  
    
    \param In
           CV - A Petra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - A Petra_BlockMap (or Petra_Map or Petra_LocalMap).
    \param In 
           ColMap - A Petra_BlockMap (or Petra_Map or Petra_LocalMap).
    \param In
           NumIndicesPerRow - An integer that indicates the (approximate) number of entries in the each row.
	   Note that it is possible to use 0 for this value and let fill occur during the insertion phase.
	   
  */
  Petra_CRS_Graph(Petra_DataAccess CV, const Petra_BlockMap& RowMap, const Petra_BlockMap& ColMap, int NumIndicesPerRow);
  
  //! Copy constructor.
  Petra_CRS_Graph(const Petra_CRS_Graph & Graph);

  //! Petra_CRS_Graph Destructor
  virtual ~Petra_CRS_Graph();
  
  //! Enter a list of elements in a specified global row of the matrix.
  /*!
    \param In
           Row - Global row number of indices.
    \param In
           NumIndices - Number of Indices.
    \param In
           Indices - Global column indices to insert.

    \return Integer error code, set to 0 if successful.
  */
  int InsertGlobalIndices(int GlobalRow, int NumIndices, int *Indices);
  
  //! Remove a list of elements from a specified global row of the matrix.
  /*!
    \param In
           Row - Global row number of indices.
    \param In
           NumIndices - Number of Indices.
    \param In
           Indices - Global column indices to remove.
	   
    \return Integer error code, set to 0 if successful.
  */
  int RemoveGlobalIndices(int GlobalRow, int NumIndices, int *Indices);
  
  //! Remove all indices from a specified global row of the matrix.
  /*!
    \param In
           Row - Global row number of indices.

    \return Integer error code, set to 0 if successful.
  */
  int RemoveGlobalIndices(int Row);

  
  //! Enter a list of elements in a specified local row of the matrix.
  /*!
    \param In
           Row - Local row number of indices.
    \param In
           NumIndices - Number of Indices.
    \param In
           Indices - Local column indices to insert.

    \return Integer error code, set to 0 if successful.
  */
  int InsertMyIndices(int LocalRow, int NumIndices, int *Indices);
  
  //! Remove a list of elements from a specified local row of the matrix.
  /*!
    \param In
           Row - Local row number of indices.
    \param In
           NumIndices - Number of Indices.
    \param In
           Indices - Local column indices to remove.
	   
    \return Integer error code, set to 0 if successful.
  */
  int RemoveMyIndices(int LocalRow, int NumIndices, int *Indices);
  
  //! Remove all indices from a specified local row of the matrix.
  /*!
    \param In
           Row - Local row number of indices.

    \return Integer error code, set to 0 if successful.
  */
  int RemoveMyIndices(int Row);

  //! Tranform matrix representation to local index space.  Perform other operations to allow optimal matrix operations.
  int TransformToLocal();

  //! Tranform to local index space using specified Domain/Range maps.  Perform other operations to allow optimal matrix operations.
  int TransformToLocal(Petra_BlockMap *DomainMap, Petra_BlockMap *RangeMap);
  
  //! Sort column indices, row-by-row, in ascending order.
  int SortIndices();

  //! Extract a list of elements in a specified global row of the matrix. Put into storage allocated by calling routine.
  /*!
    \param In
           Row - Global row number to get indices.
    \param In
           LenOfIndices - Length of Indices array.
    \param Out
           NumIndices - Number of Indices.
    \param Out
           Indices - Global column indices corresponding to values.

    \return Integer error code, set to 0 if successful.
  */
  int ExtractGlobalRowCopy(int GlobalRow, int LenOfIndices, int & NumIndices, int * Indices) const;

  //! Extract a list of elements in a specified local row of the matrix. Put into storage allocated by calling routine.
  /*!
    \param In
           Row - Local row number to get indices.
    \param In
           LenOfIndices - Length of Indices array.
    \param Out
           NumIndices - Number of Indices.
    \param Out
           Indices - Local column indices corresponding to values.

    \return Integer error code, set to 0 if successful.
  */
  int ExtractMyRowCopy(int LocalRow, int LenOfIndices, int & NumIndices, int * Indices) const;

  //! Get a view of the elements in a specified global row of the matrix.
  /*!
    This function requires that the graph not be completed (TransformToLocal() was \e not called).
    \param In
           Row - Local row number to get indices.
    \param Out
           NumIndices - Number of Indices.
    \param Out
           Indices - Column indices corresponding to values.

    \return Integer error code, set to 0 if successful. Returns -1 if invalid row.  Returns -2 if graph is completed.
  */
    int ExtractGlobalRowView(int GlobalRow, int & NumIndices, int *& Indices) const;

  //! Get a view of the elements in a specified local row of the matrix.
  /*!
    This function requires that the graph be completed TransformToLocal() was called).
    \param In
           Row - Local row number to get indices.
    \param Out
           NumIndices - Number of Indices.
    \param Out
           Indices - Column indices corresponding to values.

    \return Integer error code, set to 0 if successful. Returns -1 if invalid row.  Returns -2 if graph is not completed.
  */
    int ExtractMyRowView(int LocalRow, int & NumIndices, int *& Indices) const;

    //! If TransformToLocal() has been called, this query returns true, otherwise it returns false.
    bool Filled() const {return(Filled_);};

    //! If SortIndices() has been called, this query returns true, otherwise it returns false.
    bool Sorted() const {return(Sorted_);};

    //! Eliminates memory that is used for construction.  Make consecutive row index sections contiguous.
    int OptimizeStorage();

    //! If OptimizeStorage() has been called, this query returns true, otherwise it returns false.
    bool StorageOptimized() const {return(StorageOptimized_);};

    //! If column indices are in global range, this query returns true, otherwise it returns false.
    bool IndicesAreGlobal() const {return(IndicesAreGlobal_);};

    //! If column indices are in local range, this query returns true, otherwise it returns false.
    bool IndicesAreLocal() const {return(IndicesAreLocal_);};

    //! If graph is lower triangular, this query returns true, otherwise it returns false.
    bool LowerTriangular() const {return(LowerTriangular_);};

    //! If graph is upper triangular, this query returns true, otherwise it returns false.
    bool UpperTriangular() const {return(UpperTriangular_);};

    //! If graph is lower triangular, this query returns true, otherwise it returns false.
    bool NoDiagonal() const {return(NoDiagonal_);};

    //! Returns true of GID is owned by the calling processor, otherwise it returns false.
    bool MyGlobalRow(int GID) const {return(RowMap_.MyGID(GID));};

    // Atribute access functions
    
    //! Returns the number of matrix rows on this processor.
    int NumMyRows() const {return(NumMyRows_);};

    //! Returns the number of matrix rows in global matrix.
    int NumGlobalRows() const {return(NumGlobalRows_);};

    //! Returns the number of matrix columns on this processor.
    int NumMyCols() const {return(NumMyCols_);};

    //! Returns the number of matrix columns in global matrix.
    int NumGlobalCols() const {return(NumGlobalCols_);};

    //! Returns the number of indices in the global graph.
    int NumGlobalNonzeros() const {return(NumGlobalNonzeros_);};

    //! Returns the number of diagonal entries in the global graph.
    int NumGlobalDiagonals() const {return(NumGlobalDiagonals_);};

    //! Returns the number of diagonal entries in the local graph.
    int NumMyDiagonals() const {return(NumMyDiagonals_);};
    
    //! Returns the number of block matrix rows on this processor.
    int NumMyBlockRows() const {return(NumMyBlockRows_);};

    //! Returns the number of Block matrix rows in global matrix.
    int NumGlobalBlockRows() const {return(NumGlobalBlockRows_);};

    //! Returns the number of Block matrix columns on this processor.
    int NumMyBlockCols() const {return(NumMyBlockCols_);};

    //! Returns the number of Block matrix columns in global matrix.
    int NumGlobalBlockCols() const {return(NumGlobalBlockCols_);};

    //! Returns the number of Block diagonal entries in the local graph.
    int NumMyBlockDiagonals() const {return(NumMyBlockDiagonals_);};

    //! Returns the number of Block diagonal entries in the global graph.
    int NumGlobalBlockDiagonals() const {return(NumGlobalBlockDiagonals_);};

    //! Returns the number of entries in the global graph.
    int NumGlobalEntries() const {return(NumGlobalEntries_);};

    //! Returns the number of entries on this processor.
    int NumMyEntries() const {return(NumMyEntries_);};
    //! Returns the max row dimension of block entries on the processor.
    int MaxRowDim() const {return(MaxRowDim_);};

    //! Returns the max row dimension of block entries across all processors.
    int GlobalMaxRowDim() const {return(GlobalMaxRowDim_);};

    //! Returns the max column dimension of block entries on the processor.
    int MaxColDim() const {return(MaxColDim_);};

    //! Returns the max column dimension of block entries across all processors.
    int GlobalMaxColDim() const {return(GlobalMaxColDim_);};

    //! Returns the local row index for given global row index, returns -1 if no local row for this global row.
    int LRID( int GRID) const {return(RowMap_.LID(GRID));};

    //! Returns the global row index for give local row index, returns IndexBase-1 if we don't have this local row.
    int GRID( int LRID) const {return(RowMap_.GID(LRID));};

    //! Returns the local column index for given global column index, returns -1 if no local column for this global column.
    int LCID( int GCID) const;

    //! Returns the global column index for give local column index, returns IndexBase-1 if we don't have this local column.
    int GCID( int LCID) const;
 
    //! Returns true if the GRID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyGRID(int GRID) const {return(LRID(GRID)!=-1);};
   
    //! Returns true if the LRID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyLRID(int LRID) const {return(GRID(LRID)!=IndexBase_-1);};

    //! Returns true if the GCID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyGCID(int GCID) const {return(LCID(GCID)!=-1);};
   
    //! Returns true if the LRID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyLCID(int LCID) const {return(GCID(LCID)!=IndexBase_-1);};

    //! Returns the number of indices in the local graph.
    int NumMyNonzeros() const {return(NumMyNonzeros_);};

    //! Returns the current number of nonzero entries in specified global row on this processor.
    int NumGlobalIndices(int Row) const;

    //! Returns the allocated number of nonzero entries in specified global row on this processor.
    int NumAllocatedGlobalIndices(int Row) const;

    //! Returns the maximum number of nonzero entries across all rows on this processor.
    int MaxNumIndices() const {return(MaxNumIndices_);};

    //! Returns the maximun number of nonzero entries across all rows across all processors.
    int GlobalMaxNumIndices() const {return(GlobalMaxNumIndices_);};

    //! Returns the maximum number of nonzero entries across all rows on this processor.
    int MaxNumNonzeros() const {return(MaxNumNonzeros_);};

    //! Returns the maximun number of nonzero entries across all rows across all processors.
    int GlobalMaxNumNonzeros() const {return(GlobalMaxNumNonzeros_);};

    //! Returns the current number of nonzero entries in specified local row on this processor.
    int NumMyIndices(int Row) const {return(NumIndicesPerRow_[Row]);};

    //! Returns the allocated number of nonzero entries in specified local row on this processor.
    int NumAllocatedMyIndices(int Row) const {return(NumAllocatedIndicesPerRow_[Row]);};

    //! Returns the index base for row and column indices for this graph.
    int IndexBase() const {return(IndexBase_);};
    
    //! Returns the RowMap associated with this matrix.
    const Petra_BlockMap& RowMap() const {return(RowMap_);};
    
    //! Returns the ColMap associated with this matrix.
    const Petra_BlockMap& ColMap() const {return(ColMap_);};
    
    //! Returns the DomainMap associated with this matrix.
    const Petra_BlockMap& DomainMap() const {return(*DomainMap_);};
    
    //! Returns the RangeMap associated with this matrix.
    const Petra_BlockMap& RangeMap() const {return(*RangeMap_);};
    
    //! Returns the Petra_Comm object associated with this matrix.
    const Petra_Comm& Comm() const {return(RowMap_.Comm());};

    //! Returns the ImportMap associated with this matrix.
    /*! ImportMap is the a map of global indices that are needed for column access by this graph.
        This information will be used by the matrix classes for distributed matrix operations.
     */
    const Petra_BlockMap& ImportMap() const 
      {if (ImportMap_==0) return (RowMap_); else return(*ImportMap_);};

    //! Returns the Importer associated with this matrix.
    const Petra_Import * Importer() const {return(Importer_);};

    //! Returns the ExportMap associated with this matrix.
    /*! ExportMap is the a map of global indices that are needed for updating elements in the
        range of this graph that are computed on this processor but not owned by this processor.
        This information will be used by the matrix classes for distributed matrix operations.
     */
    const Petra_BlockMap& ExportMap() const {return(*ExportMap_);};

    //! Returns the Exporter associated with this matrix.
    const Petra_Export * Exporter() const {return(Exporter_);};

    //! Removes any redundant column indices in the rows of the graph.
    int RemoveRedundantIndices();

    //! Fills a graph with rows from a source graph based on the specified importer.
  /*!
    \param In
           SourceGraph - Graph from which values are imported into the "\e this" graph.
    \param In
           Importer - A Petra_Import object specifying the communication required.

    \param In
           CombineMode - A Petra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
    int Import(const Petra_CRS_Graph& SourceGraph, const Petra_Import & Importer, 
	       Petra_CombineMode CombineMode);

    //! Fills a graph with rows from a source graph based on the specified exporter.
  /*!
    \param In
           SourceGraph - Graph from which values are imported into the "\e this" graph.
    \param In
           Exporter - A Petra_Export object specifying the communication required.

    \param In
           CombineMode - A Petra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
    int Export(const Petra_CRS_Graph& SourceGraph, 
	       const Petra_Export & Exporter, Petra_CombineMode CombineMode);

    //! Fills a graph with rows from a source graph based on the specified exporter.
  /*!
    \param In
           SourceGraph - Graph from which values are imported into the "\e this" graph.
    \param In
           Exporter - A Petra_Export object specifying the communication required.  Communication is done
	   in reverse of an export.

    \param In
           CombineMode - A Petra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
    int Import(const Petra_CRS_Graph& SourceGraph, 
	       const Petra_Export & Exporter, Petra_CombineMode CombineMode);

    //! Fills a graph with rows from a source graph based on the specified importer.
  /*!
    \param In
           SourceGraph - Graph from which values are imported into the "\e this" graph.
    \param In
           Importer - A Petra_Import object specifying the communication required.Communication is done
	   in reverse of an import.

    \param In
           CombineMode - A Petra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
    int Export(const Petra_CRS_Graph& SourceGraph, 
	       const Petra_Import & Importer, Petra_CombineMode CombineMode);

 protected:

    // If column indices are stored in one long array (via a call to OptimizeStorage), this query returns true, 
    // otherwise it returns false.
    bool IndicesAreContiguous() const {return(IndicesAreContiguous_);};
    int ** Indices() const {return(Indices_);};
    int * NumIndicesPerRow() const {return(NumIndicesPerRow_);};
    int * NumAllocatedIndicesPerRow() const {return(NumAllocatedIndicesPerRow_);};
    void SetSorted(bool Flag) {Sorted_ = Flag;};
    void SetIndicesAreGlobal(bool Flag) {IndicesAreGlobal_ = Flag;};
    void SetIndicesAreLocal(bool Flag) {IndicesAreLocal_ = Flag;};
    void SetIndicesAreContiguous(bool Flag) {IndicesAreContiguous_ = Flag;};
    void SetNoRedundancies(bool Flag) {NoRedundancies_ = Flag;};
    int InsertIndices(int Row, int NumIndices, int *Indices);
    bool FindGlobalIndexLoc(int LocalRow, int Index, int Start, int & Loc);
    bool FindMyIndexLoc(int LocalRow, int Index, int Start, int &
 Loc);
    int MakeIndicesLocal(const Petra_BlockMap & DomainMap, const Petra_BlockMap & RangeMap);

    friend class Petra_RDP_CRS_Matrix;
    friend class Petra_RDP_VBR_Matrix;

 private:
    void InitializeDefaults();
    int Allocate(int * NumIndicesPerRow, int Inc );
    int ComputeGlobalConstants();
    int ReAllocate();
    void SetFilled(bool Flag) {Filled_ = Flag;};
    bool Allocated() const {return(Allocated_);};
    bool NoRedundancies() const {return(NoRedundancies_);};
    void SetAllocated(bool Flag) {Allocated_ = Flag;};
  int DoTransfer(const Petra_CRS_Graph& SourceGraph, 
				       Petra_CombineMode CombineMode,
				       int NumSameIDs, int NumPermuteIDs, int NumRemoteIDs, 
				       int NumExportIDs, 
				       int *PermuteToLIDs, int *PermuteFromLIDs, int *RemoteLIDs, 
				       int * ExportLIDs,
				       int Nsend, int Nrecv, int SizeOfPacket,
				       int & LenExports, int * & Exports,
				       int & LenImports, int * & Imports,
#ifdef PETRA_MPI
				     GSComm_Plan & Plan, 
#endif
				     bool DoReverse);

  int CopyAndPermute(Petra_CRS_Graph & Target, 
					   const Petra_CRS_Graph & Source,
					   int NumSameIDs, 
					   int NumPermuteIDs, int * PermuteToLIDs,
					   int *PermuteFromLIDs);

  int Pack(const Petra_CRS_Graph & Source,
				 int NumSendIDs, int * SendLIDs, int * Sends);

  int UnpackAndCombine(Petra_CRS_Graph & Target, int SizeOfPacket,
					     int NumRecvIDs, int * RecvLIDs, 
					     int * Recvs, Petra_CombineMode CombineMode);

  // Non-trivially initialized by all constructors
  const Petra_BlockMap& RowMap_;
  const Petra_BlockMap& ColMap_;

  // Defined by TransformToLocal and related
  Petra_BlockMap * DomainMap_;
  Petra_BlockMap * RangeMap_;
  Petra_BlockMap * ImportMap_;
  Petra_Import * Importer_;
  Petra_BlockMap * ExportMap_;
  Petra_Export * Exporter_;

  bool Filled_;
  bool Allocated_;
  bool Sorted_;
  bool StorageOptimized_;
  bool NoRedundancies_;
  bool IndicesAreGlobal_;
  bool IndicesAreLocal_;
  bool IndicesAreContiguous_;
  bool LowerTriangular_;
  bool UpperTriangular_;
  bool NoDiagonal_;

  bool GlobalConstantsComputed_;
  
  int IndexBase_;
  
  int NumGlobalEntries_;
  int NumGlobalBlockRows_;
  int NumGlobalBlockCols_;
  int NumGlobalBlockDiagonals_;
  int NumMyEntries_;
  int NumMyBlockRows_;
  int NumMyBlockCols_;
  int NumMyBlockDiagonals_;
  
  int MaxRowDim_;
  int MaxColDim_;
  int GlobalMaxRowDim_;
  int GlobalMaxColDim_;
  int MaxNumNonzeros_;
  int GlobalMaxNumNonzeros_;
  
  int NumGlobalNonzeros_;
  int NumGlobalRows_;
  int NumGlobalCols_;
  int NumGlobalDiagonals_;
  int NumMyNonzeros_;
  int NumMyRows_;
  int NumMyCols_;
  int NumMyDiagonals_;
  
  int **Indices_;
  int * NumAllocatedIndicesPerRow_;
  int * NumIndicesPerRow_;
  int MaxNumIndices_;
  int GlobalMaxNumIndices_;
  int * All_Indices_;
  mutable int LenImports_;
  mutable int LenExports_;
  mutable int * Imports_;
  mutable int * Exports_;
  Petra_DataAccess CV_;
  
  
#ifdef PETRA_LEVELSCHEDULING
 public:
  //! Build level scheduling information for triangular matrix, upper or lower.
  /*! Computes level scheduling information for the current triangular graph.  
    
    \param In
           NumThreads - The number of threads intended for parallel execution.  Each level will be
	                partitioned so that each thread will get roughly the same number of nonzero
			terms at each level and thus perform approximately the same amount of work.
  */
  int ComputeLevels(int NumThreads);
  
 protected:
  

  int NumThreads_;
  int ** ThreadStartRows_;
  int NumLevels_;
  int * LevelIndices_;
  int * LevelOrder_;

  
 private:
  void listfill(int n, int head, int * link, int * vector);

#endif
  
};

//! << operator will work for Petra_CRS_Graph.
ostream& operator << (ostream& os, const Petra_CRS_Graph& A);

#endif /* _PETRA_CRS_GRAPH_H_ */
