
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _EPETRA_CRSGRAPH_H_
#define _EPETRA_CRSGRAPH_H_

#include "Epetra_DistObject.h" 
#include "Epetra_BlockMap.h" 
class Epetra_Util;
class Epetra_Time;
class Epetra_Import;
class Epetra_Export;
class Epetra_Distributor; 

//! Epetra_CrsGraph: A class for constructing and using sparse compressed row graphs.

/*! The Epetra_CrsGraph enable the piecewise construction and use of sparse matrix graphs (the integer structure without
    values) where entries are intended for row access.  

    Epetra_CrsGraph is a base class for all row-based matrix classes.

<b>Constructing Epetra_CrsGraph objects</b>

Constructing Epetra_CrsGraph objects is a multi-step process.  The basic steps are as follows:
<ol>
  <li> Create Epetra_CrsGraph instance, including storage,  via constructor.
  <li> Enter row and column entry information via calls to the InsertIndices function.
  <li> Complete construction via TransformToLocal call.
  <li> (Optional) Optimize the graph storage via a call to StorageOptimize.
</ol>

Note that, even after a matrix is constructed, it is possible to add or remove entries.  However, 
TransformToLocal must be
called again before the graph is used for subsequent operations.



*/    

class Epetra_CrsGraph: public Epetra_DistObject {
      
 public:

  //@{ \name Constructors/Destructor.
  //! Epetra_CrsGraph constuctor with variable number of indices per row.
  /*! Creates a Epetra_CrsGraph object and allocates storage.  
    
    \param In
           CV - A Epetra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - An Epetra_BlockMap (or Epetra_Map or Epetra_LocalMap) listing the rows that this 
	   processor will contribute to.
    \param In
           NumIndicesPerRow - An integer array of length NumMyRows
	   such that NumIndicesPerRow[i] indicates the (approximate) number of entries in the ith row.
  */
  Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int *NumIndicesPerRow);
  
  //! Epetra_CrsGraph constuctor with fixed number of indices per row.
  /*! Creates a Epetra_CrsGraph object and allocates storage.  
    
    \param In
           CV - A Epetra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - An Epetra_BlockMap (or Epetra_Map or Epetra_LocalMap) listing the rows that this 
	   processor will contribute to.
    \param In
           NumIndicesPerRow - An integer that indicates the (approximate) number of entries in the each row.
	   Note that it is possible to use 0 for this value and let fill occur during the insertion phase.
	   
  */
  Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int NumIndicesPerRow);
  
  //! Epetra_CrsGraph constuctor with variable number of indices per row.
  /*! Creates a Epetra_CrsGraph object and allocates storage.  
    
    \param In
           CV - A Epetra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - An Epetra_BlockMap (or Epetra_Map or Epetra_LocalMap) listing the rows that this 
	   processor will contribute to.
    \param In 
           ColMap - An Epetra_BlockMap (or Epetra_Map or Epetra_LocalMap) listing the columns that this 
	   processor will contribute to.
    \param In
           NumIndicesPerRow - An integer array of length NumMyRows
	   such that NumIndicesPerRow[i] indicates the (approximate) number of entries in the ith row.
  */
  Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, 
		  const Epetra_BlockMap& ColMap, int *NumIndicesPerRow);
  
  //! Epetra_CrsGraph constuctor with fixed number of indices per row.
  /*! Creates a Epetra_CrsGraph object and allocates storage.  
    
    \param In
           CV - A Epetra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - An Epetra_BlockMap (or Epetra_Map or Epetra_LocalMap) listing the rows that this 
	   processor will contribute to.
    \param In 
           ColMap - An Epetra_BlockMap (or Epetra_Map or Epetra_LocalMap) listing the columns that this 
	   processor will contribute to.
    \param In
           NumIndicesPerRow - An integer that indicates the (approximate) number of entries in the each row.
	   Note that it is possible to use 0 for this value and let fill occur during the insertion phase.
	   
  */
  Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, 
		  const Epetra_BlockMap& ColMap, int NumIndicesPerRow);
  
  //! Copy constructor.
  Epetra_CrsGraph(const Epetra_CrsGraph & Graph);

  //! Epetra_CrsGraph Destructor
  virtual ~Epetra_CrsGraph();
  //@}
  
  //@{ \name Insertion/Removal methods.
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
  //@}

  //@{ \name Transformation methods
  
  //! Tranform matrix representation to local index space.  Perform other operations to allow optimal matrix operations.
  int TransformToLocal();

  //! Tranform to local index space using specified Domain/Range maps.  Perform other operations to allow optimal matrix operations.
  int TransformToLocal(Epetra_BlockMap *DomainMap, Epetra_BlockMap *RangeMap);

    //! Eliminates memory that is used for construction.  Make consecutive row index sections contiguous.
    int OptimizeStorage();


  //! Sort column indices, row-by-row, in ascending order.
  int SortIndices();


    //! Removes any redundant column indices in the rows of the graph.
    int RemoveRedundantIndices();
    //@}

  //@{ \name Extraction methods.

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
  //@}

  //@{ \name Graph Properties Query Methods.
    //! If TransformToLocal() has been called, this query returns true, otherwise it returns false.
    bool Filled() const {return(Filled_);};

    //! If SortIndices() has been called, this query returns true, otherwise it returns false.
    bool Sorted() const {return(Sorted_);};
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

    //! If graph has no diagonal entries, this query returns true, otherwise it returns false.
    bool NoDiagonal() const {return(NoDiagonal_);};

    //! Returns true of GID is owned by the calling processor, otherwise it returns false.
    bool MyGlobalRow(int GID) const {return(RowMap().MyGID(GID));};
  //@}
  
  //@{ \name Atribute access functions
    
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
    const Epetra_BlockMap& RowMap() const {return(Epetra_DistObject::Map());};
    
    //! Returns the Column Map associated with this matrix.
    const Epetra_BlockMap& ColMap() const {return(*ColMap_);};
    
    //! Returns the DomainMap associated with this matrix.
    const Epetra_BlockMap& DomainMap() const {return(*DomainMap_);};
    
    //! Returns the RangeMap associated with this matrix.
    const Epetra_BlockMap& RangeMap() const {return(*RangeMap_);};

    //! Returns the Importer associated with this matrix.
    const Epetra_Import * Importer() const {return(Importer_);};

    //! Returns the Exporter associated with this matrix.
    const Epetra_Export * Exporter() const {return(Exporter_);};

    //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
    const Epetra_Comm & Comm() const {return(Epetra_DistObject::Comm());};
  //@}
  
  //@{ \name Local/Global ID methods

    //! Returns the local row index for given global row index, returns -1 if no local row for this global row.
    int LRID( int GRID) const {return(RowMap().LID(GRID));};

    //! Returns the global row index for give local row index, returns IndexBase-1 if we don't have this local row.
    int GRID( int LRID) const {return(RowMap().GID(LRID));};

    //! Returns the local column index for given global column index, returns -1 if no local column for this global column.
    int LCID( int GCID) const {if (ColMap_==0) return(-1); else return(ColMap().LID(GCID));};

    //! Returns the global column index for give local column index, returns IndexBase-1 if we don't have this local column.
    int GCID( int LCID) const {if (ColMap_==0) return(-1); else return(ColMap().GID(LCID));};
 
    //! Returns true if the GRID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyGRID(int GRID) const {return(LRID(GRID)!=-1);};
   
    //! Returns true if the LRID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyLRID(int LRID) const {return(GRID(LRID)!=IndexBase_-1);};

    //! Returns true if the GCID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyGCID(int GCID) const {return(LCID(GCID)!=-1);};
   
    //! Returns true if the LRID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyLCID(int LCID) const {return(GCID(LCID)!=IndexBase_-1);};
  //@}
  
  //@{ \name Inlined Operator Methods.

    //! Inlined bracket operator for fast access to data. (Const and Non-const versions)
    /*! No error checking and dangerous for optimization purposes.
    \param Loc (In) - Local row.
	  
    \return reference to pointer to locally indexed Loc row in matrix.
  */
    inline int *& operator[]( int Loc ) { return Indices_[Loc]; }
    inline int * const & operator[]( int Loc ) const { return Indices_[Loc]; }
  //@}

  //@{ \name I/O Methods.

  //! Print method
  virtual void Print(ostream & os) const;
  //@}
  //@{ \name Deprecated methods:  These methods still work, but will be removed in a future version.


    //! Use ColMap() instead. 
    const Epetra_BlockMap& ImportMap() const {return(*ColMap_);};
  //@}

    
 protected:

    // If column indices are stored in one long array (via a call to OptimizeStorage), this query returns true, 
    // otherwise it returns false.
    bool IndicesAreContiguous() const {return(IndicesAreContiguous_);};
    int ** Indices() const {return(Indices_);};
    int * NumIndicesPerRow() const {return(NumIndicesPerRow_);};
    int * NumAllocatedIndicesPerRow() const {return(NumAllocatedIndicesPerRow_);};
    void SetSorted(bool Flag) {Sorted_ = Flag;};
    void SetGlobalConstantsComputed(bool Flag) {GlobalConstantsComputed_ = Flag;};
    void SetIndicesAreGlobal(bool Flag) {IndicesAreGlobal_ = Flag;};
    void SetIndicesAreLocal(bool Flag) {IndicesAreLocal_ = Flag;};
    void SetIndicesAreContiguous(bool Flag) {IndicesAreContiguous_ = Flag;};
    void SetNoRedundancies(bool Flag) {NoRedundancies_ = Flag;};
    int InsertIndices(int Row, int NumIndices, int *Indices);
    bool FindGlobalIndexLoc(int LocalRow, int Index, int Start, int & Loc);
    bool FindMyIndexLoc(int LocalRow, int Index, int Start, int &
 Loc);
    int MakeIndicesLocal(const Epetra_BlockMap & DomainMap, const Epetra_BlockMap & RangeMap);
    bool GlobalConstantsComputed() const;
    void ComputeIndexState();

    friend class Epetra_CrsMatrix;
    friend class Epetra_VbrMatrix;

 private:
    int MakeColMap(const Epetra_BlockMap & DomainMap, const Epetra_BlockMap & RangeMap);
    int MakeImportExport();
    void InitializeDefaults();
    int Allocate(int * NumIndicesPerRow, int Inc );
    int ComputeGlobalConstants();
    int ReAllocate();
    void SetFilled(bool Flag) {Filled_ = Flag;};
    bool Allocated() const {return(Allocated_);};
    bool NoRedundancies() const {return(NoRedundancies_);};
    void SetAllocated(bool Flag) {Allocated_ = Flag;};

    int CheckSizes(const Epetra_DistObject& A){return(0);};
    int CopyAndPermute(const Epetra_DistObject & Source,
		     int NumSameIDs, 
		     int NumPermuteIDs, int * PermuteToLIDs,
		     int *PermuteFromLIDs);
  
  int PackAndPrepare(const Epetra_DistObject & Source,
				     int NumExportIDs, int * ExportLIDs,
				     int Nsend, int Nrecv,
				     int & LenExports, char * & Exports, int & LenImports, 
				     char * & Imports, 
				     int & SizeOfPacket, Epetra_Distributor & Distor);
  
  int UnpackAndCombine(const Epetra_DistObject & Source,
		       int NumImportIDs, int * ImportLIDs, 
		       char * Imports, int & SizeOfPacket, 
		       Epetra_Distributor & Distor, Epetra_CombineMode CombineMode);


  // Defined by TransformToLocal and related
  Epetra_BlockMap * DomainMap_;
  Epetra_BlockMap * RangeMap_;
  Epetra_BlockMap * ColMap_;
  Epetra_Import * Importer_;
  Epetra_Export * Exporter_;

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
  Epetra_DataAccess CV_;
  
  

};
#endif /* _EPETRA_CRSGRAPH_H_ */
