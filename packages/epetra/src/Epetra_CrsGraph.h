
//@HEADER
/*
************************************************************************

              Epetra: Linear Algebra Services Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#ifndef EPETRA_CRSGRAPH_H
#define EPETRA_CRSGRAPH_H

#include "Epetra_DistObject.h" 
#include "Epetra_CrsGraphData.h"
class Epetra_BlockMap;
class Epetra_Util;
class Epetra_Time;
class Epetra_Import;
class Epetra_Export;
class Epetra_Distributor;
class Epetra_RowMatrix;

//! Epetra_CrsGraph: A class for constructing and using sparse compressed row graphs.

/*! The Epetra_CrsGraph enable the piecewise construction and use of sparse matrix graphs (the integer structure without
    values) where entries are intended for row access.  

    Epetra_CrsGraph is a base class for all row-based matrix classes.

<b>Constructing Epetra_CrsGraph objects</b>

Constructing Epetra_CrsGraph objects is a multi-step process.  The basic steps are as follows:
<ol>
  <li> Create Epetra_CrsGraph instance, including storage,  via constructor.
  <li> Enter row and column entry information via calls to the InsertIndices function.
  <li> Complete construction via FillComplete call.
  <li> (Optional) Optimize the graph storage via a call to OptimizeStorage.
</ol>

Note that, even after a matrix is constructed, it is possible to add or remove entries.  However, 
FillComplete must be
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
  Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int* NumIndicesPerRow);
  
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
									const Epetra_BlockMap& ColMap, int* NumIndicesPerRow);
  
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
	/*! This will create a Level 1 deep copy. This Graph will share ownership
		  of the CrsGraphData object with the right hand side Graph.
	*/
  Epetra_CrsGraph(const Epetra_CrsGraph& Graph);

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

    \return Integer error code, set to 0 if successful. Returns 1 if data is shared.
  */
  int InsertGlobalIndices(int GlobalRow, int NumIndices, int* Indices);
  
  //! Remove a list of elements from a specified global row of the matrix.
  /*!
    \param In
           Row - Global row number of indices.
    \param In
           NumIndices - Number of Indices.
    \param In
           Indices - Global column indices to remove.
	   
    \return Integer error code, set to 0 if successful. Returns 1 if data is shared.
  */
  int RemoveGlobalIndices(int GlobalRow, int NumIndices, int* Indices);
  
  //! Remove all indices from a specified global row of the matrix.
  /*!
    \param In
           Row - Global row number of indices.

    \return Integer error code, set to 0 if successful. Returns 1 if data is shared.
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

    \return Integer error code, set to 0 if successful. Returns 1 if data is shared.
  */
  int InsertMyIndices(int LocalRow, int NumIndices, int* Indices);
  
  //! Remove a list of elements from a specified local row of the matrix.
  /*!
    \param In
           Row - Local row number of indices.
    \param In
           NumIndices - Number of Indices.
    \param In
           Indices - Local column indices to remove.
	   
    \return Integer error code, set to 0 if successful. Returns 1 if data is shared.
  */
  int RemoveMyIndices(int LocalRow, int NumIndices, int* Indices);
  
  //! Remove all indices from a specified local row of the matrix.
  /*!
    \param In
           Row - Local row number of indices.

    \return Integer error code, set to 0 if successful. Returns 1 if data is shared.
  */
  int RemoveMyIndices(int Row);
  //@}

  //@{ \name Transformation methods
  
  //! Tranform matrix representation to local index space.  Perform other operations to allow optimal matrix operations.
  /*!
    \return Integer error code, set to 0 if successful. Returns 1 if data is shared.
  */
  int FillComplete();

  //! Tranform to local index space using specified Domain/Range maps.  Perform other operations to allow optimal matrix operations.
  /*!
    \return Integer error code, set to 0 if successful. Returns 1 if data is shared.
  */
  int FillComplete(const Epetra_BlockMap& DomainMap, const Epetra_BlockMap& RangeMap);

	//! Eliminates memory that is used for construction.  Make consecutive row index sections contiguous.
  /*!
    \return Integer error code, set to 0 if successful. Returns 1 if data is shared.
  */
	int OptimizeStorage();

	
  //! Sort column indices, row-by-row, in ascending order.
  /*!
    \return Integer error code, set to 0 if successful. Returns 1 if data is shared.
  */
  int SortIndices();


	//! Removes any redundant column indices in the rows of the graph.
  /*!
    \return Integer error code, set to 0 if successful. Returns 1 if data is shared.
  */
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
  int ExtractGlobalRowCopy(int GlobalRow, int LenOfIndices, int& NumIndices, int* Indices) const;

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
  int ExtractMyRowCopy(int LocalRow, int LenOfIndices, int& NumIndices, int* Indices) const;

  //! Get a view of the elements in a specified global row of the matrix.
  /*!
    This function requires that the graph not be completed (FillComplete() was \e not called).
    \param In
           Row - Local row number to get indices.
    \param Out
           NumIndices - Number of Indices.
    \param Out
           Indices - Column indices corresponding to values.

    \return Integer error code, set to 0 if successful. Returns -1 if invalid row.  Returns -2 if graph is completed.
  */
	int ExtractGlobalRowView(int GlobalRow, int& NumIndices, int*& Indices) const;

  //! Get a view of the elements in a specified local row of the matrix.
  /*!
    This function requires that the graph be completed FillComplete() was called).
    \param In
           Row - Local row number to get indices.
    \param Out
           NumIndices - Number of Indices.
    \param Out
           Indices - Column indices corresponding to values.

    \return Integer error code, set to 0 if successful. Returns -1 if invalid row.  Returns -2 if graph is not completed.
  */
	int ExtractMyRowView(int LocalRow, int& NumIndices, int*& Indices) const;
  //@}

  //@{ \name Graph Properties Query Methods.
  //! If FillComplete() has been called, this query returns true, otherwise it returns false.
	bool Filled() const {return(CrsGraphData_->Filled_);};
	
	//! If SortIndices() has been called, this query returns true, otherwise it returns false.
	bool Sorted() const {return(CrsGraphData_->Sorted_);};
	
	//! If OptimizeStorage() has been called, this query returns true, otherwise it returns false.
	bool StorageOptimized() const {return(CrsGraphData_->StorageOptimized_);};
	
	//! If column indices are in global range, this query returns true, otherwise it returns false.
	bool IndicesAreGlobal() const {return(CrsGraphData_->IndicesAreGlobal_);};
	
	//! If column indices are in local range, this query returns true, otherwise it returns false.
	bool IndicesAreLocal() const {return(CrsGraphData_->IndicesAreLocal_);};
	
	//! If graph is lower triangular in local index space, this query returns true, otherwise it returns false.
	bool LowerTriangular() const {return(CrsGraphData_->LowerTriangular_);};
	
	//! If graph is upper triangular in local index space, this query returns true, otherwise it returns false.
	bool UpperTriangular() const {return(CrsGraphData_->UpperTriangular_);};
	
	//! If graph has no diagonal entries in global index space, this query returns true, otherwise it returns false.
	bool NoDiagonal() const {return(CrsGraphData_->NoDiagonal_);};
	
	//! Returns true of GID is owned by the calling processor, otherwise it returns false.
	bool MyGlobalRow(int GID) const {return(RowMap().MyGID(GID));};

	//! Returns true if we have a well-defined ColMap, and returns false otherwise.
	/*! We have a well-defined ColMap if a) a ColMap was passed in at construction, 
		or b) the MakeColMap function has been called. (Calling either of the FillComplete functions
		will result in MakeColMap being called.) 
	*/
	bool HaveColMap() const {return(CrsGraphData_->HaveColMap_);};
  //@}
  
  //@{ \name Atribute access functions
    
	//! Returns the number of matrix rows on this processor.
	int NumMyRows() const {return(CrsGraphData_->NumMyRows_);};
	
	//! Returns the number of matrix rows in global matrix.
	int NumGlobalRows() const {return(CrsGraphData_->NumGlobalRows_);};
	
	//! Returns the number of matrix columns on this processor.
	int NumMyCols() const {return(CrsGraphData_->NumMyCols_);};
	
	//! Returns the number of matrix columns in global matrix.
	int NumGlobalCols() const {return(CrsGraphData_->NumGlobalCols_);};
	
	//! Returns the number of indices in the global graph.
	int NumGlobalNonzeros() const {return(CrsGraphData_->NumGlobalNonzeros_);};
	
	//! Returns the number of diagonal entries in the global graph, based on global row/column index comparisons.
	int NumGlobalDiagonals() const {return(CrsGraphData_->NumGlobalDiagonals_);};
	
	//! Returns the number of diagonal entries in the local graph, based on global row/column index comparisons.
	int NumMyDiagonals() const {return(CrsGraphData_->NumMyDiagonals_);};
	
	//! Returns the number of block matrix rows on this processor.
	int NumMyBlockRows() const {return(CrsGraphData_->NumMyBlockRows_);};
	
	//! Returns the number of Block matrix rows in global matrix.
	int NumGlobalBlockRows() const {return(CrsGraphData_->NumGlobalBlockRows_);};
	
	//! Returns the number of Block matrix columns on this processor.
	int NumMyBlockCols() const {return(CrsGraphData_->NumMyBlockCols_);};
	
	//! Returns the number of Block matrix columns in global matrix.
	int NumGlobalBlockCols() const {return(CrsGraphData_->NumGlobalBlockCols_);};
	
	//! Returns the number of Block diagonal entries in the local graph, based on global row/column index comparisons.
	int NumMyBlockDiagonals() const {return(CrsGraphData_->NumMyBlockDiagonals_);};
	
	//! Returns the number of Block diagonal entries in the global graph, based on global row/column index comparisons.
	int NumGlobalBlockDiagonals() const {return(CrsGraphData_->NumGlobalBlockDiagonals_);};
	
	//! Returns the number of entries in the global graph.
	int NumGlobalEntries() const {return(CrsGraphData_->NumGlobalEntries_);};
	
	//! Returns the number of entries on this processor.
	int NumMyEntries() const {return(CrsGraphData_->NumMyEntries_);};
	//! Returns the max row dimension of block entries on the processor.
	int MaxRowDim() const {return(CrsGraphData_->MaxRowDim_);};
	
	//! Returns the max row dimension of block entries across all processors.
	int GlobalMaxRowDim() const {return(CrsGraphData_->GlobalMaxRowDim_);};
	
	//! Returns the max column dimension of block entries on the processor.
	int MaxColDim() const {return(CrsGraphData_->MaxColDim_);};
	
	//! Returns the max column dimension of block entries across all processors.
	int GlobalMaxColDim() const {return(CrsGraphData_->GlobalMaxColDim_);};
	
	//! Returns the number of indices in the local graph.
	int NumMyNonzeros() const {return(CrsGraphData_->NumMyNonzeros_);};
	
	//! Returns the current number of nonzero entries in specified global row on this processor.
	int NumGlobalIndices(int Row) const;
	
	//! Returns the allocated number of nonzero entries in specified global row on this processor.
	int NumAllocatedGlobalIndices(int Row) const;
	
	//! Returns the maximum number of nonzero entries across all rows on this processor.
	int MaxNumIndices() const {return(CrsGraphData_->MaxNumIndices_);};
	
	//! Returns the maximun number of nonzero entries across all rows across all processors.
	int GlobalMaxNumIndices() const {return(CrsGraphData_->GlobalMaxNumIndices_);};
	
	//! Returns the maximum number of nonzero points across all rows on this processor.
	/*! For each entry in the graph, let i = the GRID of the entry and j = the CGID of the entry.  Then
	    the entry size is the product of the rowmap elementsize of i and the colmap elementsize of i.
	    Let ki = sum of all entry sizes for the entries in the ith row. 
	    For example,
            if the ith block row had 5 block entries and the element size of each entry was 4-by-4, ki would be 80.
	    Then this function return the max over all ki for all row on this processor.
	*/
	int MaxNumNonzeros() const {return(CrsGraphData_->MaxNumNonzeros_);};
	
	//! Returns the maximun number of nonzero points across all rows across all processors.
	/*! This function returns the max over all processor of MaxNumNonzeros().
	 */
	int GlobalMaxNumNonzeros() const {return(CrsGraphData_->GlobalMaxNumNonzeros_);};
	
	//! Returns the current number of nonzero entries in specified local row on this processor.
	int NumMyIndices(int Row) const {return(CrsGraphData_->NumIndicesPerRow_[Row]);};
	
	//! Returns the allocated number of nonzero entries in specified local row on this processor.
	int NumAllocatedMyIndices(int Row) const {return(CrsGraphData_->NumAllocatedIndicesPerRow_[Row]);};
	
	//! Returns the index base for row and column indices for this graph.
	int IndexBase() const {return(CrsGraphData_->IndexBase_);};
	
	//! Returns the RowMap associated with this matrix.
	const Epetra_BlockMap& RowMap() const {return(Epetra_DistObject::Map());};
    
	/** Replaces the current RowMap with the user-specified map object, but only
	    if currentmap->PointSameAs(newmap) is true. This is a collective function.
	    Returns 0 if map is replaced, -1 if not.
	*/
	int ReplaceRowMap(const Epetra_BlockMap& newmap);

	//! Returns the Column Map associated with this matrix.
	const Epetra_BlockMap& ColMap() const {return(CrsGraphData_->ColMap_);};
	
	//! Returns the DomainMap associated with this matrix.
	const Epetra_BlockMap& DomainMap() const {return(CrsGraphData_->DomainMap_);};
	
	//! Returns the RangeMap associated with this matrix.
	const Epetra_BlockMap& RangeMap() const {return(CrsGraphData_->RangeMap_);};
	
	//! Returns the Importer associated with this matrix.
	const Epetra_Import* Importer() const {return(CrsGraphData_->Importer_);};
	
	//! Returns the Exporter associated with this matrix.
	const Epetra_Export* Exporter() const {return(CrsGraphData_->Exporter_);};
	
	//! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
	const Epetra_Comm& Comm() const {return(Epetra_DistObject::Comm());};
  //@}
  
  //@{ \name Local/Global ID methods

	//! Returns the local row index for given global row index, returns -1 if no local row for this global row.
	int LRID(int GRID) const {return(RowMap().LID(GRID));};
	
	//! Returns the global row index for give local row index, returns IndexBase-1 if we don't have this local row.
	int GRID(int LRID) const {return(RowMap().GID(LRID));};
	
	//! Returns the local column index for given global column index, returns -1 if no local column for this global column.
	int LCID(int GCID) const
	  {
	    return( CrsGraphData_->HaveColMap_ ? ColMap().LID(GCID) : -1 );
	  }
	
	//! Returns the global column index for give local column index, returns IndexBase-1 if we don't have this local column.
	int GCID(int LCID) const
	  {
	    return( CrsGraphData_->HaveColMap_ ? ColMap().GID(LCID) : -1 );
	  }
	
	//! Returns true if the GRID passed in belongs to the calling processor in this map, otherwise returns false.
	bool MyGRID(int GRID) const {return(LRID(GRID) != -1);};
	
	//! Returns true if the LRID passed in belongs to the calling processor in this map, otherwise returns false.
	bool MyLRID(int LRID) const {return(GRID(LRID) != IndexBase() - 1);};
	
	//! Returns true if the GCID passed in belongs to the calling processor in this map, otherwise returns false.
	bool MyGCID(int GCID) const {return(LCID(GCID) != -1);};
	
	//! Returns true if the LRID passed in belongs to the calling processor in this map, otherwise returns false.
	bool MyLCID(int LCID) const {return(GCID(LCID) != IndexBase() - 1);};
  //@}
  
  //@{ \name Inlined Operator Methods.
	
	//! Inlined bracket operator for fast access to data. (Const and Non-const versions)
	/*! No error checking and dangerous for optimization purposes.
    \param Loc (In) - Local row.
	  
    \return reference to pointer to locally indexed Loc row in matrix.
  */

	//inline int* & operator[]( int Loc ) { return(CrsGraphData_->Indices_[Loc]); }
	//inline int* const & operator[]( int Loc ) const { return(CrsGraphData_->Indices_[Loc]); }
	/* NOTE: the commented out code above is the previous versions of operator[]. 
		 They returned a reference to int*, or reference to int* const (pointer is const,
		 data is not). This was deemed unnecessary, dangerous, (and now with ISDVs, 
		 sometimes not possible). All Epetra code compiles and works with the new versions 
		 (listed below). It is possible that some user code depends on the old functionality; 
		 those users should contact an Epetra developer.
	*/
	inline int* operator [] (int Loc) {return(CrsGraphData_->Sidekick()[Loc]);}
	inline int* operator [] (int Loc) const {return(CrsGraphData_->Sidekick()[Loc]);}

  //@}

	//! Assignment operator
	/*! This will do a Level 1 deep copy. It will share ownership of the CrsGraphData
		  with the right hand side Graph.
	*/
	Epetra_CrsGraph& operator = (const Epetra_CrsGraph& Source);

  //@{ \name I/O Methods.

  //! Print method
  virtual void Print(ostream& os) const;

	void PrintGraphData(ostream& os) const {CrsGraphData_->Print(os);};
	void PrintGraphData(ostream& os, int level) const {CrsGraphData_->Print(os, level);};
  //@}

  //@{ \name Deprecated methods:  These methods still work, but will be removed in a future version.

	//! Use ColMap() instead. 
	const Epetra_BlockMap& ImportMap() const {return(CrsGraphData_->ColMap_);};
	
	//! Use FillComplete() instead.
	int TransformToLocal();
	
	//! Use FillComplete(const Epetra_BlockMap& DomainMap, const Epetra_BlockMap& RangeMap) instead.
	int TransformToLocal(const Epetra_BlockMap* DomainMap, const Epetra_BlockMap* RangeMap);
	
  //@}

  //@{ \name Expert Users and Developers Only

	//! Returns the reference count of CrsGraphData.
	/*! (Intended for testing purposes.) */
	int ReferenceCount() const {return(CrsGraphData_->ReferenceCount());}

	//! Returns a pointer to the CrsGraphData instance this CrsGraph uses. 
	/*! (Intended for developer use only for testing purposes.) */
	const Epetra_CrsGraphData* DataPtr() const {return(CrsGraphData_);}

  //@}	

	// functions listed in protected are the ones used by CrsMatrix and VbrMatrix.
	// functions listed in private are the ones that are really private.
	// (just pretend CrsMatrix and VbrMatrix derive from CrsGraph to understand the distinction.)
	friend class Epetra_CrsMatrix;
	friend class Epetra_VbrMatrix;
	friend class Epetra_FECrsGraph;
	friend class Epetra_FECrsMatrix;
	friend class Epetra_FEVbrMatrix;
	friend class Epetra_OffsetIndex;

 protected:
	int* NumIndicesPerRow() const {return(CrsGraphData_->NumIndicesPerRow_.Values());}
	int* NumAllocatedIndicesPerRow() const {return(CrsGraphData_->NumAllocatedIndicesPerRow_.Values());}
	int** Indices() const {return(CrsGraphData_->Sidekick());}
	int* Indices(int LocalRow) {return(CrsGraphData_->Indices_[LocalRow].Values());}
	// If column indices are stored in one long array (via a call to OptimizeStorage), 
	// IndicesAreContiguous returns true, otherwise it returns false.
	bool IndicesAreContiguous() const {return(CrsGraphData_->IndicesAreContiguous_);}
	bool NoRedundancies() const {return(CrsGraphData_->NoRedundancies_);}
	bool GlobalConstantsComputed() const;
	bool FindGlobalIndexLoc(int LocalRow, int Index, int Start, int& Loc) const;
	bool FindGlobalIndexLoc(int NumIndices, const int* Indices, int Index, int Start, int& Loc) const;
	bool FindMyIndexLoc(int LocalRow, int Index, int Start, int& Loc) const;
	bool FindMyIndexLoc(int NumIndices, const int* Indices, int Index, int Start, int& Loc) const;
	int InsertIndices(int Row, int NumIndices, int* Indices);
	int MakeIndicesLocal(const Epetra_BlockMap& DomainMap, const Epetra_BlockMap& RangeMap);
	void SetIndicesAreLocal(bool Flag) {CrsGraphData_->IndicesAreLocal_ = Flag;}
	void SetIndicesAreGlobal(bool Flag) {CrsGraphData_->IndicesAreGlobal_ = Flag;}
	void SetSorted(bool Flag) {CrsGraphData_->Sorted_ = Flag;}

 private:
	void SetGlobalConstantsComputed(bool Flag) {CrsGraphData_->GlobalConstantsComputed_ = Flag;}
	void SetIndicesAreContiguous(bool Flag) {CrsGraphData_->IndicesAreContiguous_ = Flag;}
	void SetNoRedundancies(bool Flag) {CrsGraphData_->NoRedundancies_ = Flag;}
	void ComputeIndexState();
	int MakeColMap(const Epetra_BlockMap& DomainMap, const Epetra_BlockMap& RangeMap);
	int Allocate(int* NumIndicesPerRow, int Inc);
	int ReAllocate();
	int ComputeGlobalConstants();
	void SetFilled(bool Flag) {CrsGraphData_->Filled_ = Flag;}
	bool Allocated() const {return(CrsGraphData_->Allocated_);}
	void SetAllocated(bool Flag) {CrsGraphData_->Allocated_ = Flag;}
	
	int CheckSizes(const Epetra_SrcDistObject& A);

	int CopyAndPermute(const Epetra_SrcDistObject& Source,
                           int NumSameIDs, 
                           int NumPermuteIDs,
                           int* PermuteToLIDs,
                           int* PermuteFromLIDs,
                           const Epetra_OffsetIndex * Indexor);
	int CopyAndPermuteRowMatrix(const Epetra_RowMatrix& A,
                                    int NumSameIDs, 
                                    int NumPermuteIDs,
                                    int* PermuteToLIDs,
                                    int* PermuteFromLIDs,
                                    const Epetra_OffsetIndex * Indexor);
	int CopyAndPermuteCrsGraph(const Epetra_CrsGraph& A,
                                   int NumSameIDs, 
                                   int NumPermuteIDs,
                                   int* PermuteToLIDs,
                                   int* PermuteFromLIDs,
                                   const Epetra_OffsetIndex * Indexor);
  
        int PackAndPrepare(const Epetra_SrcDistObject& Source,
                           int NumExportIDs,
                           int* ExportLIDs,
                           int& LenExports,
                           char*& Exports,
                           int& SizeOfPacket,
                           int * Sizes,
                           bool & VarSizes,
                           Epetra_Distributor& Distor);
        int PackAndPrepareCrsGraph(const Epetra_CrsGraph& A,
                                   int NumExportIDs,
                                   int* ExportLIDs,
                                   int& LenExports,
                                   char*& Exports,
                                   int& SizeOfPacket,
                                   int* Sizes,
                                   bool& VarSizes,
                                   Epetra_Distributor& Distor);
        int PackAndPrepareRowMatrix(const Epetra_RowMatrix& A,
                                    int NumExportIDs,
                                    int* ExportLIDs,
                                    int& LenExports,
                                    char*& Exports,
                                    int& SizeOfPacket,
                                    int* Sizes,
                                    bool& VarSizes,
                                    Epetra_Distributor& Distor);

        int UnpackAndCombine(const Epetra_SrcDistObject& Source,
                             int NumImportIDs,
                             int* ImportLIDs, 
                             int LenImports, 
                             char* Imports,
                             int& SizeOfPacket, 
                             Epetra_Distributor& Distor,
                             Epetra_CombineMode CombineMode,
                             const Epetra_OffsetIndex * Indexor);

	void CleanupData();

	Epetra_CrsGraphData* CrsGraphData_;	

};
#endif /* EPETRA_CRSGRAPH_H */
