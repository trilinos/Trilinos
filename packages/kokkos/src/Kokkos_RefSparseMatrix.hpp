
//@HEADER
// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#ifndef KOKKOS_CRSMATRIX_H
#define KOKKOS_CRSMATRIX_H

#include "Kokkos_DistObject.hpp" 
#include "Kokkos_CompObject.hpp" 
#include "Kokkos_BLAS.hpp"
#include "Kokkos_RowMatrix.hpp"
#include "Kokkos_Operator.hpp"
#include "Kokkos_CrsGraph.hpp"
class Kokkos_Map;
class Kokkos_Import;
class Kokkos_Export;
class Kokkos_Vector;
class Kokkos_MultiVector;

//! Kokkos_CrsMatrix: A class for constructing and using real-valued double-precision sparse compressed row matrices.

/*! The Kokkos_CrsMatrix enables the piecewise construction and use of real-valued double-precision sparse matrices
    where matrix entries are intended for row access.

    At this time, the primary function provided by Kokkos_CrsMatrix is matrix times vector and matrix 
    times multi-vector multiplication.  It is also possible to extract matrix rows from a constructed matrix.

<b>Constructing Kokkos_CrsMatrix objects</b>

Constructing Kokkos_CrsMatrix objects is a multi-step process.  The basic steps are as follows:
<ol>
  <li> Create Kokkos_CrsMatrix instance, including storage,  via constructor.
  <li> Enter values via one or more Put or SumInto functions.
  <li> Complete construction via FillComplete call.
</ol>

Note that, even after a matrix is constructed, it is possible to update existing matrix entries.  It is \e not possible to
create new entries.

<b> Counting Floating Point Operations </b>

Each Kokkos_CrsMatrix object keeps track of the number
of \e serial floating point operations performed using the specified object as the \e this argument
to the function.  The Flops() function returns this number as a double precision number.  Using this 
information, in conjunction with the Kokkos_Time class, one can get accurate parallel performance
numbers.  The ResetFlops() function resets the floating point counter.

\warning A Kokkos_Map is required for the Kokkos_CrsMatrix constructor.

*/    

class Kokkos_CrsMatrix: public Kokkos_DistObject, public Kokkos_CompObject, public Kokkos_BLAS, public virtual Kokkos_RowMatrix {
 public:

  //@{ \name Constructors/Destructor.
  //! Kokkos_CrsMatrix constuctor with variable number of indices per row.
  /*! Creates a Kokkos_CrsMatrix object and allocates storage.  
    
	\param In
	CV - A Kokkos_DataAccess enumerated type set to Copy or View.
	\param In 
	RowMap - An Kokkos_Map listing the rows that this processor will contribute to.
	\param In
	NumEntriesPerRow - An integer array of length NumRows
	such that NumEntriesPerRow[i] indicates the (approximate) number of entries in the ith row.
  */
  Kokkos_CrsMatrix(Kokkos_DataAccess CV, const Kokkos_Map& RowMap, int* NumEntriesPerRow);
  
  //! Kokkos_CrsMatrix constuctor with fixed number of indices per row.
  /*! Creates a Kokkos_CrsMatrix object and allocates storage.  
    
	\param In
	CV - A Kokkos_DataAccess enumerated type set to Copy or View.
	\param In 
	RowMap - An Kokkos_Map listing the rows that this processor will contribute to.
	\param In
	NumEntriesPerRow - An integer that indicates the (approximate) number of entries in the each row.
	Note that it is possible to use 0 for this value and let fill occur during the insertion phase.
	
  */
  Kokkos_CrsMatrix(Kokkos_DataAccess CV, const Kokkos_Map& RowMap, int NumEntriesPerRow);

  //! Kokkos_CrsMatrix constuctor with variable number of indices per row.
  /*! Creates a Kokkos_CrsMatrix object and allocates storage.  
    
	\param In
	CV - A Kokkos_DataAccess enumerated type set to Copy or View.
	\param In 
	RowMap - An Kokkos_Map listing the rows that this processor will contribute to.
	\param In 
	ColMap - An Kokkos_Map listing the columns that this processor will contribute to.
	\param In
	NumEntriesPerRow - An integer array of length NumRows
	such that NumEntriesPerRow[i] indicates the (approximate) number of entries in the ith row.
  */
  Kokkos_CrsMatrix(Kokkos_DataAccess CV, const Kokkos_Map& RowMap, const Kokkos_Map& ColMap, int* NumEntriesPerRow);
  
  //! Kokkos_CrsMatrix constuctor with fixed number of indices per row.
  /*! Creates a Kokkos_CrsMatrix object and allocates storage.  
    
	\param In
	CV - A Kokkos_DataAccess enumerated type set to Copy or View.
	\param In 
	RowMap - An Kokkos_Map listing the rows that this processor will contribute to.
	\param In 
	ColMap - An Kokkos_Map listing the columns that this processor will contribute to.
	\param In
	NumEntriesPerRow - An integer that indicates the (approximate) number of entries in the each row.
	Note that it is possible to use 0 for this value and let fill occur during the insertion phase.
	
  */
  Kokkos_CrsMatrix(Kokkos_DataAccess CV, const Kokkos_Map& RowMap, const Kokkos_Map& ColMap, int NumEntriesPerRow);

  //! Construct a matrix using an existing Kokkos_CrsGraph object.
  /*! Allows the nonzero structure from another matrix, or a structure that was
		constructed independently, to be used for this matrix.
    \param In
		CV - A Kokkos_DataAccess enumerated type set to Copy or View.
    \param In
		Graph - A Kokkos_CrsGraph object, extracted from another Kokkos matrix object or constructed directly from
		using the Kokkos_CrsGraph constructors.
  */
	
  Kokkos_CrsMatrix(Kokkos_DataAccess CV, const Kokkos_CrsGraph& Graph);
	
  //! Copy constructor.
  Kokkos_CrsMatrix(const Kokkos_CrsMatrix& Matrix);
	
  //! Kokkos_CrsMatrix Destructor
  virtual ~Kokkos_CrsMatrix();
  //@}
  
  //@{ \name Insertion/Replace/SumInto methods.
	
  //! Initialize all values in the matrix with constant value.
  /*!
    \param In
		ScalarConstant - Value to use.
		
    \return Integer error code, set to 0 if successful.
  */
	int PutScalar(double ScalarConstant);
	
  //! Multiply all values in the matrix by a constant value (in place: A <- ScalarConstant * A).
  /*!
    \param In
		ScalarConstant - Value to use.
		
    \return Integer error code, set to 0 if successful.
  */
	int Scale(double ScalarConstant);
	
  //! Insert a list of elements in a given global row of the matrix.
  /*!
    \param In
		GlobalRow - Row number (in global coordinates) to put elements.
    \param In
		NumEntries - Number of entries.
    \param In
		Values - Values to enter.
    \param In
		Indices - Global column indices corresponding to values.
		
    \return Integer error code, set to 0 if successful. Note that if the
    allocated length of the row has to be expanded, a positive warning code
    will be returned.
  */
	int InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
	
  //! Replace current values with this list of entries for a given global row of the matrix.
  /*!
    \param In
		GlobalRow - Row number (in global coordinates) to put elements.
    \param In
		NumEntries - Number of entries.
    \param In
		Values - Values to enter.
    \param In
		Indices - Global column indices corresponding to values.
		
    \return Integer error code, set to 0 if successful. Note that if a value
    is not already present for the specified location in the matrix, the
    input value will be ignored and a positive warning code will be returned.
  */
	int ReplaceGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
	
  //! Add this list of entries to existing values for a given global row of the matrix.
  /*!
    \param In
		GlobalRow - Row number (in global coordinates) to put elements.
    \param In
		NumEntries - Number of entries.
    \param In
		Values - Values to enter.
    \param In
		Indices - Global column indices corresponding to values.
		
    \return Integer error code, set to 0 if successful. Note that if a value
    is not already present for the specified location in the matrix, the
    input value will be ignored and a positive warning code will be returned.
  */
	int SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);

  //! Insert a list of elements in a given local row of the matrix.
  /*!
    \param In
		MyRow - Row number (in local coordinates) to put elements.
    \param In
		NumEntries - Number of entries.
    \param In
		Values - Values to enter.
    \param In
		Indices - Local column indices corresponding to values.
		
    \return Integer error code, set to 0 if successful. Note that if the
    allocated length of the row has to be expanded, a positive warning code
    will be returned.
  */
	int InsertMyValues(int MyRow, int NumEntries, double* Values, int* Indices);

  //! Replace current values with this list of entries for a given local row of the matrix.
  /*!
    \param In
		MyRow - Row number (in local coordinates) to put elements.
    \param In
		NumEntries - Number of entries.
    \param In
		Values - Values to enter.
    \param In
		Indices - Local column indices corresponding to values.
		
    \return Integer error code, set to 0 if successful. Note that if a value
    is not already present for the specified location in the matrix, the
    input value will be ignored and a positive warning code will be returned.
  */
	int ReplaceMyValues(int MyRow, int NumEntries, double* Values, int* Indices);

  //! Add this list of entries to existing values for a given local row of the matrix.
  /*!
    \param In
		MyRow - Row number (in local coordinates) to put elements.
    \param In
		NumEntries - Number of entries.
    \param In
		Values - Values to enter.
    \param In
		Indices - Local column indices corresponding to values.
		
    \return Integer error code, set to 0 if successful. Note that if the
    allocated length of the row has to be expanded, a positive warning code
    will be returned.
  */
	int SumIntoMyValues(int MyRow, int NumEntries, double* Values, int* Indices);

	//! Replaces diagonal values of the with those in the user-provided vector.
	/*! This routine is meant to allow replacement of {\bf existing} diagonal values.
		If a diagonal value does not exist for a given row, the corresponding value in
		the input Kokkos_Vector will be ignored and the return code will be set to 1.
		
		The Kokkos_Map associated with the input Kokkos_Vector must be compatible with
		the RowMap of the matrix.
		
    \param Diagonal (In) - New values to be placed in the main diagonal.
		
    \return Integer error code, set to 0 if successful, 1 of one or more diagonal entries not present in matrix.
  */
	int ReplaceDiagonalValues(const Kokkos_Vector& Diagonal);
	
  //@}

  //@{ \name Transformation methods
  

	//! Signal that data entry is complete.  Perform transformations to local index space.
	/* This version of FillComplete assumes that the domain and range distributions are
		 identical to the matrix row distributions.
	*/
	int FillComplete();

	//! Signal that data entry is complete.  Perform transformations to local index space.
	/* This version of FillComplete requires the explicit specification of the domain
		 and range distribution maps.  These maps are used for importing and exporting vector
		 and multi-vector elements that are needed for distributed matrix computations.  For
		 example, to compute y = Ax in parallel, we would specify the DomainMap as the distribution
		 of the vector x and the RangeMap as the distribution of the vector y.
    \param In
		DomainMap - Map that describes the distribution of vector and multi-vectors in the
		matrix domain.
    \param In
		RangeMap - Map that describes the distribution of vector and multi-vectors in the
		matrix range.
    */
	int FillComplete(const Kokkos_Map& DomainMap, const Kokkos_Map& RangeMap);
    
	//! Sort column entries, row-by-row, in ascending order.
	int SortEntries();

	//! Add entries that have the same column index. Remove redundant entries from list.
	int MergeRedundantEntries();
	
	//! Analyzes matrix and attempts to optimize storage for matrix operations.
	int OptimizeStorage();
	
	//! Eliminates memory that is used for construction.  Make consecutive row index sections contiguous.
	int MakeDataContiguous() {KOKKOS_CHK_ERR(OptimizeStorage()); return(0);};
	//@}
	
  //@{ \name Extraction methods.
	
	//! Returns a copy of the specified global row in user-provided arrays.
	/*! 
    \param In
		GlobalRow - Global row to extract.
    \param In
		Length - Length of Values and Indices.
    \param Out
		NumEntries - Number of nonzero entries extracted.
    \param Out
		Values - Extracted values for this row.
    \param Out
		Indices - Extracted global column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful.
  */
	int ExtractGlobalRowCopy(int GlobalRow, int Length, int& NumEntries, double* Values, int* Indices) const;

	//! Returns a copy of the specified local row in user-provided arrays.
	/*! 
    \param In
		MyRow - Local row to extract.
    \param In
		Length - Length of Values and Indices.
    \param Out
		NumEntries - Number of nonzero entries extracted.
    \param Out
		Values - Extracted values for this row.
    \param Out
	   	Indices - Extracted local column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful.
  */
	int ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, double* Values, int* Indices) const;

	//! Returns a copy of the specified global row values in user-provided array.
	/*! 
    \param In
		GlobalRow - Global row to extract.
    \param In
		Length - Length of Values and Indices.
    \param Out
		NumEntries - Number of nonzero entries extracted.
    \param Out
		Values - Extracted values for this row.
	  
    \return Integer error code, set to 0 if successful.
  */
	int ExtractGlobalRowCopy(int GlobalRow, int Length, int& NumEntries, double* Values) const;

	//! Returns a copy of the specified local row values in user-provided array.
	/*! 
    \param In
		MyRow - Local row to extract.
    \param In
		Length - Length of Values and Indices.
    \param Out
		NumEntries - Number of nonzero entries extracted.
    \param Out
		Values - Extracted values for this row.
	  
    \return Integer error code, set to 0 if successful.
  */
	int ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, double* Values) const;

	//! Returns a copy of the main diagonal in a user-provided vector.
	/*! 
    \param Out
		Diagonal - Extracted main diagonal.
		
    \return Integer error code, set to 0 if successful.
  */
	int ExtractDiagonalCopy(Kokkos_Vector& Diagonal) const;
	
	//! Returns a view of the specified global row values via pointers to internal data.
	/*! 
    \param In
		GlobalRow - Global row to view.
    \param Out
		NumEntries - Number of nonzero entries extracted.
    \param Out
		Values - Extracted values for this row.
    \param Out
		Indices - Extracted global column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful. Returns -1 of row not on this processor.  
		Returns -2 if matrix is not in global form (if FillComplete() has already been called).
  */
	int ExtractGlobalRowView(int GlobalRow, int& NumEntries, double*& Values, int*& Indices) const;
	
	//! Returns a view of the specified local row values via pointers to internal data.
	/*! 
    \param In
		MyRow - Local row to view.
    \param Out
		NumEntries - Number of nonzero entries extracted.
    \param Out
		Values - Extracted values for this row.
    \param Out
		Indices - Extracted local column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful. Returns -1 of row not on this processor.  
		Returns -2 if matrix is not in local form (if FillComplete() has \e not been called).
  */
	int ExtractMyRowView(int MyRow, int& NumEntries, double*& Values, int*& Indices) const;
	
	//! Returns a view of the specified global row values via pointers to internal data.
	/*! 
    \param In
		GlobalRow - Global row to extract.
    \param Out
		NumEntries - Number of nonzero entries extracted.
    \param Out
		Values - Extracted values for this row.
	  
    \return Integer error code, set to 0 if successful.
  */
	int ExtractGlobalRowView(int GlobalRow, int& NumEntries, double*& Values) const;

	//! Returns a view of the specified local row values via pointers to internal data.
	/*! 
    \param In
		MyRow - Local row to extract.
    \param Out
		NumEntries - Number of nonzero entries extracted.
    \param Out
		Values - Extracted values for this row.
	  
    \return Integer error code, set to 0 if successful.
  */
	int ExtractMyRowView(int MyRow, int& NumEntries, double*& Values) const;
	//@}
	
  //@{ \name Computational methods.
	
	//! Returns the result of a Kokkos_CrsMatrix multiplied by a Kokkos_Vector x in y.
	/*! 
    \param In
		TransA -If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param In
		x -A Kokkos_Vector to multiply by.
    \param Out
		y -A Kokkos_Vector containing result.
		
    \return Integer error code, set to 0 if successful.
  */
	int Multiply(bool TransA, const Kokkos_Vector& x, Kokkos_Vector& y) const;

	//! Returns the result of a Kokkos_CrsMatrix multiplied by a Kokkos_MultiVector X in Y.
	/*! 
    \param In
		TransA -If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param In
		X - A Kokkos_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
		Y -A Kokkos_MultiVector of dimension NumVectorscontaining result.
		
    \return Integer error code, set to 0 if successful.
  */
	int Multiply(bool TransA, const Kokkos_MultiVector& X, Kokkos_MultiVector& Y) const;

	//! Returns the result of a local solve using the Kokkos_CrsMatrix on a Kokkos_Vector x in y.
	/*! This method solves a triangular system of equations asynchronously on each processor.
    \param In
		Upper -If true, solve Uy = x, otherwise solve Ly = x.
    \param In
		Trans -If true, solve transpose problem.
    \param In
		UnitDiagonal -If true, assume diagonal is unit (whether it's stored or not).
    \param In
		x -A Kokkos_Vector to solve for.
    \param Out
		y -A Kokkos_Vector containing result.
		
    \return Integer error code, set to 0 if successful.
  */
	int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Kokkos_Vector& x, Kokkos_Vector& y) const;
	
	//! Returns the result of a local solve using the Kokkos_CrsMatrix a Kokkos_MultiVector X in Y.
	/*! This method solves a triangular system of equations asynchronously on each processor.
    \param In
		Upper -If true, solve Uy = x, otherwise solve Ly = x.
    \param In
		Trans -If true, solve transpose problem.
    \param In
		UnitDiagonal -If true, assume diagonal is unit (whether it's stored or not).
    \param In
		X - A Kokkos_MultiVector of dimension NumVectors to solve for.
    \param Out
		Y -A Kokkos_MultiVector of dimension NumVectors containing result.
		
    \return Integer error code, set to 0 if successful.
  */
	int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Kokkos_MultiVector& X, Kokkos_MultiVector& Y) const;
	
	//! Computes the sum of absolute values of the rows of the Kokkos_CrsMatrix, results returned in x.
	/*! The vector x will return such that x[i] will contain the inverse of sum of the absolute values of the 
		\e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the global row number of A
		and j denotes the global column number of A.  Using the resulting vector from this function as input to LeftScale()
		will make the infinity norm of the resulting matrix exactly 1.
    \param Out
		x -A Kokkos_Vector containing the row sums of the \e this matrix. 
		\warning It is assumed that the distribution of x is the same as the rows of \e this.
		
    \return Integer error code, set to 0 if successful.
  */
	int InvRowSums(Kokkos_Vector& x) const;
	
	//! Scales the Kokkos_CrsMatrix on the left with a Kokkos_Vector x.
	/*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
		and j denotes the column number of A.
    \param In
		x -A Kokkos_Vector to solve for.
		
    \return Integer error code, set to 0 if successful.
  */
	int LeftScale(const Kokkos_Vector& x);
	
	//! Computes the sum of absolute values of the columns of the Kokkos_CrsMatrix, results returned in x.
	/*! The vector x will return such that x[j] will contain the inverse of sum of the absolute values of the 
		\e this matrix will be sca such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
		and j denotes the global column number of A.  Using the resulting vector from this function as input to 
		RighttScale() will make the one norm of the resulting matrix exactly 1.
    \param Out
		x -A Kokkos_Vector containing the column sums of the \e this matrix. 
		\warning It is assumed that the distribution of x is described by the DomainMap() of \e this.
		
    \return Integer error code, set to 0 if successful.
  */
	int InvColSums(Kokkos_Vector& x) const;

	//! Scales the Kokkos_CrsMatrix on the right with a Kokkos_Vector x.
	/*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
		and j denotes the global column number of A.
    \param In
		x -The Kokkos_Vector used for scaling \e this.
		
    \return Integer error code, set to 0 if successful.
  */
	int RightScale(const Kokkos_Vector& x);
  //@}
	
  //@{ \name Matrix Properties Query Methods.
	
	
	//! If MergeRedundantEntries() has been called, this query returns true, otherwise it returns false.
	bool NoRedundancies() const {return(Graph_->NoRedundancies());};
	
	//! If SortEntries() has been called, this query returns true, otherwise it returns false.
	bool Sorted() const {return(Graph_->Sorted());};
	
	//! If FillComplete() has been called, this query returns true, otherwise it returns false.
	bool Filled() const {return(Graph_->Filled());};
	
	//! If OptimizeStorage() has been called, this query returns true, otherwise it returns false.
	bool StorageOptimized() const {return(Graph_->StorageOptimized());};
	
	//! If matrix indices has not been transformed to local, this query returns true, otherwise it returns false.
	bool IndicesAreGlobal() const {return(Graph_->IndicesAreGlobal());};
	
	//! If matrix indices has been transformed to local, this query returns true, otherwise it returns false.
	bool IndicesAreLocal() const {return(Graph_->IndicesAreLocal());};
	
	//! If matrix indices are packed into single array (done in OptimizeStorage()) return true, otherwise false.
	bool IndicesAreContiguous() const {return(Graph_->IndicesAreContiguous());};
	
	//! If matrix is lower triangular, this query returns true, otherwise it returns false.
	bool LowerTriangular() const {return(Graph_->LowerTriangular());};
	
	//! If matrix is upper triangular, this query returns true, otherwise it returns false.
	bool UpperTriangular() const {return(Graph_->UpperTriangular());};
	
	//! If matrix has no diagonal entries, this query returns true, otherwise it returns false.
	bool NoDiagonal() const {return(Graph_->NoDiagonal());};
	
  //@}
  
  //@{ \name Atribute access functions
	
	//! Returns the infinity norm of the global matrix.
	/* Returns the quantity \f$ \| A \|_\infty\f$ such that
		 \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].
	*/ 
	double NormInf() const;
	
	//! Returns the one norm of the global matrix.
	/* Returns the quantity \f$ \| A \|_1\f$ such that
		 \f[\| A \|_1= \max_{1\lej\len} \sum_{i=1}^m |a_{ij}| \f].
	*/ 
	double NormOne() const;
	
	//! Returns the number of nonzero entries in the global matrix.
	int NumGlobalNonzeros() const {return(Graph_->NumGlobalNonzeros());};
	
	//! Returns the number of global matrix rows.
	int NumGlobalRows() const {return(Graph_->NumGlobalRows());};
	
	//! Returns the number of global matrix columns.
	int NumGlobalCols() const {return(Graph_->NumGlobalCols());};
	
	//! Returns the number of global nonzero diagonal entries.
	int NumGlobalDiagonals() const {return(Graph_->NumGlobalDiagonals());};
	
	//! Returns the number of nonzero entries in the calling processor's portion of the matrix.
	int NumMyNonzeros() const {return(Graph_->NumMyNonzeros());};
	
	//! Returns the number of matrix rows owned by the calling processor.
	int NumMyRows() const {return(Graph_->NumMyRows());};
	
	//! Returns the number of matrix columns owned by the calling processor.
	int NumMyCols() const {return(Graph_->NumMyCols());};
	
	//! Returns the number of local nonzero diagonal entries.
	int NumMyDiagonals() const {return(Graph_->NumMyDiagonals());};
	
	//! Returns the current number of nonzero entries in specified global row on this processor.
	int NumGlobalEntries(int Row) const {return(Graph_->NumGlobalIndices(Row));};
	
	//! Returns the allocated number of nonzero entries in specified global row on this processor.
	int NumAllocatedGlobalEntries(int Row) const{return(Graph_->NumAllocatedGlobalIndices(Row));};
	
	//! Returns the maximum number of nonzero entries across all rows on this processor.
	int MaxNumEntries() const {return(Graph_->MaxNumIndices());};

	//! Returns the maximum number of nonzero entries across all rows on this processor.
	int GlobalMaxNumEntries() const {return(Graph_->GlobalMaxNumIndices());};
	
	//! Returns the current number of nonzero entries in specified local row on this processor.
	int NumMyEntries(int Row) const {return(Graph_->NumMyIndices(Row));};
	
	//! Returns the allocated number of nonzero entries in specified local row on this processor.
	int NumAllocatedMyEntries(int Row) const {return(Graph_->NumAllocatedMyIndices(Row));};
	
	//! Returns the index base for row and column indices for this graph.
	int IndexBase() const {return(Graph_->IndexBase());};
	
	
	//! Returns true if the graph associated with this matrix was pre-constructed and therefore not changeable.
	bool StaticGraph() {return(StaticGraph_);};
	//! Returns a pointer to the Kokkos_CrsGraph object associated with this matrix.
	const Kokkos_CrsGraph& Graph() const {return(*Graph_);};
	
	//! Returns the Kokkos_Map object associated with the rows of this matrix.
	const Kokkos_Map& RowMap() const {return((Kokkos_Map &)Graph_->RowMap());};

	//! Replaces the current RowMap with the user-specified map object.
	int ReplaceRowMap(const Kokkos_BlockMap& newmap)
	  {return( Graph_->ReplaceRowMap(newmap) ); }

	//! Returns the Kokkos_Map object that describes the column distribution across processors.
	const Kokkos_Map& ColMap() const {return((Kokkos_Map &) Graph_->ColMap());};
	
	//! Returns the Kokkos_Map object associated with the domain of this matrix operator.
	const Kokkos_Map& DomainMap() const {return((Kokkos_Map &)Graph_->DomainMap());}
	
	//! Returns the Kokkos_Map object associated with the range of this matrix operator.
	const Kokkos_Map& RangeMap() const  {return((Kokkos_Map &)Graph_->RangeMap());}
	
	//! Returns the Kokkos_Import object that contains the import operations for distributed operations.
	const Kokkos_Import* Importer() const {return(Graph_->Importer());};
	
	//! Returns the Kokkos_Export object that contains the export operations for distributed operations.
	const Kokkos_Export* Exporter() const {return(Graph_->Exporter());};
	
	//! Returns a pointer to the Kokkos_Comm communicator associated with this matrix.
	const Kokkos_Comm& Comm() const {return(Kokkos_DistObject::Comm());};
  //@}
  
  //@{ \name Local/Global ID methods
	//! Returns the local row index for given global row index, returns -1 if no local row for this global row.
	int LRID( int GRID) const {return(Graph_->LRID(GRID));};
	
	//! Returns the global row index for give local row index, returns IndexBase-1 if we don't have this local row.
	int GRID( int LRID) const {return(Graph_->GRID(LRID));};
	
	//! Returns the local column index for given global column index, returns -1 if no local column for this global column.
	int LCID( int GCID) const {return(Graph_->LCID(GCID));};
	
	//! Returns the global column index for give local column index, returns IndexBase-1 if we don't have this local column.
	int GCID( int LCID) const {return(Graph_->GCID(LCID));};
	
	//! Returns true if the GRID passed in belongs to the calling processor in this map, otherwise returns false.
	bool MyGRID(int GRID) const {return(Graph_->MyGRID(GRID));};
	
	//! Returns true if the LRID passed in belongs to the calling processor in this map, otherwise returns false.
	bool MyLRID(int LRID) const {return(Graph_->MyLRID(LRID));};
	
	//! Returns true if the GCID passed in belongs to the calling processor in this map, otherwise returns false.
	bool MyGCID(int GCID) const {return(Graph_->MyGCID(GCID));};
   
	//! Returns true if the LRID passed in belongs to the calling processor in this map, otherwise returns false.
	bool MyLCID(int LCID) const {return(Graph_->MyLCID(LCID));};

	//! Returns true of GID is owned by the calling processor, otherwise it returns false.
	bool MyGlobalRow(int GID) const {return(Graph_->MyGlobalRow(GID));};
  //@}
  
  
  //@{ \name I/O Methods.

  //! Print method
  virtual void Print(ostream& os) const;
  //@}

  //@{ \name Additional methods required to support the Kokkos_Operator interface.

	//! Returns a character string describing the operator
	char* Label() const {return(Kokkos_Object::Label());};
	
	//! If set true, transpose of this operator will be applied.
	/*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
		affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
		does not support transpose use, this method should return a value of -1.
		
    \param In
		UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.
		
    \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

	//! Returns the result of a Kokkos_Operator applied to a Kokkos_MultiVector X in Y.
	/*! 
    \param In
		X - A Kokkos_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
		Y -A Kokkos_MultiVector of dimension NumVectors containing result.
		
    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Kokkos_MultiVector& X, Kokkos_MultiVector& Y) const {
    return(Kokkos_CrsMatrix::Multiply(Kokkos_CrsMatrix::UseTranspose(), X, Y));};
	
	//! Returns the result of a Kokkos_Operator inverse applied to an Kokkos_MultiVector X in Y.
	/*! In this implementation, we use several existing attributes to determine how virtual
		method ApplyInverse() should call the concrete method Solve().  We pass in the UpperTriangular(), 
		the Kokkos_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
		if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
		be accounted for when doing a triangular solve.
		
		\param In
		X - A Kokkos_MultiVector of dimension NumVectors to solve for.
    \param Out
		Y -A Kokkos_MultiVector of dimension NumVectors containing result.
		
    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Kokkos_MultiVector& X, Kokkos_MultiVector& Y) const {
    return(Solve(UpperTriangular(), Kokkos_CrsMatrix::UseTranspose(), NoDiagonal(), X, Y));};

	//! Returns true because this class can compute an Inf-norm.
	bool HasNormInf() const {return(true);};
	
	//! Returns the current UseTranspose setting.
	bool UseTranspose() const {return(UseTranspose_);};
	
	//! Returns the Kokkos_Map object associated with the domain of this matrix operator.
	const Kokkos_Map& OperatorDomainMap() const {return(DomainMap());};
	
	//! Returns the Kokkos_Map object associated with the range of this matrix operator.
	const Kokkos_Map& OperatorRangeMap() const {return(RangeMap());};

  //@}
  //@{ \name Additional methods required to implement RowMatrix interface.

	//! Return the current number of values stored for the specified local row.
	/*! Similar to NumMyEntries() except NumEntries is returned as an argument
		and error checking is done on the input value MyRow.
    \param In
		MyRow - Local row.
    \param Out
		NumEntries - Number of nonzero values.
	  
    \return Integer error code, set to 0 if successful.
  */
	int NumMyRowEntries(int MyRow, int& NumEntries) const;

	//! Returns the Kokkos_Map object associated with the rows of this matrix.
	const Kokkos_Map& RowMatrixRowMap() const {return(RowMap());};
	
	//! Returns the Kokkos_Map object associated with columns of this matrix.
	const Kokkos_Map& RowMatrixColMap() const {return(ColMap());};
	
	//! Returns the Kokkos_Import object that contains the import operations for distributed operations.
	const Kokkos_Import* RowMatrixImporter() const {return(Importer());};
	
  //@}
	
  //@{ \name Inlined Operator Methods.
	
	//! Inlined bracket operator for fast access to data. (Const and Non-const versions)
	/*! No error checking and dangerous for optimization purposes.
    \param In
		Loc - Local row.
	  
    \return reference to pointer to locally indexed Loc row in matrix.
  */
	inline double*& operator[] (int Loc) {return Values_[Loc];}
	inline double* const & operator[] (int Loc) const {return Values_[Loc];}
  //@}
	
  //@{ \name Deprecated methods:  These methods still work, but will be removed in a future version.
	
	//! Use ColMap() instead. 
	const Kokkos_Map& ImportMap() const {return((Kokkos_Map&) Graph_->ImportMap());};

	//! Use FillComplete() instead.
	int TransformToLocal();

	//! Use FillComplete(const Kokkos_Map& DomainMap, const Kokkos_Map& RangeMap) instead.
	int TransformToLocal(const Kokkos_Map* DomainMap, const Kokkos_Map* RangeMap);

  //@}
	
	
 protected:
	bool Allocated() const {return(Allocated_);};
	int SetAllocated(bool Flag) {Allocated_ = Flag; return(0);};
	double** Values() const {return(Values_);};
	
  void InitializeDefaults();
  int Allocate();

  int InsertValues(int LocalRow, int NumEntries, double* Values, int* Indices);
  void SetStaticGraph(bool Flag) {StaticGraph_ = Flag;};
  int CheckSizes(const Kokkos_SrcDistObject& A);
  int CopyAndPermute(const Kokkos_SrcDistObject& Source,
										 int NumSameIDs, 
										 int NumPermuteIDs, int* PermuteToLIDs,
										 int* PermuteFromLIDs);
  int CopyAndPermuteCrsMatrix(const Kokkos_CrsMatrix& A,
															int NumSameIDs, 
															int NumPermuteIDs, int* PermuteToLIDs,
															int* PermuteFromLIDs);
  int CopyAndPermuteRowMatrix(const Kokkos_RowMatrix& A,
															int NumSameIDs, 
															int NumPermuteIDs, int* PermuteToLIDs,
															int* PermuteFromLIDs);
  
  int PackAndPrepare(const Kokkos_SrcDistObject& Source, int NumExportIDs, int* ExportLIDs,
										 int Nsend, int Nrecv,
										 int& LenExports, char*& Exports, int& LenImports, 
										 char*& Imports, 
										 int& SizeOfPacket, Kokkos_Distributor& Distor);
  
  int UnpackAndCombine(const Kokkos_SrcDistObject& Source, 
											 int NumImportIDs, int* ImportLIDs, 
											 char* Imports, int& SizeOfPacket, 
											 Kokkos_Distributor& Distor, Kokkos_CombineMode CombineMode);
  
  Kokkos_CrsGraph* Graph_;
  bool Allocated_;
  bool StaticGraph_;
  bool UseTranspose_;
  
  double** Values_;
  double* All_Values_;
  mutable double NormInf_;
  mutable double NormOne_;

  int NumMyRows_;
  int* NumEntriesPerRow_;
  int* NumAllocatedEntriesPerRow_;
  int** Indices_;
  mutable Kokkos_MultiVector* ImportVector_;
  mutable Kokkos_MultiVector* ExportVector_;

  Kokkos_DataAccess CV_;

};
#endif /* KOKKOS_CRSMATRIX_H */
