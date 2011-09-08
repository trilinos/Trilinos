//@HEADER
// ***********************************************************************
// 
//                Komplex: Complex Linear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#ifndef KOMPLEX_ROWMATRIX_H
#define KOMPLEX_ROWMATRIX_H

#include "Epetra_RowMatrix.h"

/*! Komplex_RowMatrix: A class for the construction and use of complex-valued matrices
    stored in the equivalent real Komplex format.  This class implements the pure virtual
    Epetra_RowMatrix class.

    <p>
    <b>Constructing Komplex_RowMatrix objects</b>

    Constructing the Komplex_RowMatrix:
    The user defines the complex-valued matrix C = (c0r+i*c0i)*A0 +(c1r+i*c1i)*A1.
    Using this general expression for the complex matrix allows easy formulation of a variety of common
    complex problems.  A0 and A1 are stored in Epetra_VbrMatrix or Epetra_CrsMatrix format.

    The different K forms (K1, K2, K3, K4, K14, and K23) can easily convert back and forth by
    going from one K form to the canonical form to another K form.  The Komplex_Ordering that each 
    Komplex_RowMatrix object has is what determines the conversions.  Let Kanon stand for the canonical
    form of a complex matrix in equivalent real formulation.  Then any K form is equivalent to:
    P_l * Kanon * P_r * D_r,
    where P_l, P_r are specific permutation matrices and D_r is a specific right diagonal scaling matrix.

    This is helpful because certain forms are advantageous in certain conditions.  To be able to convert
    back and forth during preconditioning and solving should allow for faster, more accurate solutions.
*/    

class Komplex_RowMatrix : public virtual Epetra_RowMatrix {

 public:

  //@{ \name Constructors/Destructor.
  //! Default Komplex_RowMatrix constuctor with fixed number of indices per row.
  /*! Creates a Komplex_RowMatrix object and allocates storage. 
    \param DataAccess (In) Copy or view. 
    \param RowMap (In) A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param NumEntriesPerRow (In) Integer giving the approximate number of entries per row.
    \param KForm (In) The Komplex_KForms to use for this RowMatrix; by default, it is set to K1.

    \warning Note that, because Epetra_LocalMap derives from
    Epetra_Map, and Epetra_Map derives from Epetra_BlockMap,
    this constructor works for all three types of Epetra map classes.
  */
  Komplex_RowMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int NumEntriesPerRow, Komplex_KForms KForm = K1);

  //! Default Komplex_RowMatrix constructor with variable number of indices per row.
  /*! Creates a Komplex_RowMatrix object and allocates storage.
    \param DataAccess (In) Copy or view. 
    \param RowMap (In) A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param NumEntriesPerRow (In) Integer array giving the approximate number of entries per row.
    \param KForm (In) The Komplex_KForms to use for this RowMatrix; by default, it is set to K1.

    \warning Note that, because Epetra_LocalMap derives from
    Epetra_Map, and Epetra_Map derives from Epetra_BlockMap,
    this constructor works for all three types of Epetra map classes.
  */
  Komplex_RowMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int* NumEntriesPerRow, Komplex_KForms KForm = K1);

  //! Constructor with fixed number of indices per row and both row and column maps.
  /*! Creates a Komplex_RowMatrix object and allocates storage.
    \param DataAccess (In) Copy or view.
    \param RowMap (In) A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap describing the layout of the rows.
    \param ColMap (In) A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap describing the layout of the columns.
    \param NumEntriesPerRow (In) Integer giving the approximate number of entries per row.
    \param KForm (In) The Komplex_KForms to use for this RowMatrix; by default, it is set to K1.

    \warning Note that, because Epetra_LocalMap derives from
    Epetra_Map, and Epetra_Map derives from Epetra_BlockMap,
    this constructor works for all three types of Epetra map classes.
  */
  Komplex_RowMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, const Epetra_BlockMap& ColMap, 
                    int NumEntriesPerRow, Komplex_KForms KForm = K1);

  //! Constructor with variable number of indices per row and both row and column maps.
  /*! Creates a Komplex_RowMatrix object and allocates storage.
    \param DataAccess (In) Copy or view.
    \param RowMap (In) A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap describing the layout of the rows.
    \param ColMap (In) A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap describing the layout of the columns.
    \param NumEntriesPerRow (In) Integer array giving the approximate number of entries per row.
    \param KForm (In) The Komplex_KForms to use for this RowMatrix; by default, it is set to K1.

    \warning Note that, because Epetra_LocalMap derives from
    Epetra_Map, and Epetra_Map derives from Epetra_BlockMap,
    this constructor works for all three types of Epetra map classes.
  */
  Komplex_RowMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, const Epetra_BlockMap& ColMap, 
			  int* NumEntriesPerRow, Komplex_KForms KForm = K1);

  //! Constructor that uses an existing Epetra_CrsGraph object.
  /*! Creates a Komplex_RowMatrix object and allocates storage.
    \param DataAccess (In) Copy or view.
    \param Graph (In) A Epetra_CrsGraph.
    \param KForm (In) The Komplex_KForms to use for this RowMatrix; by default, it is set to K1.
  */
  Komplex_RowMatrix(Epetra_DataAccess CV, const Epetra_CrsGraph &Graph, Komplex_KForms KForm = K1);

  //! General Komplex_RowMatrix constructor taking one Epetra_RowMatrix object.
  /*! Creates a Komplex_RowMatrix object and fills it.
    \param cr (In) The real part of the complex coefficient multiplying A.
    \param ci (In) The imag part of the complex coefficient multiplying A.
    \param A (In) An Epetra_RowMatrix that is used to define the true complex matrix, real and imaginary values are interleaved.
    \param KForm (In) The Komplex_KForms to use for this RowMatrix; by default, it is set to K1.
  */
  Komplex_RowMatrix(double cr, double ci, Epetra_RowMatrix* A, Komplex_KForms KForm = K1);

  //! General Komplex_RowMatrix constructor.
  /*! Creates a Komplex_RowMatrix object and fills it.
    \param c0r (In) The real part of the complex coefficient multiplying A0.
    \param c0i (In) The imag part of the complex coefficient multiplying A0.
    \param A0 (In) An Epetra_RowMatrix that is one of the matrices used to define the true complex matrix.
    \param c1r (In) The real part of the complex coefficient multiplying A1.
    \param c1i (In) The imag part of the complex coefficient multiplying A1.
    \param A1 (In) An Epetra_RowMatrix that is the second of the matrices used to define the true complex matrix.
    \param KForm (In) The Komplex_KForms to use for this operator; by default, it is set to K1.
  */
  Komplex_RowMatrix(double c0r, double c0i, Epetra_RowMatrix* A0, double c1r, double c1i, 
		 	  Epetra_RowMatrix* A1, Komplex_KForms KForm = K1);

  //! 
  //! Copy constructor.
  Komplex_RowMatrix(const Komplex_RowMatrix& Matrix);

  //! Komplex_RowMatrix Destructor
  virtual ~Komplex_RowMatrix();
  //@}
  
  //@{ \name Equivalence and querying methods.

  Komplex_RowMatrix & operator=(const Komplex_RowMatrix& src);

  //! If this matrix has been filled, this query returns true; otherwise it returns false.
  bool Filled() const;

  //! If OptimizeStorage() has been called, this query returns true; otherwise it returns false.
  bool StorageOptimized() const;
	
  //! If matrix indices have not been transformed to local, this query returns true; otherwise it returns false.
  bool IndicesAreGlobal() const;
	
  //! If matrix indices have been transformed to local, this query returns true; otherwise it returns false.
  bool IndicesAreLocal() const;
	
  //! If matrix indices are packed into single array (done in OptimizeStorage()) return true, otherwise false.
  bool IndicesAreContiguous() const;

  //! If matrix has no diagonal entries in global index space, this query returns true; otherwise it returns false.
  bool NoDiagonal() const;

  //! Returns the index base for row and column indices for this graph.
  int IndexBase() const;
  
  //! Returns a reference to the Epetra_CrsGraph object associated with this matrix.
  const Epetra_CrsGraph & Graph() const;
	
  //! Returns the Epetra_Map object associated with the rows of this matrix.
  const Epetra_Map & RowMap() const;

  //! Returns the Epetra_Map object that describes the set of column-indices that appear in each processor's locally owned matrix rows.
  /*! Note that if the matrix was constructed with only a row-map, then until FillComplete() is called, this method returns
      a column-map that is a copy of the row-map. That 'initial' column-map is replaced with a computed column-map (that
      contains the set of column-indices appearing in each processor's local portion of the matrix) when FillComplete() is
      called.
    \pre HaveColMap()==true
  */
  const Epetra_Map & ColMap() const;
	
  //! Returns the Epetra_Map object associated with the domain of this matrix operator.
  /*!
      \pre Filled()==true
  */
  const Epetra_Map & DomainMap() const;
	
  //! Returns the Epetra_Map object associated with the range of this matrix operator.
  /*!
      \pre Filled()==true
  */
  const Epetra_Map & RangeMap() const;
	
  //! Returns the Epetra_Import object that contains the import operations for distributed operations.
  const Epetra_Import * Importer() const;
	
  //! Returns the Epetra_Export object that contains the export operations for distributed operations.
  const Epetra_Export * Exporter() const;
	
  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm & Comm() const;

  //@}
 
  //@{ \name Insertion, replace and sum into methods. ##NEED THE REPLACE AND SUM METHODS!!!!!!!!!!!!!

  //! Initialize all values in the matrix with constant value.
  /*!
    \param ScalarConstant (In) Value to use.
		
    \return Integer error code, set to 0 if successful.
    \pre None.
    \post All values in \e this set to ScalarConstant.
  */
  int PutScalar(double ScalarConstant);
	
  //! Multiply all values in the matrix by a constant value (in place: A <- ScalarConstant * A).
  /*!
    \param ScalarConstant (In) Value to use.
		
    \return Integer error code, set to 0 if successful.
    \pre None.
    \post All values of \e this have been multiplied by ScalarConstant.
  */
  int Scale(double ScalarConstant);

  //! Replaces diagonal values of the matrix with those in the user-provided vector.
  /*! This routine is meant to allow replacement of {\bf existing} diagonal values.
      If a diagonal value does not exist for a given row, the corresponding value in
      the input Epetra_Vector will be ignored and the return code will be set to 1.
		
      The Epetra_Map associated with the input Epetra_Vector must be compatible with
      the RowMap of the matrix.
		
      \param Diagonal (In) New values to be placed in the main diagonal.
		
      \return Integer error code, set to 0 if successful, 1 if one or more diagonal entries not present in matrix.
      \pre Filled()==true
      \post Diagonal values have been replaced with the values of Diagonal.
  */
  int ReplaceDiagonalValues(const Epetra_Vector& Diagonal);

  //@}

  //@{ \name Transformation methods
  
  //! Signal that data entry is complete.  Perform transformations to local index space.
  /*  This version of FillComplete assumes that the domain and range distributions are
      identical to the matrix row distributions.
  */
  int FillComplete();

  //! Signal that data entry is complete.  Perform transformations to local index space.
  /*  This version of FillComplete requires the explicit specification of the domain
      and range distribution maps.  These maps are used for importing and exporting vector
      and multi-vector elements that are needed for distributed matrix computations.  For
      example, to compute y = Ax in parallel, we would specify the DomainMap as the distribution
      of the vector x and the RangeMap as the distribution of the vector y.

    \param DomainMap (In) Map that describes the distribution of vector and multivectors in the matrix domain.
    \param RangeMap (In) Map that describes the distribution of vector and multivectors in the matrix range.

    \return Error code, set to 0 if successful. Positive warning code of 2 if it is detected that the
            matrix-graph got out of sync since this matrix was constructed (for instance, if
            graph.FillComplete() was called by another matrix that shares the graph).

    \post IndicesAreLocal()==true
  */
  int FillComplete(const Epetra_Map& DomainMap, const Epetra_Map& RangeMap);

  //@}

  //@{ \name Computational methods.

  //! Returns the result of a Komplex_RowMatrix multiplied by a Epetra_MultiVector X in Y.
  /*! 
    \param TransA (In) If true, multiply by the transpose of the matrix, otherwise just use the matrix.
    \param X (In) A Epetra_MultiVector of dimension NumVectors to multiply with the matrix.
    \param Y (Out) A Epetra_MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.
  */
  int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of a solve using the Komplex_RowMatrix on a Epetra_Vector x in y.
  /*! 
    \param Upper (In) If true, solve Ux = y, otherwise solve Lx = y.
    \param Trans (In) If true, solve the transpose problem.
    \param UnitDiagonal (In) If true, assume diagonal is unit (whether it's stored or not).
    \param x (In) A Epetra_Vector to solve for.
    \param y (Out) A Epetra_Vector containing the result.

    \return Integer error code, set to 0 if successful.
  */
  int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_Vector& x, Epetra_Vector& y) const;

  //! Returns result of a local-only solve using a triangular Komplex_RowMatrix with Epetra_MultiVectors X and Y.
  /*! This method will perform a triangular solve independently on each processor of the parallel machine.
      No communication is performed.
  \param Upper (In) If true, solve Ux = y, otherwise solve Lx = y.
  \param Trans (In) If true, solve transpose problem.
  \param UnitDiagonal (In) If true, assume diagonal is unit (whether it's stored or not).
  \param X (In) A Epetra_MultiVector of dimension NumVectors to solve for.
  \param Y (Out) A Epetra_MultiVector of dimension NumVectors containing result.

  \return Integer error code, set to 0 if successful.
  */
  int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Computes the sum of absolute values of the rows of the Komplex_RowMatrix, results returned in x.
  /*! The vector x will return such that x[i] will contain the inverse of sum of the absolute values of the 
      \e this matrix and will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the global row number of A
      and j denotes the global column number of A.  Using the resulting vector from this function as input to LeftScale()
	will make the infinity norm of the resulting matrix exactly 1.
    \param x (Out) A Epetra_Vector containing the row sums of the \e this matrix.
 
    \return Integer error code, set to 0 if successful.
    \warning It is assumed that the distribution of x is the same as the rows of \e this.
  */
  int InvRowSums(Epetra_Vector& x) const;

  //! Scales the Komplex_RowMatrix on the left with a Epetra_Vector x.
  /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
      and j denotes the column number of A.
    \param x (In) A Epetra_Vector to use to scale.

    \return Integer error code, set to 0 if successful.
  */
  int LeftScale(const Epetra_Vector& x);

  //! Computes the sum of absolute values of the columns of the Komplex_RowMatrix, results returned in x.
  /*! The vector x will return such that x[j] will contain the inverse of sum of the absolute values of the 
      \e this matrix and will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
      and j denotes the global column number of A.  Using the resulting vector from this function as input to 
	RightScale() will make the one norm of the resulting matrix exactly 1.
    \param x (Out) A Epetra_Vector containing the column sums of the \e this matrix. 

    \warning It is assumed that the distribution of x is the same as the rows of \e this.
    \return Integer error code, set to 0 if successful.
  */
  int InvColSums(Epetra_Vector& x) const;

  //! Scales the Komplex_RowMatrix on the right with a Epetra_Vector x.
  /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
      and j denotes the global column number of A.
    \param x (In) The Epetra_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
  */
  int RightScale(const Epetra_Vector& x);
  //@}

  //@{ \name Matrix Properties Query Methods.

  //! If matrix is lower triangular in local index space, this query returns true; otherwise it returns false.
  bool LowerTriangular() const;

  //! If matrix is upper triangular in local index space, this query returns true; otherwise it returns false.
  bool UpperTriangular() const;
  //@}
  
  //@{ \name Atribute access functions

  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].
  */ 
  double NormInf() const;

  //! Returns the one norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_1\f$ such that
     \f[\| A \|_1 = \max_{1\lej\len} \sum_{i=1}^m |a_{ij}| \f].
  */ 
  double NormOne() const;

  //! Returns the number of nonzero entries owned by the calling processor.
  int NumMyNonzeros() const;

  //! Returns the number of matrix rows owned by the calling processor.
  int NumMyRows() const;

  //! Returns the number of matrix columns owned by the calling processor.
  int NumMyCols() const;

  //! Returns the number of global matrix rows.
  int NumGlobalRows() const;

  //! Returns the number of global matrix columns.
  int NumGlobalCols() const;

  //! Returns the number of nonzero entries in the global matrix.
  int NumGlobalNonzeros() const;

  //! Returns the number of local nonzero diagonal entries, based on global row/column index comparisons.
  int NumMyDiagonals() const;
    
  //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
  int NumGlobalDiagonals() const;
  //@}
    
  //@{ \name Additional methods required to implement the Epetra_RowMatrix interface.
	
  //! Returns a copy of the specified global row in user-provided arrays.
  /*! 
    \param GlobalRow (In) Global row to extract.
    \param ILength (In) Length of Values and Indices.
    \param NumEntries (Out) Number of nonzero entries extracted.
    \param Values (Out) Extracted values for this row.
    \param Indices (Out) Extracted global column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful, non-zero if global row is not owned by calling process
     or if the number of entries in this row exceed the Length parameter.
  */
  int ExtractGlobalRowCopy(int GlobalRow, int Length, int& NumEntries, double* Values, int* Indices) const;

  //! Returns a copy of the specified global row values in user-provided array.
  /*! 
    \param GlobalRow (In) Global row to extract.
    \param Length (In) Length of Values.
    \param NumEntries (Out) Number of nonzero entries extracted.
    \param Values (Out) Extracted values for this row.
	  
    \return Integer error code, set to 0 if successful.
  */
  int ExtractGlobalRowCopy(int GlobalRow, int Length, int& NumEntries, double* Values) const;

  //! Returns a copy of the specified local row in user-provided arrays.
  /*! 
  \param MyRow (In) Local row to extract.
  \param Length (In) Length of Values and Indices.
  \param NumEntries (Out) Number of nonzero entries extracted.
  \param Values (Out) Extracted values for this row.
  \param Indices (Out) Extracted local column indices for the corresponding values. 

  \return Integer error code, set to 0 if successful.
  */
  int ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, double* Values, int* Indices) const;

  //! Return the current number of values stored for the specified local row.
  /*! 
  \param MyRow (In) Local row.
  \param NumEntries (Out) Number of nonzero values.

  \return Integer error code, set to 0 if successful.
  */
  int NumMyRowEntries(int MyRow, int& NumEntries) const;

  //! Returns a copy of the main diagonal in a user-provided vector.
  /*! 
    \param Diagonal (Out) Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
  */
  int ExtractDiagonalCopy(Epetra_Vector& Diagonal) const;

  //! Returns the maximum of NumMyRowEntries() over all rows.
  int MaxNumEntries() const;

  //! Returns the Epetra_Map object associated with the rows of this matrix.
  const Epetra_Map & RowMatrixRowMap() const; 

  //! Returns the Epetra_Map object associated with columns of this matrix.
  const Epetra_Map & RowMatrixColMap() const;

  //! Returns the Epetra_Import object that contains the import operations for distributed operations.
  const Epetra_Import * RowMatrixImporter() const;
  //@}

  //@{ \name I/O Methods.

  //! Print method
  void Print(ostream& os) const;
  //@}

  //@{ \name Additional methods required to support the Epetra_Operator interface.

  //! Returns a character string describing the operator
  const char * Label() const;
 	
  //! If set true, transpose of this operator will be applied.
  /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
	affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	does not support transpose use, this method should return a value of -1.
    \param UseTranspose (In) If true, multiply by the transpose of operator, otherwise just use operator.
		
    \return Always returns 0 unless the implementation of this interface does not support transpose use.
  */
  int SetUseTranspose(bool UseTranspose);

  //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
  /*! 
    \param X (In) An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Y (Out) An Epetra_MultiVector of dimension NumVectors containing result.
	
    \return Integer error code, set to 0 if successful.
    \pre Filled()==true
    \post Unchanged.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
	
  //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
  /*! In this implementation, we use several existing attributes to determine how virtual
      method ApplyInverse() should call the concrete method Solve().  We pass in the UpperTriangular(), 
	the Epetra_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
	if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
	be accounted for when doing a triangular solve.
	  
    \param X (In) An Epetra_MultiVector of dimension NumVectors to solve for.
    \param Y (Out) An Epetra_MultiVector of dimension NumVectors containing result.
		
    \return Integer error code, set to 0 if successful.
    \pre Filled()==true
    \post Unchanged.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns true because this class can compute an Inf-norm.
  bool HasNormInf() const;
	
  //! Returns the current UseTranspose setting.
  bool UseTranspose() const;
	
  //! Returns the Epetra_Map object associated with the domain of this matrix operator.
  const Epetra_Map & OperatorDomainMap() const;
	
  //! Returns the Epetra_Map object associated with the range of this matrix operator.
  const Epetra_Map & OperatorRangeMap() const;

  //@}

  protected:

  private:
  Komplex_Ordering* Ordering_;
  Epetra_RowMatrix* Real_;
  Epetra_RowMatrix* Imag_;
  bool IsOneObject_;
  };

#endif /* KOMPLEX_ROWMATRIX_H */