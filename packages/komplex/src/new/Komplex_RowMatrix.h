//@HEADER
/*
************************************************************************

              Komplex: Complex Linear Solver Package 
                Copyright (2002) Sandia Corporation

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
/*! Komplex_RowMatrix: A class for the construction and use of complex-valued matrices
stored in the equivalent real Komplex format.  This class implements the pure virtual
Epetra_RowMatrix class.

<p>
<b>Constructing Komplex_RowMatrix objects</b>

Constructing the Komplex_RowMatrix:
      The user defines the complex-valued matrix C = (c0r+i*c0i)*A0 +(c1r+i*c1i)*A1.
      Using this general expression for the complex matrix allows easy formulation of a variety of common
      complex problems.  A0 and A1 are stored in Epetra_VbrMatrix format.

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
  //! Default Komplex_RowMatrix constuctor.
  /*! Creates a Komplex_RowMatrix object and allocates storage.  
  */
  Komplex_RowMatrix(void);
  
  //! General Komplex_RowMatrix constructor.
  /*! Creates a Komplex_RowMatrix object and fills it.
 	\param In 
		 c0r - The real part of the complex coefficient multiplying A0.
      \param In
		 c0i - The imag part of the complex coefficient multiplying A0.
      \param In 
		 A0 - An Epetra_RowMatrix that is one of the matrices used to define the true complex matrix.
      \param In
		 c1r - The real part of the complex coefficient multiplying A1.
      \param In
		 c1i - The imag part of the complex coefficient multiplying A1.
      \param In
		 A1 - An Epetra_RowMatrix that is the second of the matrices used to define the true complex matrix.
  */
  Komplex_RowMatrix(double c0r, double c0i, Epetra_VbrMatrix A0, double c1r, double c1i, Epetra_VbrMatrix A1);

  //! Copy constructor.
  Komplex_RowMatrix(const Komplex_RowMatrix & Matrix);

  //! Komplex_RowMatrix Destructor
  virtual ~Komplex_RowMatrix();
  //@}
  
  //@{ \name Equivalence and querying methods.

  Komplex_RowMatrix& operator=(const Komplex_RowMatrix & src);

  //! If this matrix has been filled, this query returns true; otherwise it returns false.
  bool Filled() const;
  //@}

  //@{ \name Computational methods.

  //! Returns the result of a Komplex_RowMatrix multiplied by a Epetra_MultiVector X in Y.
  /*! 
  \param In
 	   TransA - If true, multiply by the transpose of the matrix, otherwise just use the matrix.
  \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with the matrix.
  \param Out
	   Y - A Epetra_MultiVector of dimension NumVectors containing the result.
  \return Integer error code, set to 0 if successful.
  */
  int Multiply(bool TransA, const Epetra_MultiVector & X, Epetra_MultiVector & Y) const;

  //! Returns the result of a solve using the Komplex_RowMatrix on a Epetra_Vector x in y.
  /*! 
  \param In
	   Upper - If true, solve Ux = y, otherwise solve Lx = y.
  \param In
	   Trans - If true, solve the transpose problem.
  \param In
	   UnitDiagonal - If true, assume diagonal is unit (whether it's stored or not).
  \param In
	   x - A Epetra_Vector to solve for.
  \param Out
	   y - A Epetra_Vector containing the  result.
  \return Integer error code, set to 0 if successful.
  */
  int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_Vector& x, Epetra_Vector& y) const;

  //! Returns result of a local-only solve using a triangular Komplex_RowMatrix with Epetra_MultiVectors X and Y.
  /*! This method will perform a triangular solve independently on each processor of the parallel machine.
      No communication is performed.

  \param In
	   Upper - If true, solve Ux = y, otherwise solve Lx = y.
  \param In
	   Trans - If true, solve transpose problem.
  \param In
	   UnitDiagonal - If true, assume diagonal is unit (whether it's stored or not).
  \param In
	   X - A Epetra_MultiVector of dimension NumVectors to solve for.
  \param Out
	   Y - A Epetra_MultiVector of dimension NumVectors containing result.
  \return Integer error code, set to 0 if successful.
  */
  int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector & X, Epetra_MultiVector & Y) const;

  //! Computes the sum of absolute values of the rows of the Komplex_RowMatrix, results returned in x.
  /*! The vector x will return such that x[i] will contain the inverse of sum of the absolute values of the 
      \e this matrix and will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the global row number of A
      and j denotes the global column number of A.  Using the resulting vector from this function as input to LeftScale()
	will make the infinity norm of the resulting matrix exactly 1.
  \param Out
         x - A Epetra_Vector containing the row sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.
  \return Integer error code, set to 0 if successful.
  */
  int InvRowSums(Epetra_Vector & x) const;

  //! Scales the Komplex_RowMatrix on the left with a Epetra_Vector x.
  /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
      and j denotes the column number of A.
  \param In
	   x - A Epetra_Vector to solve for.
  \return Integer error code, set to 0 if successful.
  */
  int LeftScale(const Epetra_Vector & x);

  //! Computes the sum of absolute values of the columns of the Komplex_RowMatrix, results returned in x.
  /*! The vector x will return such that x[j] will contain the inverse of sum of the absolute values of the 
      \e this matrix and will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
      and j denotes the global column number of A.  Using the resulting vector from this function as input to 
	RightScale() will make the one norm of the resulting matrix exactly 1.
  \param Out
	   x - A Epetra_Vector containing the column sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.
  \return Integer error code, set to 0 if successful.
  */
  int InvColSums(Epetra_Vector & x) const ;

  //! Scales the Komplex_RowMatrix on the right with a Epetra_Vector x.
  /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
      and j denotes the global column number of A.
  \param In
	   x - The Epetra_Vector used for scaling \e this.
  \return Integer error code, set to 0 if successful.
  */
  int RightScale(const Epetra_Vector & x);
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
  int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double * Values, int * Indices) const;

  //! Return the current number of values stored for the specified local row.
  /*! 
  \param In
         MyRow - Local row.
  \param Out
    	   NumEntries - Number of nonzero values.
  \return Integer error code, set to 0 if successful.
  */
  int NumMyRowEntries(int MyRow, int & NumEntries) const;

  //! Returns a copy of the main diagonal in a user-provided vector.
  /*! 
  \param Out
	   Diagonal - Extracted main diagonal.
  \return Integer error code, set to 0 if successful.
  */
  int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const;

  //! Returns the maximum of NumMyRowEntries() over all rows.
  int MaxNumEntries() const;

  //! Returns the Epetra_Map object associated with the rows of this matrix.
  const Epetra_Map & RowMatrixRowMap() const; 

  //! Returns the Epetra_Map object associated with columns of this matrix.
  const Epetra_Map & RowMatrixColMap() const;

  //! Returns the Epetra_Import object that contains the import operations for distributed operations.
  const Epetra_Import * RowMatrixImporter() const;
  //@}

  //@{ \name Non-interface Required Methods

  //! Returns the \e this matrix as a Epetra_VbrMatrix object.
  int GenerateVbrMatrix(Epetra_VbrMatrix & Matrix);
  //@}

  protected:

  private:
  Komplex_Ordering Ordering_;
  bool Filled_;
  bool LowerTriangular_;
  bool UpperTriangular_;
  double NormInf_;
  double NormOne_;
  };