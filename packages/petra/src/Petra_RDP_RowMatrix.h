#ifndef _PETRA_RDP_ROWMATRIX_H_
#define _PETRA_RDP_ROWMATRIX_H_
//! Petra_RDP_RowMatrix: A pure virtual class for using real-valued double-precision row matrices.

/*! The Petra_RDP_RowMatrix class is a pure virtual class (specifies interface only) that 
    enable the use of real-valued double-precision sparse matrices
    where matrix entries are intended for row access.  It is currently implemented by both the
    Petra_RDP_CRS_Matrix and Petra_RDP_VBR_Matrix classes.

   
*/    

#include "Petra_Import.h"
#include "Petra_Export.h"
#include "Petra_RDP_Vector.h"
#include "Petra_RDP_MultiVector.h"


class Petra_RDP_RowMatrix {
      
 public:

  //! Destructor
  virtual ~Petra_RDP_RowMatrix() {};
  //! If FillComplete() has been called, this query returns true, otherwise it returns false.
  virtual bool Filled() const = 0;
  
  // Matrix data extraction routines
  
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
    Indices - Extracted global column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful.
  */
  virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const = 0;

  //! Returns a copy of the main diagonal in a user-provided vector.
  /*! 
    \param Out
    Diagonal - Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
  */
  virtual int ExtractDiagonalCopy(Petra_RDP_Vector & Diagonal) const = 0;

  // Mathematical functions.

  //! Returns the result of a Petra_RDP_RowMatrix multiplied by a Petra_RDP_MultiVector X in Y.
  /*! 
    \param In
    TransA -If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param In
    X - A Petra_RDP_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
    Y -A Petra_RDP_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Multiply(bool TransA, const Petra_RDP_MultiVector& X, Petra_RDP_MultiVector& Y) const = 0;

  //! Returns the result of a Petra_RDP_RowMatrix multiplied by a Petra_RDP_MultiVector X in Y.
  /*! 
    \param In
    Upper -If true, solve Ux = y, otherwise solve Lx = y.
    \param In
    Trans -If true, solve transpose problem.
    \param In
    UnitDiagonal -If true, assume diagonal is unit (whether it's stored or not).
    \param In
    X - A Petra_RDP_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y -A Petra_RDP_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Petra_RDP_MultiVector& X, 
		    Petra_RDP_MultiVector& Y) const = 0;

  //! Computes the sum of absolute values of the rows of the Petra_RDP_RowMatrix, results returned in x.
  /*! The vector x will return such that x[i] will contain the inverse of sum of the absolute values of the 
    \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the global row number of A
    and j denotes the global column number of A.  Using the resulting vector from this function as input to LeftScale()
    will make the infinity norm of the resulting matrix exactly 1.
    \param Out
    x -A Petra_RDP_Vector containing the row sums of the \e this matrix. 
    \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
  */
  virtual int InvRowSums(Petra_RDP_Vector& x) const = 0;

  //! Scales the Petra_RDP_RowMatrix on the left with a Petra_RDP_Vector x.
  /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
    and j denotes the column number of A.
    \param In
    x -A Petra_RDP_Vector to solve for.

    \return Integer error code, set to 0 if successful.
  */
  virtual int LeftScale(const Petra_RDP_Vector& x) = 0;

  //! Computes the sum of absolute values of the columns of the Petra_RDP_RowMatrix, results returned in x.
  /*! The vector x will return such that x[j] will contain the inverse of sum of the absolute values of the 
    \e this matrix will be sca such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
    and j denotes the global column number of A.  Using the resulting vector from this function as input to 
    RighttScale() will make the one norm of the resulting matrix exactly 1.
    \param Out
    x -A Petra_RDP_Vector containing the column sums of the \e this matrix. 
    \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
  */
  virtual int InvColSums(Petra_RDP_Vector& x) const = 0;

  //! Scales the Petra_RDP_RowMatrix on the right with a Petra_RDP_Vector x.
  /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
    and j denotes the global column number of A.
    \param In
    x -The Petra_RDP_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
  */
  virtual int RightScale(const Petra_RDP_Vector& x) = 0;

  // Atribute access functions

  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\len} \sum_{i=1}^m |a_{ij}| \f].
  */ 
  virtual double NormInf() const = 0;

  //! Returns the one norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_1\f$ such that
     \f[\| A \|_1= \max_{1\lej\len} \sum_{j=1}^n |a_{ij}| \f].
  */ 
  virtual double NormOne() const = 0;

  //! Returns the number of nonzero entries in the global matrix.
  virtual int NumGlobalNonzeros() const = 0;

  //! Returns the number of global matrix rows.
  virtual int NumGlobalRows() const = 0;

  //! Returns the number of global matrix columns.
  virtual int NumGlobalCols() const= 0;

  //! Returns the number of global nonzero diagonal entries.
  virtual int NumGlobalDiagonals() const = 0;
    
  //! Returns the number of nonzero entries in the calling processor's portion of the matrix.
  virtual int NumMyNonzeros() const = 0;

  //! Returns the number of matrix rows owned by the calling processor.
  virtual int NumMyRows() const = 0;

  //! Returns the number of matrix columns owned by the calling processor.
  virtual int NumMyCols() const = 0;

  //! Returns the number of local nonzero diagonal entries.
  virtual int NumMyDiagonals() const = 0;

  //! If matrix is lower triangular, this query returns true, otherwise it returns false.
  virtual bool LowerTriangular() const = 0;

  //! If matrix is upper triangular, this query returns true, otherwise it returns false.
  virtual bool UpperTriangular() const = 0;

  //! Returns a pointer to the Petra_Comm communicator associated with this matrix.
  virtual const Petra_Comm & Comm() const = 0;

  //! Returns the Petra_Map object associated with the rows of this matrix.
  virtual const Petra_BlockMap & BlockRowMap() const = 0;

  //! Returns the Petra_Map object that describes the import vector for distributed operations.
  virtual const Petra_BlockMap & BlockImportMap() const = 0;

  //! Returns the Petra_Import object that contains the import operations for distributed operations.
  virtual const Petra_Import * Importer() const = 0;

};

#endif /* _PETRA_RDP_ROWMATRIX_H_ */
