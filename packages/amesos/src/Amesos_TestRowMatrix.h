#ifndef IFPACK_ROWMATRIX_H
#define IFPACK_ROWMATRIX_H

#include "Amesos_ConfigDefs.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_RowMatrix.h"
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Import;
class Epetra_BlockMap;

//! Amesos_TestRowMatrix: a class to test Epetra_RowMatrix based codes.
/*!
Class Amesos_TestRowMatrix enables the creation of a Epetra_RowMatrix
derived class for testing purposed. This class requires another Epetra_RowMatrix
as input, and minimic the behavior of this matrix. However, as it 
\e this object is not derived from Epetra_CrsMatrix or 
Epetra_VbrMatrix, a dynamic_cast
will not result in any Epetra_CrsMatrix or Epetra_VrbMatrix object. 

\author Marzio Sala, SNL 9214

\date Sep-04
 
 */ 
class Amesos_TestRowMatrix : public virtual Epetra_RowMatrix {

public:
  //@{ \name Constructor.
  //! Constructor
  Amesos_TestRowMatrix(Epetra_RowMatrix* Matrix) :
    Matrix_(Matrix)
  {}

  //@}
  //@{ \name Destructor.
  //! Destructor
  virtual ~Amesos_TestRowMatrix()
  {}

  //@}

  //@{ \name Matrix data extraction routines

  //! Returns the number of nonzero entries in MyRow.
  /*! 
    \param 
    MyRow - (In) Local row.
    \param 
    NumEntries - (Out) Number of nonzero values present.

    \return Integer error code, set to 0 if successful.
    */
  virtual int NumMyRowEntries(int MyRow, int & NumEntries) const
  {
    return(Matrix().NumMyRowEntries(MyRow,NumEntries));
  }

  //! Returns the maximum of NumMyRowEntries() over all rows.
  virtual int MaxNumEntries() const
  {
    return(Matrix().MaxNumEntries());
  }

  //! Returns a copy of the specified local row in user-provided arrays.
  /*! 
    \param
    MyRow - (In) Local row to extract.
    \param
    Length - (In) Length of Values and Indices.
    \param
    NumEntries - (Out) Number of nonzero entries extracted.
    \param
    Values - (Out) Extracted values for this row.
    \param 
    Indices - (Out) Extracted global column indices for the corresponding values.

    \return Integer error code, set to 0 if successful.
    */
  virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const
  {
    return(Matrix().ExtractMyRowCopy(MyRow,Length, NumEntries,
				     Values, Indices));
  }

  //! Returns a copy of the main diagonal in a user-provided vector.
  /*! 
    \param
    Diagonal - (Out) Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
    */
  virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
  {
    return(Matrix().ExtractDiagonalCopy(Diagonal));
  }
  //@}

  //@{ \name Mathematical functions.

  //! Returns the result of a Epetra_RowMatrix multiplied by a Epetra_MultiVector X in Y.
  /*! 
    \param 
    TransA -(In) If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param 
    X - (In) A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param 
    Y -(Out) A Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
    */
  virtual int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
  {
    return(Matrix().Multiply(TransA,X,Y));
  }

  //! Returns result of a local-only solve using a triangular Epetra_RowMatrix with Epetra_MultiVectors X and Y (NOT IMPLEMENTED).
  virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, 
		    Epetra_MultiVector& Y) const
  {
    return(Matrix().Solve(Upper,Trans,UnitDiagonal,X,Y));
  }

  virtual int Apply(const Epetra_MultiVector& X,
		    Epetra_MultiVector& Y) const
  {
    return(Matrix().Apply(X,Y));
  }

  virtual int ApplyInverse(const Epetra_MultiVector& X,
			   Epetra_MultiVector& Y) const
  {
    return(Matrix().ApplyInverse(X,Y));
  }
  //! Computes the sum of absolute values of the rows of the Epetra_RowMatrix, results returned in x (NOT IMPLEMENTED).
  virtual int InvRowSums(Epetra_Vector& x) const
  {
    return(Matrix().InvRowSums(x));
  }

  //! Scales the Epetra_RowMatrix on the left with a Epetra_Vector x (NOT IMPLEMENTED).
  virtual int LeftScale(const Epetra_Vector& x)
  {
    return(Matrix().LeftScale(x));
  }

  //! Computes the sum of absolute values of the columns of the Epetra_RowMatrix, results returned in x (NOT IMPLEMENTED).
  virtual int InvColSums(Epetra_Vector& x) const
  {
    return(Matrix().InvColSums(x));
  }


  //! Scales the Epetra_RowMatrix on the right with a Epetra_Vector x (NOT IMPLEMENTED).
  virtual int RightScale(const Epetra_Vector& x) 
  {
    return(Matrix().RightScale(x));
  }

  //@}

  //@{ \name Atribute access functions

  //! If FillComplete() has been called, this query returns true, otherwise it returns false.
  virtual bool Filled() const
  {
    return(Matrix().Filled());
  }

  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\len} \sum_{i=1}^m |a_{ij}| \f].
     */ 
  virtual double NormInf() const
  {
    return(Matrix().NormInf());
  }

  //! Returns the one norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_1\f$ such that
     \f[\| A \|_1= \max_{1\lej\len} \sum_{j=1}^n |a_{ij}| \f].
     */ 
  virtual double NormOne() const
  {
    return(Matrix().NormOne());
  }

  //! Returns the number of nonzero entries in the global matrix.
  virtual int NumGlobalNonzeros() const
  {
    return(Matrix().NumGlobalNonzeros());
  }

  //! Returns the number of global matrix rows.
  virtual int NumGlobalRows() const
  {
    return(Matrix().NumGlobalRows());
  }

  //! Returns the number of global matrix columns.
  virtual int NumGlobalCols() const
  {
    return(Matrix().NumGlobalCols());
  }

  //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
  virtual int NumGlobalDiagonals() const
  {
    return(Matrix().NumGlobalDiagonals());
  }

  //! Returns the number of nonzero entries in the calling processor's portion of the matrix.
  virtual int NumMyNonzeros() const
  {
    return(Matrix().NumMyNonzeros());
  }

  //! Returns the number of matrix rows owned by the calling processor.
  virtual int NumMyRows() const
  {
    return(Matrix().NumMyRows());
  }

  //! Returns the number of matrix columns owned by the calling processor.
  virtual int NumMyCols() const
  {
    return(Matrix().NumMyCols());
  }

  //! Returns the number of local nonzero diagonal entries, based on global row/column index comparisons.
  virtual int NumMyDiagonals() const
  {
    return(Matrix().NumMyDiagonals());
  }

  //! If matrix is lower triangular in local index space, this query returns true, otherwise it returns false.
  virtual bool LowerTriangular() const
  {
    return(Matrix_->LowerTriangular());
  }

  //! If matrix is upper triangular in local index space, this query returns true, otherwise it returns false.
  virtual bool UpperTriangular() const
  {
    return(Matrix_->UpperTriangular());
  }

  //! Returns the Epetra_Map object associated with the rows of this matrix.
  virtual const Epetra_Map & RowMatrixRowMap() const
  {
    return(Matrix().RowMatrixRowMap());
  }
  //! Returns the Epetra_Map object associated with the columns of this matrix.
  virtual const Epetra_Map & RowMatrixColMap() const
  {
    return(Matrix().RowMatrixColMap());
  }

  //! Returns the Epetra_Import object that contains the import operations for distributed operations.
  virtual const Epetra_Import * RowMatrixImporter() const
  {
    return(Matrix().RowMatrixImporter());
  }
  //@}

  // following functions are required to derive Epetra_RowMatrix objects.

#ifdef FIXME
  //! Sets ownership.
  int SetOwnership(bool ownership)
  {
    return(Matrix().SetOwnership(ownership));
  }
#endif

  //! Sets use transpose (not implemented).
  int SetUseTranspose(bool UseTranspose)
  {
    return(Matrix().SetUseTranspose(UseTranspose));
  }

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const 
  {
    return(Matrix().UseTranspose());
  }

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const
  {
    return(Matrix().HasNormInf());
  }

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & Comm() const
  {
    return(Matrix().Comm());
  }

  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const 
  {
    return(Matrix().OperatorDomainMap());
  }

  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const 
  {
    return(Matrix().OperatorRangeMap());
  }
  //@}

  const Epetra_BlockMap& Map() const
  {
    return(Matrix().Map());
  }

  char* Label() const
  {
    return(Matrix().Label());
  }


private:

  Epetra_RowMatrix& Matrix()
  {
    return(*Matrix_);
  }

  const Epetra_RowMatrix& Matrix() const
  {
    return(*Matrix_);
  }

  //! Pointer to the matrix to be preconditioned.
  Epetra_RowMatrix* Matrix_;
  
};

#endif /* IFPACK_ROWMATRIX_H */

