#ifndef AMESOS_ROWMATRIX_H
#define AMESOS_ROWMATRIX_H

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

class Amesos_LocalRowMatrix : public virtual Epetra_RowMatrix {
      
 public:
  //@{ \name Constructor.
    //! Constructor
    Amesos_LocalRowMatrix(Epetra_RowMatrix* Matrix);

  //@}
  //@{ \name Destructor.
    //! Destructor
    virtual ~Amesos_LocalRowMatrix();

  //@}
  
  //@{ \name Matrix data extraction routines

    //! Returns the number of nonzero entries in MyRow.
    /*! 
    \param In
           MyRow - Local row.
    \param Out
	   NumEntries - Number of nonzero values present.
	  
    \return Integer error code, set to 0 if successful.
  */
    virtual int NumMyRowEntries(int MyRow, int & NumEntries) const
    {
      return(NumEntries_[MyRow]);
    }

    //! Returns the maximum of NumMyRowEntries() over all rows.
    virtual int MaxNumEntries() const
    {
      return(MaxNumLocalEntries_);
    }

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
    virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;

    //! Returns a copy of the main diagonal in a user-provided vector.
    /*! 
    \param Out
	   Diagonal - Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
  */
    virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const;
  //@}
  
  //@{ \name Mathematical functions.

    //! Returns the result of a Epetra_RowMatrix multiplied by a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   TransA -If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
  */
    virtual int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
      if (TransA == true) {
	AMESOS_CHK_ERR(-1);
      }

      AMESOS_CHK_ERR(Apply(X,Y));
    }

    //! Returns result of a local-only solve using a triangular Epetra_RowMatrix with Epetra_MultiVectors X and Y.
    /*! This method will perform a triangular solve independently on each processor of the parallel machine.
        No communication is performed.
    \param In
	   Upper -If true, solve Ux = y, otherwise solve Lx = y.
    \param In
	   Trans -If true, solve transpose problem.
    \param In
	   UnitDiagonal -If true, assume diagonal is unit (whether it's stored or not).
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
    virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, 
		      Epetra_MultiVector& Y) const
    {
      AMESOS_RETURN(-1); // not implemented 
    }

    virtual int Apply(const Epetra_MultiVector& X,
		      Epetra_MultiVector& Y) const;

    virtual int ApplyInverse(const Epetra_MultiVector& X,
			     Epetra_MultiVector& Y) const;
    //! Computes the sum of absolute values of the rows of the Epetra_RowMatrix, results returned in x.
    /*! The vector x will return such that x[i] will contain the inverse of sum of the absolute values of the 
        \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.  Using the resulting vector from this function as input to LeftScale()
	will make the infinity norm of the resulting matrix exactly 1.
    \param Out
	   x -A Epetra_Vector containing the row sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
  */
    virtual int InvRowSums(Epetra_Vector& x) const
    {
      AMESOS_RETURN(-1); // not implemented
    }

    //! Scales the Epetra_RowMatrix on the left with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
        and j denotes the column number of A.
    \param In
	   x -A Epetra_Vector to solve for.

    \return Integer error code, set to 0 if successful.
  */
    virtual int LeftScale(const Epetra_Vector& x)
    {
      AMESOS_RETURN(-1); // not implemented
    }

    //! Computes the sum of absolute values of the columns of the Epetra_RowMatrix, results returned in x.
    /*! The vector x will return such that x[j] will contain the inverse of sum of the absolute values of the 
        \e this matrix will be sca such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.  Using the resulting vector from this function as input to 
	RighttScale() will make the one norm of the resulting matrix exactly 1.
    \param Out
	   x -A Epetra_Vector containing the column sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
  */
    virtual int InvColSums(Epetra_Vector& x) const
    {
      AMESOS_RETURN(-1); // not implemented
    }
    

    //! Scales the Epetra_RowMatrix on the right with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.
    \param In
	   x -The Epetra_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
  */
    virtual int RightScale(const Epetra_Vector& x) 
    {
      AMESOS_RETURN(-1); // not implemented
    }

  //@}
  
  //@{ \name Atribute access functions

    //! If FillComplete() has been called, this query returns true, otherwise it returns false.
    virtual bool Filled() const
    {
      return true;
    }

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\len} \sum_{i=1}^m |a_{ij}| \f].
    */ 
    virtual double NormInf() const
    {
      return(-1.0);
    }

    //! Returns the one norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_1\f$ such that
       \f[\| A \|_1= \max_{1\lej\len} \sum_{j=1}^n |a_{ij}| \f].
    */ 
    virtual double NormOne() const
    {
      AMESOS_RETURN(-1.0);
    }

    //! Returns the number of nonzero entries in the global matrix.
    virtual int NumGlobalNonzeros() const
    {
      return(NumNonzeros_);
    }

    //! Returns the number of global matrix rows.
    virtual int NumGlobalRows() const
    {
      return(NumRows_);
    }

    //! Returns the number of global matrix columns.
    virtual int NumGlobalCols() const
    {
      return(NumCols_);
    }

    //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
    virtual int NumGlobalDiagonals() const
    {
      return(NumDiagonals_);
    }
    
    //! Returns the number of nonzero entries in the calling processor's portion of the matrix.
    virtual int NumMyNonzeros() const
    {
      return(NumNonzeros_);
    }

    //! Returns the number of matrix rows owned by the calling processor.
    virtual int NumMyRows() const
    {
      return(NumRows_);
    }

    //! Returns the number of matrix columns owned by the calling processor.
    virtual int NumMyCols() const
    {
      return(NumCols_);
    }

    //! Returns the number of local nonzero diagonal entries, based on global row/column index comparisons.
    virtual int NumMyDiagonals() const
    {
      return(NumDiagonals_);
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
      return(*Map_);
    }

    //! Returns the Epetra_Map object associated with the columns of this matrix.
    virtual const Epetra_Map & RowMatrixColMap() const
    {
      return(*Map_);
    }

    //! Returns the Epetra_Import object that contains the import operations for distributed operations.
    virtual const Epetra_Import * RowMatrixImporter() const
    {
      return(0);
      }
  //@}
  
  // following functions are required to derive Epetra_RowMatrix objects.

  //! Sets ownership.
  int SetOwnership(bool ownership)
  {
    AMESOS_RETURN(-1);
  }

  //! Sets use transpose (not implemented).
  int SetUseTranspose(bool UseTranspose)
  {
    AMESOS_RETURN(-1);
  }

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const 
  {
    AMESOS_RETURN(false);
  }
  
  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const
  {
    return(false);
  }
  
  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & Comm() const
  {
    return(*SerialComm_);
  }
  
  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const 
  {
    return(*Map_);
  }
  
  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const 
  {
    return(*Map_);
  }
  //@}

  const Epetra_BlockMap& Map() const 
  {
    return(*Map_);
  }

  char* Label() const{
    return(Label_);
  };

private:

  Epetra_RowMatrix* Matrix_;
  Epetra_Map* Map_;
  int NumRows_;
  int NumCols_;
  int NumNonzeros_;
  int MaxNumEntries_;
  vector<int> NumEntries_;
  int* Indices_;
  double* Values_;
  int MaxNumLocalEntries_;
  int NumDiagonals_;

  char* Label_;
#ifdef HAVE_MPI
  Epetra_MpiComm* SerialComm_;
#else
  Epetra_SerialComm* SerialComm_;
#endif

};

#endif /* AMESOS_ROWMATRIX_H */

