/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef ML_ROWMATRIX_H
#define ML_ROWMATRIX_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*! \file ml_RowMatrix.h
 *  \brief Wrapper from ML_Operator to Epetra_RowMatrix
 */

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA)

#include <vector>
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "ml_epetra.h"
#include "Epetra_Operator.h"
class Epetra_MultiVector;
#include "Epetra_RowMatrix.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"

class Epetra_Vector;
class Epetra_Importer;

namespace ML_Epetra {

/*!
 * \class RowMatrix
 *
 * \brief Basic wrapper from ML_Operator to Epetra_RowMatrix.
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \date Last updated on 15-Mar-05.
 */

class RowMatrix : public virtual Epetra_RowMatrix {

 public:
  //@{ \name Constructor.
    //! Constructor, constructs Comm object if not provided
    RowMatrix(ML_Operator* Op, const Epetra_Comm* Comm = 0,
              const bool cheap = false, const USR_COMM =
#ifdef HAVE_MPI
              MPI_COMM_WORLD
#else
              0
#endif
              );

  //@}
  //@{ \name Destructor.
    //! Destructor
    virtual ~RowMatrix();

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
    virtual int NumMyRowEntries(int MyRow, int & NumEntries) const;


    //! Returns the maximum of NumMyRowEntries() over all rows.
    virtual int MaxNumEntries() const;

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
    virtual int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

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
    virtual int Solve(bool /* Upper */, bool /* Trans */, bool /* UnitDiagonal */, const Epetra_MultiVector& /* X */,
		      Epetra_MultiVector& /* Y */) const
    {
      ML_RETURN(-1); // not implemented
    }

    virtual int Apply(const Epetra_MultiVector& X,
		      Epetra_MultiVector& Y) const
    {
      ML_RETURN(Multiply(false,X,Y));
    }

    virtual int ApplyInverse(const Epetra_MultiVector& /* X */,
			     Epetra_MultiVector& /* Y */) const
    {
      ML_RETURN(-1);
    }
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
    virtual int InvRowSums(Epetra_Vector& /* x */) const
    {
      ML_RETURN(-1); // not implemented
    }

    //! Scales the Epetra_RowMatrix on the left with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
        and j denotes the column number of A.
    \param In
	   x -A Epetra_Vector to solve for.

    \return Integer error code, set to 0 if successful.
  */
    virtual int LeftScale(const Epetra_Vector& /* x */)
    {
      ML_RETURN(-1); // not implemented
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
    virtual int InvColSums(Epetra_Vector& /* x */) const
    {
      ML_RETURN(-1); // not implemented
    }


    //! Scales the Epetra_RowMatrix on the right with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.
    \param In
	   x -The Epetra_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
  */
    virtual int RightScale(const Epetra_Vector& /* x */)
    {
      ML_RETURN(-1); // not implemented
    }

  //@}

  //@{ \name Attribute access functions

    //! If FillComplete() has been called, this query returns true, otherwise it returns false.
    virtual bool Filled() const
    {
      return true;
    }

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\len} \sum_{i=1}^m |a_{ij}| \f].
    */
    virtual double NormInf() const;

    //! Returns the one norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_1\f$ such that
       \f[\| A \|_1= \max_{1\lej\len} \sum_{j=1}^n |a_{ij}| \f].
    */
    virtual double NormOne() const
    {
      return(-1.0);
    }

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    //! Returns the number of nonzero entries in the global matrix.
    virtual int NumGlobalNonzeros() const;

    //! Returns the number of global matrix rows.
    virtual int NumGlobalRows() const;

    //! Returns the number of global matrix columns.
    virtual int NumGlobalCols() const;

    //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
    virtual int NumGlobalDiagonals() const;
#endif

     //! Returns the number of nonzero entries in the global matrix.
    virtual long long NumGlobalNonzeros64() const;

    //! Returns the number of global matrix rows.
    virtual long long NumGlobalRows64() const;

    //! Returns the number of global matrix columns.
    virtual long long NumGlobalCols64() const;

    //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
    virtual long long NumGlobalDiagonals64() const;

    //! Returns the number of nonzero entries in the calling processor's portion of the matrix.
    virtual int NumMyNonzeros() const;

    //! Returns the number of matrix rows owned by the calling processor.
    virtual int NumMyRows() const;

    //! Returns the number of matrix columns owned by the calling processor.
    virtual int NumMyCols() const;

    //! Returns the number of local nonzero diagonal entries, based on global row/column index comparisons.
    virtual int NumMyDiagonals() const;

    //! If matrix is lower triangular in local index space, this query returns true, otherwise it returns false.
    virtual bool LowerTriangular() const;

    //! If matrix is upper triangular in local index space, this query returns true, otherwise it returns false.
    virtual bool UpperTriangular() const;

    //! Returns the Epetra_Map object associated with the rows of this matrix.
    virtual const Epetra_Map & RowMatrixRowMap() const;

    //! Returns the Epetra_Map object associated with the columns of this matrix.
    virtual const Epetra_Map & RowMatrixColMap() const;

    //! Returns the Epetra_Import object that contains the import operations for distributed operations.
    virtual const Epetra_Import * RowMatrixImporter() const;
  //@}

  // following functions are required to derive Epetra_RowMatrix objects.

  //! Sets ownership.
  int SetOwnership(bool /* ownership */){return(-1);};

  //! Sets use transpose (not implemented).
  int SetUseTranspose(bool /* UseTransposeFlag */){return(-1);}

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(false);};

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const{return(false);};

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & Comm() const{return(*Comm_);};

  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const {return(*DomainMap_);};

  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const {return(*RangeMap_);};
  //@}

  void SetLabel(const char* label)
  {
    strcpy(Label_,label);
  };

  const char* Label() const{
    return(Label_);
  };

  //!  Returns a reference to RowMatrix->Map().
  const Epetra_BlockMap & Map() const
  {
    return(*DomainMap_);
  }

  //!  Print the global matrix.  This uses a potentially artificial numbering.
  int Print() const;

private:

  //! Pointer to the ML_Operator structure that is wrapped.
  ML_Operator* Op_;
  //! Communicator object, given by the user or allocated here.
  const Epetra_Comm* Comm_;
  //! If \c true, the dtor will destroy the communicator.
  bool FreeCommObject_;
  //! Number of local rows.
  int NumMyRows_;
  //! Number of global rows.
  long long NumGlobalRows_;
  //! Number of local columns.
  int NumMyCols_;
  //! Number of global columns.j
  long long NumGlobalCols_;
  //! Map for row distribution.
  Epetra_Map* DomainMap_;
  //! Map for row distribution.
  Epetra_Map* RangeMap_;
  //! Map for column distribution.
  Epetra_Map* ColMap_;
  //! Maximum number of elements in a row.
  int MaxNumEntries_;
  //! Diagonal elements of the matrix.
  std::vector<double> Diagonal_;
  //! Contains the nonzero elements in a local row.
  std::vector<int> NumMyRowEntries_;
  //! Work vector for getrow().
  mutable int Allocated_;
  //! Work vector for getrow().
  mutable std::vector<int> Indices_;
  //! Work vector for getrow().
  mutable std::vector<double> Values_;
  //! Contains the infinity norm of the matrix.
  double NormInf_;
  //! Number of local nonzeros.
  int NumMyNonzeros_;
  //! Number of global nonzeros.
  long long NumGlobalNonzeros_;
  //! Number of nonzero local diagonal elements.
  int NumMyDiagonals_;
  //! Number of nonzero global diagonal elements.
  long long NumGlobalDiagonals_;
  //! Importer.
  mutable Epetra_Import* Importer_;
  //! Label of \c this object.
  char* Label_;

}; // class RowMatrix

} // namespace ML_Epetra

#endif /* HAVE_ML_EPETRA */
#endif /* ML_ROWMATRIX_H */
