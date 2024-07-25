// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_NODEFILTER_HPP
#define IFPACK2_NODEFILTER_HPP

#ifdef IFPACK2_NODE_AWARE_CODE

#include "Ifpack2_ConfigDefs.hpp"
#ifdef HAVE_MPI
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialComm.hpp"
#endif
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_IntVector.hpp"
#include "Teuchos_RCP.hpp"
#include "Ifpack2_OverlappingRowMatrix.hpp"
class Tpetra_Map;
class Tpetra_MultiVector;
class Tpetra_Vector;
class Tpetra_Import;
class Tpetra_BlockMap;

//! Ifpack2_NodeFilter a class for light-weight extraction of the submatrix corresponding to local rows and columns.

/*! Class Ifpack2_NodeFilter enables a light-weight contruction of an
 Tpetra_RowMatrix-derived object, containing only the elements of the original, 
 distributed matrix with local row and column ID. The local
 submatrix is based on a communicator containing the local process only. 
 Each process will have its local object, corresponding to the local submatrix.
 Submatrices may or may not overlap.
 
 The following instructions can be used to create "localized" matrices:
 \code
 #include "Ifpack2_NodeFilter.hpp"
 ...
 Teuchos::RCP<Tpetra_RowMatrix> A;             // fill the elements of A,
 A->FillComplete();

 Ifpack2_NodeFilter LocalA(A);
 \endcode

 Once created, \c LocalA defined, on each process, the submatrix 
 corresponding to local rows and columns only. The creation 
 and use of
 \c LocalA is "cheap", as the elements of the local matrix are
 obtained through calls to ExtractMyRowCopy on the original, distributed
 matrix, say A. This means that \c A must remain in scope every time 
 \c LocalA is accessed.

 A very convenient use of this class is to use Ifpack2 solvers to
 compute the LU factorizations of local blocks. If applied to
 a localized matrix, several Ifpack2 objects can operator in the same
 phase in a safe way, without non-required data exchange.

 \author Michael Heroux, SNL 9214

 \date Sep-04
 
 */ 
class Ifpack2_NodeFilter : public virtual Tpetra_RowMatrix {

public:
  //@{ \name Constructor.
  //! Constructor
  Ifpack2_NodeFilter(const Teuchos::RCP<const Tpetra_RowMatrix>& Matrix, int nodeID);
  //Ifpack2_NodeFilter(const Teuchos::RCP<const Tpetra_RowMatrix>& Matrix,const Tpetra_Comm *);

  //@}
  //@{ \name Destructor.
  //! Destructor
  ~Ifpack2_NodeFilter(){
    if(Ac_LIDMap_) delete [] Ac_LIDMap_;
    if(Bc_LIDMap_) delete [] Bc_LIDMap_;
  };
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
    NumEntries = NumEntries_[MyRow];
    return(0);
  }

  //! Returns the maximum of NumMyRowEntries() over all rows.
  virtual int MaxNumEntries() const
  {
    return(MaxNumEntries_);
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
  virtual inline int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;

  //! Returns a copy of the main diagonal in a user-provided vector.
  /*! 
    \param
    Diagonal - (Out) Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
    */
  virtual int ExtractDiagonalCopy(Tpetra_Vector & Diagonal) const;
  //@}

  //@{ \name Mathematical functions.

  //! Returns the result of a Tpetra_RowMatrix multiplied by a Tpetra_MultiVector X in Y.
  /*! 
    \param 
    TransA -(In) If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param 
    X - (In) A Tpetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param 
    Y -(Out) A Tpetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
    */
  virtual int Multiply(bool TransA, const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
  {
    if (TransA == true) {
      IFPACK2_CHK_ERR(-1);
    }

    IFPACK2_CHK_ERR(Apply(X,Y));
    return(0);
  }

  //! Returns result of a local-only solve using a triangular Tpetra_RowMatrix with Tpetra_MultiVectors X and Y (NOT IMPLEMENTED).
  virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
		    Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
  {
    IFPACK2_RETURN(-1); // not implemented 
  }

  virtual int Apply(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
		    Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual int ApplyInverse(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
			   Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;
  //! Computes the sum of absolute values of the rows of the Tpetra_RowMatrix, results returned in x (NOT IMPLEMENTED).
  virtual int InvRowSums(Tpetra_Vector& x) const
  {
    IFPACK2_RETURN(-1); // not implemented
  }

  //! Scales the Tpetra_RowMatrix on the left with a Tpetra_Vector x (NOT IMPLEMENTED).
  virtual int LeftScale(const Tpetra_Vector& x)
  {
    IFPACK2_RETURN(-1); // not implemented
  }

  //! Computes the sum of absolute values of the columns of the Tpetra_RowMatrix, results returned in x (NOT IMPLEMENTED).
  virtual int InvColSums(Tpetra_Vector& x) const
  {
    IFPACK2_RETURN(-1); // not implemented
  }


  //! Scales the Tpetra_RowMatrix on the right with a Tpetra_Vector x (NOT IMPLEMENTED).
  virtual int RightScale(const Tpetra_Vector& x) 
  {
    IFPACK2_RETURN(-1); // not implemented
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
    IFPACK2_RETURN(-1.0);
  }

  //! Returns the number of nonzero entries in the global matrix.
  virtual int NumGlobalNonzeros() const
  {
    return(NumGlobalNonzeros_);
  }

  //! Returns the number of global matrix rows.
  virtual int NumGlobalRows() const
  {
    return(NumGlobalRows_);
  }

  //! Returns the number of global matrix columns.
  virtual int NumGlobalCols() const
  {
    return(NumGlobalRows_);
  }

  //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
  virtual int NumGlobalDiagonals() const
  {
    return(NumGlobalRows_);
  }

  //! Returns the number of nonzero entries in the calling processor's portion of the matrix.
  virtual int NumMyNonzeros() const
  {
    return(NumMyNonzeros_);
  }

  //! Returns the number of matrix rows owned by the calling processor.
  virtual int NumMyRows() const
  {
    return(NumMyRows_);
  }

  //! Returns the number of matrix columns owned by the calling processor.
  virtual int NumMyCols() const
  {
    return(NumMyCols_);
  }

  //! Returns the number of local nonzero diagonal entries, based on global row/column index comparisons.
  virtual int NumMyDiagonals() const
  {
    return(NumMyRows_);
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

  //! Returns the Tpetra_Map object associated with the rows of this matrix.
  virtual const Tpetra_Map & RowMatrixRowMap() const
  {
    return(*Map_);
  }

  //! Returns the Tpetra_Map object associated with the columns of this matrix.
  virtual const Tpetra_Map & RowMatrixColMap() const
  {
    return(*colMap_);
  }

  //! Returns the Tpetra_Import object that contains the import operations for distributed operations.
  virtual const Tpetra_Import * RowMatrixImporter() const
  {
    return(&*Importer_);
  }
  //@}

  virtual const Tpetra_Import* Importer() const {return(&*Importer_);}

  virtual const Tpetra_Export* Exporter() const {return(&*Exporter_);}

  // following functions are required to derive Tpetra_RowMatrix objects.

  //! Sets ownership.
  int SetOwnership(bool ownership)
  {
    IFPACK2_RETURN(-1);
  }

  //! Sets use transpose (not implemented).
  int SetUseTranspose(bool UseTranspose_in)
  {
    UseTranspose_ = UseTranspose_in;
    return(0);
  }

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const 
  {
    return(UseTranspose_);
  }

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const
  {
    return(false);
  }

  //! Returns a pointer to the Tpetra_Comm communicator associated with this operator.
  const Tpetra_Comm & Comm() const
  {
    return(*SubComm_);
  }

  //! Returns the Tpetra_Map object associated with the domain of this operator.
  const Tpetra_Map & OperatorDomainMap() const 
  {
    return(*Map_);
   }

  //! Returns the Tpetra_Map object associated with the range of this operator.
  const Tpetra_Map & OperatorRangeMap() const 
  {
    return(*Map_);
  }
  //@}

const Tpetra_BlockMap& Map() const;

const char* Label() const{
  return(Label_);
};

private:
  void UpdateImportVector(int NumVectors) const;
  void UpdateExportVector(int NumVectors) const;

  //! Pointer to the matrix to be preconditioned.
  Teuchos::RCP<const Tpetra_RowMatrix> Matrix_;
  const Ifpack2_OverlappingRowMatrix* ovA_;
#ifdef HAVE_MPI
  //! Communicator containing this process only.
  Teuchos::RCP<Tpetra_MpiComm> SubComm_;
  MPI_Comm nodeMPIComm_;
#else
  //! Communicator containing this process only.
  Teuchos::RCP<Tpetra_SerialComm> SubComm_;
#endif
  //! Row, domain, and range map, based on SubComm_, containing the local rows only.
  Teuchos::RCP<Tpetra_Map> Map_;
  //! Column map based on SubComm_, containing the local rows only.
  Teuchos::RCP<Tpetra_Map> colMap_;
  //! Number of rows in the local matrix.
  int NumMyRows_;
  //! Number of cols in the local matrix.
  int NumMyCols_;
  //! Number of nonzeros in the local matrix.
  int NumMyNonzeros_;
  //! Number of rows in the global filtered matrix.
  int NumGlobalRows_;
  //! Number of nonzeros in the global filtered matrix.
  int NumGlobalNonzeros_;
  //! Maximum number of nonzero entries in a row for the filtered matrix.
  int MaxNumEntries_;
  //! Maximum number of nonzero entries in a row for Matrix_.
  int MaxNumEntriesA_;
  //! NumEntries_[i] contains the nonzero entries in row `i'.
  std::vector<int> NumEntries_;
  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable std::vector<int> Indices_;
  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable std::vector<double> Values_;
  //! If true, the tranpose of the local matrix will be used.
  bool UseTranspose_;
  //! Label for \c this object.
  char Label_[80];
  Teuchos::RCP<Tpetra_Vector> Diagonal_;
  double NormOne_;
  double NormInf_;

  //! Maps to speed LID-LID conversions
  int* Ac_LIDMap_;
  int* Bc_LIDMap_;
  int* Ar_LIDMap_;
  int* Br_LIDMap_;

  //! CrsMatrix pointer, if needed
  const Tpetra_CrsMatrix* Acrs_;

  int NumMyRowsA_;
  int NumMyRowsB_;
  
  //mutable Teuchos::RCP<Tpetra_MultiVector> ImportVector_;
  //mutable Teuchos::RCP<Tpetra_MultiVector> ExportVector_;
  mutable Tpetra_MultiVector* ExportVector_;
  mutable Tpetra_MultiVector* ImportVector_;
  Teuchos::RCP<Tpetra_Import> Importer_;
  Teuchos::RCP<Tpetra_Export> Exporter_;

};
#endif //ifdef IFPACK2_NODE_AWARE_CODE
#endif /* IFPACK2_NODEFILTER_HPP */
