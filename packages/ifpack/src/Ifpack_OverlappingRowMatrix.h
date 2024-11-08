
/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
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
*/

#ifndef IFPACK_OVERLAPPINGROWMATRIX_H
#define IFPACK_OVERLAPPINGROWMATRIX_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CombineMode.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Import.h"
#include "Epetra_Map.h"
#ifdef HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS
#include "Epetra_IntVector.h"
#else
# ifdef IFPACK_NODE_AWARE_CODE
# include "Epetra_IntVector.h"
# endif
#endif

class Epetra_Map;
class Epetra_BlockMap;
class Epetra_CrsMatrix;
class Epetra_Comm;

//! Ifpack_OverlappingRowMatrix: matrix with ghost rows, based on Epetra_RowMatrix
//
class Ifpack_OverlappingRowMatrix : public virtual Epetra_RowMatrix {

public:

  //@{ Constructors/Destructors
#ifdef HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS
  Ifpack_OverlappingRowMatrix(const Teuchos::RefCountPtr<const Epetra_RowMatrix>& Matrix_in,
                              int OverlapLevel_in, int subdomainID);
#else
# ifdef IFPACK_NODE_AWARE_CODE
  Ifpack_OverlappingRowMatrix(const Teuchos::RefCountPtr<const Epetra_RowMatrix>& Matrix_in,
                              int OverlapLevel_in, int myNodeID);
# endif
#endif
  Ifpack_OverlappingRowMatrix(const Teuchos::RefCountPtr<const Epetra_RowMatrix>& Matrix_in,
                              int OverlapLevel_in);

#ifdef HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS
  ~Ifpack_OverlappingRowMatrix() {};
#else
# ifdef IFPACK_NODE_AWARE_CODE
  ~Ifpack_OverlappingRowMatrix();
# else
  ~Ifpack_OverlappingRowMatrix() {};
# endif
#endif
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
  virtual int NumMyRowEntries(int MyRow, int & NumEntries) const;

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
  virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;
#ifdef HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS
  virtual int ExtractGlobalRowCopy(int MyRow, int Length, int & NumEntries, double* Values, int* Indices) const;
#else
# ifdef IFPACK_NODE_AWARE_CODE
  virtual int ExtractGlobalRowCopy(int MyRow, int Length, int & NumEntries, double* Values, int* Indices) const;
# endif
#endif

  //! Returns a copy of the main diagonal in a user-provided vector.
  /*!
    \param
    Diagonal - (Out) Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
    */
  virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const;
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
  virtual int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns result of a local-only solve using a triangular Epetra_RowMatrix with Epetra_MultiVectors X and Y (NOT IMPLEMENTED).
  virtual int Solve(bool /* Upper */, bool /* Trans */, bool /* UnitDiagonal */, const Epetra_MultiVector& /* X */,
                    Epetra_MultiVector& /* Y */) const
  {
    IFPACK_RETURN(-1); // not implemented
  }

  virtual int Apply(const Epetra_MultiVector& X,
                    Epetra_MultiVector& Y) const;

  virtual int ApplyInverse(const Epetra_MultiVector& X,
                           Epetra_MultiVector& Y) const;
  //! Computes the sum of absolute values of the rows of the Epetra_RowMatrix, results returned in x (NOT IMPLEMENTED).
  virtual int InvRowSums(Epetra_Vector& /* x */) const
  {
    IFPACK_RETURN(-1); // not implemented
  }

  //! Scales the Epetra_RowMatrix on the left with a Epetra_Vector x (NOT IMPLEMENTED).
  virtual int LeftScale(const Epetra_Vector& /* x */)
  {
    IFPACK_RETURN(-1); // not implemented
  }

  //! Computes the sum of absolute values of the columns of the Epetra_RowMatrix, results returned in x (NOT IMPLEMENTED).
  virtual int InvColSums(Epetra_Vector& /* x */) const
  {
    IFPACK_RETURN(-1); // not implemented
  }


  //! Scales the Epetra_RowMatrix on the right with a Epetra_Vector x (NOT IMPLEMENTED).
  virtual int RightScale(const Epetra_Vector& /* x */)
  {
    IFPACK_RETURN(-1); // not implemented
  }

  //@}

  //@{ \name Attribute access functions

  //! If FillComplete() has been called, this query returns true, otherwise it returns false.
  virtual bool Filled() const
  {
    return(true);
  }

  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\len} \sum_{i=1}^m |a_{ij}| \f].
     */
  virtual double NormInf() const
  {
    return(A().NormInf());
  }

  //! Returns the one norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_1\f$ such that
     \f[\| A \|_1= \max_{1\lej\len} \sum_{j=1}^n |a_{ij}| \f].
     */
  virtual double NormOne() const
  {
    return(A().NormOne());
  }

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  //! Returns the number of nonzero entries in the global matrix.
  virtual int NumGlobalNonzeros() const
  {
    if(A().RowMatrixRowMap().GlobalIndicesInt())
       return (int) NumGlobalNonzeros_;
    else
       throw "Ifpack_OverlappingRowMatrix::NumGlobalNonzeros: Global indices not int";
  }

  //! Returns the number of global matrix rows.
  virtual int NumGlobalRows() const
  {
    return(A().NumGlobalRows());
  }

  //! Returns the number of global matrix columns.
  virtual int NumGlobalCols() const
  {
    return(A().NumGlobalCols());
  }

  //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
  virtual int NumGlobalDiagonals() const
  {
    return(A().NumGlobalDiagonals());
  }
#endif
  //! Returns the number of nonzero entries in the global matrix.
  virtual long long NumGlobalNonzeros64() const
  {
    return(NumGlobalNonzeros_);
  }

  //! Returns the number of global matrix rows.
  virtual long long NumGlobalRows64() const
  {
    return(A().NumGlobalRows64());
  }

  //! Returns the number of global matrix columns.
  virtual long long NumGlobalCols64() const
  {
    return(A().NumGlobalCols64());
  }

  //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
  virtual long long NumGlobalDiagonals64() const
  {
    return(A().NumGlobalDiagonals64());
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
    return(NumMyDiagonals_);
  }

  //! If matrix is lower triangular in local index space, this query returns true, otherwise it returns false.
  virtual bool LowerTriangular() const
  {
    return(A().LowerTriangular());
  }

  //! If matrix is upper triangular in local index space, this query returns true, otherwise it returns false.
  virtual bool UpperTriangular() const
  {
    return(A().UpperTriangular());
  }

  //! Returns the Epetra_Map object associated with the rows of this matrix.
  virtual const Epetra_Map & RowMatrixRowMap() const
  {
    return(*Map_);
  }

  //! Returns the Epetra_Map object associated with the columns of this matrix.
  virtual const Epetra_Map & RowMatrixColMap() const
  {
#ifdef HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS
    return(*colMap_);
#else
#   ifdef IFPACK_NODE_AWARE_CODE
    return(*colMap_);
#   else
    return(*Map_);
#   endif
#endif
  }

  //! Returns the Epetra_Import object that contains the import operations for distributed operations.
  virtual const Epetra_Import * RowMatrixImporter() const
  {
    return(&*Importer_);
  }
  //@}

  // following functions are required to derive Epetra_RowMatrix objects.

  //! Sets ownership.
  int SetOwnership(bool /* ownership */)
  {
    IFPACK_RETURN(-1);
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
    return(A().HasNormInf());
  }

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & Comm() const
  {
    return(A().Comm());
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

const Epetra_BlockMap& Map() const;

const char* Label() const{
  return(Label_.c_str());
};

int OverlapLevel() const
{
  return(OverlapLevel_);
}

int ImportMultiVector(const Epetra_MultiVector& X,
                      Epetra_MultiVector& OvX,
                      Epetra_CombineMode CM = Insert);

int ExportMultiVector(const Epetra_MultiVector& OvX,
                      Epetra_MultiVector& X,
                      Epetra_CombineMode CM = Add);
#ifdef HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS
  inline const Epetra_RowMatrix& A() const
  {
    return(*Matrix_);
  }

  inline Epetra_CrsMatrix& B() const
  {
    return(*ExtMatrix_);
  }
#else
# ifdef IFPACK_NODE_AWARE_CODE
  inline const Epetra_RowMatrix& A() const
  {
    return(*Matrix_);
  }

  inline Epetra_CrsMatrix& B() const
  {
    return(*ExtMatrix_);
  }
# endif
#endif

private:
#ifndef HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS
# ifndef IFPACK_NODE_AWARE_CODE
  inline const Epetra_RowMatrix& A() const
  {
    return(*Matrix_);
  }

  inline Epetra_RowMatrix& B() const;
# endif
#endif

  int NumMyRows_;
  int NumMyCols_;
  int NumMyDiagonals_;
  int NumMyNonzeros_;

  long long NumGlobalNonzeros_;
  int MaxNumEntries_;

  int NumMyRowsA_;
  int NumMyRowsB_;

  bool UseTranspose_;

  Teuchos::RefCountPtr<const Epetra_Map> Map_;
#ifdef HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS
  const Epetra_Map *colMap_;
#else
# ifdef IFPACK_NODE_AWARE_CODE
  const Epetra_Map *colMap_;
# endif
#endif
  Teuchos::RefCountPtr<const Epetra_Import> Importer_;

  Teuchos::RefCountPtr<const Epetra_RowMatrix> Matrix_;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> ExtMatrix_;
  Teuchos::RefCountPtr<Epetra_Map> ExtMap_;
  Teuchos::RefCountPtr<Epetra_Import> ExtImporter_;

  int OverlapLevel_;
  std::string Label_;

  template<typename int_type>
  void BuildMap(int OverlapLevel_in);

}; // class Ifpack_OverlappingRowMatrix

#endif // IFPACK_OVERLAPPINGROWMATRIX_H
