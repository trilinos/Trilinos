/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_ROWMATRIX_H
#define EPETRA_ROWMATRIX_H

class Epetra_Comm;
class Epetra_Import;
class Epetra_Export;
class Epetra_Vector;
class Epetra_MultiVector;
#include "Epetra_ConfigDefs.h"
#include "Epetra_Operator.h"
#include "Epetra_SrcDistObject.h"


//! Epetra_RowMatrix: A pure virtual class for using real-valued double-precision row matrices.

/*! The Epetra_RowMatrix class is a pure virtual class (specifies interface only) that 
    enable the use of real-valued double-precision sparse matrices
    where matrix entries are intended for row access.  It is currently implemented by both the
    Epetra_CrsMatrix and Epetra_VbrMatrix classes.

   
*/    


class EPETRA_LIB_DLL_EXPORT Epetra_RowMatrix: public virtual Epetra_Operator, public virtual Epetra_SrcDistObject {
      
 public:
   //! @name Destructor
  //@{ 
    //! Destructor
    virtual ~Epetra_RowMatrix() {};

  //@}
  
  //! @name Matrix data extraction routines
  //@{ 

    //! Returns the number of nonzero entries in MyRow.
    /*! 
    \param In
           MyRow - Local row.
    \param Out
     NumEntries - Number of nonzero values present.
    
    \return Integer error code, set to 0 if successful.
  */
    virtual int NumMyRowEntries(int MyRow, int & NumEntries) const = 0;


    //! Returns the maximum of NumMyRowEntries() over all rows.
    virtual int MaxNumEntries() const = 0;

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
    virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const = 0;

    //! Returns a copy of the main diagonal in a user-provided vector.
    /*! 
    \param Out
     Diagonal - Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
  */
    virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const = 0;
  //@}
  
  //! @name Mathematical functions
  //@{ 

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
    virtual int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

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
          Epetra_MultiVector& Y) const = 0;

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
    virtual int InvRowSums(Epetra_Vector& x) const = 0;

    //! Scales the Epetra_RowMatrix on the left with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
        and j denotes the column number of A.
    \param In
     x -A Epetra_Vector to solve for.

    \return Integer error code, set to 0 if successful.
  */
    virtual int LeftScale(const Epetra_Vector& x) = 0;

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
    virtual int InvColSums(Epetra_Vector& x) const = 0;

    //! Scales the Epetra_RowMatrix on the right with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.
    \param In
     x -The Epetra_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
  */
    virtual int RightScale(const Epetra_Vector& x) = 0;
  //@}
  
  //! @name Attribute access functions
  //@{ 

    //! If FillComplete() has been called, this query returns true, otherwise it returns false.
    virtual bool Filled() const = 0;

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
    /*
      Note that depending on the matrix implementation, it is sometimes
      possible to have some nonzeros that appear on multiple processors.
      In that case, those nonzeros may be counted multiple times (also
      depending on the matrix implementation).
    */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    virtual int NumGlobalNonzeros() const = 0;
#endif
    virtual long long NumGlobalNonzeros64() const = 0;

    //! Returns the number of global matrix rows.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    virtual int NumGlobalRows() const = 0;
#endif
    virtual long long NumGlobalRows64() const = 0;

    //! Returns the number of global matrix columns.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    virtual int NumGlobalCols() const= 0;
#endif
    virtual long long NumGlobalCols64() const= 0;

    //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    virtual int NumGlobalDiagonals() const = 0;
#endif
    virtual long long NumGlobalDiagonals64() const = 0;
    
    //! Returns the number of nonzero entries in the calling processor's portion of the matrix.
    virtual int NumMyNonzeros() const = 0;

    //! Returns the number of matrix rows owned by the calling processor.
    virtual int NumMyRows() const = 0;

    //! Returns the number of matrix columns owned by the calling processor.
    virtual int NumMyCols() const = 0;

    //! Returns the number of local nonzero diagonal entries, based on global row/column index comparisons.
    virtual int NumMyDiagonals() const = 0;

    //! If matrix is lower triangular in local index space, this query returns true, otherwise it returns false.
    virtual bool LowerTriangular() const = 0;

    //! If matrix is upper triangular in local index space, this query returns true, otherwise it returns false.
    virtual bool UpperTriangular() const = 0;

    //! Returns the Epetra_Map object associated with the rows of this matrix.
    virtual const Epetra_Map & RowMatrixRowMap() const = 0;

    //! Returns the Epetra_Map object associated with the columns of this matrix.
    virtual const Epetra_Map & RowMatrixColMap() const = 0;

    //! Returns the Epetra_Import object that contains the import operations for distributed operations.
    virtual const Epetra_Import * RowMatrixImporter() const = 0;
  //@}
};

#endif /* EPETRA_ROWMATRIX_H */
