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

#ifndef EPETRA_VBRROWMATRIX_H
#define EPETRA_VBRROWMATRIX_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif



#include "Epetra_BasicRowMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

//! Epetra_VbrRowMatrix: A class for using an existing Epetra_VbrMatrix object as an Epetra_RowMatrix object.

/*! The Epetra_VbrRowMatrix class takes an existing Epetra_VbrMatrix object and allows its
    use as an Epetra_RowMatrix without allocating additional storage.  Although the Epetra_VbrMatrix itself
    inherits from Epetra_RowMatrix, a design flaw in the inheritance structure of Epetra prohibits the use of
    an Epetra_VbrMatrix object as an Epetra_RowMatrix in some important situations.  Therefore we recommend the
    use of this class to wrap an Epetra_VbrMatrix object.

    \warning This class takes a pointer to an existing Epetra_VbrMatrix object.  It is assumed that the user
    will pass in a pointer to a valid Epetra_VbrMatrix object, and will retain it throughout the life of the
    Epetra_VbrRowMatrix object.

*/

class Epetra_VbrRowMatrix: public Epetra_BasicRowMatrix {

 public:

   //! @name Constructors/Destructor
  //@{
  //! Epetra_VbrRowMatrix constuctor.
  /* The constructor for this class requires a pointer to a fully constructed instance of an Epetra_VbrMatrix
     object.
     \param Matrix (In) Pointer to an existing Epetra_VbrMatrix. The input matrix must be retained by the user
     throughout the existance of the dependent Epetra_VbrRowmatrix object.
     \pre Matrix must have Matrix->Filled()==true.
  */
  Epetra_VbrRowMatrix(Epetra_VbrMatrix * Matrix): Epetra_BasicRowMatrix(Matrix->Comm()), matrix_(Matrix) {
  if (Matrix==0) throw Matrix->RowMatrixRowMap().ReportError("Input matrix must have called FillComplete()", -1);
  SetMaps(Matrix->RowMatrixRowMap(), Matrix->RowMatrixColMap(), Matrix->OperatorDomainMap(), Matrix->OperatorRangeMap());
  if (!Matrix->Filled()) throw Matrix->RowMatrixRowMap().ReportError("Input matrix must have called FillComplete()", -1);
  SetLabel("Epetra::VbrRowMatrix");
}

  //! Epetra_VbrRowMatrix Destructor
    virtual ~Epetra_VbrRowMatrix(){}
  //@}

  //! @name Post-construction modifications
  //@{
  //! Update the matrix to which this object points.
  /* Updates the matrix that the Epetra_VbrRowMatrix will use to satisfy the Epetra_RowMatrix functionality.
     \param Matrix (In) A pointer to an existing, fully constructed Epetra_VbrMatrix.
     \pre Matrix must have Matrix->Filled()==true.
  */
  int UpdateMatrix(Epetra_VbrMatrix * Matrix){
    if (Matrix ==0) {
      EPETRA_CHK_ERR(-1);
    }
    else matrix_ = Matrix;
    return(0);
  }
  //@}

  //! @name Methods required for implementing Epetra_BasicRowMatrix
  //@{

    //! Returns a copy of the specified local row in user-provided arrays.
    /*!
    \param MyRow (In) - Local row to extract.
    \param Length (In) - Length of Values and Indices.
    \param NumEntries (Out) - Number of nonzero entries extracted.
    \param Values (Out) - Extracted values for this row.
    \param Indices (Out) - Extracted global column indices for the corresponding values.
	
    \return Integer error code, set to 0 if successful, set to -1 if MyRow not valid, -2 if Length is too short (NumEntries will have required length).
  */
  int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const {

    EPETRA_CHK_ERR(matrix_->ExtractMyRowCopy(MyRow, Length, NumEntries, Values, Indices));
    return(0);
  }

    //! Returns a reference to the ith entry in the matrix, along with its row and column index
    /*!
    \param CurEntry (In) - Local entry to extract.
    \param Value (Out) - Extracted reference to current values.
    \param RowIndex (Out) - Row index for current entry.
    \param ColIndex (Out) - Column index for current entry.
	
    \return Integer error code, set to 0 if successful, set to -1 if CurEntry not valid.
  */
    int ExtractMyEntryView(int CurEntry, double * &Value, int & RowIndex, int & ColIndex) {
      return(-1);
    }

    //! Returns a const reference to the ith entry in the matrix, along with its row and column index.
    /*!
    \param CurEntry (In) - Local entry to extract.
    \param Value (Out) - Extracted reference to current values.
    \param RowIndex (Out) - Row index for current entry.
    \param ColIndex (Out) - Column index for current entry.
	
    \return Integer error code, set to 0 if successful, set to -1 if CurEntry not valid.
  */
    int ExtractMyEntryView(int CurEntry, double const * & Value, int & RowIndex, int & ColIndex) const {
      return(-1);
    }

    //! Return the current number of values stored for the specified local row.
    /*! Similar to NumMyEntries() except NumEntries is returned as an argument
      and error checking is done on the input value MyRow.
      \param MyRow - (In) Local row.
      \param NumEntries - (Out) Number of nonzero values.

      \return Integer error code, set to 0 if successful, set to -1 if MyRow not valid.
      \pre None.
      \post Unchanged.
    */
    int NumMyRowEntries(int MyRow, int & NumEntries) const {
      EPETRA_CHK_ERR(matrix_->NumMyRowEntries(MyRow, NumEntries));
      return(0);
    }

    //@}

    //! @name Computational methods
  //@{

    //! Scales the Epetra_VbrMatrix on the right with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.
    \param In
	   x -The Epetra_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
  */
    int RightScale(const Epetra_Vector& x){
      HaveNumericConstants_ = false;
      UpdateFlops(NumGlobalNonzeros());
      EPETRA_CHK_ERR(matrix_->RightScale(x));
      return(0);
    }

    //! Scales the Epetra_VbrMatrix on the left with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
        and j denotes the column number of A.
    \param In
	   x -A Epetra_Vector to solve for.

    \return Integer error code, set to 0 if successful.
  */
    int LeftScale(const Epetra_Vector& x){
      HaveNumericConstants_ = false;
      UpdateFlops(NumGlobalNonzeros());
      EPETRA_CHK_ERR(matrix_->LeftScale(x));
      return(0);
    }

    //! Returns the result of a Epetra_VbrRowMatrix multiplied by a Epetra_MultiVector X in Y.
    /*!
    \param In
	   TransA -If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
  */
    int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const{
      EPETRA_CHK_ERR(matrix_->Multiply(TransA, X, Y));
      return(0);
    }

    //! Returns the result of a Epetra_VbrRowMatrix solve with a Epetra_MultiVector X in Y (not implemented).
    /*!
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
    int Solve(bool Upper, bool Trans, bool UnitDiagonal,
              const Epetra_MultiVector& X,
              Epetra_MultiVector& Y) const {
      EPETRA_CHK_ERR(matrix_->Solve(Upper, Trans, UnitDiagonal, X, Y));
      return(0);
    }  //@}



 private:

    Epetra_VbrMatrix * matrix_;

};
#endif /* EPETRA_VBRROWMATRIX_H */
