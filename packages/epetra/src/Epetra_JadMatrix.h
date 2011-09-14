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

#ifndef EPETRA_JADMATRIX_H
#define EPETRA_JADMATRIX_H

#include "Epetra_BasicRowMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"


class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_Import;
class Epetra_Export;

//! Epetra_JadMatrix: A class for constructing matrix objects optimized for common kernels.

/*! The Epetra_JadMatrix class takes an existing Epetra_RowMatrix ojbect, analyzes it and 
    builds a jagged diagonal equivalent of it. Once constructed, it is also possible to 
    update the values of the matrix with values from another Epetra_RowMatrix that has
    the identical structure.
    
*/    

class EPETRA_LIB_DLL_EXPORT Epetra_JadMatrix: public Epetra_BasicRowMatrix {
      
 public:

   //! @name Constructors/Destructor
  //@{ 
  //! Epetra_JadMatrix constuctor.
  /* The constructor for this class requires a fully constructed instance of an Epetra_RowMatrix
     object.
     \param Matrix (In) An existing Epetra_RowMatrix.
     \pre Matrix must have Matrix.Filled()==true.
  */
  Epetra_JadMatrix(const Epetra_RowMatrix & Matrix);

  //! Epetra_JadMatrix Destructor
  virtual ~Epetra_JadMatrix();
  //@}
  
  //! @name Post-construction modifications
  //@{ 
  //! Update values using a matrix with identical structure.
  /* Updates the values only using a matrix that has exactly the same structure as
     the matrix used to construct this Epetra_JadMatrix object.  Once the constructor
     is called, the Matrix argument is no longer needed.
     \param Matrix (In) An existing Epetra_RowMatrix with \e identical structure to 
            the matrix used to create this Epetra_JadMatrix.
     \param CheckStructure (In) Optional argument, by default is false.  If set to true, 
            the method will check to see if the structure of Matrix is compatible with
	    the structure of matrix used to create this Epetra_JadMatrix.  Performing
	    this check has signficant overhead, so it should only be turned on for debugging.
     \pre Matrix must have Matrix.Filled()==true.
  */ 
  int UpdateValues(const Epetra_RowMatrix & Matrix, bool CheckStructure = false);
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
    int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;

    //! Returns a reference to the ith entry in the matrix, along with its row and column index.
    /*! 
    \param CurEntry (In) - Local entry to extract.
    \param Value (Out) - Extracted reference to current values.
    \param RowIndex (Out) - Row index for current entry.
    \param ColIndex (Out) - Column index for current entry.
	  
    \return Integer error code, set to 0 if successful, set to -1 if CurEntry not valid.
  */
    int ExtractMyEntryView(int CurEntry, double * &Value, int & RowIndex, int & ColIndex) { 
      if (CurEntry>=NumMyNonzeros_) EPETRA_CHK_ERR(-1); 
      Value = &Values_[CurEntry];
      ColIndex = Indices_[CurEntry];
      for (int j=0; j<NumJaggedDiagonals_; j++) if (CurEntry<IndexOffset_[j+1]) {RowIndex = RowPerm_[CurEntry-IndexOffset_[j]]; break;}
      return(0);
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
      if (CurEntry>=NumMyNonzeros_) EPETRA_CHK_ERR(-1); 
      Value = &Values_[CurEntry];
      ColIndex = Indices_[CurEntry];
      for (int j=0; j<NumJaggedDiagonals_; j++) if (CurEntry<IndexOffset_[j+1]) {RowIndex = RowPerm_[CurEntry-IndexOffset_[j]]; break;}
      return(0);
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
    int NumMyRowEntries(int MyRow, int & NumEntries) const;

    //@}

    //! @name Computational methods
  //@{ 

    //! Returns the result of a Epetra_JadMatrix multiplied by a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   TransA -If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
  */
    int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Returns the result of a Epetra_JadMatrix solve with a Epetra_MultiVector X in Y (not implemented).
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
              Epetra_MultiVector& Y) const
    {
      (void)Upper;
      (void)Trans;
      (void)UnitDiagonal;
      (void)X;
      (void)Y;
      return(-1);
    }
  //@}



 protected:

  void GeneralMV(bool TransA, double * x, double * y) const;
  void GeneralMM(bool TransA, double ** X, int LDX, double ** Y, int LDY, int NumVectors) const;
  void GeneralMM3RHS(bool TransA, double ** X, int LDX, double ** Y, int LDY, int NumVectors) const;
  void GeneralMM2RHS(bool TransA, double * x, int ldx, double * y, int ldy) const;
  void Allocate(const Epetra_RowMatrix & Matrix);

  Epetra_SerialDenseVector Values_;
  Epetra_IntSerialDenseVector Indices_;
  Epetra_IntSerialDenseVector IndexOffset_;
  Epetra_IntSerialDenseVector Profile_;
  Epetra_IntSerialDenseVector RowPerm_;
  Epetra_IntSerialDenseVector InvRowPerm_;
  int NumJaggedDiagonals_;

};
#endif /* EPETRA_JADMATRIX_H */
