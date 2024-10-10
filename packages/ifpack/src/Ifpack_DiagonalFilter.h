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

#ifndef IFPACK_DIAGONALFILTER_H
#define IFPACK_DIAGONALFILTER_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_ConfigDefs.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Time.h"
#include "Teuchos_RefCountPtr.hpp"

class Epetra_Comm;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Import;
class Epetra_BlockMap;

//! Ifpack_DiagonalFilter: Filter to modify the diagonal entries of a given Epetra_RowMatrix.
/*!

Ifpack_DiagonalFilter modifies the elements on the diagonal.

A typical use is as follows:
\code
Teuchos::RefCountPtr<Epetra_RowMatrix> A;
// creates a matrix B such that
// B(i,i) = AbsoluteThreshold * sgn(B(i,i)) + 
//          RelativeThreshold * B(i,i)
double AbsoluteThreshold = 1e-3;
double RelativeThreshold = 1.01;

Ifpack_DiagonalFilter B(A, AbsoluteThreshold, RelativeThreshold);
\endcode

\author Marzio Sala, SNL 9214.

\data Last modified on 24-Jan-05.

*/
class Ifpack_DiagonalFilter : public virtual Epetra_RowMatrix {

public:
  //! Constructor.
  Ifpack_DiagonalFilter(const Teuchos::RefCountPtr<Epetra_RowMatrix>& Matrix,
                        double AbsoluteThreshold,
                        double RelativeThreshold);
  
  //! Destructor.
  virtual ~Ifpack_DiagonalFilter() {};

  //! Returns the number of entries in MyRow.
  virtual int NumMyRowEntries(int MyRow, int& NumEntries) const
  {
    return(A_->NumMyRowEntries(MyRow, NumEntries));
  }

  //! Returns the maximum number of entries.
  virtual int MaxNumEntries() const
  {
    return(A_->MaxNumEntries());
  }

  inline virtual int ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, 
                               double* Values, int* Indices) const;

  virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
  {
    int ierr = A_->ExtractDiagonalCopy(Diagonal);
    IFPACK_RETURN(ierr);
  }

  virtual int Multiply(bool TransA, const Epetra_MultiVector& X, 
		       Epetra_MultiVector& Y) const;

  virtual int Solve(bool /* Upper */, bool /* Trans */, bool /* UnitDiagonal */, 
		    const Epetra_MultiVector& /* X */,
		    Epetra_MultiVector& /* Y */) const
  {
    IFPACK_CHK_ERR(-1);
  }

  virtual int Apply(const Epetra_MultiVector& X,
		    Epetra_MultiVector& Y) const
  {
    int ierr = Multiply(UseTranspose(),X,Y);
    IFPACK_RETURN(ierr);
  }

  virtual int ApplyInverse(const Epetra_MultiVector& /* X */,
			   Epetra_MultiVector& /* Y */) const
  {
    IFPACK_CHK_ERR(-1);
  }

  virtual int InvRowSums(Epetra_Vector& /* x */) const
  {
    IFPACK_CHK_ERR(-1);
  }

  virtual int LeftScale(const Epetra_Vector& x)
  {
    return(A_->LeftScale(x));
  }

  virtual int InvColSums(Epetra_Vector& /* x */) const
  {
    IFPACK_CHK_ERR(-1);;
  }

  virtual int RightScale(const Epetra_Vector& x) 
  {
    return(A_->RightScale(x));
  }

  virtual bool Filled() const
  {
    return(A_->Filled());
  }

  //! Not implemented for efficiency reasons.
  virtual double NormInf() const
  {
    return(-1.0);
  }

  //! Not implemented for efficiency reasons.
  virtual double NormOne() const
  {
    return(-1.0);
  }

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  virtual int NumGlobalNonzeros() const
  {
    return(A_->NumGlobalNonzeros());
  }

  virtual int NumGlobalRows() const
  {
    return(A_->NumGlobalRows());
  }

  virtual int NumGlobalCols() const
  {
    return(A_->NumGlobalCols());
  }

  virtual int NumGlobalDiagonals() const
  {
    return(A_->NumGlobalDiagonals());
  }
#endif

  virtual long long NumGlobalNonzeros64() const
  {
    return(A_->NumGlobalNonzeros64());
  }

  virtual long long NumGlobalRows64() const
  {
    return(A_->NumGlobalRows64());
  }

  virtual long long NumGlobalCols64() const
  {
    return(A_->NumGlobalCols64());
  }

  virtual long long NumGlobalDiagonals64() const
  {
    return(A_->NumGlobalDiagonals64());
  }

  virtual int NumMyNonzeros() const
  {
    return(A_->NumMyNonzeros());
  }

  virtual int NumMyRows() const
  {
    return(A_->NumMyRows());
  }

  virtual int NumMyCols() const
  {
    return(A_->NumMyCols());
  }

  virtual int NumMyDiagonals() const
  {
    return(A_->NumMyDiagonals());
  }

  virtual bool LowerTriangular() const
  {
    return(A_->LowerTriangular());
  }

  virtual bool UpperTriangular() const
  {
    return(A_->UpperTriangular());
  }

  virtual const Epetra_Map& RowMatrixRowMap() const
  {
    return(A_->RowMatrixRowMap());
  }

  virtual const Epetra_Map& RowMatrixColMap() const
  {
    return(A_->RowMatrixColMap());
  }

  virtual const Epetra_Import* RowMatrixImporter() const
  {
    return(A_->RowMatrixImporter());
  }

  int SetUseTranspose(bool UseTranspose_in)
  {
    return(A_->SetUseTranspose(UseTranspose_in));
  }

  bool UseTranspose() const 
  {
    return(A_->UseTranspose());
  }

  //! Not implemented for efficiency reasons.
  bool HasNormInf() const
  {
    return(false);
  }

  const Epetra_Comm& Comm() const
  {
    return(A_->Comm());
  }

  const Epetra_Map& OperatorDomainMap() const 
  {
    return(A_->OperatorDomainMap());
  }

  const Epetra_Map& OperatorRangeMap() const 
  {
    return(A_->OperatorRangeMap());
  }

  const Epetra_BlockMap& Map() const 
  {
    return(A_->Map());
  }

  const char* Label() const{
    return(A_->Label());
  }

private:

  //! Pointer to the matrix to be filtered
  Teuchos::RefCountPtr<Epetra_RowMatrix> A_;
  //! This value (times the sgn(A(i,i)) is added to the diagonal elements
  double AbsoluteThreshold_;
  //! Multiplies A(i,i) by this value.
  double RelativeThreshold_;
  //! Stores the position of the diagonal element, or -1 if not present.
  std::vector<int> pos_;
  //! Stores as additional diagonal contribution due to the filter.
  std::vector<double> val_;

};


#endif /* IFPACK_DIAGONALFILTER_H */
