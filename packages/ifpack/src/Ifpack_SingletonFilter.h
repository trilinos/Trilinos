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

#ifndef IFPACK_SINGLETONFILTER_H
#define IFPACK_SINGLETONFILTER_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_RefCountPtr.hpp"

class Epetra_Comm;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Import;
class Epetra_BlockMap;

//! Ifpack_SingletonFilter: Filter based on matrix entries.
//
class Ifpack_SingletonFilter : public virtual Epetra_RowMatrix {

public:
  //! Constructor.
  Ifpack_SingletonFilter(const Teuchos::RefCountPtr<Epetra_RowMatrix>& Matrix);

  //! Destructor.
  virtual ~Ifpack_SingletonFilter() {};

  //! Returns the number of entries in MyRow.
  virtual inline int NumMyRowEntries(int MyRow, int & /* NumEntries */) const
  {
    return(NumEntries_[MyRow]);
  }

  //! Returns the maximum number of entries.
  virtual int MaxNumEntries() const
  {
    return(MaxNumEntries_);
  }

  virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;

  virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const;

  virtual int Multiply(bool TransA, const Epetra_MultiVector& X, 
		       Epetra_MultiVector& Y) const;

  virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, 
		    const Epetra_MultiVector& X,
		    Epetra_MultiVector& Y) const;

  virtual int Apply(const Epetra_MultiVector& X,
		    Epetra_MultiVector& Y) const;

  virtual int ApplyInverse(const Epetra_MultiVector& X,
			   Epetra_MultiVector& Y) const;

  virtual int InvRowSums(Epetra_Vector& /* x */) const
  {
    return(-98); // NOT IMPLEMENTED
  }

  virtual int LeftScale(const Epetra_Vector& /* x */)
  {
    return(-98); // NOT IMPLEMENTED
  }

  virtual int InvColSums(Epetra_Vector& /* x */) const
  {
    return(-98); // NOT IMPLEMENTED
  }

  virtual int RightScale(const Epetra_Vector& /* x */) 
  {
    return(-98); // NOT IMPLEMENTED
  }

  virtual bool Filled() const
  {
    return(A_->Filled());
  }

  virtual double NormInf() const
  {
    return(-1.0);
  }

  virtual double NormOne() const
  {
    return(-1.0);
  }

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  virtual int NumGlobalNonzeros() const
  {
    return(NumNonzeros_);
  }

  virtual int NumGlobalRows() const
  {
    return(NumRows_);
  }

  virtual int NumGlobalCols() const
  {
    return(NumRows_);
  }

  virtual int NumGlobalDiagonals() const
  {
    return(NumRows_);
  }
#endif

  virtual long long NumGlobalNonzeros64() const
  {
    return(NumNonzeros_);
  }

  virtual long long NumGlobalRows64() const
  {
    return(NumRows_);
  }

  virtual long long NumGlobalCols64() const
  {
    return(NumRows_);
  }

  virtual long long NumGlobalDiagonals64() const
  {
    return(NumRows_);
  }

  virtual int NumMyNonzeros() const
  {
    return(NumNonzeros_);
  }

  virtual int NumMyRows() const
  {
    return(NumRows_);
  }

  virtual int NumMyCols() const
  {
    return(NumRows_);
  }

  virtual int NumMyDiagonals() const
  {
    return(NumRows_);
  }

  virtual bool LowerTriangular() const
  {
    return(false);
  }

  virtual bool UpperTriangular() const
  {
    return(false);
  }

  virtual const Epetra_Map & RowMatrixRowMap() const
  {
    return(*Map_);
  }

  virtual const Epetra_Map & RowMatrixColMap() const
  {
    return(*Map_);
  }

  virtual const Epetra_Import * RowMatrixImporter() const
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

  bool HasNormInf() const
  {
    return(false);
  }

  const Epetra_Comm & Comm() const
  {
    return(A_->Comm());
  }

  const Epetra_Map & OperatorDomainMap() const 
  {
    return(*Map_);
  }

  const Epetra_Map & OperatorRangeMap() const 
  {
    return(*Map_);
  }

  const Epetra_BlockMap& Map() const 
  {
    return(*(const Epetra_BlockMap*)(&*Map_));
  }

  const char* Label() const{
    return(Label_);
  }

  int SolveSingletons(const Epetra_MultiVector& RHS, 
		      Epetra_MultiVector& LHS);

  int CreateReducedRHS(const Epetra_MultiVector& LHS,
		       const Epetra_MultiVector& RHS, 
		       Epetra_MultiVector& ReducedRHS);

  int UpdateLHS(const Epetra_MultiVector& ReducedLHS,
		Epetra_MultiVector& LHS);

private:

  //! Pointer to the matrix to be preconditioned.
  Teuchos::RefCountPtr<Epetra_RowMatrix> A_;

  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable std::vector<int> Indices_;
  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable std::vector<double> Values_;
  //! Label for \c this object.
  char Label_[80];
  
  int NumSingletons_;
  std::vector<int> SingletonIndex_;

  std::vector<int> Reorder_;
  std::vector<int> InvReorder_;

  std::vector<int> NumEntries_;

  int NumRows_;
  int NumRowsA_;
  int MaxNumEntries_;
  int MaxNumEntriesA_;
  int NumNonzeros_;
  Teuchos::RefCountPtr<Epetra_Map> Map_;
  
  Teuchos::RefCountPtr<Epetra_Vector> Diagonal_;

};

#endif /* IFPACK_SINGLETONFILTER_H */
