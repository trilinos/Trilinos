/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef TIFPACK_SINGLETONFILTER_HPP
#define TIFPACK_SINGLETONFILTER_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Teuchos_RefCountPtr.hpp"

class Tpetra_Comm;
class Tpetra_Map;
class Tpetra_MultiVector;
class Tpetra_Import;
class Tpetra_BlockMap;

//! Tifpack_SingletonFilter: Filter based on matrix entries.
//
class Tifpack_SingletonFilter : public virtual Tpetra_RowMatrix {

public:
  //! Constructor.
  Tifpack_SingletonFilter(const Teuchos::RCP<Tpetra_RowMatrix>& Matrix);

  //! Destructor.
  virtual ~Tifpack_SingletonFilter() {};

  //! Returns the number of entries in MyRow.
  virtual inline int NumMyRowEntries(int MyRow, int & NumEntries) const
  {
    return(NumEntries_[MyRow]);
  }

  //! Returns the maximum number of entries.
  virtual int MaxNumEntries() const
  {
    return(MaxNumEntries_);
  }

  virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;

  virtual int ExtractDiagonalCopy(Tpetra_Vector & Diagonal) const;

  virtual int Multiply(bool TransA, const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
		       Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, 
		    const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
		    Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual int Apply(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
		    Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual int ApplyInverse(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
			   Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual int InvRowSums(Tpetra_Vector& x) const
  {
    return(-98); // NOT IMPLEMENTED
  }

  virtual int LeftScale(const Tpetra_Vector& x)
  {
    return(-98); // NOT IMPLEMENTED
  }

  virtual int InvColSums(Tpetra_Vector& x) const
  {
    return(-98); // NOT IMPLEMENTED
  }

  virtual int RightScale(const Tpetra_Vector& x) 
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

  virtual const Tpetra_Map & RowMatrixRowMap() const
  {
    return(*Map_);
  }

  virtual const Tpetra_Map & RowMatrixColMap() const
  {
    return(*Map_);
  }

  virtual const Tpetra_Import * RowMatrixImporter() const
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

  const Tpetra_Comm & Comm() const
  {
    return(A_->Comm());
  }

  const Tpetra_Map & OperatorDomainMap() const 
  {
    return(*Map_);
  }

  const Tpetra_Map & OperatorRangeMap() const 
  {
    return(*Map_);
  }

  const Tpetra_BlockMap& Map() const 
  {
    return(*(const Tpetra_BlockMap*)(&*Map_));
  }

  const char* Label() const{
    return(Label_);
  }

  int SolveSingletons(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS, 
		      Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& LHS);

  int CreateReducedRHS(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& LHS,
		       const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS, 
		       Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& ReducedRHS);

  int UpdateLHS(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& ReducedLHS,
		Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& LHS);

private:

  //! Pointer to the matrix to be preconditioned.
  Teuchos::RCP<Tpetra_RowMatrix> A_;

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
  Teuchos::RCP<Tpetra_Map> Map_;
  
  Teuchos::RCP<Tpetra_Vector> Diagonal_;

};

#endif /* TIFPACK_SINGLETONFILTER_HPP */
