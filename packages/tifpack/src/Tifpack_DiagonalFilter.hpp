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

#ifndef TIFPACK_DIAGONALFILTER_HPP
#define TIFPACK_DIAGONALFILTER_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Time.hpp"
#include "Teuchos_RefCountPtr.hpp"

class Tpetra_Comm;
class Tpetra_Map;
class Tpetra_MultiVector;
class Tpetra_Import;
class Tpetra_BlockMap;

//! Tifpack_DiagonalFilter: Filter to modify the diagonal entries of a given Tpetra_RowMatrix.
/*!

Tifpack_DiagonalFilter modifies the elements on the diagonal.

A typical use is as follows:
\code
Teuchos::RCP<Tpetra_RowMatrix> A;
// creates a matrix B such that
// B(i,i) = AbsoluteThreshold * sgn(B(i,i)) + 
//          RelativeThreshold * B(i,i)
double AbsoluteThreshold = 1e-3;
double RelativeThreshold = 1.01;

Tifpack_DiagonalFilter B(A, AbsoluteThreshold, RelativeThreshold);
\endcode

\author Michael Heroux, SNL 9214.

\data Last modified on 24-Jan-05.

*/
class Tifpack_DiagonalFilter : public virtual Tpetra_RowMatrix {

public:
  //! Constructor.
  Tifpack_DiagonalFilter(const Teuchos::RCP<Tpetra_RowMatrix>& Matrix,
                        double AbsoluteThreshold,
                        double RelativeThreshold);
  
  //! Destructor.
  virtual ~Tifpack_DiagonalFilter() {};

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

  virtual int ExtractDiagonalCopy(Tpetra_Vector & Diagonal) const
  {
    TIFPACK_RETURN(A_->ExtractDiagonalCopy(Diagonal));
  }

  virtual int Multiply(bool TransA, const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
		       Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, 
		    const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
		    Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
  {
    TIFPACK_CHK_ERR(-1);
  }

  virtual int Apply(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
		    Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
  {
    TIFPACK_RETURN(Multiply(UseTranspose(),X,Y));
  }

  virtual int ApplyInverse(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
			   Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
  {
    TIFPACK_CHK_ERR(-1);
  }

  virtual int InvRowSums(Tpetra_Vector& x) const
  {
    TIFPACK_CHK_ERR(-1);
  }

  virtual int LeftScale(const Tpetra_Vector& x)
  {
    return(A_->LeftScale(x));
  }

  virtual int InvColSums(Tpetra_Vector& x) const
  {
    TIFPACK_CHK_ERR(-1);;
  }

  virtual int RightScale(const Tpetra_Vector& x) 
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

  virtual const Tpetra_Map& RowMatrixRowMap() const
  {
    return(A_->RowMatrixRowMap());
  }

  virtual const Tpetra_Map& RowMatrixColMap() const
  {
    return(A_->RowMatrixColMap());
  }

  virtual const Tpetra_Import* RowMatrixImporter() const
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

  const Tpetra_Comm& Comm() const
  {
    return(A_->Comm());
  }

  const Tpetra_Map& OperatorDomainMap() const 
  {
    return(A_->OperatorDomainMap());
  }

  const Tpetra_Map& OperatorRangeMap() const 
  {
    return(A_->OperatorRangeMap());
  }

  const Tpetra_BlockMap& Map() const 
  {
    return(A_->Map());
  }

  const char* Label() const{
    return(A_->Label());
  }

private:

  //! Pointer to the matrix to be filtered
  Teuchos::RCP<Tpetra_RowMatrix> A_;
  //! This value (times the sgn(A(i,i)) is added to the diagonal elements
  double AbsoluteThreshold_;
  //! Multiplies A(i,i) by this value.
  double RelativeThreshold_;
  //! Stores the position of the diagonal element, or -1 if not present.
  std::vector<int> pos_;
  //! Stores as additional diagonal contribution due to the filter.
  std::vector<double> val_;

};


#endif /* TIFPACK_DIAGONALFILTER_HPP */
