#ifndef IFPACK_SPARSITYFILTER_H
#define IFPACK_SPARSITYFILTER_H

#include "Ifpack_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_RefCountPtr.hpp"

class Epetra_Comm;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Import;
class Epetra_BlockMap;

//! Ifpack_SparsityFilter: a class to drop based on sparsity.

class Ifpack_SparsityFilter : public virtual Epetra_RowMatrix {

public:
  Ifpack_SparsityFilter(const Teuchos::RefCountPtr<Epetra_RowMatrix>& Matrix,
			int AllowedNumEntries,
			int AllowedBandwidth = -1);

  virtual ~Ifpack_SparsityFilter() {};

  virtual inline int NumMyRowEntries(int MyRow, int & NumEntries) const
  {
    NumEntries = NumEntries_[MyRow];
    return(0);
  }

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

  virtual int InvRowSums(Epetra_Vector& x) const
  {
    return(-98);
  }

  virtual int LeftScale(const Epetra_Vector& x)
  {
    return(-98);
  }

  virtual int InvColSums(Epetra_Vector& x) const
  {
    return(-98);
  }

  virtual int RightScale(const Epetra_Vector& x) 
  {
    return(A_->RightScale(x));
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

  virtual const Epetra_Map & RowMatrixRowMap() const
  {
    return(A_->RowMatrixRowMap());
  }

  virtual const Epetra_Map & RowMatrixColMap() const
  {
    return(A_->RowMatrixColMap());
  }

  virtual const Epetra_Import * RowMatrixImporter() const
  {
    return(A_->RowMatrixImporter());
  }

  int SetUseTranspose(bool UseTranspose)
  {
    return(A_->SetUseTranspose(UseTranspose));
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
    return(A_->OperatorDomainMap());
  }

  const Epetra_Map & OperatorRangeMap() const 
  {
    return(A_->OperatorRangeMap());
  }

  const Epetra_BlockMap& Map() const 
  {
    return(A_->Map());
  }

  const char* Label() const{
    return(Label_);
  }

private:

  //! Pointer to the matrix to be preconditioned.
  Teuchos::RefCountPtr<Epetra_RowMatrix> A_;
  //! Maximum entries in each row.
  int MaxNumEntries_;
  int MaxNumEntriesA_;

  //! Maximum allowed bandwidth.
  int AllowedBandwidth_;
  //! Maximum allowed entries per row.
  int AllowedEntries_;
  
  //! Number of nonzeros for the dropped matrix.
  int NumNonzeros_;

  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable vector<int> Indices_;
  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable vector<double> Values_;
  //! Label for \c this object.
  char Label_[80];

  int NumRows_;
  vector<int> NumEntries_;

};


#endif /* IFPACK_SPARSITYFILTER_H */
