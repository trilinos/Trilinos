#ifndef IFPACK_DROPFILTER_H
#define IFPACK_DROPFILTER_H

#include "Ifpack_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
#include <vector>
class Epetra_Comm;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Import;
class Epetra_BlockMap;
using namespace std;

//! Ifpack_DropFilter: Filter based on matrix entries.
//
class Ifpack_DropFilter : public virtual Epetra_RowMatrix {

public:
  //! Constructor.
  Ifpack_DropFilter(Epetra_RowMatrix* Matrix,
		    double DropTop, double MinDiagValue = 0.0,
		    double AddToDiag = 0.0);

  //! Destructor.
  virtual ~Ifpack_DropFilter();

  //! Returns the number of entries in MyRow.
  virtual inline int NumMyRowEntries(int MyRow, int & NumEntries) const
  {
    return(A_.NumMyRowEntries(MyRow, NumEntries));
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

  virtual int InvRowSums(Epetra_Vector& x) const
  {
    return(A_.InvRowSums(x));
  }

  virtual int LeftScale(const Epetra_Vector& x)
  {
    return(A_.LeftScale(x));
  }

  virtual int InvColSums(Epetra_Vector& x) const
  {
    return(A_.InvColSums(x));
  }

  virtual int RightScale(const Epetra_Vector& x) 
  {
    return(A_.RightScale(x));
  }

  virtual bool Filled() const
  {
    return(A_.Filled());
  }

  virtual double NormInf() const
  {
    return(NormInf_);
  }

  virtual double NormOne() const
  {
    return(NormOne_);
  }

  virtual int NumGlobalNonzeros() const
  {
    return(A_.NumGlobalNonzeros());
  }

  virtual int NumGlobalRows() const
  {
    return(A_.NumGlobalRows());
  }

  virtual int NumGlobalCols() const
  {
    return(A_.NumGlobalCols());
  }

  virtual int NumGlobalDiagonals() const
  {
    return(A_.NumGlobalDiagonals());
  }

  virtual int NumMyNonzeros() const
  {
    return(A_.NumMyNonzeros());
  }

  virtual int NumMyRows() const
  {
    return(A_.NumMyRows());
  }

  virtual int NumMyCols() const
  {
    return(A_.NumMyCols());
  }

  virtual int NumMyDiagonals() const
  {
    return(A_.NumMyDiagonals());
  }

  virtual bool LowerTriangular() const
  {
    return(LowerTriangular_);
  }

  virtual bool UpperTriangular() const
  {
    return(UpperTriangular_);
  }

  virtual const Epetra_Map & RowMatrixRowMap() const
  {
    return(A_.RowMatrixRowMap());
  }

  virtual const Epetra_Map & RowMatrixColMap() const
  {
    return(A_.RowMatrixColMap());
  }

  virtual const Epetra_Import * RowMatrixImporter() const
  {
    return(A_.RowMatrixImporter());
  }

#ifdef FIXME
  int SetOwnership(bool ownership)
  {
    return(A_.SetOwnership(ownership));
  }
#endif

  int SetUseTranspose(bool UseTranspose)
  {
    return(A_.SetUseTranspose(UseTranspose));
  }

  bool UseTranspose() const 
  {
    return(A_.UseTranspose());
  }

  bool HasNormInf() const
  {
    return(true);
  }

  const Epetra_Comm & Comm() const
  {
    return(A_.Comm());
  }

  const Epetra_Map & OperatorDomainMap() const 
  {
    return(A_.OperatorDomainMap());
  }

  const Epetra_Map & OperatorRangeMap() const 
  {
    return(A_.OperatorRangeMap());
  }

  const Epetra_BlockMap& Map() const 
  {
    return(A_.Map());
  }

  char* Label() const{
    return((char*)Label_);
  }

private:

  //! Pointer to the matrix to be preconditioned.
  Epetra_RowMatrix& A_;
  //! Drop tolerance.
  double DropTol_;
  //! This value will be added to all diagonal entries.
  double AddToDiag_;
  //! Minimum abs of diagonal value.
  double MinDiagValue_;
  //! Maximum entries in each row.
  int MaxNumEntries_;
  
  //! If \c true, the dropped matrix is upper triangular.
  bool UpperTriangular_;
  //! If \c true, the dropped matrix is lower triangular.
  bool LowerTriangular_;
  //! Inserse of row sum for the dropped matrix.
  Epetra_Vector* InvRowSum_;
  //! Inserse of colum sum for the dropped matrix.
  Epetra_Vector* InvColSum_;
  //! Number of local nonzeros for the dropped matrix.
  int NumMyNonzeros_;
  //! Number of global nonzeros for the dropped matrix.
  int NumGlobalNonzeros_;
  //! Infinite norm of the dropped matrix.
  double NormInf_;
  //! 1-norm of the dropped matrix.
  double NormOne_;

  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable vector<int> Indices_;
  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable vector<double> Values_;
  //! Label for \c this object.
  char Label_[80];

};


#endif /* IFPACK_DROPFILTER_H */
