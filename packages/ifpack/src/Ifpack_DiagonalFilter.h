#ifndef IFPACK_DIAGONALFILTER_H
#define IFPACK_DIAGONALFILTER_H

#include "Ifpack_ConfigDefs.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
class Epetra_Comm;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Import;
class Epetra_BlockMap;
#include "Epetra_Time.h"

//! Ifpack_DiagonalFilter: Filter to modify the diagonal entries of a given Epetra_RowMatrix.
/*!

Ifpack_DiagonalFilter modifies the elements on the diagonal.

A typical use is as follows:
\code
Epetra_RowMatrix* A;
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
  Ifpack_DiagonalFilter(Epetra_RowMatrix* Matrix,
                        double AbsoluteThreshold,
                        double RelativeThreshold) :
    A_(*Matrix),
    AbsoluteThreshold_(AbsoluteThreshold),
    RelativeThreshold_(RelativeThreshold)
  {
    Epetra_Time Time(Comm());

    pos_.resize(NumMyRows());
    val_.resize(NumMyRows());

    vector<int> Indices(MaxNumEntries());
    vector<double> Values(MaxNumEntries());
    int NumEntries;

    for (int MyRow = 0 ; MyRow < NumMyRows() ; ++MyRow) {

      pos_[MyRow] = -1;
      val_[MyRow] = 0.0;
      int ierr = A_.ExtractMyRowCopy(MyRow, MaxNumEntries(), NumEntries, 
                                     &Values[0], &Indices[0]);
      assert (ierr == 0);

      for (int i = 0 ; i < NumEntries ; ++i) {
        if (Indices[i] == MyRow) {
          pos_[MyRow] = i;
          val_[MyRow] = Values[i] * (RelativeThreshold_ - 1) +
            AbsoluteThreshold_ * EPETRA_SGN(Values[i]);
        }
        break;
      }
    }
    cout << "TIME = " << Time.ElapsedTime() << endl;
  }
  
  //! Destructor.
  virtual ~Ifpack_DiagonalFilter() {};

  //! Returns the number of entries in MyRow.
  virtual int NumMyRowEntries(int MyRow, int& NumEntries) const
  {
    return(A_.NumMyRowEntries(MyRow, NumEntries));
  }

  //! Returns the maximum number of entries.
  virtual int MaxNumEntries() const
  {
    return(A_.MaxNumEntries());
  }

  inline virtual int ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, 
                               double* Values, int* Indices) const;

  virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
  {
    IFPACK_RETURN(A_.ExtractDiagonalCopy(Diagonal));
  }

  virtual int Multiply(bool TransA, const Epetra_MultiVector& X, 
		       Epetra_MultiVector& Y) const;

  virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, 
		    const Epetra_MultiVector& X,
		    Epetra_MultiVector& Y) const
  {
    IFPACK_CHK_ERR(-1);
  }

  virtual int Apply(const Epetra_MultiVector& X,
		    Epetra_MultiVector& Y) const
  {
    IFPACK_RETURN(Multiply(UseTranspose(),X,Y));
  }

  virtual int ApplyInverse(const Epetra_MultiVector& X,
			   Epetra_MultiVector& Y) const
  {
    IFPACK_CHK_ERR(-1);
  }

  virtual int InvRowSums(Epetra_Vector& x) const
  {
    IFPACK_CHK_ERR(-1);
  }

  virtual int LeftScale(const Epetra_Vector& x)
  {
    return(A_.LeftScale(x));
  }

  virtual int InvColSums(Epetra_Vector& x) const
  {
    IFPACK_CHK_ERR(-1);;
  }

  virtual int RightScale(const Epetra_Vector& x) 
  {
    return(A_.RightScale(x));
  }

  virtual bool Filled() const
  {
    return(A_.Filled());
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
    return(A_.LowerTriangular());
  }

  virtual bool UpperTriangular() const
  {
    return(A_.UpperTriangular());
  }

  virtual const Epetra_Map& RowMatrixRowMap() const
  {
    return(A_.RowMatrixRowMap());
  }

  virtual const Epetra_Map& RowMatrixColMap() const
  {
    return(A_.RowMatrixColMap());
  }

  virtual const Epetra_Import* RowMatrixImporter() const
  {
    return(A_.RowMatrixImporter());
  }

  int SetUseTranspose(bool UseTranspose)
  {
    return(A_.SetUseTranspose(UseTranspose));
  }

  bool UseTranspose() const 
  {
    return(A_.UseTranspose());
  }

  //! Not implemented for efficiency reasons.
  bool HasNormInf() const
  {
    return(false);
  }

  const Epetra_Comm& Comm() const
  {
    return(A_.Comm());
  }

  const Epetra_Map& OperatorDomainMap() const 
  {
    return(A_.OperatorDomainMap());
  }

  const Epetra_Map& OperatorRangeMap() const 
  {
    return(A_.OperatorRangeMap());
  }

  const Epetra_BlockMap& Map() const 
  {
    return(A_.Map());
  }

  const char* Label() const{
    return(A_.Label());
  }

private:

  //! Pointer to the matrix to be filtered
  Epetra_RowMatrix& A_;
  //! This value (times the sgn(A(i,i)) is added to the diagonal elements
  double AbsoluteThreshold_;
  //! Multiplies A(i,i) by this value.
  double RelativeThreshold_;
  //! Stores the position of the diagonal element, or -1 if not present.
  vector<int> pos_;
  //! Stores as additional diagonal contribution due to the filter.
  vector<double> val_;

};


#endif /* IFPACK_DIAGONALFILTER_H */
