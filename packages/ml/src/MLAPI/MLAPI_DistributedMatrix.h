#ifndef MLAPI_DISTRIBUTEDMATRIX_H
#define MLAPI_DISTRIBUTEDMATRIX_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_DistributedMatrix.h

\brief MLAPI wrapper for Epetra_FECrsMatrix, which allows MATLAB-like syntax.

\author Marzio Sala, D-INFK/ETHZ.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_common.h"

#include "MLAPI_Error.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Space.h"
#include "MLAPI_BaseLinearCombination.h"
#include "Epetra_Map.h"
#include "Epetra_FECrsMatrix.h"

namespace MLAPI {

class DistributedMatrix : public Epetra_RowMatrix, public Operator {

public:

  DistributedMatrix(const Space& RowSpace, const Space& ColSpace)
  {
    FillCompleted_ = false;

    RowSpace_ = RowSpace;
    ColSpace_ = ColSpace;

    int locNumMyRows = RowSpace_.GetNumMyElements();
    int locNumMyCols = ColSpace_.GetNumMyElements();

    // FIXME: add MyGlobalElements()
    RangeMap_ = new Epetra_Map(-1, locNumMyRows, 0, GetEpetra_Comm());
    DomainMap_ = new Epetra_Map(-1, locNumMyCols, 0, GetEpetra_Comm());

    Matrix_ = new Epetra_FECrsMatrix(Copy, *RangeMap_, 0);
  }

  virtual ~DistributedMatrix() {
    //    delete Matrix_;
    delete RangeMap_;
    delete DomainMap_;
  }


  virtual int NumMyRowEntries(int MyRow, int& NumEntries) const
  {
    ML_RETURN(Matrix_->NumMyRowEntries(MyRow, NumEntries));
  }

  virtual int MaxNumEntries() const
  {
    return(Matrix_->MaxNumEntries());
  }

  virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries,
                               double* Values, int* Indices) const
  {
    ML_RETURN(Matrix_->ExtractMyRowCopy(MyRow, Length, NumEntries,
                                        Values, Indices));
  }

  virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
  {
    ML_RETURN(Matrix_->ExtractDiagonalCopy(Diagonal));
  }

  virtual int Multiply(bool TransA, const Epetra_MultiVector& X,
                       Epetra_MultiVector& Y) const
  {
    ML_RETURN(Matrix_->Multiply(TransA, X, Y));
  }

  virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X,
                    Epetra_MultiVector& Y) const
  {
    ML_RETURN(Matrix_->Solve(Upper, Trans, UnitDiagonal, X, Y));
  }

  virtual int InvRowSums(Epetra_Vector& x) const
  {
    ML_RETURN(Matrix_->InvRowSums(x));
  }

  virtual int LeftScale(const Epetra_Vector& x)
  {
    ML_RETURN(Matrix_->LeftScale(x));
  }

  virtual int InvColSums(Epetra_Vector& x) const
  {
    ML_RETURN(Matrix_->InvColSums(x));
  }

  virtual int RightScale(const Epetra_Vector& x)
  {
    ML_RETURN(Matrix_->RightScale(x));
  }

  virtual bool Filled() const
  {
    return(Matrix_->Filled());
  }

  virtual double NormInf() const
  {
    return(Matrix_->NormInf());
  }

  virtual double NormOne() const
  {
    return(Matrix_->NormOne());
  }

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  virtual int NumGlobalNonzeros() const
  {
    return(Matrix_->NumGlobalNonzeros());
  }

  virtual int NumGlobalRows() const
  {
    return(Matrix_->NumGlobalRows());
  }

  virtual int NumGlobalCols() const
  {
    return(Matrix_->NumGlobalCols());
  }

  virtual int NumGlobalDiagonals() const
  {
    return(Matrix_->NumGlobalDiagonals());
  }
#endif

  virtual long long NumGlobalNonzeros64() const
  {
    return(Matrix_->NumGlobalNonzeros64());
  }

  virtual long long NumGlobalRows64() const
  {
    return(Matrix_->NumGlobalRows64());
  }

  virtual long long NumGlobalCols64() const
  {
    return(Matrix_->NumGlobalCols64());
  }

  virtual long long NumGlobalDiagonals64() const
  {
    return(Matrix_->NumGlobalDiagonals64());
  }

  virtual int NumMyNonzeros() const
  {
    return(Matrix_->NumMyNonzeros());
  }

  virtual int NumMyRows() const
  {
    return(Matrix_->NumMyRows());
  }

  virtual int NumMyCols() const
  {
    return(Matrix_->NumMyCols());
  }

  virtual int NumMyDiagonals() const
  {
    return(Matrix_->NumMyDiagonals());
  }

  virtual bool LowerTriangular() const
  {
    return(Matrix_->LowerTriangular());
  }

  virtual bool UpperTriangular() const
  {
    return(Matrix_->UpperTriangular());
  }

  virtual const Epetra_Map & RowMatrixRowMap() const
  {
    return(Matrix_->RowMatrixRowMap());
  }

  virtual const Epetra_Map & RowMatrixColMap() const
  {
    return(Matrix_->RowMatrixColMap());
  }

  virtual const Epetra_Import * RowMatrixImporter() const
  {
    return(Matrix_->RowMatrixImporter());
  }

  virtual const Epetra_Map& OperatorDomainMap() const
  {
    return(Matrix_->OperatorDomainMap());
  }

  virtual const Epetra_Map& OperatorRangeMap() const
  {
    return(Matrix_->OperatorRangeMap());
  }

  virtual const Epetra_Map& Map() const
  {
    return(Matrix_->RowMatrixRowMap());
  }

  //@}

  virtual int SetUseTranspose(bool what)
  {
    return(Matrix_->SetUseTranspose(what));
  }

  int Apply(const MultiVector& X, MultiVector& Y) const
  {
    return(Operator::Apply(X, Y));
  }

  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
  {
    if (!IsFillCompleted())
      ML_THROW("Matrix not yet FillComplete()'d", -1);

    return(Matrix_->Apply(X, Y));
  }

  virtual int ApplyInverse(const Epetra_MultiVector& X,
                           Epetra_MultiVector& Y) const
  {
    if (!IsFillCompleted())
      ML_THROW("Matrix not yet FillComplete()'d", -1);

    return(Matrix_->ApplyInverse(X, Y));
  }

  virtual const char* Label() const
  {
    return(Matrix_->Label());
  }

  virtual bool UseTranspose() const
  {
    return(Matrix_->UseTranspose());
  }

  virtual bool HasNormInf() const
  {
    return(Matrix_->HasNormInf());
  }

  virtual const Epetra_Comm& Comm() const
  {
    return(Matrix_->Comm());
  }

  std::ostream& Print(std::ostream& os, const bool verbose = true) const;

  Space GetDomainSpace() const
  {
    return(ColSpace_);
  }

  Space GetRangeSpace() const
  {
    return(RowSpace_);
  }

  inline double& operator()(const int GRID, const int GCID)
  {
    if (IsFillCompleted())
    {
      ML_THROW("Matrix already FillCompleted()'d", -1);
    }
    else
    {
      rows_.push_back(GRID);
      cols_.push_back(GCID);
      vals_.push_back(0.0);
      return(vals_[vals_.size() - 1]);
    }
  }

  inline void ReplaceElement(const int GRID, const int GCID,
                             const double value)
  {
    if (!IsFillCompleted())
      ML_THROW("Matrix not FillCompleted()'d yet", -1);

    int LRID = RangeMap_->LID(GRID);
    int LCID = Matrix_->ColMap().LID(GCID);
    if (Matrix_->ReplaceMyValues((int)LRID, 1, (double*)&value,
                                 (int*)&LCID) < 0)
      ML_THROW("Can only replace locally owned elements", -1);
  }

  void FillComplete()
  {
    // populate the matrix here
    for (unsigned int i = 0 ; i < vals_.size() ; ++i)
    {
      int    GRID  = rows_[i];
      int    GCID  = cols_[i];
      double value = vals_[i];
      if (Matrix_->ReplaceGlobalValues(1, &GRID, 1, &GCID, &value) > 0)
        Matrix_->InsertGlobalValues(1, &GRID, 1, &GCID, &value);
    }
    rows_.resize(0);
    cols_.resize(0);
    vals_.resize(0);

    // freeze the matrix
    if (Matrix_->GlobalAssemble())
      ML_THROW("Error in GlobalAssemble()", -1);

    if (Matrix_->FillComplete(*DomainMap_, *RangeMap_))
      ML_THROW("Error in FillComplete()", -1);

    FillCompleted_ = true;

    Reshape(ColSpace_, RowSpace_, Matrix_, true);
  }

  bool IsFillCompleted() const
  {
    return(FillCompleted_);
  }

private:

  DistributedMatrix(const DistributedMatrix& /* rhs */)
  {
  }

  DistributedMatrix& operator=(const DistributedMatrix& /* rhs */)
  {
    return(*this);
  }

  mutable std::vector<int>    rows_;
  mutable std::vector<int>    cols_;
  mutable std::vector<double> vals_;

  Epetra_FECrsMatrix* Matrix_;
  Space ColSpace_;
  Space RowSpace_;

  Epetra_Map* DomainMap_;
  Epetra_Map* RangeMap_;

  bool FillCompleted_;

}; // class DistributedMatrix

} // namespace MLAPI

#endif // ifndef MLAPI_DISTRIBUTEDMATRIX_H
