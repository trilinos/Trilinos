#ifndef MLAPI_DISTRIBUTEDMATRIX_H
#define MLAPI_DISTRIBUTEDMATRIX_H

#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

#include "ml_include.h"
#include "ml_lapack.h"
#include "ml_comm.h"
#include "MLAPI_Error.h"
#include "MLAPI_BaseObject.h"
#include "MLAPI_Space.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include <iomanip>

namespace MLAPI {

class DistributedMatrix : public Epetra_RowMatrix, public BaseObject {

public:

  DistributedMatrix(const Space& RowSpace, const Space& ColSpace)
  {
    FillCompleted_ = false;

    RowSpace_ = RowSpace;
    ColSpace_ = ColSpace;

    int NumMyRows = RowSpace_.GetNumMyElements();
    int NumMyCols = ColSpace_.GetNumMyElements();
    
    // FIXME: add MyGlobalElements()
    RangeMap_ = new Epetra_Map(-1, NumMyRows, 0, GetEpetra_Comm());
    DomainMap_ = new Epetra_Map(-1, NumMyCols, 0, GetEpetra_Comm());

    Matrix_ = new Epetra_FECrsMatrix(Copy, *RangeMap_, 0);
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

  ostream& Print(ostream& os, bool verbose = true) const
  {
    int Length = MaxNumEntries();
    vector<double> Values(Length);
    vector<int>    Indices(Length);

    if (GetMyPID() == 0) {

      os << endl;
      os << "*** MLAPI::DistributedMatrix ***" << endl;
      os << "Label = " << GetLabel() << endl;
      os << "Number of rows    = " << GetRangeSpace().GetNumGlobalElements() << endl;
      os << "Number of columns = " << GetDomainSpace().GetNumGlobalElements() << endl;
      os << endl;
      os.width(10); os << "row ID";
      os.width(10); os << "col ID";
      os.width(30); os << "value";
      os << endl;
      os << endl;
    }

    for (int iproc = 0 ; iproc < GetNumProcs() ; ++iproc) {

      if (GetMyPID() == iproc) {

        if (IsFillCompleted()) {
          for (int i = 0 ; i < NumMyRows() ; ++i) {
            int GRID = RangeMap_->GID(i);
            double* Values;
            int* Indices;
            int NumEntries;
            Matrix_->ExtractMyRowView(i, NumEntries, Values, Indices);
            for (int j = 0 ; j < NumEntries ; ++j) {
              os.width(10); os << GRID;
              os.width(10); os << Matrix_->RowMatrixColMap().GID(Indices[j]);
              os.width(30); os << Values[j];
              os << endl;
            }
          }
        }
        else {
          for (int i = 0 ; i < NumMyRows() ; ++i) {
            int GRID = RangeMap_->GID(i);
            double* Values;
            int* Indices;
            int NumEntries;
            Matrix_->ExtractGlobalRowView(GRID, NumEntries, Values, Indices);
            for (int j = 0 ; j < NumEntries ; ++j) {
              os.width(10); os << GRID;
              os.width(10); os << Indices[j];
              os.width(30); os << Values[j];
              os << endl;
            }
          }
        }
      }
      Barrier();
    }

    if (GetMyPID() == 0)
      os << endl;

    Barrier();

    return(os);
  }

  Space GetDomainSpace() const
  {
    return(ColSpace_);
  }

  Space GetRangeSpace() const
  {
    return(RowSpace_);
  }

  void SetElement(int row, int col, double value)
  {
    if (IsFillCompleted())
      ML_THROW("Matrix already FillComplete()'d", -1);

    if (Matrix_->ReplaceGlobalValues(1, &row, 1, &col, &value) > 0)
      Matrix_->InsertGlobalValues(1, &row, 1, &col, &value);
  }
  
  void FillComplete()
  {
    
    if (Matrix_->GlobalAssemble())
      ML_THROW("Error in GlobalAssemble()", -1);

    if (Matrix_->FillComplete(*DomainMap_, *RangeMap_))
      ML_THROW("Error in FillComplete()", -1);

    FillCompleted_ = true;
  }

  bool IsFillCompleted() const 
  {
    return(FillCompleted_);
  }

private:

  DistributedMatrix(const DistributedMatrix& rhs)
  {
  }

  DistributedMatrix& operator=(const DistributedMatrix& rhs)
  {
    return(*this);
  }

  vector<vector<double> > Values_;
  vector<vector<int> >    Indices_;

  Epetra_FECrsMatrix* Matrix_;
  Space ColSpace_;
  Space RowSpace_;

  Epetra_Map* DomainMap_;
  Epetra_Map* RangeMap_;

  bool FillCompleted_;

}; // class DistributedMatrix

} // namespace MLAPI

#endif // ifdef HAVE_ML_MLAPI
#endif // ifndef MLAPI_DISTRIBUTEDMATRIX_H
