#ifndef ML_PYMATRIX_H
#define ML_PYMATRIX_H

#include "ml_common.h"
#include "MLAPI_Operator.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Map.h"
#include "Epetra_FECrsMatrix.h"

using namespace std;
using namespace Teuchos;
using namespace MLAPI;

class PyMatrix : public Operator
{

public:
  //@{ \name Constructors and destructors.

  //! Constructor with given already computed ML_Operator pointer.
  PyMatrix(const Space& RowSpace, const Space& ColSpace)
  {
    ColSpace_ = ColSpace;
    RowSpace_ = RowSpace;

    // domain map
    if (ColSpace.IsLinear())
    {
      ColMap_ = rcp(new Epetra_Map(-1, ColSpace.GetNumMyElements(), 0,
                                      GetEpetra_Comm()));
    }
    else
    {
      RefCountPtr<Epetra_IntSerialDenseVector> IMap = ColSpace.GetRCPMyGlobalElements();
      ColMap_ = rcp(new Epetra_Map(-1, IMap->Length(), IMap->Values(), 0, 
                                      GetEpetra_Comm()));
    }

    // range map
    if (RowSpace.IsLinear())
    {
      RowMap_ = rcp(new Epetra_Map(-1, RowSpace.GetNumMyElements(), 0,
                                      GetEpetra_Comm()));
    }
    else
    {
      RefCountPtr<Epetra_IntSerialDenseVector> IMap = RowSpace.GetRCPMyGlobalElements();
      RowMap_ = rcp(new Epetra_Map(-1, IMap->Length(), IMap->Values(), 0, 
                                      GetEpetra_Comm()));
    }

    // I suppose that RowMap == RowMatrixRowMap
    Matrix_ = rcp(new Epetra_FECrsMatrix(Copy, *(RowMap_.get()), 0));
    
    // It seems safer to me to add zero diagonal elements
    for (int i = 0 ; i < ColMap_->NumMyElements() ; ++i)
    {
      int GlobalCol = ColMap_->GID(i);
      SetElement(GlobalCol, GlobalCol, 0.0);
    }
  }

  //! Destructor.
  ~PyMatrix() { }

  void SetElement(int row, int col, double value)
  {
    if (Matrix_->Filled())
    {
      int MyRow = RowMap_->LID(row);
      int MyCol = ColMap_->LID(col);
      Matrix_->ReplaceMyValues(MyRow, 1, &value, &MyCol);
    }
    else
    {
    if (Matrix_->ReplaceGlobalValues(1, &row, 1, &col, &value) > 0)
      Matrix_->InsertGlobalValues(1, &row, 1, &col, &value);
    }
  }
  
  void FillComplete()
  {
    if (Matrix_->GlobalAssemble(false))
      ML_THROW("Error in GlobalAssemble()", -1);

    if (Matrix_->FillComplete(*(ColMap_.get()), *(RowMap_.get())))
      ML_THROW("Error in FillComplete()", -1);

    Reshape(ColSpace_,RowSpace_, Matrix_.get(), false);
  }

  Epetra_FECrsMatrix* GetMatrix()
  {
    return(Matrix_.get());
  }

private:

  //! Col space.
  Space ColSpace_;
  //! Row space.
  Space RowSpace_;
  RefCountPtr<Epetra_Map> ColMap_;
  RefCountPtr<Epetra_Map> RowMap_;
  RefCountPtr<Epetra_FECrsMatrix> Matrix_;

}; // PyMatrix

#endif 
