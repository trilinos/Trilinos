// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef ML_PYMATRIX_H
#define ML_PYMATRIX_H

#include "ml_common.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_Operator.h"
#include "Teuchos_RCP.hpp"
#include "Epetra_Map.h"
#include "Epetra_FECrsMatrix.h"

using namespace MLAPI;

namespace PyTrilinos
{

class PyMatrix : public MLAPI::Operator
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
      ColMap_ = Teuchos::rcp(new Epetra_Map(-1, ColSpace.GetNumMyElements(), 0,
					    GetEpetra_Comm()));
    }
    else
    {
      Teuchos::RCP<Epetra_IntSerialDenseVector> IMap = ColSpace.GetRCPMyGlobalElements();
      ColMap_ = Teuchos::rcp(new Epetra_Map(-1, IMap->Length(), IMap->Values(), 0, 
					    GetEpetra_Comm()));
    }

    // range map
    if (RowSpace.IsLinear())
    {
      RowMap_ = Teuchos::rcp(new Epetra_Map(-1, RowSpace.GetNumMyElements(), 0,
					    GetEpetra_Comm()));
    }
    else
    {
      Teuchos::RCP<Epetra_IntSerialDenseVector> IMap = RowSpace.GetRCPMyGlobalElements();
      RowMap_ = Teuchos::rcp(new Epetra_Map(-1, IMap->Length(), IMap->Values(), 0, 
					    GetEpetra_Comm()));
    }

    // I suppose that RowMap == RowMatrixRowMap
    Matrix_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, *(RowMap_.get()), 0));
    
    // It seems safer to me to add zero diagonal elements
    for (int i = 0 ; i < ColMap_->NumMyElements() ; ++i)
    {
      int GlobalCol = ColMap_->GID(i);
      SetElement(GlobalCol, GlobalCol, 0.0);
    }
  }

  //! Destructor.
  ~PyMatrix() { }

  const Space & GetRowSpace() const
  {
    return RowSpace_;
  }

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
  Teuchos::RCP<Epetra_Map> ColMap_;
  Teuchos::RCP<Epetra_Map> RowMap_;
  Teuchos::RCP<Epetra_FECrsMatrix> Matrix_;

};  // PyMatrix

}  // Namespace PyTrilinos

#endif 
