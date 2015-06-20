// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef ML_PYMATRIX_HPP
#define ML_PYMATRIX_HPP

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
