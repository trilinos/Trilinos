// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "Epetra_Time.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Amesos_EpetraInterface.h"
#include "Amesos_EpetraRedistributor.h"

Amesos_EpetraRedistributor::Amesos_EpetraRedistributor(const Epetra_LinearProblem & prob,
						       const AMESOS::Parameter::List &ParameterList) :
  Amesos_EpetraInterface(prob,ParameterList),
  TargetMap_(0),
  TargetBlockMap_(0),
  TargetRHS_(0),
  IsTargetMapOK_(false),
  IsImportAndExportOK_(false),
  ImportToTarget_(0),
  ImportFromTarget_(0),
  NumTargetProcs_(-1)
{
  // do nothing
}

int Amesos_EpetraRedistributor::SetRedistributor(const int NumTargetProcs) 
{
  NumTargetProcs_ = NumTargetProcs;

  if( TargetMap_ ) {
    delete TargetMap_ ; TargetMap_ = 0;
  }  
  if( TargetBlockMap_ ) {
    delete TargetBlockMap_ ; TargetBlockMap_ = 0;
  }  
  if( ImportToTarget_ ) {
    delete ImportToTarget_; ImportToTarget_ = 0;
  }  
  if( ImportFromTarget_ ) {
    delete ImportFromTarget_; ImportFromTarget_ = 0;
  }  
  if( TargetRHS_ ) {
    delete TargetRHS_; TargetRHS_ = 0;
  }
    
  EPETRA_CHK_ERR(CreateTargetMap());
  EPETRA_CHK_ERR(CreateImportAndExport());
  return 0;
}  

int Amesos_EpetraRedistributor::CreateTargetMap() 
{
  
  if( RowA() == NULL ) return -1;
  
  if( NumTargetProcs_<1 || NumTargetProcs_ >  RowA()->Comm().NumProc() ) {
    NumTargetProcs_ = 1;
  }
    
  //  if( IsTargetMapOK_ ) return 0;
  
  // MS // need to create a linear map. It will be used for:
  // MS // - bring rhs and lhs to proc 0 is KeepMatrixDistributed_ = true
  // MS // - bring rhs, lhs and matrix to proc 0 otherwise.
  // MS // NOTE: one proc case brings to IsLocal_ == 1 (The use can fix this
  // MS //       value is required). If IsLocal_ == 1, don't build TargetMap_
  
  if( IsLocal() == false ) {
    
    int mod = ( NumGlobalRows()) % (NumTargetProcs_);
    int NumMyRowsTarget = 0 ;
    if(  RowA()->Comm().MyPID() < NumTargetProcs_ ) {
      NumMyRowsTarget =  NumGlobalRows()/NumTargetProcs_;
    }
    if(  RowA()->Comm().MyPID() == 0 ) NumMyRowsTarget += mod;
    
    if(  MatrixType() != AMESOS_VBR_MATRIX ) {
      TargetMap_ = new Epetra_Map( -1, NumMyRowsTarget, 0,  RowA()->Comm() );
      assert( TargetMap_ != NULL );
    } else {
      TargetBlockMap_ = new Epetra_BlockMap(-1,NumMyRowsTarget, NumPDEEqns(),0,RowA()->Comm());
      assert( TargetBlockMap_ );
    }
  }
  
  IsTargetMapOK_ = true;
  
  return 0;
}

int Amesos_EpetraRedistributor::CreateImportAndExport() 
{
  
  //  if( IsImportAndExportOK_ ) return 0;
  //  if( IsTargetMapOK_ == false ) CreateTargetMap();
  
  if(  IsLocal()  == false ) {
    if( MatrixType() != AMESOS_VBR_MATRIX ) {
      // kludge, don't understant very well
      ImportToTarget_   = new Epetra_Import(*TargetMap_,  RowA()->RowMatrixRowMap());
      ImportFromTarget_ = new Epetra_Import( RowA()->RowMatrixRowMap(), *TargetMap_);
    } else {
      ImportToTarget_   = new Epetra_Import(*TargetBlockMap_,  VbrA()->RowMap());
      ImportFromTarget_ = new Epetra_Import( VbrA()->RowMap(), *TargetBlockMap_);
    }
    
  }

  bool IsImportAndExportOK_ = true;
  
  return 0;
  
}

Amesos_EpetraRedistributor::~Amesos_EpetraRedistributor() 
{
  if( TargetMap_ ) {
    delete TargetMap_ ; TargetMap_ = 0;
  }  
  if( TargetBlockMap_ ) {
    delete TargetBlockMap_ ; TargetBlockMap_ = 0;
  }  
  if( ImportToTarget_ ) {
    delete ImportToTarget_; ImportToTarget_ = 0;
  }  
  if( ImportFromTarget_ ) {
    delete ImportFromTarget_; ImportFromTarget_ = 0;
  }  
  if( TargetRHS_ ) {
    delete TargetRHS_; TargetRHS_ = 0;
  }
  
}

int Amesos_EpetraRedistributor::CreateTargetRHS(int nrhs)
{

  if( IsLocal() ) return 0;
  
  if( TargetRHS_ != NULL ) 
    if( nrhs != TargetRHS_->NumVectors() ) {
      delete TargetRHS_; TargetRHS_ = NULL;
    }
  
  if( TargetRHS_ == NULL ) {
    if( MatrixType() == AMESOS_VBR_MATRIX )
      TargetRHS_ = new Epetra_MultiVector( *TargetBlockMap_, nrhs ) ;
    else
      TargetRHS_ = new Epetra_MultiVector( *TargetMap_, nrhs ) ;
    
    assert( TargetRHS_ != 0 ) ;
  }

  return 0;
}
