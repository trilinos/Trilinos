/* Copyright (2001) Sandia Corportation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 *
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

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

int Amesos_EpetraRedistributor::CreateSerialMap() 
{

  if( IsSerialMapOK_ ) return 0;
  
  // MS // need to create a linear map. It will be used for:
  // MS // - bring rhs and lhs to proc 0 is KeepMatrixDistributed_ = true
  // MS // - bring rhs, lhs and matrix to proc 0 otherwise.
  // MS // NOTE: one proc case brings to IsLocal_ == 1 (The use can fix this
  // MS //       value is required). If IsLocal_ == 1, don't build SerialMap_
  
  if( IsLocal_ == false ) {
    
    int NumMyRowsSerial = 0 ;
    if (RowA()->Comm().MyPID()==0) NumMyRowsSerial = NumGlobalRows();
    if( MatrixType() != AMESOS_VBR_MATRIX ) {
      SerialMap_ = new Epetra_Map( -1, NumMyRowsSerial, 0, RowA()->Comm() );
      assert( SerialMap_ != NULL );
    } else {
      SerialBlockMap_ = new Epetra_BlockMap(-1,NumMyRowsSerial,NumPDEEqns(),0,RowA()->Comm());
      assert( SerialBlockMap_ );
    }
  }
  
  IsSerialMapOK_ = true;
  
  return 0;
}

int Amesos_EpetraRedistributor::CreateImportAndExport() 
{
  
  if( IsImportAndExportOK_ ) return 0;
  if( IsSerialMapOK_ == false ) CreateSerialMap();
  
  if( IsLocal_ == false ) {
    if( MatrixType() != AMESOS_VBR_MATRIX ) {
      // kludge, don't understant very well
      ImportToProcZero_   = new Epetra_Import(*SerialMap_, RowA()->RowMatrixRowMap());
      ImportFromProcZero_ = new Epetra_Import(RowA()->RowMatrixRowMap(), *SerialMap_);
    } else {
      ImportToProcZero_   = new Epetra_Import(*SerialBlockMap_, VbrA()->RowMap());
      ImportFromProcZero_ = new Epetra_Import(VbrA()->RowMap(), *SerialBlockMap_);
    }
    
  }

  bool IsImportAndExportOK_ = true;
  
  return 0;
  
}

Amesos_EpetraRedistributor::Amesos_EpetraRedistributor(const Epetra_LinearProblem * Problem) :
  SerialMap_(0),
  SerialBlockMap_(0),
  SerialRHS_(0),
  IsLocal_(false),
  IsSerialMapOK_(false),
  IsImportAndExportOK_(false),
  ImportToProcZero_(0),
  ImportFromProcZero_(0),
  Amesos_EpetraInterface(Problem)
{
  if( Comm().NumProc() == 1 )  SetIsLocal(true);
}

Amesos_EpetraRedistributor::~Amesos_EpetraRedistributor() 
{
  
  if( SerialMap_ ) delete SerialMap_ ;
  if( SerialBlockMap_ ) delete SerialBlockMap_ ; 
  if( ImportToProcZero_ ) delete ImportToProcZero_;
  if( ImportFromProcZero_ ) delete ImportFromProcZero_;
  if( SerialRHS_ ) delete SerialRHS_;
  
}

int Amesos_EpetraRedistributor::CreateSerialRHS(int nrhs)
{

  if( IsLocal_ ) return 0;
  
  if( SerialRHS_ != NULL ) 
    if( nrhs != SerialRHS_->NumVectors() ) {
      delete SerialRHS_; SerialRHS_ = NULL;
    }
  
  if( SerialRHS_ == NULL ) {
    if( MatrixType() == AMESOS_VBR_MATRIX )
      SerialRHS_ = new Epetra_MultiVector( *SerialBlockMap_, nrhs ) ;
    else
      SerialRHS_ = new Epetra_MultiVector( *SerialMap_, nrhs ) ;
    
    assert( SerialRHS_ != 0 ) ;
  }

  return 0;
}
