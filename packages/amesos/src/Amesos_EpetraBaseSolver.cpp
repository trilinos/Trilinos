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
#include "Amesos_EpetraBaseSolver.h"

//=============================================================================

Amesos_EpetraBaseSolver::Amesos_EpetraBaseSolver(const Epetra_LinearProblem & Problem):
  RowIndices_(0),
  ColIndices_(0),
  Values_(0),
  MaxNumEntries_(0),
  NumMyRows_(0),
  NumMyBlockRows_(0),
  NumGlobalRows_(0),
  NumMyNonzeros_(0),
  NumGlobalNonzeros_(0),
  MyGlobalElements_(0),
  IsSetInterfaceOK_(false),
  IsLocal_(false),
  Entries_(0),
  BlockIndices_(0),
  NumPDEEqns_(1),
  Problem_(Problem),
  MatrixProperty_(AMESOS_UNSYM)
{
}
 
//=============================================================================

Amesos_EpetraBaseSolver::~Amesos_EpetraBaseSolver()
{
  if( RowIndices_ ) delete RowIndices_;
  if( ColIndices_ ) delete ColIndices_;
  if( Values_ ) delete Values_;
}

//=============================================================================

int Amesos_EpetraBaseSolver::SetInterface(Epetra_RowMatrix * RowA)
{

  //  if( IsSetInterfaceOK_ ) return 0;

  // we could do something smarter, like:
  // if all the nodes are on process zero, solve using this
  // process only ????
  
  if( RowA == 0 ) return -3;
  
  if( Comm().NumProc() == 1 )  SetIsLocal(true);
  
  RowA_ = RowA;       MatrixType_ = AMESOS_ROW_MATRIX;
  assert( RowA_ != NULL );
  
  CrsA_ = dynamic_cast<Epetra_CrsMatrix*>(RowA_) ; 
  if( CrsA_ != NULL ) MatrixType_ = AMESOS_CRS_MATRIX;

  VbrA_ = dynamic_cast<Epetra_VbrMatrix*>(RowA_) ; 
  if( VbrA_ != NULL )  MatrixType_ = AMESOS_VBR_MATRIX;

  // MS // retrive information about the distributed matrix
  // MS // Some of the information are not really needed, but will
  // MS // improve performances in allocating memory for Col, Row and Val.

  // MS // Notation: I use `rows' for matrix, and `elements' for map

  NumMyRows_ =  RowA_->NumMyRows();
  NumGlobalRows_ =  RowA_->NumGlobalRows();
  NumMyNonzeros_ = RowA_->NumMyNonzeros();
  
  NumGlobalNonzeros_ = RowA_->NumGlobalNonzeros();
 
  // matrix must be square
  assert( NumGlobalRows_ == RowA_->NumGlobalCols() );  
  
  MaxNumEntries_ = RowA_->MaxNumEntries();
  
  assert(MaxNumEntries_>0);
  
  if( MatrixType_ != AMESOS_VBR_MATRIX ) {
    
    MyGlobalElements_ = RowA_->RowMatrixRowMap().MyGlobalElements();
    NumMyBlockRows_ = NumMyRows_;
    
  } else {
    
    if( VbrA_->RowMap().ConstantElementSize () == false ) {
      // non-constant element size can indeed be supported,
      // but in a less efficient way.
      if( RowA_->Comm().MyPID() == 0 ) {
	cerr << "Amesos ERROR : VBR matrices with non-constant element size\n"
	     << "Amesos ERROR : are not supported\n";
	exit( EXIT_FAILURE );
      }
    }

    NumMyBlockRows_ = VbrA_->NumMyBlockRows();
    
    NumPDEEqns_ = VbrA_->RowMap().NumMyPoints()/NumMyBlockRows_;
    MyGlobalElements_ = VbrA_->RowMap().MyGlobalElements();
  }

  if( RowIndices_ ) delete RowIndices_;
  if( ColIndices_ ) delete ColIndices_;
  if( Values_ ) delete Values_;
  
  RowIndices_ = new Epetra_IntSerialDenseVector(MaxNumEntries_);
  ColIndices_ = new Epetra_IntSerialDenseVector(MaxNumEntries_);
  Values_     = new Epetra_SerialDenseVector(MaxNumEntries_);
  
  IsSetInterfaceOK_ = true;

  return 0;
}

//=============================================================================

int Amesos_EpetraBaseSolver::GetRow(int BlockRow, int & NumIndices,
				    int * & RowIndices,
				    int * & ColIndices, double * & MatrixValues)
{

  int * ColIndices1;
  
  if( IsSetInterfaceOK_ == false ) return -1;
  
  switch( MatrixType_ ) {

  case AMESOS_CRS_MATRIX: {
    
    CrsA_->ExtractMyRowView(BlockRow, NumIndices, MatrixValues, ColIndices1);
    for( int i=0 ; i<NumIndices ; ++i ) {
      RowIndices[i] = MyGlobalElements_[BlockRow];
      ColIndices[i] = CrsA_->ColMap().GID(ColIndices1[i]);
    }
    
  } break;

  case AMESOS_VBR_MATRIX: {

    int ierr;
    int NumBlockEntries;
    int * BlockIndices;
    int count = 0;
    int GlobalRow = MyGlobalElements_[BlockRow]*NumPDEEqns_;
      
    int RowDim;
    ierr = VbrA_->ExtractMyBlockRowView(BlockRow,RowDim,
					NumBlockEntries, BlockIndices, Entries_);
    assert(ierr==0);
    assert(RowDim==NumPDEEqns_);
      
    // insert all the nonzeros in the block
    for( int j=0 ; j<NumBlockEntries ; ++j ) {
	
      int GlobalCol = VbrA_->ColMap().GID(BlockIndices[j])*NumPDEEqns_;
      for( int k=0 ; k<NumPDEEqns_ ; ++k ) {
	for( int kk=0 ; kk<NumPDEEqns_ ; ++kk ) {
	  if( (*Entries_[j])(k,kk) ) {
	    RowIndices[count] = GlobalRow + k;
	    ColIndices[count] = GlobalCol + kk;
	    MatrixValues[count] = (*Entries_[j])(k,kk);
	    count++;
	  }
	}
      }    
    }

    NumIndices = count;

    } break;

  case AMESOS_ROW_MATRIX:  {

    // in this case BlockRow is a Row
    RowA_->ExtractMyRowCopy(BlockRow,MaxNumEntries_,NumIndices,MatrixValues,ColIndices1);
    for( int i=0 ; i<NumIndices ; ++i ) {
      RowIndices[i] = MyGlobalElements_[BlockRow];
      ColIndices[i] = RowA_->RowMatrixColMap().GID(ColIndices1[i]);
    }
    
  } break;
  
  default:
    return -1;
    
  }
  
  return 0;
  
} 

//=============================================================================

int Amesos_EpetraBaseSolver::IndexBase() const
{
  switch( MatrixType_ ) {

  case AMESOS_CRS_MATRIX:
    return( CrsA()->IndexBase() );
    break;
  case AMESOS_VBR_MATRIX:
    return( VbrA()->IndexBase() );
    break;
  }
  return 0;  
}      

//=============================================================================

bool Amesos_EpetraBaseSolver::MatrixShapeOK() const
{ 
  bool OK = true;

  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}
