// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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

#include <EpetraExt_BlockJacobi_LinearProblem.h>

#include <Epetra_VbrMatrix.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>
#include <Epetra_MultiVector.h>
#include <Epetra_LinearProblem.h>

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseSVD.h>

#include <set>

using std::vector;

namespace EpetraExt {

LinearProblem_BlockJacobi::
~LinearProblem_BlockJacobi()
{
  for( int i = 0; i < NumBlocks_; ++i )
  {
    if( SVDs_[i] ) delete SVDs_[i];
    else if( Inverses_[i] ) delete Inverses_[i];

    if( RHSBlocks_[i] ) delete RHSBlocks_[i];
  }

  if( NewProblem_ ) delete NewProblem_;
  if( NewMatrix_ ) delete NewMatrix_;
}

LinearProblem_BlockJacobi::NewTypeRef
LinearProblem_BlockJacobi::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  Epetra_VbrMatrix * OrigMatrix = dynamic_cast<Epetra_VbrMatrix*>( orig.GetMatrix() );

  if( OrigMatrix->RowMap().DistributedGlobal() )
  { cout << "FAIL for Global!\n"; abort(); }
  if( OrigMatrix->IndicesAreGlobal() )
  { cout << "FAIL for Global Indices!\n"; abort(); }

  NumBlocks_ = OrigMatrix->NumMyBlockRows();

  //extract serial dense matrices from vbr matrix
  VbrBlocks_.resize(NumBlocks_);
  VbrBlockCnt_.resize(NumBlocks_);
  VbrBlockDim_.resize(NumBlocks_);
  VbrBlockIndices_.resize(NumBlocks_);
  for( int i = 0; i < NumBlocks_; ++i )
  {
    OrigMatrix->ExtractMyBlockRowView( i, VbrBlockDim_[i], VbrBlockCnt_[i], VbrBlockIndices_[i], VbrBlocks_[i] );
  }

  SVDs_.resize(NumBlocks_);
  Inverses_.resize(NumBlocks_);
  for( int i = 0; i < NumBlocks_; ++i )
  {
    if( VbrBlockDim_[i] > 1 )
    {
      SVDs_[i] = new Epetra_SerialDenseSVD();
      SVDs_[i]->SetMatrix( *(VbrBlocks_[i][ VbrBlockCnt_[i]-1 ]) );
      SVDs_[i]->Factor();
      SVDs_[i]->Invert( rthresh_, athresh_ );
      Inverses_[i] = SVDs_[i]->InvertedMatrix();
    }
    else
    {
      SVDs_[i] = 0;
      double inv = 1. / (*(VbrBlocks_[i][ VbrBlockCnt_[i]-1 ]))(0,0);
      Inverses_[i] = new Epetra_SerialDenseMatrix( Copy, &inv, 1, 1, 1 );
    }
  }

  if( verbose_ > 2 )
  {
    cout << "SVDs and Inverses!\n";
    for( int i = 0; i < NumBlocks_; ++i )
    {
      cout << "Block: " << i << " Size: " << VbrBlockDim_[i] << endl;
      if( SVDs_[i] ) SVDs_[i]->Print(cout);
      cout << *(Inverses_[i]) << endl;
    }
  }

  Epetra_MultiVector * RHS = orig.GetRHS();
  double * A;
  int LDA;
  RHS->ExtractView( &A, &LDA );
  double * currLoc = A;
  RHSBlocks_.resize(NumBlocks_);
  for( int i = 0; i < NumBlocks_; ++i )
  {
    RHSBlocks_[i] = new Epetra_SerialDenseVector( View, currLoc, VbrBlockDim_[i] );
    currLoc += VbrBlockDim_[i];
  }

  newObj_ = &orig;

  return *newObj_;
}

bool
LinearProblem_BlockJacobi::
fwd()
{
  if( verbose_ > 2 )
  {
    cout << "-------------------\n";
    cout << "BlockJacobi\n";
    cout << "-------------------\n";
  }
  
  double MinSV =  1e15;
  double MaxSV =  0.0;

  multiset<double> SVs;

  for( int i = 0; i < NumBlocks_; ++i )
  {
    if( VbrBlockDim_[i] > 1 )
    {
      SVDs_[i]->Factor();
      if( SVDs_[i]->S()[0] > MaxSV ) MaxSV = SVDs_[i]->S()[0];
      if( SVDs_[i]->S()[VbrBlockDim_[i]-1] < MinSV ) MinSV = SVDs_[i]->S()[VbrBlockDim_[i]-1];
      for( int j = 0; j < VbrBlockDim_[i]; ++j ) SVs.insert( SVDs_[i]->S()[j] );
    }
    else
    {
      SVs.insert(1.0);
      MaxSV = max( MaxSV, 1.0 );
    }
  }

  multiset<double>::iterator iterSI = SVs.begin();
  multiset<double>::iterator endSI = SVs.end();
  int i = 0;
  if( verbose_ > 2 )
  {
    cout << endl;
    cout << "Singular Values\n";
    for( ; iterSI != endSI; ++iterSI, i++ ) cout << i << "\t" << *iterSI << endl;
    cout << endl;
  }

  Epetra_VbrMatrix * OrigMatrix = dynamic_cast<Epetra_VbrMatrix*>( origObj_->GetMatrix() );

  double abs_thresh = athresh_;
  double rel_thresh = rthresh_;
  if( thresholding_ == 1 )
  {
    abs_thresh = MaxSV * rel_thresh;
    rel_thresh = 0.0;
  }
      
  for( int i = 0; i < NumBlocks_; ++i )
  {
    if( VbrBlockDim_[i] > 1 )
      SVDs_[i]->Invert( rel_thresh, abs_thresh );
    else
      (*Inverses_[i])(0,0) = 1./(*(VbrBlocks_[i][ VbrBlockCnt_[i]-1 ]))(0,0);
  }

  for( int i = 0; i < NumBlocks_; ++i )
  {
    for( int j = 0; j < (VbrBlockCnt_[i]-1); ++j )
    {
      Epetra_SerialDenseMatrix tempMat( *(VbrBlocks_[i][j]) );
      VbrBlocks_[i][j]->Multiply( false, false, 1.0, *(Inverses_[i]), tempMat, 0.0 );
    }

    Epetra_SerialDenseMatrix tempMat2( *(RHSBlocks_[i]) );
    RHSBlocks_[i]->Multiply( false, false, 1.0, *(Inverses_[i]), tempMat2, 0.0 );

    if( verbose_ > 2 )
    {
      cout << "DiagBlock: " << i << endl;
      cout << *(VbrBlocks_[i][VbrBlockCnt_[i]-1]);
      cout << "RHSBlock: " << i << endl;
      cout << *(RHSBlocks_[i]);
    }
  }

  if( verbose_ > 2 )
  {
    cout << "Block Jacobi'd Matrix!\n";
    if( removeDiag_ ) cout << *NewMatrix_ << endl;
    else              cout << *(dynamic_cast<Epetra_VbrMatrix*>(origObj_->GetMatrix())) << endl;
    cout << "Block Jacobi'd RHS!\n";
    cout << *(origObj_->GetRHS());
    cout << endl;
  }

  if( verbose_ > 0 )
  {
    cout << "Min Singular Value: " << MinSV << endl;
    cout << "Max Singular Value: " << MaxSV << endl;
    cout << "--------------------\n";
  }

  return true;
}

bool
LinearProblem_BlockJacobi::
rvs()
{
  return true;
}

} //namespace EpetraExt
