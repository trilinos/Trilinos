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

#include <EpetraExt_Scale_LinearProblem.h>

#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

namespace EpetraExt {

LinearProblem_Scale::
~LinearProblem_Scale()
{
  int lsize = (int) lScaleVecs_.size();
  for( int i = 0; i < lsize; ++i )
    delete lScaleVecs_[i];
  int rsize = (int) rScaleVecs_.size();
  for( int i = 0; i < rsize; ++i )
    delete rScaleVecs_[i];
}

bool
LinearProblem_Scale::
fwd()
{
  Epetra_CrsMatrix & Matrix = *(dynamic_cast<Epetra_CrsMatrix*>(origObj_->GetMatrix()));

  const Epetra_BlockMap & RHSMap = origObj_->GetRHS()->Map();
  const Epetra_BlockMap & LHSMap = origObj_->GetLHS()->Map();

  if( iters_ > 0 )
  {
    if( lScale_ != None && !lScaleVecs_.size() )
    {
      lScaleVecs_.resize(iters_);
      for( int i = 0; i < iters_; ++i )
        lScaleVecs_[i] = new Epetra_Vector( RHSMap );
    }
    if( rScale_ != None && !rScaleVecs_.size() )
    {
      rScaleVecs_.resize(iters_);
      for( int i = 0; i < iters_; ++i )
        rScaleVecs_[i] = new Epetra_Vector( LHSMap );
    }

    for( int i = 0; i < iters_; ++i )
    {
      if( lScale_ != None )
      {
        switch( lScale_ )
        {
          case Max: Matrix.InvRowMaxs( *(lScaleVecs_[i]) );
                    break;
          case Sum: Matrix.InvRowSums( *(lScaleVecs_[i]) );
                    break;
          case Diag: Matrix.ExtractDiagonalCopy( *(lScaleVecs_[i]) );
		     lScaleVecs_[i]->Reciprocal( *(lScaleVecs_[i]) );
		     break;
          default:  break;
        };
        if( expFac_ != 1.0 )
        {
          int numVals = RHSMap.NumMyPoints();
          for( int j = 0; j < numVals; ++j ) (*(lScaleVecs_[i]))[j] = pow( (*(lScaleVecs_[i]))[j], expFac_ );
        }
        newObj_->LeftScale( *lScaleVecs_[i] );
      }
      if( rScale_ != None )
      {
        switch( rScale_ )
        {
          case Max: Matrix.InvColMaxs( *(rScaleVecs_[i]) );
                    break;
          case Sum: Matrix.InvColSums( *(rScaleVecs_[i]) );
                    break;
          case Diag: Matrix.ExtractDiagonalCopy( *(rScaleVecs_[i]) );
		     rScaleVecs_[i]->Reciprocal( *(rScaleVecs_[i]) );
		     break;
          default:  break;
        };
        if( expFac_ != 1.0 )
        {
          int numVals = LHSMap.NumMyPoints();
          for( int j = 0; j < numVals; ++j ) (*(rScaleVecs_[i]))[j] = pow( (*(rScaleVecs_[i]))[j], expFac_ );
        }
        newObj_->RightScale( *rScaleVecs_[i] );
      }
    }
  }

  scaled_ = true;

  return true;
}

bool
LinearProblem_Scale::
rvs()
{
  if( !scaled_ ) std::cout << "EpetraExt::LinearProblem_Scale::rvs() : Problem Not Scaled!\n";

  if( iters_ > 0 )
  {
    for( int i = 0; i < iters_; ++i )
    {
      int loc = iters_-i-1;
      if( rScale_ != None )
      {
        rScaleVecs_[loc]->Reciprocal( *(rScaleVecs_[loc]) );
        newObj_->RightScale( *(rScaleVecs_[loc]) );
      }
      if( lScale_ != None )
      {
        lScaleVecs_[loc]->Reciprocal( *(lScaleVecs_[loc]) );
        newObj_->LeftScale( *(lScaleVecs_[loc]) );
      }
    }
  }
  return true;
}

} //namespace EpetraExt

