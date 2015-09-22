//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

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

