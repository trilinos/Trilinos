/*
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
*/

#include "GenSQP_YUEpetraVector.hpp"

namespace GenSQP {

YUEpetraVector::YUEpetraVector( const Teuchos::RefCountPtr<Epetra_MultiVector> &y_epetra_vec,
                                const Teuchos::RefCountPtr<Epetra_MultiVector> &u_epetra_vec )
  :y_epetra_vec_(y_epetra_vec), u_epetra_vec_(u_epetra_vec)
    {}

// Overridden from Vector

double YUEpetraVector::innerProd( const Vector &x ) const
{
  double ydot[1];
  double udot[1];
  YUEpetraVector &ex = Teuchos::dyn_cast<YUEpetraVector>(const_cast <Vector&>(x));
  y_epetra_vec_->Dot( *ex.y_epetra_vec_, ydot );
  if (u_epetra_vec_.get() == 0)
    udot[0] = 0.0;
  else
    u_epetra_vec_->Dot( *ex.u_epetra_vec_, udot );
  return (ydot[0] + udot[0]);
}

void YUEpetraVector::linComb( const double &alpha, const Vector &x, const double &beta )
{
  YUEpetraVector &ex = Teuchos::dyn_cast<YUEpetraVector>(const_cast <Vector&>(x));
  y_epetra_vec_->Update( alpha, *ex.y_epetra_vec_, beta );
  if (u_epetra_vec_.get() != 0)
    u_epetra_vec_->Update( alpha, *ex.u_epetra_vec_, beta );
}

void YUEpetraVector::Scale( const double &alpha )
{
  y_epetra_vec_->Scale( alpha );
  if (u_epetra_vec_.get() != 0)
    u_epetra_vec_->Scale( alpha );
}

void YUEpetraVector::Set( const double &alpha )
{
  y_epetra_vec_->PutScalar( alpha );
  if (u_epetra_vec_.get() != 0)
    u_epetra_vec_->PutScalar( alpha );
}

void YUEpetraVector::Set( const double &alpha, const Vector &x )
{
  YUEpetraVector &ex = Teuchos::dyn_cast<YUEpetraVector>(const_cast <Vector&>(x));
  y_epetra_vec_->Scale( alpha, *ex.y_epetra_vec_ );
  if (u_epetra_vec_.get() != 0)
    u_epetra_vec_->Scale( alpha, *ex.u_epetra_vec_ );
}

Teuchos::RefCountPtr<Vector> YUEpetraVector::createVector() const
{
  Teuchos::RefCountPtr<Epetra_MultiVector> yptr =
      Teuchos::rcp(new Epetra_MultiVector(y_epetra_vec_->Map(),1,false));
  Teuchos::RefCountPtr<Epetra_MultiVector> uptr = Teuchos::null;
  if (u_epetra_vec_.get() != 0)
    uptr = Teuchos::rcp(new Epetra_MultiVector(u_epetra_vec_->Map(),1,false));
  
  return Teuchos::rcp( new YUEpetraVector( yptr, uptr ));
}

Teuchos::RefCountPtr<const Epetra_MultiVector> YUEpetraVector::getYVector() const
{
  return y_epetra_vec_;
}

Teuchos::RefCountPtr<const Epetra_MultiVector> YUEpetraVector::getUVector() const
{
  return u_epetra_vec_;
}

} // namespace GenSQP
