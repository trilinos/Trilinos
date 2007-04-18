// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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


/*! \file BelosEpetraAdapter.cpp
    \brief Implementation of the interfaces between Belos virtual classes and Epetra concrete classes.
*/

#include "BelosEpetraAdapter.hpp"

using namespace Belos;
  
  //-------------------------------------------------------------
  
  //////////////////////////////////////////////////////////////////////
  // Construction/Destruction
  //////////////////////////////////////////////////////////////////////
  
  
EpetraMultiVec::EpetraMultiVec(const Epetra_BlockMap& Map, double * array, 
			       const int numvecs, const int stride)
  : Epetra_MultiVector(Copy, Map, array, stride, numvecs) 
{
}


EpetraMultiVec::EpetraMultiVec(const Epetra_BlockMap& Map, const int numvecs, bool zeroOut)
  : Epetra_MultiVector(Map, numvecs, zeroOut) 
{
}


EpetraMultiVec::EpetraMultiVec(Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, 				
			       const std::vector<int>& index )
  : Epetra_MultiVector(CV, P_vec, &(const_cast<std::vector<int> &>(index))[0], index.size())
{
}


EpetraMultiVec::EpetraMultiVec(const Epetra_MultiVector& P_vec)
  : Epetra_MultiVector(P_vec) 
{
}


EpetraMultiVec::~EpetraMultiVec() 
{
}
//
//  member functions inherited from Belos::MultiVec
//
//
//  Simulating a virtual copy constructor. If we could rely on the co-variance
//  of virtual functions, we could return a pointer to EpetraMultiVec
//  (the derived type) instead of a pointer to the pure virtual base class.
//

MultiVec<double>* EpetraMultiVec::Clone ( const int numvecs ) const
{
  EpetraMultiVec * ptr_apv = new EpetraMultiVec(Map(), numvecs, false);
  return ptr_apv; // safe upcast.
}
//
//  the following is a virtual copy constructor returning
//  a pointer to the pure virtual class. vector values are
//  copied.
//

MultiVec<double>* EpetraMultiVec::CloneCopy() const
{
  EpetraMultiVec *ptr_apv = new EpetraMultiVec(*this);
  return ptr_apv; // safe upcast
}


MultiVec<double>* EpetraMultiVec::CloneCopy ( const std::vector<int>& index ) const
{
  EpetraMultiVec * ptr_apv = new EpetraMultiVec(Copy, *this, index);
  return ptr_apv; // safe upcast.
}


MultiVec<double>* EpetraMultiVec::CloneView ( const std::vector<int>& index ) 
{
  EpetraMultiVec * ptr_apv = new EpetraMultiVec(View, *this, index);
  return ptr_apv; // safe upcast.
}
  

void EpetraMultiVec::SetBlock( const MultiVec<double>& A, const std::vector<int>& index ) 
{	
  EpetraMultiVec temp_vec(View, *this, index);
  
  int numvecs = index.size();
  if ( A.GetNumberVecs() != numvecs ) {
    std::vector<int> index2( numvecs );
    for(int i=0; i<numvecs; i++)
      index2[i] = i;
    EpetraMultiVec *tmp_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); assert(tmp_vec!=NULL);
    EpetraMultiVec A_vec(View, *tmp_vec, index2);
    temp_vec.MvAddMv( 1.0, A_vec, 0.0, A_vec );
  }
  else {
    temp_vec.MvAddMv( 1.0, A, 0.0, A );
  }
}								

//-------------------------------------------------------------
//
// *this <- alpha * A * B + beta * (*this)
//
//-------------------------------------------------------------

void EpetraMultiVec::MvTimesMatAddMv ( const double alpha, const MultiVec<double>& A, 
				       const Teuchos::SerialDenseMatrix<int,double>& B, const double beta ) 
{
  Epetra_LocalMap LocalMap(B.numRows(), 0, Map().Comm());
  Epetra_MultiVector B_Pvec(Copy, LocalMap, B.values(), B.stride(), B.numCols());
  
  EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); assert(A_vec!=NULL);
  
  assert( Multiply( 'N', 'N', alpha, *A_vec, B_Pvec, beta ) == 0 );
}

//-------------------------------------------------------------
//
// *this <- alpha * A + beta * B
//
//-------------------------------------------------------------
  
void EpetraMultiVec::MvAddMv ( const double alpha , const MultiVec<double>& A, 
			       const double beta, const MultiVec<double>& B) 
{
  EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); assert(A_vec!=NULL);
  EpetraMultiVec *B_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(B)); assert(B_vec!=NULL);
  
  assert ( Update( alpha, *A_vec, beta, *B_vec, 0.0 ) == 0 ); 
}

//-------------------------------------------------------------
//
// dense B <- alpha * A^T * (*this)
//
//-------------------------------------------------------------

void EpetraMultiVec::MvTransMv ( const double alpha, const MultiVec<double>& A,
				 Teuchos::SerialDenseMatrix<int,double>& B) const
{    
  EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A));
  
  if (A_vec) {
    Epetra_LocalMap LocalMap(B.numRows(), 0, Map().Comm());
    Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());
    
    assert ( B_Pvec.Multiply( 'T', 'N', alpha, *A_vec, *this, 0.0 ) == 0 ); 
  }
}

//-------------------------------------------------------------
//
// b[i] = A[i]^T * this[i]
// 
//-------------------------------------------------------------

void EpetraMultiVec::MvDot ( const MultiVec<double>& A, std::vector<double>* b ) const
{
  EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); assert(A_vec!=NULL);
  if (A_vec && b && ( (int)b->size() >= A_vec->NumVectors() ) ) {
    assert( this->Dot( *A_vec, &(*b)[0] ) == 0 );
  }
}

//-------------------------------------------------------------
//
// alpha[i] = norm of i-th column of (*this)
//
//-------------------------------------------------------------

void EpetraMultiVec::MvNorm ( std::vector<double>* normvec, NormType norm_type ) const {
  if (normvec && ((int)normvec->size() >= GetNumberVecs())) {
    switch( norm_type ) {
    case ( OneNorm ) :
      assert( Norm1(&(*normvec)[0]) == 0 );
      break;
    case ( TwoNorm ) :
      assert( Norm2(&(*normvec)[0]) == 0 );
      break;
    case ( InfNorm ) :	
      assert( NormInf(&(*normvec)[0]) == 0 );
      break;
    default:
      break;
    }
  }
}

///////////////////////////////////////////////////////////////
//
// implementation of the BelosEpetraOp class.
//
////////////////////////////////////////////////////////////////////
//
// BelosOperator constructors
//

EpetraOp::EpetraOp( const Teuchos::RefCountPtr<Epetra_Operator> &Op ) 
  : Epetra_Op(Op)
{
}

//
// BelosOperator applications
//
ReturnType EpetraOp::Apply ( const MultiVec<double>& x, 
			     MultiVec<double>& y, ETrans trans ) const {
  int info=0;
  MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
  Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
  //
  assert( vec_x!=NULL && vec_y!=NULL );
  //
  // Set the operator to apply the transpose if that is requested.
  // (TO DO:  insert a throw if the application returns a nonzero )
  //
  if ( trans ) { 
    info = Epetra_Op->SetUseTranspose( true );
    if (info != 0) { return Undef; }
  }
  info = Epetra_Op->Apply( *vec_x, *vec_y );
  
  if ( trans ) { 
    info = Epetra_Op->SetUseTranspose( false );
    if (info != 0) { return Undef; }
  }
  
  if (info != 0) { return Error; }
  return Ok;	
}

///////////////////////////////////////////////////////////////
//
// implementation of the BelosEpetraPrecOp class.
//
////////////////////////////////////////////////////////////////////
//
// BelosOperator constructors
//

EpetraPrecOp::EpetraPrecOp( const Teuchos::RefCountPtr<Epetra_Operator> &Op ) 
  : Epetra_Op(Op)
{
}

//
// BelosOperator applications
//
ReturnType EpetraPrecOp::Apply ( const MultiVec<double>& x, 
				 MultiVec<double>& y, ETrans trans ) const {
  int info=0;
  MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
  Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
  //
  assert( vec_x!=NULL && vec_y!=NULL );
  //
  // Set the operator to apply the transpose if that is requested.
  // (TO DO:  insert a throw if the application returns a nonzero )
  //
  if ( trans ) { 
    info = Epetra_Op->SetUseTranspose( true );
    if (info != 0) { return Undef; }
  }
  info = Epetra_Op->ApplyInverse( *vec_x, *vec_y );
  
  if ( trans ) { 
    info = Epetra_Op->SetUseTranspose( false );
    if (info != 0) { return Undef; }
  }
  
  if (info != 0) { return Error; }
  return Ok;	
}

int EpetraPrecOp::Apply ( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const 
{
  //
  // This operator applies Y = A^{-1}*X
  //
  int info = 0;
  info = Epetra_Op->ApplyInverse( X, Y );
  return info;
}

int EpetraPrecOp::ApplyInverse( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const
{
  //
  // This operator applies Y = A*X
  //
  int info = 0;
  info = Epetra_Op->Apply( X, Y );
  return info;
}

