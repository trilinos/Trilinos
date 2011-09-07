//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
// ************************************************************************
//@HEADER


/*! \file BelosEpetraAdapter.cpp
    \brief Implementation of the interfaces between Belos virtual classes and Epetra concrete classes.
*/

#include "BelosEpetraAdapter.hpp"

namespace Belos {

  //-------------------------------------------------------------
  
  //////////////////////////////////////////////////////////////////////
  // Construction/Destruction
  //////////////////////////////////////////////////////////////////////
  
  
EpetraMultiVec::EpetraMultiVec(const Epetra_BlockMap& Map_in, double * array, 
			       const int numvecs, const int stride)
  : Epetra_MultiVector(Copy, Map_in, array, stride, numvecs) 
{
}


EpetraMultiVec::EpetraMultiVec(const Epetra_BlockMap& Map_in, const int numvecs, bool zeroOut)
  : Epetra_MultiVector(Map_in, numvecs, zeroOut) 
{
}


EpetraMultiVec::EpetraMultiVec(Epetra_DataAccess CV_in, const Epetra_MultiVector& P_vec, 				
			       const std::vector<int>& index )
  : Epetra_MultiVector(CV_in, P_vec, &(const_cast<std::vector<int> &>(index))[0], index.size())
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
//  a pointer to the pure virtual class. std::vector values are
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


MultiVec<double>* EpetraMultiVec::CloneViewNonConst ( const std::vector<int>& index ) 
{
  EpetraMultiVec * ptr_apv = new EpetraMultiVec(View, *this, index);
  return ptr_apv; // safe upcast.
}
  

const MultiVec<double>* EpetraMultiVec::CloneView ( const std::vector<int>& index ) const
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
    EpetraMultiVec *tmp_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
    TEST_FOR_EXCEPTION(tmp_vec==NULL, EpetraMultiVecFailure,
                       "Belos::EpetraMultiVec::SetBlock cast from Belos::MultiVec<> to Belos::EpetraMultiVec failed.");
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
  Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());
  
  EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
  TEST_FOR_EXCEPTION(A_vec==NULL, EpetraMultiVecFailure,
                     "Belos::EpetraMultiVec::MvTimesMatAddMv cast from Belos::MultiVec<> to Belos::EpetraMultiVec failed.");
  
  int info = Multiply( 'N', 'N', alpha, *A_vec, B_Pvec, beta );
  TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
		     "Belos::EpetraMultiVec::MvTimesMatAddMv call to Multiply() returned a nonzero value.");

}

//-------------------------------------------------------------
//
// *this <- alpha * A + beta * B
//
//-------------------------------------------------------------
  
void EpetraMultiVec::MvAddMv ( const double alpha , const MultiVec<double>& A, 
			       const double beta, const MultiVec<double>& B) 
{
  EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A));
  TEST_FOR_EXCEPTION( A_vec==NULL, EpetraMultiVecFailure,
                     "Belos::EpetraMultiVec::MvAddMv cast from Belos::MultiVec<> to Belos::EpetraMultiVec failed.");
  EpetraMultiVec *B_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(B));
  TEST_FOR_EXCEPTION( B_vec==NULL, EpetraMultiVecFailure,
                     "Belos::EpetraMultiVec::MvAddMv cast from Belos::MultiVec<> to Belos::EpetraMultiVec failed.");
  
  int info = Update( alpha, *A_vec, beta, *B_vec, 0.0 );
  TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
		     "Belos::EpetraMultiVec::MvAddMv call to Update() returned a nonzero value.");
}

//-------------------------------------------------------------
//
// this[i] = alpha[i] * this[i]
//
//-------------------------------------------------------------
void EpetraMultiVec::MvScale ( const std::vector<double>& alpha )
{
  // Check to make sure the vector is as long as the multivector has columns.
  int numvecs = this->NumVectors();
  TEST_FOR_EXCEPTION((int)alpha.size() != numvecs, EpetraMultiVecFailure, 
			 "Belos::MultiVecTraits<double,Epetra_MultiVec>::MvScale scaling vector (alpha) not same size as number of input vectors (mv).");
  int ret = 0;
  std::vector<int> tmp_index( 1, 0 );
  for (int i=0; i<numvecs; i++) {
    Epetra_MultiVector temp_vec(View, *this, &tmp_index[0], 1);
    ret = temp_vec.Scale( alpha[i] );
    TEST_FOR_EXCEPTION(ret!=0, EpetraMultiVecFailure, 
                      "Belos::MultiVecTraits<double,Epetra_MultiVec>::MvScale call to Scale() returned a nonzero value.");
    tmp_index[0]++;
  }
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
    
    int info = B_Pvec.Multiply( 'T', 'N', alpha, *A_vec, *this, 0.0 );
    TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
		       "Belos::EpetraMultiVec::MvTransMv call to Multiply() returned a nonzero value.");
  }
}

//-------------------------------------------------------------
//
// b[i] = A[i]^T * this[i]
// 
//-------------------------------------------------------------

void EpetraMultiVec::MvDot ( const MultiVec<double>& A, std::vector<double>& b ) const
{
  EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
  TEST_FOR_EXCEPTION(A_vec==NULL, EpetraMultiVecFailure,
                     "Belos::EpetraMultiVec::MvDot cast from Belos::MultiVec<> to Belos::EpetraMultiVec failed.");
  if (A_vec && ( (int)b.size() >= A_vec->NumVectors() ) ) {
     int info = this->Dot( *A_vec, &b[0] );
     TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
			"Belos::EpetraMultiVec::MvDot call to Dot() returned a nonzero value.");   
  }
}

//-------------------------------------------------------------
//
// alpha[i] = norm of i-th column of (*this)
//
//-------------------------------------------------------------

void EpetraMultiVec::MvNorm ( std::vector<double>& normvec, NormType norm_type ) const {
  if ((int)normvec.size() >= GetNumberVecs()) {
    int info = 0;
    switch( norm_type ) {
    case ( OneNorm ) :
      info = Norm1(&normvec[0]);
      break;
    case ( TwoNorm ) :
      info = Norm2(&normvec[0]);
      break;
    case ( InfNorm ) :	
      info = NormInf(&normvec[0]);
      break;
    default:
      break;
    }
    TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
		       "Belos::EpetraMultiVec::MvNorm call to Norm() returned a nonzero value.");
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

EpetraOp::EpetraOp( const Teuchos::RCP<Epetra_Operator> &Op ) 
  : Epetra_Op(Op)
{
}

//
// BelosOperator applications
//
void EpetraOp::Apply ( const MultiVec<double>& x, 
                       MultiVec<double>& y, ETrans trans ) const {
  int info=0;
  MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
  Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
  //
  TEST_FOR_EXCEPTION( vec_x==NULL || vec_y==NULL, EpetraOpFailure, 
		      "Belos::EpetraOp::Apply, x and/or y cannot be dynamic cast to an Epetra_MultiVector.");
  //
  // Set the operator to apply the transpose if that is requested.
  //
  if ( trans ) { 
    info = Epetra_Op->SetUseTranspose( true );
  }
  info = Epetra_Op->Apply( *vec_x, *vec_y );
  
  if ( trans ) { 
    info = Epetra_Op->SetUseTranspose( false );
  }
  
  TEST_FOR_EXCEPTION(info!=0, EpetraOpFailure, 
		     "Belos::EpetraOp::Apply call to Apply() returned a nonzero value."); 

}

// ///////////////////////////////////////////////////////////////////
//
// Implementation of the Belos::EpetraPrecOp class.
//
// ///////////////////////////////////////////////////////////////////

EpetraPrecOp::EpetraPrecOp (const Teuchos::RCP<Epetra_Operator> &Op) 
  : Epetra_Op(Op)
{}

// The version of Apply() that takes an optional 'trans' argument and
// returns void implements the Belos::Operator interface.
void 
EpetraPrecOp::Apply (const MultiVec<double>& x, 
		     MultiVec<double>& y, 
		     ETrans trans) const 
{
  int info = 0;
  MultiVec<double>&  temp_x = const_cast<MultiVec<double> &>(x);
  Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  TEST_FOR_EXCEPTION(vec_x == NULL, EpetraOpFailure, 
		     "Belos::EpetraPrecOp::Apply: The MultiVec<double> input x "
		     "cannot be dynamic cast to an Epetra_MultiVector.");
  Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
  TEST_FOR_EXCEPTION(vec_x == NULL, EpetraOpFailure, 
		     "Belos::EpetraPrecOp::Apply: The MultiVec<double> input y "
		     "cannot be dynamic cast to an Epetra_MultiVector.");

  // Remember the original value of the transpose flag, so that we can
  // restore it on exit.
  //
  // Do NOT change this behavior just because you don't want an
  // exception to be thrown.  The semantics of Belos::Operator demand
  // that this version of the Apply() method be STATELESS, and that
  // the 'trans' argument have its literal meaning.  If the underlying
  // operator doesn't implement applying the transpose, then it is
  // RIGHT for this method to throw an exception immediately, rather
  // than to silently compute the wrong thing.  However, if you don't
  // try to apply the transpose, then no exception will be thrown.
  const bool originalTransposeFlag = Epetra_Op->UseTranspose ();

  // Whether we want to apply the transpose of the operator.
  // Recall that Epetra operators are always real-valued, never
  // complex-valued, so the conjugate transpose means the same
  // thing as the transpose.
  const bool newTransposeFlag = (trans != NOTRANS);

  // If necessary, set (or unset) the transpose flag to the value
  // corresponding to 'trans'.
  if (newTransposeFlag != originalTransposeFlag) {
    // Toggle the transpose flag temporarily.  We will restore the
    // original flag's value before this method exits.
    info = Epetra_Op->SetUseTranspose (newTransposeFlag);
    TEST_FOR_EXCEPTION(info != 0, EpetraOpFailure,
		       "Belos::EpetraPrecOp::Apply: Toggling the transpose flag "
		       "of the operator failed, returning a nonzero error code " 
		       "of " << info << ".  That probably means the underlying "
		       "Epetra_Operator object doesn't know how to apply its "
		       "transpose.  Are you perhaps using a preconditioner for "
		       "which applying the transpose is not implemented?");
  }

  // EpetraPrecOp's Apply() methods apply the inverse of the
  // underlying operator.  This may not succeed for all
  // implementations of Epetra_Operator, so we have to check.
  info = Epetra_Op->ApplyInverse (*vec_x, *vec_y);
  TEST_FOR_EXCEPTION(info != 0, EpetraOpFailure, 
		     "Belos::EpetraPrecOp::Apply: Calling ApplyInverse() on the "
		     "underlying Epetra_Operator object failed, returning a "
		     "nonzero error code of " << info << ".  This probably means"
		     " that the underlying Epetra_Operator object doesn't know "
		     "how to apply its inverse.");

  // SetUseTranspose() changes the state of the operator, so if
  // applicable, we have to change the state back.
  if (newTransposeFlag != originalTransposeFlag) {
    info = Epetra_Op->SetUseTranspose (originalTransposeFlag);
    TEST_FOR_EXCEPTION(info != 0, EpetraOpFailure,
		       "Belos::EpetraPrecOp::Apply: Resetting the original "
		       "transpose flag value of the underlying Epetra_Operator "
		       "failed, returning a nonzero error code of " << info
		       << ".  This is an unusual error, since we were able to "
		       "call its SetUseTranspose method successfully before.  "
		       "This suggests a bug in the subclass of Epetra_Operator "
		       "which you are currently using.  This is probably not a "
		       "Belos bug.");
  }

  // Make sure that the transpose flag has its original value.  If
  // not, we throw std::logic_error instead of EpetraOpFailure, since
  // that's definitely a bug.  It's safe to do this check whether or
  // not we actually had to toggle the transpose flag.  Any reasonable
  // implementation of Epetra_Operator should make calling
  // UseTranspose() cheap.
  //
  // Note to code maintainers: The reason we capture the value of
  // UseTranspose() instead of just calling it twice, is that if the
  // exception test really does trigger, then there is something
  // fundamentally wrong.  If something is that wrong, then we want to
  // avoid further calls to the operator's methods.  For example, it
  // could be that the UseTranspose() method is erroneously changing
  // the value of the flag, so if we call that method twice, it might
  // have the right value on the second call.  This would make the
  // resulting exception message confusing.
  const bool finalTransposeFlag = Epetra_Op->UseTranspose ();
  TEST_FOR_EXCEPTION(originalTransposeFlag != finalTransposeFlag,
		     std::logic_error,
		     "Belos::EpetraPrecOp::Apply: The value of the underlying "
		     "Epetra_Operator's transpose flag changed unexpectedly! "
		     "The original value at the top of this method was " 
		     << originalTransposeFlag << ", and its new value is " 
		     << finalTransposeFlag << ".  This suggests either a bug in "
		     "Belos (in the implementation of this routine), or a bug in"
		     " the underlying operator.");
}

// The version of Apply() that takes two arguments and returns int
// implements the Epetra_Operator interface.
int 
EpetraPrecOp::Apply (const Epetra_MultiVector &X, 
		     Epetra_MultiVector &Y) const 
{
  //
  // This operator applies Y = A^{-1}*X
  //
  int info = Epetra_Op->ApplyInverse( X, Y );
  return info;
}

// This implements Epetra_Operator::ApplyInverse().
int 
EpetraPrecOp::ApplyInverse (const Epetra_MultiVector &X, 
			    Epetra_MultiVector &Y) const
{
  //
  // This operator applies Y = A*X
  //
  int info = Epetra_Op->Apply( X, Y );
  return info;
}

}  // end namespace Belos
