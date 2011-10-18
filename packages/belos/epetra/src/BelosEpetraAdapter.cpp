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

  // An anonymous namespace restricts the scope of its definitions to
  // this file.
  namespace {

    /// \fn implementsApplyTranspose
    /// \brief Whether Op implements applying the transpose.
    ///
    /// Epetra_Operator instances are not required to implement
    /// applying the transpose (or conjugate transpose, if
    /// applicable).  This function lets you ask.  
    ///
    /// We need this function because Epetra_Operator doesn't have a
    /// direct query method.  We query by asking it to set using the
    /// transpose, and examine the returned error code.  Epetra
    /// operators have a persistent "use the transpose" state that you
    /// can set or unset using Epetra_Operator::SetUseTranspose().
    ///
    /// \note This function is not thread-safe, because it makes no
    ///   effort to protect Op against simultaneous calls to its
    ///   Apply() or SetUseTranspose() methods.
    ///
    /// We use this function multiple times in this file, so it's
    /// worthwhile breaking it out into a separate function, rather
    /// than copying and pasting.
    bool 
    implementsApplyTranspose (const Epetra_Operator& Op)
    {
      // We must first check Op's current transpose state, so that we
      // leave Op in the same state as it was given to us.
      const bool transposed = Op.UseTranspose ();

      // If Op is already transposed, Apply() doesn't need to do anything
      // in order to apply the transpose.  (Apply() will just apply the
      // already transposed operator in that case.)
      if (transposed)
	return true;

      // SetUseTranspose() is a nonconst method, so we have to cast away
      // const-ness of Op in order to call it.  Epetra_Operator gives us
      // no other way to query whether it knows how to apply the
      // transpose.
      const int errcode = const_cast<Epetra_Operator&>(Op).SetUseTranspose (true);

      // SetUseTranspose() follows the POSIX convention that returning
      // zero means success.  If SetUseTranspose() failed, the best we can
      // do is assume that the state of Op was not changed.  Any
      // reasonable implementation of SetUseTranspose() should not change
      // the state of the operator if setting the transpose failed.
      if (errcode != 0)
	return false;

      // Now we need to restore the original transpose state.
      if (! transposed)
	// We assume that this call succeeds.
	(void) const_cast<Epetra_Operator&>(Op).SetUseTranspose (false);

      return errcode == 0;
    }

    /// \class EpetraOperatorTransposeScopeGuard
    /// \brief Safely sets and unsets an Epetra_Operator's transpose state.
    ///
    /// This class uses RAII (Resource Acquisition Is Instantiation)
    /// techniques to set the transpose state of an Epetra_Operator
    /// instance on construction, and restore its original transpose
    /// state on destruction.
    class EpetraOperatorTransposeScopeGuard {
    private:
      //! Reference to the Epetra_Operator instance to protect.
      const Epetra_Operator& Op_;
      
      //! The original transpose state of the operator.
      const bool originalTransposeFlag_;

      /// \brief Whether we want to apply the transpose of the operator.
      ///
      /// Recall that Epetra operators are always real-valued, never
      /// complex-valued, so the conjugate transpose means the same
      /// thing as the transpose.
      const bool newTransposeFlag_;

    public:
      EpetraOperatorTransposeScopeGuard (const Epetra_Operator& Op, 
					 const ETrans trans)
	: Op_ (Op),
	  originalTransposeFlag_ (Op.UseTranspose ()),
	  newTransposeFlag_ (trans != NOTRANS)
      {
	// If necessary, set (or unset) the transpose flag to the value
	// corresponding to 'trans'.
	if (newTransposeFlag_ != originalTransposeFlag_) {
	  // Toggle the transpose flag.  The destructor will restore
	  // its original value.
	  const int info = 
	    const_cast<Epetra_Operator &>(Op_).SetUseTranspose (newTransposeFlag_);
	  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, EpetraOpFailure,
			     "Toggling the transpose flag of the operator failed, "
			     "returning a nonzero error code of " << info 
			     << ".  That probably means that the Epetra_Operator "
			     "subclass instance doesn't know how to apply its "
			     "transpose.  Are you perhaps using a subclass of "
			     "Epetra_Operator for which applying the transpose "
			     "is not implemented?  Anyway, just to help with "
			     "debugging, you specified trans=" 
			     << (trans == NOTRANS ? "NOTRANS" : (trans == TRANS ? "TRANS" : "CONJTRANS"))
			     << ", and the original operator was set "
			     << (originalTransposeFlag_ ? "" : "NOT ")
			     << "to apply the transpose.");
	}
      }

      ~EpetraOperatorTransposeScopeGuard ()
      {
	// SetUseTranspose() changes the state of the operator, so if
	// applicable, we have to change the state back.
	if (newTransposeFlag_ != originalTransposeFlag_) {
	  const int info = 
	    const_cast<Epetra_Operator &>(Op_).SetUseTranspose (originalTransposeFlag_);
	  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, EpetraOpFailure,
			     "Resetting the original "
			     "transpose flag value of the Epetra_Operator failed, "
			     "returning a nonzero error code of " << info << ".  "
			     "This is an unusual error, since we were able to call "
			     "its SetUseTranspose() method successfully before.  "
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
	const bool finalTransposeFlag = Op_.UseTranspose ();
	TEUCHOS_TEST_FOR_EXCEPTION(originalTransposeFlag_ != finalTransposeFlag,
			   std::logic_error,
			   "Belos::OperatorTraits::Apply: The value of the "
			   "Epetra_Operator's transpose flag changed unexpectedly!"
			   "  The original value at the top of this method was "
			   << originalTransposeFlag_ << ", and its new value is "
			   << finalTransposeFlag << ".  This suggests either a "
			   "bug in Belos (in the implementation of this routine), "
			   "or a bug in the operator.");
      }
    };

  } // namespace (anonymous)

  // ///////////////////////////////////////////////////////////////////
  // Construction/Destruction
  // ///////////////////////////////////////////////////////////////////
  
  
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
    TEUCHOS_TEST_FOR_EXCEPTION(tmp_vec==NULL, EpetraMultiVecFailure,
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
  TEUCHOS_TEST_FOR_EXCEPTION(A_vec==NULL, EpetraMultiVecFailure,
                     "Belos::EpetraMultiVec::MvTimesMatAddMv cast from Belos::MultiVec<> to Belos::EpetraMultiVec failed.");
  
  int info = Multiply( 'N', 'N', alpha, *A_vec, B_Pvec, beta );
  TEUCHOS_TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
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
  TEUCHOS_TEST_FOR_EXCEPTION( A_vec==NULL, EpetraMultiVecFailure,
                     "Belos::EpetraMultiVec::MvAddMv cast from Belos::MultiVec<> to Belos::EpetraMultiVec failed.");
  EpetraMultiVec *B_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(B));
  TEUCHOS_TEST_FOR_EXCEPTION( B_vec==NULL, EpetraMultiVecFailure,
                     "Belos::EpetraMultiVec::MvAddMv cast from Belos::MultiVec<> to Belos::EpetraMultiVec failed.");
  
  int info = Update( alpha, *A_vec, beta, *B_vec, 0.0 );
  TEUCHOS_TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
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
  TEUCHOS_TEST_FOR_EXCEPTION((int)alpha.size() != numvecs, EpetraMultiVecFailure, 
			 "Belos::MultiVecTraits<double,Epetra_MultiVec>::MvScale scaling vector (alpha) not same size as number of input vectors (mv).");
  int ret = 0;
  std::vector<int> tmp_index( 1, 0 );
  for (int i=0; i<numvecs; i++) {
    Epetra_MultiVector temp_vec(View, *this, &tmp_index[0], 1);
    ret = temp_vec.Scale( alpha[i] );
    TEUCHOS_TEST_FOR_EXCEPTION(ret!=0, EpetraMultiVecFailure, 
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
    TEUCHOS_TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
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
  TEUCHOS_TEST_FOR_EXCEPTION(A_vec==NULL, EpetraMultiVecFailure,
                     "Belos::EpetraMultiVec::MvDot cast from Belos::MultiVec<> to Belos::EpetraMultiVec failed.");
  if (A_vec && ( (int)b.size() >= A_vec->NumVectors() ) ) {
     int info = this->Dot( *A_vec, &b[0] );
     TEUCHOS_TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
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
    TEUCHOS_TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
		       "Belos::EpetraMultiVec::MvNorm call to Norm() returned a nonzero value.");
  }
}

///////////////////////////////////////////////////////////////
//
// Implementation of the Belos::EpetraOp class.
//
///////////////////////////////////////////////////////////////

EpetraOp::EpetraOp( const Teuchos::RCP<Epetra_Operator> &Op ) 
  : Epetra_Op(Op)
{
}

void 
EpetraOp::Apply (const MultiVec<double>& x, 
		 MultiVec<double>& y, 
		 ETrans trans) const 
{
  MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
  Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);

  TEUCHOS_TEST_FOR_EXCEPTION( vec_x==NULL || vec_y==NULL, EpetraOpFailure, 
		      "Belos::EpetraOp::Apply: x and/or y could not be "
		      "dynamic cast to an Epetra_MultiVector.");

  // Temporarily set the transpose state of Op, if it's not the same
  // as trans, and restore it on exit of this scope.
  EpetraOperatorTransposeScopeGuard guard (*Epetra_Op, trans);

  // Apply the operator to x and put the result in y.
  const int info = Epetra_Op->Apply (*vec_x, *vec_y);
  TEUCHOS_TEST_FOR_EXCEPTION(info!=0, EpetraOpFailure, 
		     "Belos::EpetraOp::Apply: The Epetra_Operator's Apply() "
		     "method returned a nonzero value of " << info << ".");
}

bool 
EpetraOp::HasApplyTranspose() const
{
  return implementsApplyTranspose (*Epetra_Op);
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
  MultiVec<double>&  temp_x = const_cast<MultiVec<double> &>(x);
  Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  TEUCHOS_TEST_FOR_EXCEPTION(vec_x == NULL, EpetraOpFailure, 
		     "Belos::EpetraPrecOp::Apply: The MultiVec<double> input x "
		     "cannot be dynamic cast to an Epetra_MultiVector.");
  Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
  TEUCHOS_TEST_FOR_EXCEPTION(vec_x == NULL, EpetraOpFailure, 
		     "Belos::EpetraPrecOp::Apply: The MultiVec<double> input y "
		     "cannot be dynamic cast to an Epetra_MultiVector.");

  // Temporarily set the transpose state of Op, if it's not the same
  // as trans, and restore it on exit of this scope.
  EpetraOperatorTransposeScopeGuard guard (*Epetra_Op, trans);

  // EpetraPrecOp's Apply() methods apply the inverse of the
  // underlying operator.  This may not succeed for all
  // implementations of Epetra_Operator, so we have to check.
  const int info = Epetra_Op->ApplyInverse (*vec_x, *vec_y);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, EpetraOpFailure, 
		     "Belos::EpetraPrecOp::Apply: Calling ApplyInverse() on the "
		     "underlying Epetra_Operator object failed, returning a "
		     "nonzero error code of " << info << ".  This probably means"
		     " that the underlying Epetra_Operator object doesn't know "
		     "how to apply its inverse.");
}

// The version of Apply() that takes two arguments and returns int
// implements the Epetra_Operator interface.
int 
EpetraPrecOp::Apply (const Epetra_MultiVector &X, 
		     Epetra_MultiVector &Y) const 
{
  // This operation computes Y = A^{-1}*X.
  const int info = Epetra_Op->ApplyInverse( X, Y );
  return info;
}

// This implements Epetra_Operator::ApplyInverse().
int 
EpetraPrecOp::ApplyInverse (const Epetra_MultiVector &X, 
			    Epetra_MultiVector &Y) const
{
  // This operation computes Y = A*X.
  const int info = Epetra_Op->Apply( X, Y );
  return info;
}

bool 
EpetraPrecOp::HasApplyTranspose() const
{
  return implementsApplyTranspose (*Epetra_Op);
}


// ///////////////////////////////////////////////////////////////////
// 
// Specialization of Belos::OperatorTraits for Epetra_Operator.
//
// ///////////////////////////////////////////////////////////////////
  
void
OperatorTraits<double, Epetra_MultiVector, Epetra_Operator>::
Apply (const Epetra_Operator& Op, 
       const Epetra_MultiVector& x, 
       Epetra_MultiVector& y,
       ETrans trans)
{ 
  // Temporarily set the transpose state of Op, if it's not the same
  // as trans, and restore it on exit of this scope.
  EpetraOperatorTransposeScopeGuard guard (Op, trans);

  const int info = Op.Apply (x, y);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, EpetraOpFailure, 
		     "Belos::OperatorTraits::Apply (Epetra specialization): "
		     "Calling the Epetra_Operator object's Apply() method "
		     "failed, returning a nonzero error code of " << info
		     << ".");
}

bool
OperatorTraits<double, Epetra_MultiVector, Epetra_Operator>::
HasApplyTranspose (const Epetra_Operator& Op)
{
  return implementsApplyTranspose (Op);
}


}  // end namespace Belos
