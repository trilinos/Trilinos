// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file BelosEpetraAdapter.cpp
    \brief Implementation of the interfaces between Belos virtual classes and Epetra concrete classes.
*/

#include "BelosEpetraAdapter.hpp"

// mfh 19 Oct 2011: Enabling this compile-time option might help you
// debug tricky issues involving Epetra operators and preconditioners.
// For example, today we used this option to figure out that a user's
// wrapped preconditioner was not correctly initializing its
// UseTranspose() flag.
#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
#  undef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
#  include "Teuchos_oblackholestream.hpp"
// We have to include this file, because we're calling a member
// function of Epetra_Comm.
#  include "Epetra_Comm.h"
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING


namespace Belos {

  // An anonymous namespace restricts the scope of its definitions to
  // this file.
  namespace {

#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
    /// \fn getErrStream
    /// \brief Get a stream for debugging output that prints only on MPI Rank 0.
    ///
    /// Rank 0 is relative to the Epetra_Operator's MPI communicator,
    /// which is why this function wants an Epetra_Operator input.
    ///
    /// \warning This function is NOT reentrant, because it has static data.
    Teuchos::RCP<std::ostream>
    getErrStream (const Epetra_Operator& Op)
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcpFromRef;

      // The static object makes this function NOT reentrant.
      static RCP<std::ostream> errStream;

      if (errStream.is_null()) {
        if (Op.Comm().MyPID() == 0) {
          errStream = rcpFromRef (std::cerr);
        } else {
          errStream = rcp (new Teuchos::oblackholestream);
        }
      }
      return errStream;
    }
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

    //! Return a string representation of the given transpose enum value.
    std::string
    etransToString (const ETrans trans)
    {
      if (trans == NOTRANS) {
        return "NOTRANS";
      } else if (trans == TRANS) {
        return "TRANS";
      } else if (trans == CONJTRANS) {
        return "CONJTRANS";
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                         "Invalid ETrans value trans = " << trans << ".");
    }

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

      // If Op is already transposed, then obviously it implements
      // applying the transpose.
      if (transposed) {
        return true;
      }

      // SetUseTranspose() is a nonconst method, so we have to cast
      // away const-ness of Op in order to call it.  Epetra_Operator
      // gives us no other way to query whether it knows how to apply
      // the transpose.
      const int transposeFailed =
        const_cast<Epetra_Operator&>(Op).SetUseTranspose (true);

      // SetUseTranspose() follows the POSIX convention that returning
      // zero means success.  If SetUseTranspose() failed, the best we can
      // do is assume that the state of Op was not changed.  Any
      // reasonable implementation of SetUseTranspose() should not change
      // the state of the operator if setting the transpose failed.
      if (transposeFailed != 0) {
        // Make sure that the failed call didn't managed to change the
        // state of Op.  It shouldn't have, but if it did, then Op is
        // broken.
        const bool currentTransposedState = Op.UseTranspose();
        TEUCHOS_TEST_FOR_EXCEPTION(currentTransposedState != transposed, std::logic_error,
                           "implementsApplyTranspose: Op.SetUseTranspose(true) "
                           "failed (with error code " << transposeFailed << "), "
                           "but nevertheless managed to change the transposed "
                           "state of Op.  On input, Op.UseTranspose() = "
                           << transposed << ", but now, Op.UseTranspose() = "
                           << currentTransposedState << ".  This likely indicates"
                           " a bug in the implementation of the operator Op.");
        return false;
      }

      // Now we need to restore the original transpose state.  We can
      // do this unconditionally because if we're here, Op's state is
      // not transposed.
      {
        const int errcode =
          const_cast<Epetra_Operator&>(Op).SetUseTranspose (false);
        TEUCHOS_TEST_FOR_EXCEPTION(errcode != 0, std::logic_error,
                           "implementsApplyTranspose: Op.SetUseTranspose(true) "
                           "succeeded, but Op.SetUseTranspose(false) failed with "
                           "error code " << errcode << ".  This likely indicates "
                           "a bug in the implementation of the operator Op.");
      }
      // If we made it this far, then Op implements transpose.
      return true;
    }

    /// \class EpetraOperatorTransposeScopeGuard
    /// \brief Safely sets and unsets an Epetra_Operator's transpose state.
    ///
    /// This class uses the RAII (Resource Acquisition Is
    /// Instantiation) technique to set the transpose state of an
    /// Epetra_Operator instance on construction, and restore its
    /// original transpose state on destruction.
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
      /// \brief Constructor.
      ///
      /// \param Op [in/out] The operator whose transpose state to guard.
      ///
      /// \param trans [in] True if we want to set Op to apply the
      ///   transpose, else false.
      EpetraOperatorTransposeScopeGuard (const Epetra_Operator& Op,
                                         const ETrans trans)
        : Op_ (Op),
          originalTransposeFlag_ (Op.UseTranspose ()),
          newTransposeFlag_ (trans != NOTRANS)
      {
        // If necessary, set (or unset) the transpose flag to the value
        // corresponding to 'trans'.
        if (newTransposeFlag_ != originalTransposeFlag_) {

#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
          using std::endl;
          std::ostream& err = *(getErrStream (Op_));
          err << "---- EpetraOperatorTransposeScopeGuard constructor:" << endl
              << "------ Toggling transpose flag from " << originalTransposeFlag_
              << " to " << newTransposeFlag_ << endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

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
                             << etransToString (trans)
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

#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
          using std::endl;
          std::ostream& err = *(getErrStream (Op_));
          err << "---- EpetraOperatorTransposeScopeGuard destructor:" << endl
              << "------ Toggling transpose flag from " << newTransposeFlag_
              << " to " << originalTransposeFlag_ << endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

          const int info =
            const_cast<Epetra_Operator &>(Op_).SetUseTranspose (originalTransposeFlag_);
          TEUCHOS_TEST_FOR_TERMINATION(info != 0,
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
        TEUCHOS_TEST_FOR_TERMINATION(originalTransposeFlag_ != finalTransposeFlag,
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


EpetraMultiVec::EpetraMultiVec (const Epetra_BlockMap& Map_in,
                                double * array,
                                const int numvecs,
                                const int stride)
  : Epetra_MultiVector(BELOSEPETRACOPY, Map_in, array, stride, numvecs)
{
}


EpetraMultiVec::EpetraMultiVec (const Epetra_BlockMap& Map_in,
                                const int numvecs,
                                bool zeroOut)
  : Epetra_MultiVector(Map_in, numvecs, zeroOut)
{
}


EpetraMultiVec::EpetraMultiVec (Epetra_DataAccess CV_in,
                                const Epetra_MultiVector& P_vec,
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
  EpetraMultiVec * ptr_apv = new EpetraMultiVec(BELOSEPETRACOPY, *this, index);
  return ptr_apv; // safe upcast.
}


MultiVec<double>* EpetraMultiVec::CloneViewNonConst ( const std::vector<int>& index )
{
  EpetraMultiVec * ptr_apv = new EpetraMultiVec(BELOSEPETRAVIEW, *this, index);
  return ptr_apv; // safe upcast.
}


const MultiVec<double>*
EpetraMultiVec::CloneView (const std::vector<int>& index) const
{
  EpetraMultiVec * ptr_apv = new EpetraMultiVec(BELOSEPETRAVIEW, *this, index);
  return ptr_apv; // safe upcast.
}


void
EpetraMultiVec::SetBlock (const MultiVec<double>& A,
                          const std::vector<int>& index)
{
  EpetraMultiVec temp_vec(BELOSEPETRAVIEW, *this, index);

  int numvecs = index.size();
  if ( A.GetNumberVecs() != numvecs ) {
    std::vector<int> index2( numvecs );
    for(int i=0; i<numvecs; i++)
      index2[i] = i;
    EpetraMultiVec *tmp_vec =
      dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A));
    TEUCHOS_TEST_FOR_EXCEPTION(tmp_vec==NULL, EpetraMultiVecFailure,
                       "Belos::EpetraMultiVec::SetBlock: Dynamic cast from "
                       "Belos::MultiVec<double> to Belos::EpetraMultiVec "
                       "failed.");
    EpetraMultiVec A_vec(BELOSEPETRAVIEW, *tmp_vec, index2);
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
  Epetra_MultiVector B_Pvec(BELOSEPETRAVIEW, LocalMap, B.values(), B.stride(), B.numCols());

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
  // Check to make sure the input vector of scale factors has the same
  // number of entries as the number of columns in the multivector.
  int numvecs = this->NumVectors();
  TEUCHOS_TEST_FOR_EXCEPTION((int)alpha.size() != numvecs, EpetraMultiVecFailure,
                     "Belos::MultiVecTraits<double,Epetra_MultiVec>::MvScale: "
                     "Vector (alpha) of scaling factors has length "
                     << alpha.size() << ", which is different than the number "
                     "of columns " << numvecs << " in the input multivector "
                     "(mv).");
  int ret = 0;
  std::vector<int> tmp_index( 1, 0 );
  for (int i=0; i<numvecs; i++) {
    Epetra_MultiVector temp_vec(BELOSEPETRAVIEW, *this, &tmp_index[0], 1);
    ret = temp_vec.Scale( alpha[i] );
    TEUCHOS_TEST_FOR_EXCEPTION(ret!=0, EpetraMultiVecFailure,
                       "Belos::MultiVecTraits<double,Epetra_MultiVec>::MvScale: "
                       "Call to Epetra_MultiVector::Scale() returned a nonzero "
                       "error code " << ret << ".");
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
    Epetra_MultiVector B_Pvec(BELOSEPETRAVIEW, LocalMap, B.values(), B.stride(), B.numCols());

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
                     "Belos::EpetraMultiVec::MvDot: Dynamic cast from Belos::"
                     "MultiVec<double> to Belos::EpetraMultiVec failed.");
  if (A_vec && ( (int)b.size() >= A_vec->NumVectors() ) ) {
     int info = this->Dot( *A_vec, &b[0] );
     TEUCHOS_TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure,
                        "Belos::EpetraMultiVec::MvDot: Call to Epetra_Multi"
                        "Vector::Dot() returned a nonzero error code "
                        << info << ".");
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
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                         "Belos::EpetraMultiVec::MvNorm: Invalid norm_type "
                         << norm_type << ".  The current list of valid norm "
                         "types is {OneNorm, TwoNorm, InfNorm}.");
    }
    TEUCHOS_TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure,
                       "Belos::EpetraMultiVec::MvNorm: Call to Epetra_Multi"
                       "Vector::Norm() returned a nonzero error code "
                       << info << ".");
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
#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
  using std::endl;

  std::ostream& err = *(getErrStream (*Epetra_Op));
  err << "Belos::EpetraOp constructor" << endl;
  err << "-- On input: Epetra_Op->UseTranspose() = "
      << Epetra_Op->UseTranspose() << endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING
}

void
EpetraOp::Apply (const MultiVec<double>& x,
                 MultiVec<double>& y,
                 ETrans trans) const
{
#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
  using std::endl;

  std::ostream& err = *(getErrStream (*Epetra_Op));
  err << "-- Belos::EpetraOp::Apply:" << endl
      << "---- Implements Belos::Operator<double>::Apply()" << endl
      << "---- trans input = " << etransToString (trans) << endl
      << "---- On input: Epetra_Op->UseTranspose() = "
      << Epetra_Op->UseTranspose() << endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

  MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
  Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);

  TEUCHOS_TEST_FOR_EXCEPTION(vec_x==NULL || vec_y==NULL, EpetraOpFailure,
                     "Belos::EpetraOp::Apply: x and/or y could not be "
                     "dynamic cast to an Epetra_MultiVector.");

  // Temporarily set the transpose state of Op, if it's not the same
  // as trans, and restore it on exit of this scope.
  EpetraOperatorTransposeScopeGuard guard (*Epetra_Op, trans);

#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
  err << "---- Before calling Apply: Epetra_Op->UseTranspose() = "
      << Epetra_Op->UseTranspose() << endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

  // Apply the operator to x and put the result in y.
  const int info = Epetra_Op->Apply (*vec_x, *vec_y);
  TEUCHOS_TEST_FOR_EXCEPTION(info!=0, EpetraOpFailure,
                     "Belos::EpetraOp::Apply: The Epetra_Operator's Apply() "
                     "method returned a nonzero error code of " << info << ".");
}

bool
EpetraOp::HasApplyTranspose() const
{
#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
  std::ostream& err = *(getErrStream (*Epetra_Op));
  err << "-- Belos::EpetraOp::HasApplyTranspose" << std::endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

  return implementsApplyTranspose (*Epetra_Op);
}


// ///////////////////////////////////////////////////////////////////
//
// Implementation of the Belos::EpetraPrecOp class.
//
// ///////////////////////////////////////////////////////////////////

EpetraPrecOp::EpetraPrecOp (const Teuchos::RCP<Epetra_Operator> &Op)
  : Epetra_Op(Op)
{
#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
  using std::endl;

  std::ostream& err = *(getErrStream (*Epetra_Op));
  err << "-- Belos::EpetraPrecOp constructor" << endl
      << "---- On input: Epetra_Op->UseTranspose() = "
      << Epetra_Op->UseTranspose() << endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING
}

// The version of Apply() that takes an optional 'trans' argument and
// returns void implements the Belos::Operator interface.
void
EpetraPrecOp::Apply (const MultiVec<double>& x,
                     MultiVec<double>& y,
                     ETrans trans) const
{
#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
  using std::endl;

  std::ostream& err = *(getErrStream (*Epetra_Op));
  err << "-- Belos::EpetraPrecOp::Apply:" << endl
      << "---- Implements Belos::Operator<double>::Apply()" << endl
      << "---- trans input = " << etransToString (trans) << endl
      << "---- On input: Epetra_Op->UseTranspose() = "
      << Epetra_Op->UseTranspose() << endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

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

#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
  err << "---- Before calling ApplyInverse: Epetra_Op->UseTranspose() = "
      << Epetra_Op->UseTranspose() << endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

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
#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
  using std::endl;

  std::ostream& err = *(getErrStream (*Epetra_Op));
  err << "-- Belos::EpetraPrecOp::Apply:" << endl
      << "---- Implements Epetra_Operator::Apply()" << endl
      << "---- On input: Epetra_Op->UseTranspose() = "
      << Epetra_Op->UseTranspose() << endl
      << "---- Calling Epetra_Op->ApplyInverse (X, Y)" << endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

  // This operation computes Y = A^{-1}*X.
  const int info = Epetra_Op->ApplyInverse( X, Y );

#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
  err << "---- Epetra_Op->ApplyInverse (X, Y) returned info = " << info << endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

  return info;
}

// This implements Epetra_Operator::ApplyInverse().
int
EpetraPrecOp::ApplyInverse (const Epetra_MultiVector &X,
                            Epetra_MultiVector &Y) const
{
#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
  using std::endl;

  std::ostream& err = *(getErrStream (*Epetra_Op));
  err << "-- Belos::EpetraPrecOp::ApplyInverse:" << endl
      << "---- Implements Epetra_Operator::ApplyInverse()" << endl
      << "---- On input: Epetra_Op->UseTranspose() = "
      << Epetra_Op->UseTranspose() << endl
      << "---- Calling Epetra_Op->Apply (X, Y)" << endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

  // This operation computes Y = A*X.
  const int info = Epetra_Op->Apply( X, Y );

#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
  err << "---- Epetra_Op->Apply (X, Y) returned info = " << info << endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

  return info;
}

bool
EpetraPrecOp::HasApplyTranspose() const
{
#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
  std::ostream& err = *(getErrStream (*Epetra_Op));
  err << "-- Belos::EpetraPrecOp::HasApplyTranspose" << std::endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

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
#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
  using std::endl;

  std::ostream& err = *(getErrStream (Op));
  err << "Belos::OperatorTraits<double, Epetra_MultiVector, Epetra_Operator>::Apply:"
      << endl
      << "-- trans input = " << etransToString (trans) << endl
      << "-- On input: Op.UseTranspose() = "
      << Op.UseTranspose() << endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

  // Temporarily set the transpose state of Op, if it's not the same
  // as trans, and restore it on exit of this scope.
  EpetraOperatorTransposeScopeGuard guard (Op, trans);

#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
  err << "-- Before calling Op.Apply: Op.UseTranspose() = "
      << Op.UseTranspose() << endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

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
#ifdef BELOS_EPETRA_AGGRESSIVE_DEBUGGING
  std::ostream& err = *(getErrStream (Op));
  err << "Belos::OperatorTraits<double, Epetra_MultiVector, Epetra_Operator>::"
    "HasApplyTranspose" << std::endl;
#endif // BELOS_EPETRA_AGGRESSIVE_DEBUGGING

  return implementsApplyTranspose (Op);
}


}  // end namespace Belos
