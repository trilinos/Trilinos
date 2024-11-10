/*
// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
*/

#ifndef BELOS_LSQR_STATUS_TEST_HPP
#define BELOS_LSQR_STATUS_TEST_HPP

/*!
  \file BelosLSQRStatusTest.hpp
  \brief Belos::StatusTest class defining LSQR convergence
*/

#include "BelosStatusTest.hpp"
#include "BelosLSQRIter.hpp"

/*!  \class Belos::LSQRStatusTest
  \brief A Belos::StatusTest class for specifying convergence of LSQR.  The outer status tests passes if an inner
  status passes a user specified number of times consecutively.  The inner status test depends on information
  specific to LSQR iteration.
*/

namespace Belos {


template <class ScalarType, class MV, class OP>
class LSQRStatusTest: public Belos::StatusTest<ScalarType,MV,OP> {

public:

  // Convenience typedefs
  typedef Teuchos::ScalarTraits<ScalarType> SCT;
  typedef typename SCT::magnitudeType MagnitudeType;
  typedef Belos::MultiVecTraits<ScalarType,MV>  MVT;

  //! @name Constructor/Destructor.
  //@{

  //! Constructor
  /*! The constructor has four optional arguments, specifying the maximum observed condition number of Abar,
 the number of successive convergent iterations, an estimate of the relative error in the data defining the right-hand
 side (b), and an estimate of the relative error in the data defining the coefficinet matrix (A).  The default
 termIterMax is 1, and the other three parameters default to 0.  The defaults specified in LSQRSolMgr are more realistic.
   */
  LSQRStatusTest( MagnitudeType condMax = 0.0,
                  int term_iter_max = 1,
                  MagnitudeType rel_rhs_err = 0.0,
                  MagnitudeType rel_mat_err = 0.0 );

  //! Destructor
  virtual ~LSQRStatusTest();
  //@}

  //! @name Status method
  //@{

  //! Check convergence status of the iterative solver: Unconverged, Converged, Failed.
  /*! This method checks if the convergence criteria are met using the current information from the iterative solver.
   */
  Belos::StatusType checkStatus(Belos::Iteration<ScalarType,MV,OP> *iSolver );

  //! Return the result of the most recent CheckStatus call.
  Belos::StatusType getStatus() const {return(status_);}

  //@}

  //! @name Reset methods
  //@{

  //! Resets the status test to the initial internal state.
  void reset();

  //! Set the tolerances
  int setCondLim(MagnitudeType condMax) {
    condMax_ = condMax;
    rcondMin_ = (condMax > 0) ? (Teuchos::ScalarTraits< MagnitudeType >::one() / condMax) : Teuchos::ScalarTraits< MagnitudeType >::eps();
    return(0);}

  int setTermIterMax(int term_iter_max) {
    term_iter_max_ = term_iter_max;
    if (term_iter_max_ < 1)
      term_iter_max_ = 1;
    return(0);}

  int setRelRhsErr(MagnitudeType rel_rhs_err) {
    rel_rhs_err_ = rel_rhs_err;
    return(0);}

  int setRelMatErr(MagnitudeType rel_mat_err) {
    rel_mat_err_ = rel_mat_err;
    return(0);}

  //@}

  //! @name Accessor methods
  //@{

  //! Returns the value of the upper limit of the condition number of Abar set in the constructor.
  MagnitudeType getCondMaxLim() const {return(condMax_);}

  //! Returns the number of successful convergent iterations required set in the constructor.
  int getTermIterMax() const {return(term_iter_max_);}

  //! Returns the value of the estimate of the relative error in the data defining b set in the constructor.
  MagnitudeType getRelRhsErr() const {return(rel_rhs_err_);}

  //! Returns the value of the estimate of the relative error in the data defining A set in the constructor.
  MagnitudeType getMatErr() const {return(rel_mat_err_);}

  //! Returns the value of the observed condition number of Abar
  MagnitudeType getMatCondNum() const {return(matCondNum_);}

  //! Returns the value of the observed (Frobenius) norm of A
  MagnitudeType getMatNorm() const {return(matNorm_);}

  //! !Returns the current number of successful iterations from the most recent StatusTest call.
  int getTermIter() const { return term_iter_; }

  //! Returns the value of the observed norm of the residual r = b-Ax
  MagnitudeType getResidNorm() const {return(resNorm_);}

  //! Returns the value of the observed norm of the Least Squares residual A^T r
  MagnitudeType getLSResidNorm() const {return(matResNorm_);}
  //@}


  //! @name Print methods
  //@{

  //! Output formatted description of stopping test to output stream.
  void print(std::ostream& os, int indent = 0) const;

  //! Print message for each status specific to this stopping test.
  void printStatus(std::ostream& os, Belos::StatusType type) const;

  //@}

  /** @name Misc. */
  //@{

  /// \brief Called in checkStatus exactly once, on the first call to checkStatus.
  ///
  Belos::StatusType firstCallCheckStatusSetup(Belos::Iteration<ScalarType,MV,OP>* iSolver);
  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{

  /** \brief Method to return description of the maximum iteration status test  */
  std::string description() const
  {
    std::ostringstream oss;
    oss << "LSQRStatusTest<>: [ limit of condition number = " << condMax_ << " ]";
    return oss.str();
  }
  //@}

private:

  //! @name Private data members.
  //@{

  //! Upper limit of condition number of Abar
  MagnitudeType condMax_;

  //! How many iterations in a row a passing test for convergence is required.
  int term_iter_max_;

  //! Error in data defining b
  MagnitudeType rel_rhs_err_;

  //! Error in data defining A
  MagnitudeType rel_mat_err_;

  //! One of the tolerances defining convergence, the reciprocal of condMax_ or, if that is zero, machine epsilon
  MagnitudeType rcondMin_;

  //! Status
  Belos::StatusType status_;

  // term_iter_ records the number of consecutive "successful" iterations.
  // convergence requires that term_iter_max consecutive iterates satisfy the other convergence tests
  int term_iter_;

  // condition number of the operator
  MagnitudeType matCondNum_;

  // Frobenius norm of the operator
  MagnitudeType matNorm_;

  // residual norm for the linear system
  MagnitudeType resNorm_;

  // least squares residual, operator^Transpose * residual
  MagnitudeType matResNorm_;

  //@}

};

template <class ScalarType, class MV, class OP>
LSQRStatusTest<ScalarType,MV,OP>::
LSQRStatusTest (MagnitudeType condMax /* = 0 */,
                int term_iter_max /* = 1 */,
                MagnitudeType rel_rhs_err /* = 0 */,
                MagnitudeType rel_mat_err /* = 0 */)
  : condMax_(condMax),
    term_iter_max_ (term_iter_max),
    rel_rhs_err_ (rel_rhs_err),
    rel_mat_err_ (rel_mat_err),
    rcondMin_ ( Teuchos::ScalarTraits<MagnitudeType>::zero() ),
    status_ (Belos::Undefined),
    term_iter_ (0),
    matCondNum_ ( Teuchos::ScalarTraits<MagnitudeType>::one() ),
    matNorm_ ( Teuchos::ScalarTraits<MagnitudeType>::zero() ),
    resNorm_  ( Teuchos::ScalarTraits<MagnitudeType>::zero() ),
    matResNorm_ ( Teuchos::ScalarTraits<MagnitudeType>::zero() )
{}

template <class ScalarType, class MV, class OP>
LSQRStatusTest<ScalarType,MV,OP>::~LSQRStatusTest()
{}

template <class ScalarType, class MV, class OP>
void LSQRStatusTest<ScalarType,MV,OP>::reset()
{
  status_ = Belos::Undefined;
}

template <class ScalarType, class MV, class OP>
Belos::StatusType LSQRStatusTest<ScalarType,MV,OP>::checkStatus( Belos::Iteration<ScalarType,MV,OP>* iSolver)
{
  const MagnitudeType MTzero = Teuchos::ScalarTraits<MagnitudeType>::zero();
  const MagnitudeType MTone = Teuchos::ScalarTraits<MagnitudeType>::one();
  if (condMax_ > MTzero )
    {
        rcondMin_ = MTone / condMax_;
    }
  else
    {
        rcondMin_ = Teuchos::ScalarTraits< MagnitudeType >::eps();
    }

  bool termIterFlag = false;
  LSQRIter<ScalarType,MV,OP>* solver = dynamic_cast< LSQRIter<ScalarType,MV,OP>* > (iSolver);
  TEUCHOS_ASSERT(solver != NULL);
  LSQRIterationState< ScalarType, MV > state = solver->getState();
  //
  //   LSQR solves a least squares problem.  A converged preconditioned residual norm
  // suffices for convergence, but is not necessary.  LSQR sometimes returns a larger
  // relative residual norm than what would have been returned by a linear solver.
  // This section evaluates three stopping criteria.  In the Solver Manager, this test
  // is combined with a generic number of iteration test.
  //   If the linear system includes a preconditioner, then the least squares problem
  // is solved for the preconditioned linear system.  Preconditioning changes the least
  // squares problem (in the sense of changing the norms), and the solution depends
  // on the preconditioner in this sense.
  //   In the context of Linear Least Squares problems, preconditioning refers
  // to the regularization matrix.  Here the regularization matrix is always a scalar
  // multiple of the identity (standard form least squres).
  //   The "loss of accuracy" concept is not yet implemented here, becuase it is unclear
  // what this means for linear least squares.  LSQR solves an inconsistent system
  // in a least-squares sense.  "Loss of accuracy" would correspond to
  // the difference between the preconditioned residual and the unpreconditioned residual.
  //

  std::cout << " X " << state.sol_norm
            << "  b-AX " << state.resid_norm
            << "  Atr  " << state.mat_resid_norm
            << "  A " << state.frob_mat_norm
            << "  cond  " << state.mat_cond_num
            << "  relResNorm " << state.resid_norm/state.bnorm
            << "  LS " << state.mat_resid_norm /( state.resid_norm * state.frob_mat_norm )
            << std::endl;

  const MagnitudeType zero = Teuchos::ScalarTraits<MagnitudeType>::zero();
  const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
  ScalarType stop_crit_1 = zero; // b = 0, done
  if( state.bnorm  > zero )
    {
      stop_crit_1 = state.resid_norm / state.bnorm;
    }
  ScalarType stop_crit_2 = zero;
  if( state.frob_mat_norm  > zero  && state.resid_norm > zero )
    {
      stop_crit_2 = (state.resid_norm > zero) ? state.mat_resid_norm / (state.frob_mat_norm * state.resid_norm) : zero;
    }
  else
    {
     if( state.resid_norm == zero )
       {
         stop_crit_2 = zero;
       }
     else
       {
         stop_crit_2 = one; // Initial mat_norm always vanishes
       }
    }
  ScalarType stop_crit_3 = one / state.mat_cond_num;
  ScalarType resid_tol = rel_rhs_err_ + rel_mat_err_ * state.frob_mat_norm * state.sol_norm / state.bnorm;
  ScalarType resid_tol_mach = Teuchos::ScalarTraits< MagnitudeType >::eps() + Teuchos::ScalarTraits< MagnitudeType >::eps() * state.frob_mat_norm * state.sol_norm / state.bnorm;

  // The expected use case for our users is that the linear system will almost
  // always be compatible, but occasionally may not be.  However, some users
  // may use LSQR for more general cases.  This is why we include the full
  // suite of tests, for both compatible and incompatible systems.
  //
  // Users will have to be educated that sometimes they will get an answer X
  // that does _not_ satisfy the linear system AX=B, but _does_ satisfy the
  // corresponding least-squares problem.  Perhaps the solution manager should
  // provide them with a way to find out.

  // stop_crit_1 is for compatible linear systems.
  // stop_crit_2 is for incompatible linear systems.
  // stop_crit_3 is for either compatible or incompatible linear systems.

  // Have we met any of the stopping criteria?
  if (stop_crit_1 <= resid_tol || stop_crit_2 <= rel_mat_err_ || stop_crit_3 <= rcondMin_ || stop_crit_1 <= resid_tol_mach || stop_crit_2 <= Teuchos::ScalarTraits< MagnitudeType >::eps() || stop_crit_3 <= Teuchos::ScalarTraits< MagnitudeType >::eps()) {
    termIterFlag = true;

    if (stop_crit_1 <= resid_tol )
      std::cout << "Conv: stop_crit_1  " << stop_crit_1  << " resid_tol " << resid_tol << std::endl;

    if (stop_crit_1 <=  resid_tol_mach )
      std::cout << "Conv: stop_crit_1  " << stop_crit_1  << " resid_tol_mach " << resid_tol_mach << std::endl;

    if (stop_crit_2 <= rel_mat_err_ )
      std::cout << "Conv: stop_crit_2  " << stop_crit_2  << " rel_mat_err " << rel_mat_err_ << std::endl;

    if (stop_crit_2 <=   Teuchos::ScalarTraits< MagnitudeType >::eps() )
      std::cout << "Conv: stop_crit_2  " << stop_crit_2  << " eps " <<   Teuchos::ScalarTraits< MagnitudeType >::eps()   << std::endl;

    if (stop_crit_3 <= rcondMin_ )
      std::cout << "Conv: stop_crit_3  " << stop_crit_3  << " rcondMin_ " << rcondMin_ << std::endl;

    if (stop_crit_3 <=   Teuchos::ScalarTraits< MagnitudeType >::eps() )
      std::cout << "Conv: stop_crit_3  " << stop_crit_3  << " eps " <<   Teuchos::ScalarTraits< MagnitudeType >::eps()   << std::endl;
  }

  // update number of consecutive successful iterations
  if (!termIterFlag) {
    term_iter_ = 0;
  } else {
    term_iter_++;
  }
  status_ = (term_iter_ < term_iter_max_) ? Belos::Failed : Belos::Passed;

  matCondNum_ = state.mat_cond_num; // information that defined convergence
  matNorm_ = state.frob_mat_norm;   // in accessible variables
  resNorm_  = state.resid_norm;
  matResNorm_ = state.mat_resid_norm;

  return status_;
}

template <class ScalarType, class MV, class OP>
void LSQRStatusTest<ScalarType,MV,OP>::print(std::ostream& os, int indent) const
{
  for (int j = 0; j < indent; j++)
    os << ' ';
  printStatus(os, status_);
  os << "limit of condition number = " << condMax_ << std::endl;
  os << "limit of condition number = " << condMax_ << std::endl;
}

template <class ScalarType, class MV, class OP>
void LSQRStatusTest<ScalarType,MV,OP>::printStatus(std::ostream&os, Belos::StatusType type) const
{
  os << std::left << std::setw(13) << std::setfill('.');
  switch (type) {
  case Belos::Passed:
    os << "Passed";
    break;
  case Belos::Failed:
    os << "Failed";
    break;
  case Belos::Undefined:
  default:
    os << "Undefined";
    break;
  }
  os << std::left << std::setfill(' ');
  return;
}

} // end Belos namespace


#endif /* BELOS_LSQR_STATUS_TEST_HPP */
