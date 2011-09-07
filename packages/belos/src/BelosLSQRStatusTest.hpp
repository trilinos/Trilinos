/*
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
*/

#ifndef BELOS_LSQR_STATUS_TEST_HPP
#define BELOS_LSQR_STATUS_TEST_HPP

/*!
  \file BelosLSQRStatusTest.hpp
  \brief Belos::StatusTest class for specifying convergence of LSQR.
*/

#include "BelosStatusTest.hpp"
#include "BelosLSQRIter.hpp"

/*!  \class LSQRStatusTest:
  \brief A Belos::StatusTest class for specifying convergence of LSQR.  The outer status tests passes if an inner
  status passes a user specified number of times consecutively.  The inner status test depends on information
  specificto LSQR iteration.
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
  /*! The constructor has four optional arguments, specifying the upper limit of the apparent condition number of Abar,
 the number of successful convergent iterations, an estimate of the relative error in the data defining the right-hand
 side (b), and an estimate of the relative error in the data defining the coefficinet matrix (A).  The default
 termIterMax is 1, and the other three parameters default to 0.
   */
  LSQRStatusTest( MagnitudeType cond_lim = 0.0, int term_iter_max = 1, MagnitudeType rel_rhs_err = 0.0, MagnitudeType rel_mat_err = 0.0 );

  //! Destructor
  virtual ~LSQRStatusTest();
  //@}

 //! @name Form and parameter definition methods.
  //@{ 

  //! Set the value of the tolerance
  /*! Resetting the limit of the condition number of Abar is allowed in cases where, in the process of testing
   *  convergence, the initial is found to be either too tight or too lax.
  */
  int setCondLim(MagnitudeType cond_lim) {
    cond_lim_ = cond_lim; 
    cond_tol_ = (cond_lim > 0) ? (Teuchos::ScalarTraits< MagnitudeType >::one() / cond_lim) : Teuchos::ScalarTraits< MagnitudeType >::eps(); 
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

  //@}

  //! @name Print methods
  //@{

  //! Output formatted description of stopping test to output stream.
  void print(std::ostream& os, int indent = 0) const;

  //! Print message for each status specific to this stopping test.
  void printStatus(std::ostream& os, Belos::StatusType type) const;

  //@}

  //! @name Methods to access data members.
  //@{ 

  //! Returns the value of the upper limit of the condition number of Abar set in the constructor.
  MagnitudeType getCondLim() const {return(cond_lim_);};

  //! Returns the number of successful convergent iterations required set in the constructor.
  int getTermIterMax() const {return(term_iter_max_);};

  //! Returns the value of the estimate of the relative error in the data defining b set in the constructor.
  MagnitudeType getRelRhsErr() const {return(rel_rhs_err_);};

  //! Returns the value of the estimate of the relative error in the data defining A set in the constructor.
  MagnitudeType getMatErr() const {return(rel_mat_err_);};

  //@}

  /** @name Misc. */
  //@{

  /** \brief Call to setup initialization.
   */
  Belos::StatusType firstCallCheckStatusSetup(Belos::Iteration<ScalarType,MV,OP>* iSolver);
  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{

  /** \brief Method to return description of the maximum iteration status test  */
  std::string description() const
  {
    std::ostringstream oss;
    oss << "LSQRStatusTest<>: [ limit of condition number = " << cond_lim_ << " ]";
    return oss.str();
  }
  //@}

private:

  //! @name Private data members.
  //@{

  //! Upper limit of condition number of Abar
  MagnitudeType cond_lim_;

  //! How many iterations in a row a passing test for convergence is required.
  int term_iter_max_;

  //! Error in data defining b
  MagnitudeType rel_rhs_err_;

  //! Error in data defining A
  MagnitudeType rel_mat_err_;

  //! Tolerance that defines convergence, the reciprocal of cond_lim_ or, if that is zero, machine epsilon
  MagnitudeType cond_tol_;

  //! Status
  Belos::StatusType status_;

  //! Is this the first time CheckStatus is called?
  bool firstcallCheckStatus_;

  //! How many iterations in a row a test for convergence has passed.
  int term_iter_;

  //@}

};

template <class ScalarType, class MV, class OP>
LSQRStatusTest<ScalarType,MV,OP>::LSQRStatusTest( MagnitudeType cond_lim /* = 0 */, int term_iter_max /* = 1 */, MagnitudeType rel_rhs_err /* = 0 */, MagnitudeType rel_mat_err /* = 0 */)
  : cond_lim_(cond_lim),
    term_iter_max_(term_iter_max),
    rel_rhs_err_(rel_rhs_err),
    rel_mat_err_(rel_mat_err),
    status_(Belos::Undefined),
    firstcallCheckStatus_(true)
{}

template <class ScalarType, class MV, class OP>
LSQRStatusTest<ScalarType,MV,OP>::~LSQRStatusTest()
{}

template <class ScalarType, class MV, class OP>
void LSQRStatusTest<ScalarType,MV,OP>::reset()
{
  status_ = Belos::Undefined;
  firstcallCheckStatus_ = true;
}

template <class ScalarType, class MV, class OP>
Belos::StatusType LSQRStatusTest<ScalarType,MV,OP>::firstCallCheckStatusSetup(Belos::Iteration<ScalarType,MV,OP>* iSolver)
{
  if (firstcallCheckStatus_) {
      firstcallCheckStatus_ = false;
      const MagnitudeType MTzero = Teuchos::ScalarTraits<MagnitudeType>::zero();
      const MagnitudeType MTone = Teuchos::ScalarTraits<MagnitudeType>::one();
      term_iter_ = 0;
      if (cond_lim_ > MTzero ) 
	cond_tol_ = MTone / cond_lim_;
      else
	cond_tol_ = Teuchos::ScalarTraits< MagnitudeType >::eps();
    }
    return Belos::Undefined;
}

template <class ScalarType, class MV, class OP>
Belos::StatusType LSQRStatusTest<ScalarType,MV,OP>::checkStatus( Belos::Iteration<ScalarType,MV,OP>* iSolver) 
{
  if (firstcallCheckStatus_) {
    Belos::StatusType status = firstCallCheckStatusSetup(iSolver);
    if(status==Belos::Failed) {
      status_ = Belos::Failed;
      return(status_);
    }
  }

  const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
  const MagnitudeType zero = Teuchos::ScalarTraits<MagnitudeType>::zero();

  bool flag = false;
  LSQRIter<ScalarType,MV,OP>* solver = dynamic_cast< LSQRIter<ScalarType,MV,OP>* > (iSolver);
  LSQRIterationState< ScalarType, MV > state = solver->getState();

  // This section computes the three stopping criteria
  ScalarType stop_crit_1 = state.resid_norm / state.bnorm;
  ScalarType stop_crit_2 = (state.resid_norm > zero) ? state.mat_resid_norm / (state.frob_mat_norm * state.resid_norm) : zero;
  ScalarType stop_crit_3 = one / state.mat_cond_num;
  ScalarType resid_tol = rel_rhs_err_ + rel_mat_err_ * state.mat_resid_norm * state.sol_norm / state.bnorm;
  ScalarType resid_tol_mach = Teuchos::ScalarTraits< MagnitudeType >::eps() + Teuchos::ScalarTraits< MagnitudeType >::eps() * state.mat_resid_norm * state.sol_norm / state.bnorm;

  // This section checks if any stopping criteria have been met
  if (stop_crit_1 <= resid_tol || stop_crit_2 <= rel_mat_err_ || stop_crit_3 <= cond_tol_ || stop_crit_1 <= resid_tol_mach || stop_crit_2 <= Teuchos::ScalarTraits< MagnitudeType >::eps() || stop_crit_3 <= Teuchos::ScalarTraits< MagnitudeType >::eps()) {
    flag = true;
  }


  // history requirement:
  // converged if thresholds met for user specified number of consecutive iterations
  if (!flag) {
    term_iter_ = -1;
  }
  term_iter_++;
  status_ = (term_iter_ < term_iter_max_) ? Belos::Failed : Belos::Passed;
  return status_;
}

template <class ScalarType, class MV, class OP>
void LSQRStatusTest<ScalarType,MV,OP>::print(std::ostream& os, int indent) const
{
  for (int j = 0; j < indent; j++)
    os << ' ';
  printStatus(os, status_);
  os << "limit of condition number = " << cond_lim_ << std::endl;
}

template <class ScalarType, class MV, class OP>
void LSQRStatusTest<ScalarType,MV,OP>::printStatus(std::ostream&os, Belos::StatusType type) const
{
  os << std::left << std::setw(13) << std::setfill('.');
  switch (type) {
  case Belos::Passed:
    os << "OK";
    break;
  case Belos::Failed:
    os << "Failed";
    break;
  case Belos::Undefined:
  default:
    os << "**";
    break;
  }
  os << std::left << std::setfill(' ');
  return;
}

} // end Belos namespace


#endif /* BELOS_LSQR_STATUS_TEST_HPP */
