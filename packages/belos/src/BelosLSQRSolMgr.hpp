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

#ifndef BELOS_LSQR_SOLMGR_HPP
#define BELOS_LSQR_SOLMGR_HPP

/// \file BelosLSQRSolMgr.hpp
/// \brief LSQRSolMgr: interface to the LSQR method.

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "BelosLSQRIteration.hpp"
#include "BelosLSQRIter.hpp"
#include "BelosOrthoManagerFactory.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosLSQRStatusTest.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosStatusTestOutputFactory.hpp"
#include "BelosOutputManager.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"

#ifdef BELOS_TEUCHOS_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

namespace Belos {

  
//! @name LSQRSolMgr Exceptions
//@{

/** \brief Belos::LSQRSolMgrLinearProblemFailure is thrown when the linear problem is
 * not setup (i.e. setProblem() was not called) when solve() is called.
 *
 * This std::exception is thrown from the LSQRSolMgr::solve() method.
 *
 */
class LSQRSolMgrLinearProblemFailure : public BelosError {
public:
  LSQRSolMgrLinearProblemFailure(const std::string& what_arg)
    : BelosError(what_arg)
  {}
};

/** \brief LSQRSolMgrOrthoFailure is thrown when the orthogonalization manager is
 * unable to generate orthonormal columns from the initial basis vectors.
 *
 * This std::exception is thrown from the LSQRSolMgr::solve() method.
 *
 */
class LSQRSolMgrOrthoFailure : public BelosError {
public:
  LSQRSolMgrOrthoFailure(const std::string& what_arg)
    : BelosError(what_arg)
  {}
};

/** \brief LSQRSolMgrBlockSizeFailure is thrown when the linear problem has
 * more than one RHS.  This is unique to single vector methods.
 *
 * This std::exception is thrown from the LSQRSolMgr::solve() method.
 *
 */
class LSQRSolMgrBlockSizeFailure : public BelosError {
public:
  LSQRSolMgrBlockSizeFailure(const std::string& what_arg) 
  : BelosError(what_arg)
  {}
};

/// \class LSQRSolMgr
/// \brief LSQR method (for linear systems and linear least-squares problems).
/// \author Sarah Knepper and David Day
/// 
/// \tparam ScalarType The type of entries in the right-hand side
///   vector(s) \f$b\f$ and solution vector(s) \f$x\f$.
/// \tparam MV The multivector type; the type of the solution
///   vector(s) and right-hand side vector(s).
/// \tparam OP The type of the matrix \f$A\f$ (and any preconditioner,
///   if one is provided).
///
/// \warning Our LSQR implementation currently only compiles correctly
///   for real-valued (not complex) ScalarType types.  You may check
///   whether ScalarType is complex using the following code:
///   \code
///   if (Teuchos::ScalarTraits<ScalarType>::isComplex) {
///     // ScalarType is complex valued.
///   } else {
///     // ScalarType is real valued.
////  }
///   \endcode
///   This may be fixed in future releases.  It is not a limitation of
///   the LSQR method itself, just of our current implementation.  For
///   now, you will not be able to compile any specialization of
///   LSQRSolMgr for a complex-valued ScalarType type.
///
/// \section Algorithm
///
/// LSQR (Paige and Saunders; see References) is an iterative method
/// for solving linear least-squares problems and linear systems.  It
/// can solve any of the following problems:
///
/// 1. Solve \f$Ax=b\f$ for \f$x\f$
/// 2. Find \f$x\f$ that minimizes \f$\|Ax - b\|_2^2\f$
/// 3. Find \f$x\f$ that minimizes \f$\|Ax - b\|_2^2 + \lambda^2 \|x\|_2^2\f$
///
/// The third problem above is the most general and includes the
/// previous two.  This is the problem LSQR actually solves.  Here,
/// \f$\lambda\f$ is a user-provided positive real constant (the
/// "damping parameter") which regularizes the problem so that it
/// always has a bounded solution, even if \f$A\f$ does not have full
/// rank.
///
/// In the words of Paige and Saunders: "The method is based on the
/// Golub-Kahan bidiagonalization process. It is algebraically
/// equivalent to applying MINRES to the normal equation[s] \f$(A^T A
/// + \lambda 2I)x = A^T b\f$, but has better numerical properties,
/// especially if \f$A\f$ is ill-conditioned."  
///
/// LSQR has some special algorithmic properties:
/// 
/// 1. It reduces \f$\|b - A x\|_2\f$ (the two-norm of the residual)
///    monotonically.
/// 2. LSQR also computes a monotonically increasing estimate of the
///    two-norm condition number of the matrix \f$A\f$.
/// 
/// Property #2 makes LSQR useful for mixed-precision algorithms.  If
/// the matrix \f$A\f$ has condition number greater than the inverse
/// of machine precision in the current working precision, one can
/// reconstruct the problem to solve in the next higher precision and
/// restart, possibly using the previous solution as an initial guess.
///
/// \section Preconditioning
///
/// If the linear problem to solve includes a preconditioner (in the
/// LinearProblem object), then the least-squares problem is solved
/// for the preconditioned linear system.  Preconditioning changes the
/// least-squares problem (in the sense of changing the norms), and
/// the solution depends on the preconditioner in this sense.  In the
/// context of linear least-squares problems, "preconditioning" refers
/// to the regularization matrix.  In this solver, the regularization
/// matrix is always a scalar multiple of the identity (standard form
/// least squares).
///
/// A converged preconditioned residual norm suffices for convergence,
/// but is not necessary.  LSQR sometimes returns a larger relative
/// residual norm than what would have been returned by a linear
/// solver.  For details on the stopping criteria, see the
/// documentation of \c LSQRStatusTest, which implements the
/// three-part stopping criterion recommended by Paige and Saunders.
///
/// Some Belos solvers implement detection of "loss of accuracy."
/// That refers to the difference between convergence of the original
/// linear system and convergence of the (left-)preconditioned linear
/// system.  LSQR does not implement detection of "loss of accuracy,"
/// because it is unclear what this means for linear least squares in
/// general.  This LSQR solves a possibly inconsistent system in a
/// least-squares sense.  
///
/// \section References
/// 
/// C. C. Paige and M. A. Saunders, LSQR: An algorithm for sparse
/// linear equations and sparse least squares, TOMS 8(1), 43-71
/// (1982).
///
/// C. C. Paige and M. A. Saunders, Algorithm 583; LSQR: Sparse linear
/// equations and least-squares problems, TOMS 8(2), 195-209 (1982).
///
/// See also the <a
/// href="http://www.stanford.edu/group/SOL/software/lsqr.html">LSQR
/// web page.</a>
///
template<class ScalarType, class MV, class OP>
class LSQRSolMgr : public SolverManager<ScalarType,MV,OP> {
  
private:
  typedef MultiVecTraits<ScalarType,MV> MVT;
  typedef OperatorTraits<ScalarType,MV,OP> OPT;
  typedef Teuchos::ScalarTraits<ScalarType> SCT;
  typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
  typedef Teuchos::ScalarTraits<MagnitudeType> MT;
  
public:
  
  //! @name Construct/Destroy
  //@{ 
  
  /*! \brief Empty constructor for LSQRSolMgr.
   * This constructor takes no arguments and sets the default values for the solver.
   * The linear problem must be passed in using setProblem() before solve() is called
   * on this object.
   * The solver values can be changed using setParameters().
   */
  LSQRSolMgr();
  
  /*! \brief Basic constructor for LSQRSolMgr.
   *
   * This constructor accepts the LinearProblem to be solved in addition
   * to a parameter list of options for the solver manager.  Blocks of size > 1 are not
   * implemented.  The options are otherwise the BlockGmres options.
   *   - "Maximum Iterations" - the maximum number of iterations the LSQR solver
   *                            is allowed to perform. Default: 1000
   *   - "Condition Limit" - a \c MagnitudeType specifying the upper limit of the estimate of
   *                         the norm of Abar to decide convergence. Default: 0.
   *   - "Term Iter Max" - the number of consecutive successful iterations required before
   *                       convergence is declared.  Default: 1.
   *   - "Rel RHS Err" - an estimate of the error in the data defining the RHS.
   *                     Default: 10*sqrt(eps).
   *   - "Rel Mat Err" - an estimate of the error in the data defining the matrix.
   *                     Default: 10*sqrt(eps).
   *   - "Orthogonalization" - a string specifying the desired orthogonalization
   *                           method.  Default: "DGKS".  See \c OrthoManagerFactory
   *                           for a list of the available orthogonalization methods.
   *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Belos::Errors
   *   - "Output Style" - a OutputType specifying the style of output. Default: Belos::General
   *   - "Lambda"  - a \c MagnitudeType that specifies the regularization parameter.
   */
  // This LSQR implementation only supports block size 1.  Like CG, LSQR is a short
  // recurrence method that, in finite precision arithmetic and without reorthogonalization, 
  // does not have the "n" step convergence property.  Without either blocks or 
  // reorthogonalization, there is nothing to "Orthogonalize."

  LSQRSolMgr( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
		 const Teuchos::RCP<Teuchos::ParameterList> &pl );
  
  //! Destructor.
  virtual ~LSQRSolMgr() {};
  //@}
  
  //! @name Accessor methods
  //@{ 

  /*! \brief Get current linear problem being solved for in this object.
   */
  const LinearProblem<ScalarType,MV,OP>& getProblem() const {
    return *problem_;
  }
  
  /*! \brief Get a parameter list containing the valid parameters for this object.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  
  /*! \brief Get a parameter list containing the current parameters for this object.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters() const { return params_; }
  
  /*! \brief Return the timers for this object. 
   *
   * The timers are ordered as follows:
   *   - time spent in solve() routine
   */
  Teuchos::Array<Teuchos::RCP<Teuchos::Time> > getTimers() const {
    return Teuchos::tuple(timerSolve_);
  }
  
  //! Iteration count from the last solve.
  int getNumIters() const {
    return numIters_;
  }

  /// \brief Estimated matrix condition number from the last solve.
  ///
  /// LSQR computes a running condition number estimate of the
  /// (preconditioned, if applicable) operator.
  MagnitudeType getMatCondNum () const {
    return matCondNum_; 
  }

  /// \brief Estimated matrix Frobenius norm from the last solve.
  ///
  /// LSQR computes a running Frobenius norm estimate of the
  /// (preconditioned, if applicable) operator.
  MagnitudeType getMatNorm () const {
    return matNorm_;
  }

  /// \brief Estimated residual norm from the last solve.
  ///
  /// LSQR computes the current residual norm.  LSQR can solve
  /// inconsistent linear systems in a least-squares sense, so the
  /// residual norm may not necessarily be small, even if LSQR
  /// converges.  (LSQR defines "convergence" to allow for possibly
  /// inconsistent systems.  See the documentation of \c
  /// LSQRStatusTest for details.)
  MagnitudeType getResNorm () const {
    return resNorm_;
  }

  //! Estimate of \f$A^* r\f$ (residual vector \f$r\f$) from the last solve.
  MagnitudeType getMatResNorm () const {
    return matResNorm_;
  }
  
  /// \brief Whether a loss of accuracy was detected during the last solve.
  ///
  /// The "loss of accuracy" concept is not yet implemented here,
  /// becuase it is unclear what this means for linear least squares.
  /// LSQR solves a possibly inconsistent linear system in a
  /// least-squares sense.  "Loss of accuracy" would correspond to the
  /// difference between the preconditioned residual and the
  /// unpreconditioned residual.
  bool isLOADetected() const { return false; }
  
  //@}
  
  //! @name Set methods
  //@{
   
  //! Set the linear problem that needs to be solved. 
  void setProblem( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem ) { problem_ = problem; }
  
  //! Set the parameters the solver manager should use to solve the linear problem. 
  void setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params );
    
  //@}
   
  //! @name Reset methods
  //@{
  /*! \brief reset the solver manager as specified by the \c ResetType, informs the
   *  solver manager that the solver should prepare for the next call to solve
   *  by resetting certain elements of the iterative solver strategy.
   */
  void reset( const ResetType type ) { if ((type & Belos::Problem) && !Teuchos::is_null(problem_)) problem_->setProblem(); }
  //@}
 
  //! @name Solver application methods
  //@{ 
    
  /*! \brief method that performs possibly repeated calls to the underlying linear solver's 
   *         iterate() routine until the problem has been solved (as defined by the solver
   *         manager) or the solver manager decides to quit.
   *
   * This method calls LSQRIter::iterate(), which will return either because a 
   * specially constructed status test evaluates to ::Passed or an std::exception is thrown.
   *
   * A return from LSQRIter::iterate() signifies that either
   *    - the maximum number of iterations has been exceeded ...  "return ::Unconverged". 
   *    - ... or convergence ... "solver manager will return ::Converged"
   * In either case the current solution is in the linear problem
   *
   * \returns ::ReturnType specifying:
   *     - ::Converged: the linear problem was solved to the specification required by the
   *                      solver manager.
   *     - ::Unconverged: the linear problem was not solved to the specification desired by
   *                      the solver manager.
   */
  ReturnType solve();
    
  //@}
    
  /** \name Overridden from Teuchos::Describable */
  //@{
    
  /** \brief Method to return description of the LSQR solver manager */
  std::string description() const;
    
  //@}
    
private:

  // Linear problem.
  Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > problem_;
  
  // Output manager.
  Teuchos::RCP<OutputManager<ScalarType> > printer_;
  Teuchos::RCP<std::ostream> outputStream_;
  
  // Status test.
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > sTest_;
  Teuchos::RCP<StatusTestMaxIters<ScalarType,MV,OP> > maxIterTest_;
  Teuchos::RCP<LSQRStatusTest<ScalarType,MV,OP> > convTest_;
  Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputTest_;
  
  // Orthogonalization manager.
  Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > ortho_; 
    
  // Current parameter list.
  Teuchos::RCP<Teuchos::ParameterList> params_;
  /// \brief Default parameter list.  
  ///
  /// Cached per instance, rather than per class, for more thread
  /// safety.  It's "mutable" because \c getValidParameters() has to
  /// create it if it hasn't been created yet.
  mutable Teuchos::RCP<const Teuchos::ParameterList> validParams_;

  // Current solver input parameters
  MagnitudeType lambda_;
  MagnitudeType relRhsErr_;
  MagnitudeType relMatErr_;
  MagnitudeType condMax_;
  int maxIters_, termIterMax_;
  std::string orthoType_; 
  MagnitudeType orthoKappa_;
  int verbosity_, outputStyle_, outputFreq_;

  // Terminal solver state values
  int numIters_;
  MagnitudeType matCondNum_;
  MagnitudeType matNorm_;
  MagnitudeType resNorm_;
  MagnitudeType matResNorm_;
    
  // Timers.
  std::string label_;
  Teuchos::RCP<Teuchos::Time> timerSolve_;

  // Internal state variables.
  bool isSet_;
  bool loaDetected_;
};

// Empty Constructor
template<class ScalarType, class MV, class OP>
LSQRSolMgr<ScalarType,MV,OP>::LSQRSolMgr() :
  isSet_(false),
  loaDetected_(false)
{}


// Basic Constructor
template<class ScalarType, class MV, class OP>
LSQRSolMgr<ScalarType,MV,OP>::
LSQRSolMgr (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
	    const Teuchos::RCP<Teuchos::ParameterList> &pl) : 
  problem_(problem),
  isSet_(false),
  loaDetected_(false)
{
  // The linear problem to solve is allowed to be null here.  The user
  // must then set a nonnull linear problem (by calling setProblem())
  // before calling solve().
  // 
  // Similarly, users are allowed to set a null parameter list here,
  // but they must first set a nonnull parameter list (by calling
  // setParameters()) before calling solve().
  if (! is_null (pl)) {
    setParameters (pl);  
  }
}


template<class ScalarType, class MV, class OP>
Teuchos::RCP<const Teuchos::ParameterList> 
LSQRSolMgr<ScalarType,MV,OP>::getValidParameters() const
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef Teuchos::ScalarTraits<MagnitudeType> STM;

  // Set all the valid parameters and their default values.
  if (is_null (validParams_)) {
    // We use Teuchos::as just in case MagnitudeType doesn't have a
    // constructor that takes an int.  Otherwise, we could just write
    // "MagnitudeType(10)".
    const MagnitudeType ten = Teuchos::as<MagnitudeType> (10);
    const MagnitudeType sqrtEps = STM::squareroot (STM::eps());

    const MagnitudeType lambda = STM::zero();
    RCP<std::ostream> outputStream = rcpFromRef (std::cout);
    const MagnitudeType relRhsErr = ten * sqrtEps;
    const MagnitudeType relMatErr = ten * sqrtEps;
    const MagnitudeType condMax = STM::one() / STM::eps();
    const int maxIters = 1000;
    const int termIterMax = 1;
    const std::string orthoType ("DGKS");
    const MagnitudeType orthoKappa = Teuchos::as<MagnitudeType> (-1.0);
    const int verbosity = Belos::Errors;
    const int outputStyle = Belos::General;
    const int outputFreq = -1;
    const std::string label ("Belos");

    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("Output Stream", outputStream,
      "is a reference-counted pointer to the output stream receiving\n"
      "all solver output.");
    pl->set("Lambda", lambda, "is the damping parameter.");
    pl->set("Rel RHS Err", relRhsErr,
	    "estimates the error in the data defining the right-\n"
	    "hand side.");
    pl->set("Rel Mat Err", relMatErr,
	    "estimates the error in the data defining the matrix.");
    pl->set("Condition Limit", condMax,
      "bounds the estimated condition number of Abar.");
    pl->set("Maximum Iterations", maxIters,
      "allows at most the maximum number of iterations.");
    pl->set("Term Iter Max", termIterMax,
      "consecutive iterations meeting thresholds are necessary for\n"
       "for convergence.");
    pl->set("Orthogonalization", orthoType,
      "uses orthogonalization of either DGKS, ICGS, IMGS, or TSQR.");
    {
      OrthoManagerFactory<ScalarType, MV, OP> factory;
      pl->set("Orthogonalization", orthoType,
	      "refers to the orthogonalization method to use.  Valid "
	      "options: " + factory.validNamesString());
      RCP<const ParameterList> orthoParams = 
	factory.getDefaultParameters (orthoType);
      pl->set ("Orthogonalization Parameters", *orthoParams, 
	       "Parameters specific to the type of orthogonalization used.");
    }
    pl->set("Orthogonalization Constant", orthoKappa,
      "is the threshold used by DGKS orthogonalization to determine\n"
      "whether or not to repeat classical Gram-Schmidt.  This parameter\n"
      "is ignored if \"Orthogonalization\" is not \"DGKS\".");
    pl->set("Verbosity", verbosity,
      "type(s) of solver information are outputted to the output\n"
      "stream.");
    pl->set("Output Style", outputStyle,
      "the style used for the solver information outputted to the\n"
      "output stream.");
    pl->set("Output Frequency", outputFreq,
      "is the frequency at which information is written to the\n"
      "output stream.");  
    pl->set("Timer Label", label,
      "is the string to use as a prefix for the timer labels.");
    //  pl->set("Restart Timers", restartTimers_);
    pl->set("Block Size", 1, "Block size parameter (currently, this must always be 1).");

    validParams_ = pl;
  }
  return validParams_;
}


template<class ScalarType, class MV, class OP>
void 
LSQRSolMgr<ScalarType,MV,OP>::
setParameters (const Teuchos::RCP<Teuchos::ParameterList> &params)
{
  using Teuchos::isParameterType;
  using Teuchos::getParameter;
  using Teuchos::null;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcpFromRef;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using Teuchos::Exceptions::InvalidParameter;
  using Teuchos::Exceptions::InvalidParameterName;
  using Teuchos::Exceptions::InvalidParameterType;

  TEUCHOS_TEST_FOR_EXCEPTION(params.is_null(), std::invalid_argument,
		     "Belos::LSQRSolMgr::setParameters: "
		     "the input ParameterList is null.");
  RCP<const ParameterList> defaultParams = getValidParameters ();
  params->validateParametersAndSetDefaults (*defaultParams);

  // At this point, params is a valid parameter list with defaults
  // filled in.  Now we can "commit" it to our instance's parameter
  // list.
  params_ = params;

  // Get the damping (a.k.a. regularization) parameter lambda.
  lambda_ = params->get<MagnitudeType> ("Lambda");

  // Get the maximum number of iterations.
  maxIters_ = params->get<int>("Maximum Iterations");

  // (Re)set the timer label.
  {
    const std::string newLabel = params->get<std::string>("Timer Label");

    // Update parameter in our list and solver timer
    if (newLabel != label_) {
      label_ = newLabel;
    }

#ifdef BELOS_TEUCHOS_TIME_MONITOR
    std::string newSolveLabel = label_ + ": LSQRSolMgr total solve time";
    if (timerSolve_.is_null()) {
      // Ask TimeMonitor for a new timer.
      timerSolve_ = TimeMonitor::getNewCounter (newSolveLabel);
    } else {
      // We've already created a timer, but we may have changed its
      // label.  If we did change its name, then we have to forget
      // about the old timer and create a new one with a different
      // name.  This is because Teuchos::Time doesn't give you a way
      // to change a timer's name, once you've created it.  We assume
      // that if the user changed the timer's label, then the user
      // wants to reset the timer's results.
      const std::string oldSolveLabel = timerSolve_->name ();

      if (oldSolveLabel != newSolveLabel) {
	// Tell TimeMonitor to forget about the old timer.
	// TimeMonitor lets you clear timers by name.
	TimeMonitor::clearCounter (oldSolveLabel);
	timerSolve_ = TimeMonitor::getNewCounter (newSolveLabel);
      }
    }
#endif // BELOS_TEUCHOS_TIME_MONITOR
  }

  // Check for a change in verbosity level
  {
    int newVerbosity = 0;
    // ParameterList gets confused sometimes about enums.  This
    // ensures that no matter how "Verbosity" was stored -- either an
    // an int, or as a Belos::MsgType enum, we will be able to extract
    // it.  If it was stored as some other type, we let the exception
    // through.
    try {
      newVerbosity = params->get<Belos::MsgType> ("Verbosity");
    } catch (Teuchos::Exceptions::InvalidParameterType&) {
      newVerbosity = params->get<int> ("Verbosity");
    }
    if (newVerbosity != verbosity_) {
      verbosity_ = newVerbosity;
    }
  }

  // (Re)set the output style.
  outputStyle_ = params->get<int> ("Output Style");

  // Get the output stream for the output manager.
  //
  // In case the output stream can't be read back in, we default to
  // stdout (std::cout), just to ensure reasonable behavior.
  {
    outputStream_ = params->get<RCP<std::ostream> > ("Output Stream");

    // We assume that a null output stream indicates that the user
    // doesn't want to print anything, so we replace it with a "black
    // hole" stream that prints nothing sent to it.  (We can't use a
    // null output stream, since the output manager always sends
    // things it wants to print to the output stream.)
    if (outputStream_.is_null())
      outputStream_ = rcp (new Teuchos::oblackholestream);
  }

  // Get the frequency of solver output.  (For example, -1 means
  // "never," 1 means "every iteration.")
  outputFreq_ = params->get<int> ("Output Frequency");

  // Create output manager if we need to, using the verbosity level
  // and output stream that we fetched above.  We do this here because
  // instantiating an OrthoManager using OrthoManagerFactory requires
  // a valid OutputManager.
  if (printer_.is_null()) {
    printer_ = rcp (new OutputManager<ScalarType> (verbosity_, outputStream_));
  } else {
    printer_->setVerbosity (verbosity_);
    printer_->setOStream (outputStream_);
  }

  // Check if the orthogonalization changed, or if we need to
  // initialize it.
  typedef OrthoManagerFactory<ScalarType, MV, OP> factory_type;
  factory_type factory;
  bool mustMakeOrtho = false;
  {
    std::string tempOrthoType;
    try {
      tempOrthoType = params_->get<std::string> ("Orthogonalization");
    } catch (InvalidParameterName&) {
      tempOrthoType = orthoType_;
    }
    if (ortho_.is_null() || tempOrthoType != orthoType_) {
      mustMakeOrtho = true;

      // Ensure that the specified orthogonalization type is valid.
      if (! factory.isValidName (tempOrthoType)) {
	std::ostringstream os;
	os << "Belos::LSQRSolMgr: Invalid orthogonalization name \"" 
	   << tempOrthoType << "\".  The following are valid options "
	   << "for the \"Orthogonalization\" name parameter: ";
	factory.printValidNames (os);
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
      }
      orthoType_ = tempOrthoType; // The name is valid, so accept it.
      params_->set ("Orthogonalization", orthoType_);
    }
  }

  // Get any parameters for the orthogonalization ("Orthogonalization
  // Parameters").  If not supplied, the orthogonalization manager
  // factory will supply default values.
  //
  // NOTE (mfh 21 Oct 2011) For the sake of backwards compatibility,
  // if params has an "Orthogonalization Constant" parameter and the
  // DGKS orthogonalization manager is to be used, the value of this
  // parameter will override DGKS's "depTol" parameter.
  //
  // Users must supply the orthogonalization manager parameters as a
  // sublist (supplying it as an RCP<ParameterList> would make the
  // resulting parameter list not serializable).
  RCP<ParameterList> orthoParams;
  { // The nonmember function returns an RCP<ParameterList>, 
    // which is what we want here.
    using Teuchos::sublist;
    // Abbreviation to avoid typos.
    const std::string paramName ("Orthogonalization Parameters");

    try {
      orthoParams = sublist (params_, paramName, true);
    } catch (InvalidParameter&) {
      // We didn't get the parameter list from params, so get a
      // default parameter list from the OrthoManagerFactory.
      // Modify params_ so that it has the default parameter list,
      // and set orthoParams to ensure it's a sublist of params_
      // (and not just a copy of one).
      params_->set (paramName, factory.getDefaultParameters (orthoType_));
      orthoParams = sublist (params_, paramName, true);
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(orthoParams.is_null(), std::logic_error, 
			     "Failed to get orthogonalization parameters.  "
			     "Please report this bug to the Belos developers.");

  // If we need to, instantiate a new MatOrthoManager subclass
  // instance corresponding to the desired orthogonalization method.
  // We've already fetched the orthogonalization method name
  // (orthoType_) and its parameters (orthoParams) above.
  //
  // NOTE (mfh 21 Oct 2011) We only instantiate a new MatOrthoManager
  // subclass if the orthogonalization method name is different than
  // before.  Thus, for some orthogonalization managers, changes to
  // their parameters may not get propagated, if the manager type
  // itself didn't change.  The one exception is the "depTol"
  // (a.k.a. orthoKappa or "Orthogonalization Constant") parameter of
  // DGKS; changes to that _do_ get propagated down to the DGKS
  // instance.
  //
  // The most general way to fix this issue would be to supply each
  // orthogonalization manager class with a setParameterList() method
  // that takes a parameter list input, and changes the parameters as
  // appropriate.  A less efficient but correct way would be simply to
  // reinstantiate the OrthoManager every time, whether or not the
  // orthogonalization method name or parameters have changed.
  if (mustMakeOrtho) {
    // Create orthogonalization manager.  This requires that the
    // OutputManager (printer_) already be initialized.  LSQR
    // currently only orthogonalizes with respect to the Euclidean
    // inner product, so we set the inner product matrix M to null.
    RCP<const OP> M = null;
    ortho_ = factory.makeMatOrthoManager (orthoType_, M, printer_, 
					  label_, orthoParams);
  }
  TEUCHOS_TEST_FOR_EXCEPTION(ortho_.is_null(), std::logic_error, 
			     "The MatOrthoManager is not yet initialized, but "
			     "should be by this point.  "
			     "Please report this bug to the Belos developers.");

  // Check which orthogonalization constant to use.  We only need this
  // if orthoType_ == "DGKS" (and we already fetched the orthoType_
  // parameter above).
  if (orthoType_ == "DGKS") {
    if (params->isParameter ("Orthogonalization Constant")) {
      orthoKappa_ = params_->get<MagnitudeType> ("Orthogonalization Constant");

      if (orthoKappa_ > 0 && ! ortho_.is_null()) {
	typedef DGKSOrthoManager<ScalarType,MV,OP> ortho_impl_type;
	rcp_dynamic_cast<ortho_impl_type> (ortho_)->setDepTol (orthoKappa_);
      }
    }
  }

  // Check for condition number limit, number of consecutive passed
  // iterations, relative RHS error, and relative matrix error.
  // Create the LSQR convergence test if necessary.
  {
    condMax_ = params->get<MagnitudeType> ("Condition Limit");
    termIterMax_ = params->get<int>("Term Iter Max");
    relRhsErr_ = params->get<MagnitudeType> ("Rel RHS Err");
    relMatErr_ = params->get<MagnitudeType> ("Rel Mat Err");

    // Create the LSQR convergence test if it doesn't exist yet. 
    // Otherwise, update its parameters.
    if (convTest_.is_null()) {
      convTest_ = 
	rcp (new LSQRStatusTest<ScalarType,MV,OP> (condMax_, termIterMax_, 
						   relRhsErr_, relMatErr_));
    } else {
      convTest_->setCondLim (condMax_);
      convTest_->setTermIterMax (termIterMax_);
      convTest_->setRelRhsErr (relRhsErr_);
      convTest_->setRelMatErr (relMatErr_);
    }
  }

  // Create the status test for maximum number of iterations if
  // necessary.  Otherwise, update it with the new maximum iteration
  // count.
  if (maxIterTest_.is_null()) {
    maxIterTest_ = rcp (new StatusTestMaxIters<ScalarType,MV,OP> (maxIters_));
  } else {
    maxIterTest_->setMaxIters (maxIters_);
  }

  // The stopping criterion is an OR combination of the test for
  // maximum number of iterations, and the LSQR convergence test.
  // ("OR combination" means that both tests will always be evaluated,
  // as opposed to a SEQ combination.)
  typedef StatusTestCombo<ScalarType,MV,OP>  StatusTestCombo_t;
  // If sTest_ is not null, then maxIterTest_ and convTest_ were
  // already constructed on entry to this routine, and sTest_ has
  // their pointers.  Thus, maxIterTest_ and convTest_ have gotten any
  // parameter changes, so we don't need to do anything to sTest_.
  if (sTest_.is_null()) {
    sTest_ = rcp (new StatusTestCombo_t (StatusTestCombo_t::OR, 
					 maxIterTest_, 
					 convTest_));
  }
  
  if (outputTest_.is_null()) {
    // Create the status test output class.
    // This class manages and formats the output from the status test.
    StatusTestOutputFactory<ScalarType,MV,OP> stoFactory (outputStyle_);
    outputTest_ = stoFactory.create (printer_, sTest_, outputFreq_, 
				     Passed + Failed + Undefined);
    // Set the solver string for the output test.
    const std::string solverDesc = " LSQR ";
    outputTest_->setSolverDesc (solverDesc);
  } else {
    // FIXME (mfh 18 Sep 2011) How do we change the output style of
    // outputTest_, without destroying and recreating it?
    outputTest_->setOutputManager (printer_);
    outputTest_->setChild (sTest_);
    outputTest_->setOutputFrequency (outputFreq_);
    // Since outputTest_ can only be created here, I'm assuming that
    // the fourth constructor argument ("printStates") was set
    // correctly on constrution; I don't need to reset it (and I can't
    // set it anyway, given StatusTestOutput's interface).
  }

  // Inform the solver manager that the current parameters were set.
  isSet_ = true;
}


template<class ScalarType, class MV, class OP>
Belos::ReturnType LSQRSolMgr<ScalarType,MV,OP>::solve() {
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Set the current parameters if they were not set before.  NOTE:
  // This may occur if the user generated the solver manager with the
  // default constructor, but did not set any parameters using
  // setParameters().
  if (!isSet_) {
    setParameters (Teuchos::parameterList (*getValidParameters()));
  }

  TEUCHOS_TEST_FOR_EXCEPTION(problem_.is_null(), LSQRSolMgrLinearProblemFailure,
		     "The linear problem to solve is null.");
  TEUCHOS_TEST_FOR_EXCEPTION(!problem_->isProblemSet(), LSQRSolMgrLinearProblemFailure,
                     "LSQRSolMgr::solve(): The linear problem is not ready, "
		     "as its setProblem() method has not been called.");
  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs (*(problem_->getRHS ())) != 1, 
		     LSQRSolMgrBlockSizeFailure, 
		     "LSQRSolMgr::solve(): The current implementation of LSQR "
		     "only knows how to solve problems with one right-hand "
		     "side, but the linear problem to solve has " 
		     << MVT::GetNumberVecs (*(problem_->getRHS ())) 
		     << " right-hand sides.");

  // We've validated the LinearProblem instance above.  If any of the
  // StatusTests needed to be initialized using information from the
  // LinearProblem, now would be the time to do so.  (This is true of
  // GMRES, where the residual convergence test(s) to instantiate
  // depend on knowing whether there is a left preconditioner.  This
  // is why GMRES has an "isSTSet_" Boolean member datum, which tells
  // you whether the status tests have been instantiated and are ready
  // for use.

  // test isFlexible might go here.

  // Next the right-hand sides to solve are identified.  Among other things,
  // this enables getCurrLHSVec() to get the current initial guess vector,
  // and getCurrRHSVec() to get the current right-hand side (in Iter).
  std::vector<int> currRHSIdx(1, 0); 
  problem_->setLSIndex(currRHSIdx);

  // Reset the status test.  
  outputTest_->reset();

  // Don't assume convergence unless we've verified that the
  // convergence test passed.
  bool isConverged = false;

  // FIXME: Currently we are setting the initial guess to zero, since
  // the solver doesn't yet know how to handle a nonzero initial
  // guess.  This could be fixed by rewriting the solver to work with
  // the residual and a delta.
  //
  // In a least squares problem with a nonzero initial guess, the
  // minimzation problem involves the distance (in a norm depending on
  // the preconditioner) between the solution and the the initial
  // guess.

  ////////////////////////////////////////////////////////////////////
  // Solve the linear problem using LSQR
  ////////////////////////////////////////////////////////////////////

  // Parameter list for the LSQR iteration.
  Teuchos::ParameterList plist;

  // Use the solver manager's "Lambda" parameter to set the
  // iteration's "Lambda" parameter.  We know that the solver
  // manager's parameter list (params_) does indeed contain the
  // "Lambda" parameter, because solve() always ensures that
  // setParameters() has been called.
  plist.set ("Lambda", lambda_);

  typedef LSQRIter<ScalarType,MV,OP> iter_type;
  RCP<iter_type> lsqr_iter = 
    rcp (new iter_type (problem_, printer_, outputTest_, plist));
#ifdef BELOS_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor slvtimer(*timerSolve_);
#endif

  // Reset the number of iterations.
  lsqr_iter->resetNumIters();
  // Reset the number of calls that the status test output knows about.
  outputTest_->resetNumCalls();
  // Set the new state and initialize the solver.
  LSQRIterationState<ScalarType,MV> newstate;
  lsqr_iter->initializeLSQR(newstate);
  // tell lsqr_iter to iterate
  try {
    lsqr_iter->iterate();    

    // First check for convergence.  If we didn't converge, then check
    // whether we reached the maximum number of iterations.  If
    // neither of those happened, there must have been a bug.
    if (convTest_->getStatus() == Belos::Passed) {
      isConverged = true;
    } else if (maxIterTest_->getStatus() == Belos::Passed) {
      isConverged = false;
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
			 "LSQRSolMgr::solve(): LSQRIteration::iterate() "
			 "returned without either the convergence test or "
			 "the maximum iteration count test passing."
			 "Please report this bug to the Belos developers.");
    }
  } catch (const std::exception &e) {
    printer_->stream(Belos::Errors) << "Error! Caught std::exception in LSQRIter::iterate() at iteration " 
				    << lsqr_iter->getNumIters() << std::endl 
				    << e.what() << std::endl;
    throw;
  }
      
  // identify current linear system as solved LinearProblem 
  problem_->setCurrLS();
  // print final summary
  sTest_->print (printer_->stream (Belos::FinalSummary));

  // Print timing information, if the corresponding compile-time and
  // run-time options are enabled.
#ifdef BELOS_TEUCHOS_TIME_MONITOR
  // Calling summarize() can be expensive, so don't call unless the
  // user wants to print out timing details.  summarize() will do all
  // the work even if it's passed a "black hole" output stream.
  if (verbosity_ & TimingDetails)
    Teuchos::TimeMonitor::summarize (printer_->stream (Belos::TimingDetails));
#endif // BELOS_TEUCHOS_TIME_MONITOR

  // A posteriori solve information
  numIters_ = maxIterTest_->getNumIters();
  matCondNum_ = convTest_->getMatCondNum();
  matNorm_ = convTest_->getMatNorm();
  resNorm_ = convTest_->getResidNorm();
  matResNorm_ = convTest_->getLSResidNorm();
 
  if (! isConverged) {
    return Belos::Unconverged;
  } else {
    return Belos::Converged;
  }
}

// LSQRSolMgr requires the solver manager to return an eponymous std::string.
template<class ScalarType, class MV, class OP>
std::string LSQRSolMgr<ScalarType,MV,OP>::description() const
{
  std::ostringstream oss;
  oss << "LSQRSolMgr<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
  oss << "{";
  oss << "Ortho Type='"<<orthoType_<<"'";
  oss << ", Lambda="<< lambda_;
  oss << ", condition number limit="<< condMax_;
  oss << ", relative RHS Error="<< relRhsErr_;
  oss << ", relative Matrix Error="<< relMatErr_;
  oss << ", maximum number of iterations="<< maxIters_;
  oss << ", termIterMax="<<termIterMax_;
  oss << "}";
  return oss.str();
}
  
} // end Belos namespace

#endif /* BELOS_LSQR_SOLMGR_HPP */
