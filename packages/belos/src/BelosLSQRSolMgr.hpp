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
 * \warning DO NOT USE; DEPRECATED.
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
/// \section Belos_LSQR_alg Algorithm
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
/// \section Belos_LSQR_real ScalarType must be real
///
/// This LSQR implementation currently only supports real-valued (not
/// complex-valued) ScalarType types.  You may check whether ScalarType is
/// complex using the following code:
/// \code
/// if (Teuchos::ScalarTraits<ScalarType>::isComplex) {
///  // ScalarType is complex valued.
/// } else {
///  // ScalarType is real valued.
/// }
/// \endcode
///
/// This is not a limitation of the LSQR method itself, just of the
/// current implementation.  If there is sufficient interest, we can
/// remedy this deficiency.  For now, if you attempt to invoke the
/// constructor when <tt>ScalarType</tt> is complex, the constructor
/// will throw an exception.  This is why this class inherits from
/// Details::RealSolverManager.  LSQRSolMgr can still compile if
/// <tt>ScalarType</tt> is complex, but you will not be able to
/// construct a LSQRSolMgr instance in that case, due to the
/// aforementioned run-time error that the constructor raises.  We do
/// this so that the class will still compile, whether ScalarType is
/// real or complex.  This helps make SolverFactory valid to compile,
/// whether ScalarType is real or complex.
///
/// \section Belos_LSQR_prec Preconditioning
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
/// \section Belos_LSQR_refs References
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


// Partial specialization for complex ScalarType.
// This contains a trivial implementation.
// See discussion in the class documentation above.
template<class ScalarType, class MV, class OP,
         const bool scalarTypeIsComplex = Teuchos::ScalarTraits<ScalarType>::isComplex>
class LSQRSolMgr :
    public Details::RealSolverManager<ScalarType, MV, OP,
                                      Teuchos::ScalarTraits<ScalarType>::isComplex>
{
  static const bool isComplex = Teuchos::ScalarTraits<ScalarType>::isComplex;
  typedef Details::RealSolverManager<ScalarType, MV, OP, isComplex> base_type;

public:
  LSQRSolMgr () :
    base_type ()
  {}
  LSQRSolMgr (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
              const Teuchos::RCP<Teuchos::ParameterList> &pl) :
    base_type ()
  {}
  virtual ~LSQRSolMgr () {}

  //! clone for Inverted Injection (DII)
  Teuchos::RCP<SolverManager<ScalarType, MV, OP> > clone () const override {
    return Teuchos::rcp(new LSQRSolMgr<ScalarType,MV,OP,scalarTypeIsComplex>);
  }
};


// Partial specialization for real ScalarType.
// This contains the actual working implementation of LSQR.
// See discussion in the class documentation above.
template<class ScalarType, class MV, class OP>
class LSQRSolMgr<ScalarType, MV, OP, false> :
    public Details::RealSolverManager<ScalarType, MV, OP, false> {
private:
  typedef MultiVecTraits<ScalarType,MV> MVT;
  typedef OperatorTraits<ScalarType,MV,OP> OPT;
  typedef Teuchos::ScalarTraits<ScalarType> STS;
  typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
  typedef Teuchos::ScalarTraits<MagnitudeType> STM;

public:

  //! @name Construct/Destroy
  //@{

  /*! \brief Empty constructor for LSQRSolMgr.
   * This constructor takes no arguments and sets the default values for the solver.
   * The linear problem must be passed in using setProblem() before solve() is called
   * on this object.
   * The solver values can be changed using setParameters().
   */
  LSQRSolMgr ();

  /*! \brief Basic constructor for LSQRSolMgr.
   *
   * This constructor accepts the LinearProblem to be solved in addition
   * to a parameter list of options for the solver manager.  Blocks of size > 1 are not
   * implemented.  The options are otherwise the BlockGmres options.
   *   - "Maximum Iterations" - the maximum number of iterations the LSQR solver
   *                            is allowed to perform. Default: 1000
   *   - "Condition Limit" - a \c MagnitudeType specifying the upper limit of the estimate of
   *                         the norm of Abar to decide convergence. Default: 0.
   *   - "Term Iter Max": The number of consecutive successful
   *     iterations required before LSQR considers the problem
   *     converged.  Default: 1.
   *   - "Rel RHS Err" (or "Convergence Tolerance"): an estimate of
   *     the error in the data defining the right-hand side.  Default:
   *     10*sqrt(eps).
   *   - "Rel Mat Err" - an estimate of the error in the data defining the matrix.
   *                     Default: 10*sqrt(eps).
   *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Belos::Errors
   *   - "Output Style" - a OutputType specifying the style of output. Default: Belos::General
   *   - "Lambda"  - a \c MagnitudeType that specifies the regularization parameter.
   *
   * This LSQR implementation only supports block size 1.  Like CG,
   * LSQR is a short recurrence method that, in finite precision
   * arithmetic and without reorthogonalization, does not have the "n"
   * step convergence property.  Without either blocks or
   * reorthogonalization, there is nothing to "Orthogonalize."
   */
  LSQRSolMgr (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >& problem,
              const Teuchos::RCP<Teuchos::ParameterList>& pl);

  //! Destructor (declared virtual for memory safety of base classes).
  virtual ~LSQRSolMgr () {}

  //! clone for Inverted Injection (DII)
  Teuchos::RCP<SolverManager<ScalarType, MV, OP> > clone () const override {
    return Teuchos::rcp(new LSQRSolMgr<ScalarType,MV,OP>);
  }
  //@}
  //! \name Accessor methods
  //@{

  /*! \brief Get current linear problem being solved for in this object.
   */
  const LinearProblem<ScalarType,MV,OP>& getProblem () const override {
    return *problem_;
  }

  /*! \brief Get a parameter list containing the valid parameters for this object.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

  /*! \brief Get a parameter list containing the current parameters for this object.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters() const override {
    return params_;
  }

  /*! \brief Return the timers for this object.
   *
   * The timers are ordered as follows:
   *   - time spent in solve() routine
   */
  Teuchos::Array<Teuchos::RCP<Teuchos::Time> > getTimers () const {
    return Teuchos::tuple (timerSolve_);
  }

  //! Iteration count from the last solve.
  int getNumIters () const override {
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
  bool isLOADetected () const override { return false; }

  //@}

  //! @name Set methods
  //@{

  //! Set the linear problem that needs to be solved.
  void setProblem (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >& problem) override {
    problem_ = problem;
  }

  //! Set the parameters the solver manager should use to solve the linear problem.
  void setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params) override;

  //@}

  //! @name Reset methods
  //@{
  /*! \brief reset the solver manager as specified by the \c ResetType, informs the
   *  solver manager that the solver should prepare for the next call to solve
   *  by resetting certain elements of the iterative solver strategy.
   */
  void reset (const ResetType type) override {
    if ((type & Belos::Problem) && ! problem_.is_null ()) {
      problem_->setProblem ();
    }
  }

  //@}
  //! \name Solver application methods
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
  ReturnType solve() override;

  //@}
  //! \name Overridden from Teuchos::Describable
  //@{

  //! One-line description of this solver.
  std::string description () const override;

  //@}

private:

  //! The linear problem to solve.
  Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > problem_;
  //! The output manager.
  Teuchos::RCP<OutputManager<ScalarType> > printer_;
  //! Output stream to which to write status output.
  Teuchos::RCP<std::ostream> outputStream_;

  //! The "master" status test (that includes all status tests).
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > sTest_;
  Teuchos::RCP<StatusTestMaxIters<ScalarType,MV,OP> > maxIterTest_;
  Teuchos::RCP<LSQRStatusTest<ScalarType,MV,OP> > convTest_;
  Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputTest_;

  //! Current parameter list.
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

template<class ScalarType, class MV, class OP>
LSQRSolMgr<ScalarType,MV,OP,false>::LSQRSolMgr () :
  lambda_ (STM::zero ()),
  relRhsErr_ (Teuchos::as<MagnitudeType> (10) * STM::squareroot (STM::eps ())),
  relMatErr_ (Teuchos::as<MagnitudeType> (10) * STM::squareroot (STM::eps ())),
  condMax_ (STM::one () / STM::eps ()),
  maxIters_ (1000),
  termIterMax_ (1),
  verbosity_ (Belos::Errors),
  outputStyle_ (Belos::General),
  outputFreq_ (-1),
  numIters_ (0),
  matCondNum_ (STM::zero ()),
  matNorm_ (STM::zero ()),
  resNorm_ (STM::zero ()),
  matResNorm_ (STM::zero ()),
  isSet_ (false),
  loaDetected_ (false)
{}

template<class ScalarType, class MV, class OP>
LSQRSolMgr<ScalarType,MV,OP,false>::
LSQRSolMgr (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >& problem,
            const Teuchos::RCP<Teuchos::ParameterList>& pl) :
  problem_ (problem),
  lambda_ (STM::zero ()),
  relRhsErr_ (Teuchos::as<MagnitudeType> (10) * STM::squareroot (STM::eps ())),
  relMatErr_ (Teuchos::as<MagnitudeType> (10) * STM::squareroot (STM::eps ())),
  condMax_ (STM::one () / STM::eps ()),
  maxIters_ (1000),
  termIterMax_ (1),
  verbosity_ (Belos::Errors),
  outputStyle_ (Belos::General),
  outputFreq_ (-1),
  numIters_ (0),
  matCondNum_ (STM::zero ()),
  matNorm_ (STM::zero ()),
  resNorm_ (STM::zero ()),
  matResNorm_ (STM::zero ()),
  isSet_ (false),
  loaDetected_ (false)
{
  // The linear problem to solve is allowed to be null here.  The user
  // must then set a nonnull linear problem (by calling setProblem())
  // before calling solve().
  //
  // Similarly, users are allowed to set a null parameter list here,
  // but they must first set a nonnull parameter list (by calling
  // setParameters()) before calling solve().
  if (! pl.is_null ()) {
    setParameters (pl);
  }
}


template<class ScalarType, class MV, class OP>
Teuchos::RCP<const Teuchos::ParameterList>
LSQRSolMgr<ScalarType,MV,OP,false>::getValidParameters() const
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;

  // Set all the valid parameters and their default values.
  if (validParams_.is_null ()) {
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
    const int verbosity = Belos::Errors;
    const int outputStyle = Belos::General;
    const int outputFreq = -1;
    const std::string label ("Belos");

    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set ("Output Stream", outputStream, "Teuchos::RCP<std::ostream> "
             "(reference-counted pointer to the output stream) receiving "
             "all solver output");
    pl->set ("Lambda", lambda, "Damping parameter");
    pl->set ("Rel RHS Err", relRhsErr, "Estimates the error in the data "
             "defining the right-hand side");
    pl->set ("Rel Mat Err", relMatErr, "Estimates the error in the data "
             "defining the matrix.");
    pl->set ("Condition Limit", condMax, "Bounds the estimated condition "
             "number of Abar.");
    pl->set ("Maximum Iterations", maxIters, "Maximum number of iterations");
    pl->set ("Term Iter Max", termIterMax, "The number of consecutive "
             "iterations must that satisfy all convergence criteria in order "
             "for LSQR to stop iterating");
    pl->set ("Verbosity", verbosity, "Type(s) of solver information written to "
             "the output stream");
    pl->set ("Output Style", outputStyle, "Style of solver output");
    pl->set ("Output Frequency", outputFreq, "Frequency at which information "
             "is written to the output stream (-1 means \"not at all\")");
    pl->set ("Timer Label", label, "String to use as a prefix for the timer "
             "labels");
    pl->set ("Block Size", 1, "Block size parameter (currently, LSQR requires "
             "this must always be 1)");
    validParams_ = pl;
  }
  return validParams_;
}


template<class ScalarType, class MV, class OP>
void
LSQRSolMgr<ScalarType,MV,OP,false>::
setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params)
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

  TEUCHOS_TEST_FOR_EXCEPTION
    (params.is_null (), std::invalid_argument,
     "Belos::LSQRSolMgr::setParameters: The input ParameterList is null.");
  RCP<const ParameterList> defaultParams = getValidParameters ();

  // FIXME (mfh 29 Apr 2015) Our users would like to supply one
  // ParameterList that works for both GMRES and LSQR.  Thus, we want
  // LSQR (the less-used solver) to ignore parameters it doesn't
  // recognize).  For now, therefore, it should not validate, since
  // validation cannot distinguish between misspellings and
  // unrecognized parameters.  (Perhaps Belos should have a central
  // facility for all parameters recognized by some solver in Belos,
  // so we could use that for spell checking.)
  //
  //params->validateParameters (*defaultParams);

  // mfh 29 Apr 2015: The convention in Belos is that the input
  // ParameterList is a "delta" from the current state.  Thus, we
  // don't fill in the input ParameterList with defaults, and we only
  // change the current state if the corresponding parameter was
  // explicitly set in the input ParameterList.  We set up the solver
  // with the default state on construction.

  // Get the damping (regularization) parameter lambda.
  if (params->isParameter ("Lambda")) {
    lambda_ = params->get<MagnitudeType> ("Lambda");
  } else if (params->isParameter ("lambda")) {
    lambda_ = params->get<MagnitudeType> ("lambda");
  }

  // Get the maximum number of iterations.
  if (params->isParameter ("Maximum Iterations")) {
    maxIters_ = params->get<int> ("Maximum Iterations");
  }
  TEUCHOS_TEST_FOR_EXCEPTION
    (maxIters_ < 0, std::invalid_argument, "Belos::LSQRSolMgr::setParameters: "
     "\"Maximum Iterations\" = " << maxIters_ << " < 0.");

  // (Re)set the timer label.
  {
    const std::string newLabel =
      params->isParameter ("Timer Label") ?
      params->get<std::string> ("Timer Label") :
      label_;

    // Update parameter in our list and solver timer
    if (newLabel != label_) {
      label_ = newLabel;
    }

#ifdef BELOS_TEUCHOS_TIME_MONITOR
    const std::string newSolveLabel = (newLabel != "") ?
      (newLabel + ": Belos::LSQRSolMgr total solve time") :
      std::string ("Belos::LSQRSolMgr total solve time");
    if (timerSolve_.is_null ()) {
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
  if (params->isParameter ("Verbosity")) {
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
  if (params->isParameter ("Output Style")) {
    outputStyle_ = params->get<int> ("Output Style");
  }

  // Get the output stream for the output manager.
  //
  // In case the output stream can't be read back in, we default to
  // stdout (std::cout), just to ensure reasonable behavior.
  if (params->isParameter ("Output Stream")) {
    outputStream_ = params->get<RCP<std::ostream> > ("Output Stream");
  }
  // We assume that a null output stream indicates that the user
  // doesn't want to print anything, so we replace it with a "black
  // hole" stream that prints nothing sent to it.  (We can't use a
  // null output stream, since the output manager always sends
  // things it wants to print to the output stream.)
  if (outputStream_.is_null ()) {
    outputStream_ = rcp (new Teuchos::oblackholestream ());
  }

  // Get the frequency of solver output.  (For example, -1 means
  // "never," and 1 means "every iteration.")
  if (params->isParameter ("Output Frequency")) {
    outputFreq_ = params->get<int> ("Output Frequency");
  }

  // Create output manager if we need to, using the verbosity level
  // and output stream that we fetched above.  Status tests (i.e.,
  // stopping criteria) need this.
  if (printer_.is_null ()) {
    printer_ = rcp (new OutputManager<ScalarType> (verbosity_, outputStream_));
  } else {
    printer_->setVerbosity (verbosity_);
    printer_->setOStream (outputStream_);
  }

  // Check for condition number limit, number of consecutive passed
  // iterations, relative RHS error, and relative matrix error.
  // Create the LSQR convergence test if necessary.
  {
    if (params->isParameter ("Condition Limit")) {
      condMax_ = params->get<MagnitudeType> ("Condition Limit");
    }
    if (params->isParameter ("Term Iter Max")) {
      termIterMax_ = params->get<int> ("Term Iter Max");
    }
    if (params->isParameter ("Rel RHS Err")) {
      relRhsErr_ = params->get<MagnitudeType> ("Rel RHS Err");
    }
    else if (params->isParameter ("Convergence Tolerance")) {
      // NOTE (mfh 29 Apr 2015) We accept this parameter as an alias
      // for "Rel RHS Err".
      relRhsErr_ = params->get<MagnitudeType> ("Convergence Tolerance");
    }

    if (params->isParameter ("Rel Mat Err")) {
      relMatErr_ = params->get<MagnitudeType> ("Rel Mat Err");
    }

    // Create the LSQR convergence test if it doesn't exist yet.
    // Otherwise, update its parameters.
    if (convTest_.is_null ()) {
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
  typedef StatusTestCombo<ScalarType,MV,OP> combo_type;
  // If sTest_ is not null, then maxIterTest_ and convTest_ were
  // already constructed on entry to this routine, and sTest_ has
  // their pointers.  Thus, maxIterTest_ and convTest_ have gotten any
  // parameter changes, so we don't need to do anything to sTest_.
  if (sTest_.is_null()) {
    sTest_ = rcp (new combo_type (combo_type::OR, maxIterTest_, convTest_));
  }

  if (outputTest_.is_null ()) {
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

  // At this point, params is a valid ParameterList.  Now we can
  // "commit" it to our instance's ParameterList.
  params_ = params;

  // Inform the solver manager that the current parameters were set.
  isSet_ = true;
}


template<class ScalarType, class MV, class OP>
Belos::ReturnType
LSQRSolMgr<ScalarType,MV,OP,false>::solve ()
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Set the current parameters if they were not set before.  NOTE:
  // This may occur if the user generated the solver manager with the
  // default constructor, but did not set any parameters using
  // setParameters().
  if (! isSet_) {
    this->setParameters (Teuchos::parameterList (* (getValidParameters ())));
  }

  TEUCHOS_TEST_FOR_EXCEPTION
    (problem_.is_null (), LSQRSolMgrLinearProblemFailure,
     "Belos::LSQRSolMgr::solve: The linear problem to solve is null.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (! problem_->isProblemSet (), LSQRSolMgrLinearProblemFailure,
     "Belos::LSQRSolMgr::solve: The linear problem is not ready, "
     "as its setProblem() method has not been called.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (MVT::GetNumberVecs (*(problem_->getRHS ())) != 1,
     LSQRSolMgrBlockSizeFailure, "Belos::LSQRSolMgr::solve: "
     "The current implementation of LSQR only knows how to solve problems "
     "with one right-hand side, but the linear problem to solve has "
     << MVT::GetNumberVecs (* (problem_->getRHS ()))
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
  std::vector<int> currRHSIdx (1, 0);
  problem_->setLSIndex (currRHSIdx);

  // Reset the status test.
  outputTest_->reset ();

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
  Teuchos::TimeMonitor slvtimer (*timerSolve_);
#endif

  // Reset the number of iterations.
  lsqr_iter->resetNumIters ();
  // Reset the number of calls that the status test output knows about.
  outputTest_->resetNumCalls ();
  // Set the new state and initialize the solver.
  LSQRIterationState<ScalarType, MV> newstate;
  lsqr_iter->initializeLSQR (newstate);
  // tell lsqr_iter to iterate
  try {
    lsqr_iter->iterate ();

    // First check for convergence.  If we didn't converge, then check
    // whether we reached the maximum number of iterations.  If
    // neither of those happened, there must have been a bug.
    if (convTest_->getStatus () == Belos::Passed) {
      isConverged = true;
    } else if (maxIterTest_->getStatus () == Belos::Passed) {
      isConverged = false;
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "Belos::LSQRSolMgr::solve: "
         "LSQRIteration::iterate returned without either the convergence test "
         "or the maximum iteration count test passing.  "
         "Please report this bug to the Belos developers.");
    }
  } catch (const std::exception& e) {
    printer_->stream(Belos::Errors)
      << "Error! Caught std::exception in LSQRIter::iterate at iteration "
      << lsqr_iter->getNumIters () << std::endl << e.what () << std::endl;
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
std::string LSQRSolMgr<ScalarType,MV,OP,false>::description () const
{
  std::ostringstream oss;
  oss << "LSQRSolMgr<...," << STS::name () << ">";
  oss << "{";
  oss << "Lambda: " << lambda_;
  oss << ", condition number limit: " << condMax_;
  oss << ", relative RHS Error: " << relRhsErr_;
  oss << ", relative Matrix Error: " << relMatErr_;
  oss << ", maximum number of iterations: " << maxIters_;
  oss << ", termIterMax: " << termIterMax_;
  oss << "}";
  return oss.str ();
}

} // end Belos namespace

#endif /* BELOS_LSQR_SOLMGR_HPP */
