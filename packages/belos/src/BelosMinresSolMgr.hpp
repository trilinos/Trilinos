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

#ifndef BELOS_MINRES_SOLMGR_HPP
#define BELOS_MINRES_SOLMGR_HPP

/// \file BelosMinresSolMgr.hpp
/// \brief Solver manager for the MINRES linear solver

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "BelosMinresIter.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestGenResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosStatusTestOutputFactory.hpp"
#include "BelosOutputManager.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#ifdef BELOS_TEUCHOS_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

#include "Teuchos_StandardParameterEntryValidators.hpp"
// Teuchos::ScalarTraits<int> doesn't define rmax(), alas, so we get
// INT_MAX from here.
#include <climits>

namespace Belos {

  //! @name MinresSolMgr Exceptions
  //@{

  /// \class MinresSolMgrLinearProblemFailure
  /// \brief Thrown when solve() called and problem not set up
  ///
  /// MinresSolMgrLinearProblemFailure is thrown when the linear
  /// problem has not been set up (e.g., setProblem() was not called;
  /// or the constructor was not provided with a linear problem to
  /// solve), but solve() was called anyway; or when the linear
  /// problem cannot be solved by MINRES (e.g., if there is more than
  /// one right-hand side).
  //
  /// This subclass of std::exception may be thrown from the
  /// MinresSolMgr::solve() method.
  ///
  class MinresSolMgrLinearProblemFailure : public BelosError {
  public:
    MinresSolMgrLinearProblemFailure (const std::string& what_arg) :
      BelosError(what_arg)
    {}
  };

  ///
  /// \class Belos::MinresSolMgr
  /// \brief MINRES linear solver solution manager
  /// \author Nico Schl\"omer
  /// \ingroup belos_solver_framework
  ///
  /// The Minimal Residual Method (MINRES) is a Krylov subspace method
  /// for solving symmetric (in real arithmetic, or Hermitian in complex
  /// arithmetic), nonsingular, but possibly indefinite linear systems
  /// \fn$Ax=b\fn$.  It works on one right-hand side \fn$b\fn$ at a
  /// time.
  ///
  /// References:
  ///
  /// C. Paige and M. Saunders.  "Solution of sparse indefinite systems
  /// of linear equations."  SIAM J. Numer. Anal., vol. 12, pp. 617-629,
  /// 1975.
  ///
  template<class ScalarType, class MV, class OP>
  class MinresSolMgr : public SolverManager<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits< MagnitudeType > MT;

  public:

    /// \brief List of valid MINRES parameters and their default values.
    ///
    /// One way to set up this solver manager with nondefault
    /// parameters, is to make a deep non-const copy of the default
    /// parameters, and change the parameters you want to change.
    ///
    /// \note This is a class ("static") method, so that it can be
    ///   called before you have constructed a MinresSolMgr object
    ///   (useful for you, so you can have valid parameters before
    ///   constructing one), or within the MinresSolMgr constructor
    ///   (useful for me, so I can set default parameters).
    static Teuchos::RCP<const Teuchos::ParameterList> defaultParameters();

    //! @name Constructors/Destructor
    //@{

    /// \brief Default constructor
    ///
    /// This constructor takes no arguments and sets the default
    /// values for the solver.  The linear problem must be passed in
    /// using setProblem() before solve() is called on this object,
    /// otherwise an exception is thrown.  The solver's parameters
    /// (which this constructor sets to their default values) may be
    /// changed using setParameters().
    MinresSolMgr();

    /// \brief Basic constructor for MinresSolMgr.
    ///
    /// \param problem [in/out] The LinearProblem to be solved
    /// \param params [in/out] Parameter list of options for the solver manager.
    ///   Parameters not provided will be filled in with default values.
    ///   These are the options accepted by the solver manager:
    ///   - "Block Size" - an \c int specifying the block size to be used by the
    ///     underlying MINRES solver. Default: 1 (which is the only valid value!)
    ///   - "Convergence Tolerance" - a \c MagnitudeType specifying the level that
    ///     residual norms must reach to decide convergence. Default value: 1e-8.
    ///   - "Maximum Iterations" - an \c int specifying the maximum number of
    ///     iterations the underlying solver is allowed to perform. Default: 1000
    ///   - "Verbosity" - a sum of MsgType (stored as an int) specifying the
    ///     verbosity.  Default value: Belos::Errors
    ///   - "Output Style" - a OutputType specifying the style of output.
    ///     Default value: Belos::General
    ///   - "Output Stream" - a reference-counted pointer to the output stream
    ///     where all solver output is sent.  Default value:
    ///     Teuchos::rcp(&std::cout,false)
    ///   - "Output Frequency" - an \c int specifying how often (in
    ///     terms of number of iterations) intermediate convergence
    ///     information should be written to the output stream.
    ///     Default value: -1 (which means no intermediate convergence
    ///     information is ever written to the output stream)
    ///   - "Timer Label" - an \c std::string to use as a prefix for the timer
    ///     labels.  Default value: "Belos"
    MinresSolMgr (const Teuchos::RCP<LinearProblem< ScalarType, MV, OP> > &problem,
                  const Teuchos::RCP<Teuchos::ParameterList> &params);

    //! Destructor.
    virtual ~MinresSolMgr() {};
    //@}

    //! @name Accessor methods
    //@{

    //! Return the linear problem to be solved.
    const LinearProblem<ScalarType,MV,OP>& getProblem() const {
      return *problem_;
    }

    //! Return the list of default parameters for this object.
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const {
      if (defaultParams_.is_null()) {
	defaultParams_ = defaultParameters ();
      }
      return defaultParams_;
    }

    //! Return the list of current parameters for this object.
    Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters() const {
      return params_;
    }

    /// \brief Return all timers for this object.
    ///
    /// Currently only one timer is being used, which is the total
    /// time spent in the \c solve() routine.  Thus, the returned
    /// Array currently has only one element.
    ///
    /// \warning If \c setParameters() has not yet been called, or if
    ///   you change the timer label, that invalidates the pointer
    ///   timer(s).
    Teuchos::Array<Teuchos::RCP<Teuchos::Time> > getTimers() const {
      return Teuchos::tuple (timerSolve_);
    }

    /// \brief Tolerance achieved by the last \c solve() invocation.
    ///
    /// This is the maximum over all right-hand sides' achieved
    /// convergence tolerances, and is set whether or not the solve
    /// actually managed to achieve the desired convergence tolerance.
    MagnitudeType achievedTol() const {
      return achievedTol_;
    }

    //! Get the iteration count for the most recent call to \c solve().
    int getNumIters() const {
      return numIters_;
    }

    /// Whether a loss of accuracy was detected in the solver.
    ///
    /// \warning This implementation of MINRES does not currently
    ///   attempt to detect a loss of accuracy in the solver; thus we
    ///   always return false (for now).
    bool isLOADetected() const { return false; }

    //@}

    //! @name Set methods (overridden from \c SolverManager)
    //@{

    void
    setProblem (const Teuchos::RCP<LinearProblem<ScalarType, MV, OP> > &problem)
    {
      problem_ = problem;
    }

    void
    setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params);

    //@}

    //! @name Reset methods (overridden from \c SolverManager)
    //@{

    void
    reset (const ResetType type)
    {
      if ((type & Belos::Problem) && ! problem_.is_null()) {
	problem_->setProblem ();
      }
    }
    //@}

    //! @name Solver application methods (overridden from \c SolverManager)
    //@{

    /// \brief Iterate until the status test tells us to stop.
    ///
    /// This method implements SolverManager::solve() (which see).
    ///
    /// MINRES' implementation of this method invokes \c
    /// MinresIter::iterate(), which will return either because a
    /// specially constructed status test that evaluates to ::Passed,
    /// or an std::exception is thrown.
    ///
    /// A return from MinresIter::iterate() signifies one of the
    /// following scenarios:
    /// - the maximum number of iterations has been exceeded. In this
    ///   scenario, the current solutions to the linear system will be
    ///   placed in the linear problem and return ::Unconverged.
    /// - global convergence has been met. In this case, the current
    ///   solutions to the linear system will be placed in the linear
    ///   problem and the solver manager will return ::Converged.
    ReturnType solve();

    //@}

    /** \name Overridden from Teuchos::Describable */
    //@{

    std::string description() const;

    //@}

  private:
    //! Linear problem to solve
    Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > problem_;

    //! Output manager.
    Teuchos::RCP<OutputManager<ScalarType> > printer_;
    Teuchos::RCP<std::ostream> outputStream_;

    /// \brief The full status test.
    ///
    /// sTest_ is an OR combination of maxIterTest_ and convTest_.  If
    /// you reallocate either of these, you have to give them to
    /// sTest_ again.  If you reallocate sTest_, you have to tell
    /// outputTest_.
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > sTest_;

    /// \brief The status test for maximum iteration count.
    ///
    /// If you reallocate this, sTest_ needs the new RCP.
    Teuchos::RCP<StatusTestMaxIters<ScalarType,MV,OP> > maxIterTest_;

    /// \brief The combined status test for convergence.
    ///
    /// If you reallocate this, sTest_ needs the new RCP.
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > convTest_;

    /// \brief The implicit (a.k.a. "recursive") residual norm test.
    ///
    /// If you reallocate this, convTest_ needs the new RCP.
    Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP> > impConvTest_;

    /// \brief The explicit residual norm test.
    ///
    /// If you reallocate this, convTest_ needs the new RCP.
    Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP> > expConvTest_;

    /// \brief The "status test" that handles output.
    ///
    /// This object keeps a pointer to printer_ and sTest_.  If you
    /// reallocate either of them, outputTest_ needs to know.
    Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputTest_;

    /// \brief List of default parameters.
    ///
    /// This is declared "mutable" because it is computed on demand.
    mutable Teuchos::RCP<const Teuchos::ParameterList> defaultParams_;

    //! List of current parameters
    Teuchos::RCP<Teuchos::ParameterList> params_;

    //! Current relative residual 2-norm convergence tolerance
    MagnitudeType convtol_;

    //! Tolerance achieved by the last \c solve() invocation.
    MagnitudeType achievedTol_;

    //! Maximum number of iterations before stopping
    int maxIters_;

    //! Current number of iterations
    int numIters_;

    //! Current block size (i.e., number of right-hand sides): always 1 (one).
    int blockSize_;

    //! Current output verbosity
    int verbosity_;

    //! Current output style
    int outputStyle_;

    //! Current frequency of output
    int outputFreq_;

    //! Timer label
    std::string label_;

    //! Total time to solution
    Teuchos::RCP<Teuchos::Time> timerSolve_;

    //! Whether the solver manager's parameters have been set
    bool parametersSet_;

    /// \brief Validate the given linear problem.
    ///
    /// We do this by raising std::invalid_argument (with an
    /// informative message) if the problem is null or its essential
    /// components are null.
    static void
    validateProblem (const Teuchos::RCP<LinearProblem<ScalarType, MV, OP> >& problem);
  };


  template<class ScalarType, class MV, class OP>
  Teuchos::RCP<const Teuchos::ParameterList>
  MinresSolMgr<ScalarType, MV, OP>::defaultParameters()
  {
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    using Teuchos::EnhancedNumberValidator;
    typedef MagnitudeType MT;
    typedef Teuchos::ScalarTraits<MT> MST;
    typedef Teuchos::ScalarTraits<int> IST;

    // List of parameters accepted by MINRES, and their default values.
    RCP<ParameterList> pl = parameterList ("MINRES");

    pl->set ("Convergence Tolerance", MST::squareroot (MST::eps()),
	     "Relative residual tolerance that needs to be achieved by "
	     "the iterative solver, in order for the linear system to be "
	     "declared converged.",
	     rcp (new EnhancedNumberValidator<MT> (MST::zero(), MST::rmax())));
    pl->set ("Maximum Iterations", static_cast<int>(1000),
	     "Maximum number of iterations allowed for each right-hand "
	     "side solved.",
	     rcp (new EnhancedNumberValidator<int> (0, INT_MAX)));
    pl->set ("Block Size", static_cast<int>(1),
	     "Number of vectors in each block.  WARNING: The current "
	     "implementation of MINRES only accepts a block size of 1, "
	     "since it can only solve for 1 right-hand side at a time.",
	     rcp (new EnhancedNumberValidator<int> (1, 1)));
    pl->set ("Verbosity", (int) Belos::Errors,
	     "The type(s) of solver information that should "
	     "be written to the output stream.");
    pl->set ("Output Style", (int) Belos::General,
	     "What style is used for the solver information written "
	     "to the output stream.");
    pl->set ("Output Frequency", static_cast<int>(-1),
	     "How often (in terms of number of iterations) intermediate "
	     "convergence information should be written to the output stream."
	     "  -1 means never.");
    pl->set ("Output Stream", rcpFromRef(std::cout),
	     "A reference-counted pointer to the output stream where all "
	     "solver output is sent.  The output stream defaults to stdout.");
    pl->set ("Timer Label", std::string("Belos"),
	     "The string to use as a prefix for the timer labels.");
    return pl;
  }

  //
  // Empty Constructor
  //
  template<class ScalarType, class MV, class OP>
  MinresSolMgr<ScalarType,MV,OP>::MinresSolMgr () :
    numIters_ (0),
    parametersSet_ (false)
  {}

  //
  // Primary constructor (use this one)
  //
  template<class ScalarType, class MV, class OP>
  MinresSolMgr<ScalarType, MV, OP>::
  MinresSolMgr (const Teuchos::RCP<LinearProblem<ScalarType, MV, OP> > &problem,
		const Teuchos::RCP<Teuchos::ParameterList>& params) :
    problem_ (problem),
    numIters_ (0),
    parametersSet_ (false)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(problem_.is_null(), std::invalid_argument,
			       "MinresSolMgr: The version of the constructor "
			       "that takes a LinearProblem to solve was given a "
			       "null LinearProblem.");
    setParameters (params);
  }

  template<class ScalarType, class MV, class OP>
  void
  MinresSolMgr<ScalarType, MV, OP>::
  validateProblem (const Teuchos::RCP<LinearProblem<ScalarType, MV, OP> >& problem)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(problem.is_null(),
      MinresSolMgrLinearProblemFailure,
      "MINRES requires that you have provided a nonnull LinearProblem to the "
      "solver manager, before you call the solve() method.");
    TEUCHOS_TEST_FOR_EXCEPTION(problem->getOperator().is_null(),
      MinresSolMgrLinearProblemFailure,
      "MINRES requires a LinearProblem object with a non-null operator (the "
      "matrix A).");
    TEUCHOS_TEST_FOR_EXCEPTION(problem->getRHS().is_null(),
      MinresSolMgrLinearProblemFailure,
      "MINRES requires a LinearProblem object with a non-null right-hand side.");
    TEUCHOS_TEST_FOR_EXCEPTION( ! problem->isProblemSet(),
      MinresSolMgrLinearProblemFailure,
      "MINRES requires that before you give it a LinearProblem to solve, you "
      "must first call the linear problem's setProblem() method.");
  }

  template<class ScalarType, class MV, class OP>
  void
  MinresSolMgr<ScalarType, MV, OP>::
  setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    using Teuchos::null;
    using Teuchos::is_null;
    using std::string;
    using std::ostream;
    using std::endl;

    if (params_.is_null()) {
      params_ = parameterList (*getValidParameters());
    }
    RCP<ParameterList> pl = params;
    pl->validateParametersAndSetDefaults (*params_);

    //
    // Read parameters from the parameter list.  We have already
    // populated it with defaults.
    //
    blockSize_ = pl->get<int> ("Block Size");
    verbosity_ = pl->get<int> ("Verbosity");
    outputStyle_ = pl->get<int> ("Output Style");
    outputFreq_ = pl->get<int>("Output Frequency");
    outputStream_ = pl->get<RCP<std::ostream> > ("Output Stream");
    convtol_ = pl->get<MagnitudeType> ("Convergence Tolerance");
    maxIters_ = pl->get<int> ("Maximum Iterations");
    //
    // All done reading parameters from the parameter list.
    // Now we know it's valid and we can store it.
    //
    params_ = pl;

    // Change the timer label, and create the timer if necessary.
    const string newLabel = pl->get<string> ("Timer Label");
    {
      if (newLabel != label_ || timerSolve_.is_null()) {
	label_ = newLabel;
#ifdef BELOS_TEUCHOS_TIME_MONITOR
	const string solveLabel = label_ + ": MinresSolMgr total solve time";
	// Unregister the old timer before creating a new one.
	if (! timerSolve_.is_null()) {
	  Teuchos::TimeMonitor::clearCounter (label_);
	  timerSolve_ = Teuchos::null;
	}
	timerSolve_ = Teuchos::TimeMonitor::getNewCounter (solveLabel);
#endif // BELOS_TEUCHOS_TIME_MONITOR
      }
    }

    // Create output manager, if necessary; otherwise, set its parameters.
    bool recreatedPrinter = false;
    if (printer_.is_null()) {
      printer_ = rcp (new OutputManager<ScalarType> (verbosity_, outputStream_));
      recreatedPrinter = true;
    } else {
      // Set the output stream's verbosity level.
      printer_->setVerbosity (verbosity_);
      // Tell the output manager about the new output stream.
      printer_->setOStream (outputStream_);
    }

    //
    // Set up the convergence tests
    //
    typedef StatusTestGenResNorm<ScalarType, MV, OP> res_norm_type;
    typedef StatusTestCombo<ScalarType, MV, OP> combo_type;

    // Do we need to allocate at least one of the implicit or explicit
    // residual norm convergence tests?
    const bool allocatedConvergenceTests =
      impConvTest_.is_null() || expConvTest_.is_null();

    // Allocate or set the tolerance of the implicit residual norm
    // convergence test.
    if (impConvTest_.is_null()) {
      impConvTest_ = rcp (new res_norm_type (convtol_));
      impConvTest_->defineResForm (res_norm_type::Implicit, TwoNorm);
      // TODO (mfh 03 Nov 2011) Allow users to define the type of
      // scaling (or a custom scaling factor).
      impConvTest_->defineScaleForm (NormOfInitRes, TwoNorm);
    } else {
      impConvTest_->setTolerance (convtol_);
    }

    // Allocate or set the tolerance of the explicit residual norm
    // convergence test.
    if (expConvTest_.is_null()) {
      expConvTest_ = rcp (new res_norm_type (convtol_));
      expConvTest_->defineResForm (res_norm_type::Explicit, TwoNorm);
      // TODO (mfh 03 Nov 2011) Allow users to define the type of
      // scaling (or a custom scaling factor).
      expConvTest_->defineScaleForm (NormOfInitRes, TwoNorm);
    } else {
      expConvTest_->setTolerance (convtol_);
    }

    // Whether we need to recreate the full status test.  We only need
    // to do that if at least one of convTest_ or maxIterTest_ had to
    // be reallocated.
    bool needToRecreateFullStatusTest = sTest_.is_null();

    // Residual status test is a combo of the implicit and explicit
    // convergence tests.
    if (convTest_.is_null() || allocatedConvergenceTests) {
      convTest_ = rcp (new combo_type (combo_type::SEQ, impConvTest_, expConvTest_));
      needToRecreateFullStatusTest = true;
    }

    // Maximum number of iterations status test.  It tells the solver to
    // stop iteration, if the maximum number of iterations has been
    // exceeded.  Initialize it if we haven't yet done so, otherwise
    // tell it the new maximum number of iterations.
    if (maxIterTest_.is_null()) {
      maxIterTest_ = rcp (new StatusTestMaxIters<ScalarType,MV,OP> (maxIters_));
      needToRecreateFullStatusTest = true;
    } else {
      maxIterTest_->setMaxIters (maxIters_);
    }

    // Create the full status test if we need to.
    //
    // The full status test: the maximum number of iterations have
    // been reached, OR the residual has converged.
    //
    // "If we need to" means either that the status test was never
    // created before, or that its two component tests had to be
    // reallocated.
    if (needToRecreateFullStatusTest) {
      sTest_ = rcp (new combo_type (combo_type::OR, maxIterTest_, convTest_));
    }

    // If necessary, create the status test output class.  This class
    // manages and formats the output from the status test.  We have
    // to recreate the output test if we had to (re)allocate either
    // printer_ or sTest_.
    if (outputTest_.is_null() || needToRecreateFullStatusTest || recreatedPrinter) {
      StatusTestOutputFactory<ScalarType,MV,OP> stoFactory (outputStyle_);
      outputTest_ = stoFactory.create (printer_, sTest_, outputFreq_,
				       Passed+Failed+Undefined);
    } else {
      outputTest_->setOutputFrequency (outputFreq_);
    }
    // Set the solver string for the output test.
    // StatusTestOutputFactory has no constructor argument for this.
    outputTest_->setSolverDesc (std::string (" MINRES "));

    // Inform the solver manager that the current parameters were set.
    parametersSet_ = true;

    if (verbosity_ & Debug) {
      using std::endl;

      std::ostream& dbg = printer_->stream (Debug);
      dbg << "MINRES parameters:" << endl << params_ << endl;
    }
  }


  template<class ScalarType, class MV, class OP>
  ReturnType MinresSolMgr<ScalarType,MV,OP>::solve()
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    using std::endl;

    if (! parametersSet_) {
      setParameters (params_);
    }
    std::ostream& dbg = printer_->stream (Debug);

#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor solveTimerMonitor (*timerSolve_);
#endif // BELOS_TEUCHOS_TIME_MONITOR

    // We need a problem to solve, else we can't solve it.
    validateProblem (problem_);

    // Reset the status test for this solve.
    outputTest_->reset();

    // The linear problem has this many right-hand sides to solve.
    // MINRES can solve only one at a time, so we solve for each
    // right-hand side in succession.
    const int numRHS2Solve = MVT::GetNumberVecs (*(problem_->getRHS()));

    // Create MINRES iteration object.  Pass along the solver
    // manager's parameters, which have already been validated.
    typedef MinresIter<ScalarType, MV, OP> iter_type;
    RCP<iter_type> minres_iter =
      rcp (new iter_type (problem_, printer_, outputTest_, *params_));

    // The index/indices of the right-hand sides for which MINRES did
    // _not_ converge.  Hopefully this is empty after the for loop
    // below!  If it is not empty, at least one right-hand side did
    // not converge.
    std::vector<int> notConverged;
    std::vector<int> currentIndices(1);

    numIters_ = 0;

    // Solve for each right-hand side in turn.
    for (int currentRHS = 0; currentRHS < numRHS2Solve; ++currentRHS) {
      // Inform the linear problem of the right-hand side(s) currently
      // being solved.  MINRES only knows how to solve linear problems
      // with one right-hand side, so we only include one index, which
      // is the index of the current right-hand side.
      currentIndices[0] = currentRHS;
      problem_->setLSIndex (currentIndices);

      dbg << "-- Current right-hand side index being solved: "
	  << currentRHS << endl;

      // Reset the number of iterations.
      minres_iter->resetNumIters();
      // Reset the number of calls that the status test output knows about.
      outputTest_->resetNumCalls();
      // Set the new state and initialize the solver.
      MinresIterationState<ScalarType, MV> newstate;

      // Get the residual vector for the current linear system
      // (that is, for the current right-hand side).
      newstate.Y = MVT::CloneViewNonConst (*(rcp_const_cast<MV> (problem_->getInitResVec())), currentIndices);
      minres_iter->initializeMinres (newstate);

      // Attempt to solve for the solution corresponding to the
      // current right-hand side.
      while (true) {
	try {
	  minres_iter->iterate();

	  // First check for convergence
	  if (convTest_->getStatus() == Passed) {
	    dbg << "---- Converged after " << maxIterTest_->getNumIters()
		<< " iterations" << endl;
	    break;
	  }
	  // Now check for max # of iterations
	  else if (maxIterTest_->getStatus() == Passed) {
	    dbg << "---- Did not converge after " << maxIterTest_->getNumIters()
		<< " iterations" << endl;
	    // This right-hand side didn't converge!
	    notConverged.push_back (currentRHS);
	    break;
	  } else {
	    // If we get here, we returned from iterate(), but none of
	    // our status tests Passed.  Something is wrong, and it is
	    // probably our fault.
	    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
              "Belos::MinresSolMgr::solve(): iterations neither converged, "
	      "nor reached the maximum number of iterations " << maxIters_
	      << ".  That means something went wrong.");
	  }
	} catch (const std::exception &e) {
	  printer_->stream (Errors)
	    << "Error! Caught std::exception in MinresIter::iterate() at "
	    << "iteration " << minres_iter->getNumIters() << endl
	    << e.what() << endl;
	  throw e;
	}
      }

      // Inform the linear problem that we are finished with the
      // current right-hand side.  It may or may not have converged,
      // but we don't try again if the first time didn't work.
      problem_->setCurrLS();

      // Get iteration information for this solve: total number of
      // iterations for all right-hand sides.
      numIters_ += maxIterTest_->getNumIters();
    }

    // Print final summary of the solution process
    sTest_->print (printer_->stream (FinalSummary));

    // Print timing information, if the corresponding compile-time and
    // run-time options are enabled.
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    // Calling summarize() can be expensive, so don't call unless the
    // user wants to print out timing details.  summarize() will do all
    // the work even if it's passed a "black hole" output stream.
    if (verbosity_ & TimingDetails) {
      Teuchos::TimeMonitor::summarize (printer_->stream (TimingDetails));
    }
#endif // BELOS_TEUCHOS_TIME_MONITOR

    // Save the convergence test value ("achieved tolerance") for this
    // solve.  This solver always has two residual norm status tests:
    // an explicit and an implicit test.  The master convergence test
    // convTest_ is a SEQ combo of the implicit resp. explicit tests.
    // If the implicit test never passes, then the explicit test won't
    // ever be executed.  This manifests as
    // expConvTest_->getTestValue()->size() < 1.  We deal with this
    // case by using the values returned by
    // impConvTest_->getTestValue().
    {
      const std::vector<MagnitudeType>* pTestValues = expConvTest_->getTestValue();
      if (pTestValues == NULL || pTestValues->size() < 1) {
	pTestValues = impConvTest_->getTestValue();
      }
      TEUCHOS_TEST_FOR_EXCEPTION(pTestValues == NULL, std::logic_error,
        "Belos::MinresSolMgr::solve(): The implicit convergence test's getTestValue() "
	"method returned NULL.  Please report this bug to the Belos developers.");
      TEUCHOS_TEST_FOR_EXCEPTION(pTestValues->size() < 1, std::logic_error,
        "Belos::MinresSolMgr::solve(): The implicit convergence test's getTestValue() "
        "method returned a vector of length zero.  Please report this bug to the "
        "Belos developers.");

      // FIXME (mfh 12 Dec 2011) Does pTestValues really contain the
      // achieved tolerances for all vectors in the current solve(), or
      // just for the vectors from the last deflation?
      achievedTol_ = *std::max_element (pTestValues->begin(), pTestValues->end());
    }

    if (notConverged.size() > 0) {
      return Unconverged;
    } else {
      return Converged;
    }
  }

  //  This method requires the solver manager to return a std::string that describes itself.
  template<class ScalarType, class MV, class OP>
  std::string MinresSolMgr<ScalarType,MV,OP>::description() const
  {
    std::ostringstream oss;
    oss << "Belos::MinresSolMgr< "
	<< Teuchos::ScalarTraits<ScalarType>::name()
	<<", MV, OP >";
    // oss << "{";
    // oss << "Block Size=" << blockSize_;
    // oss << "}";
    return oss.str();
  }

} // end Belos namespace

#endif /* BELOS_MINRES_SOLMGR_HPP */
