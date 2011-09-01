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

    /// \brief List of default parameters
    /// 
    /// Return a smart pointer to a list of valid parameters and their
    /// default values.  The right way to set up this solver manager
    /// with nondefault parameters, is to make a deep non-const copy
    /// of the default parameters, and change the parameters you want
    /// to change.
    ///
    /// \note This is a class ("static") method, so that it can be
    /// called before you have constructed a MinresSolMgr object
    /// (useful for you, so you can have valid parameters before
    /// constructing one), or within the MinresSolMgr constructor
    /// (useful for me, so I can set default parameters).
    static Teuchos::RCP< const Teuchos::ParameterList > defaultParameters();

    /// \brief Validate parameters
    ///
    /// Validate the given non-null list of parameters.  Raise
    /// std::invalid_argument if any provided parameters have invalid
    /// values.  If any parameters that should be there are not
    /// provided, fill them in with default values.
    ///
    /// \param params [in/out] Non-null smart pointer to parameter list
    ///
    /// \note We make this method a public class ("static") method, so
    ///   you can validate your own list of parameters before creating
    ///   a MinresSolMgr.
    ///
    static void
    validateParameters (const Teuchos::RCP< Teuchos::ParameterList >& params);
    
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
    MinresSolMgr (const Teuchos::RCP< LinearProblem< ScalarType, MV, OP > > &problem,
                  const Teuchos::RCP< Teuchos::ParameterList > &params);
    
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
    Teuchos::RCP< const Teuchos::ParameterList > getValidParameters() const {
      return defaultParameters();
    }
    
    //! Return the list of current parameters for this object.
    Teuchos::RCP< const Teuchos::ParameterList > getCurrentParameters() const { 
      return params_; 
    }
    
    /// \brief Return all timers for this object.
    ///
    /// Return all the timers for this object.  Currently only one
    /// timer is being used, which is the total time spent in the
    /// solve() routine.  Thus, the returned Array currently has only
    /// one element.
    Teuchos::Array< Teuchos::RCP< Teuchos::Time > > getTimers() const {
      return tuple(timerSolve_);
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
    
    //! @name Set methods
    //@{
   
    //! Set the linear problem to be solved.
    void 
    setProblem (const Teuchos::RCP< LinearProblem< ScalarType, MV, OP > > &problem) 
    { 
      validateProblem (problem);
      problem_ = problem; 
    }
   
    /// Set the parameters to use when solving the linear problem.
    ///
    /// \param params [in/out] List of parameters to use when solving
    ///   the linear problem.  This list will be modified as necessary
    ///   to include default parameters that need not be provided.
    ///   If params is null, then use default parameters entirely.
    ///
    /// \note If you want nondefault parameter values, the right way
    ///   to start is to make a deep copy of the ParameterList
    ///   returned by defaultParameters().  That ParameterList has all
    ///   the parameters that MinresSolMgr wants, along with
    ///   human-readable documentation and validators.
    void 
    setParameters (const Teuchos::RCP< Teuchos::ParameterList >& params);
    
    //@}
   
    //! @name Reset methods
    //@{

    /// \brief Reset the solver manager
    ///
    /// Reset the solver manager in a way specified by the \c
    /// ResetType.  This informs the solver manager that the solver
    /// should prepare for the next call to solve by resetting certain
    /// elements of the iterative solver strategy.
    void 
    reset (const ResetType type) 
    { 
      if ((type & Belos::Problem) && ! Teuchos::is_null (problem_)) 
	problem_->setProblem(); 
    }
    //@}
 
    //! @name Solver application methods
    //@{ 
    
    /// \brief Iterate until status test tells us to stop
    //
    /// This method performs possibly repeated calls to the underlying
    /// linear solver's iterate() routine, until the problem has been
    /// solved (as decided by the solver manager via the status
    /// test(s)), or the solver manager decides to quit.
    ///
    /// This method calls MinresIter::iterate(), which will return
    /// either because a specially constructed status test that
    /// evaluates to ::Passed, or an std::exception is thrown.
    ///
    /// A return from MinresIter::iterate() signifies one of the
    /// following scenarios:
    /// - the maximum number of iterations has been exceeded. In this
    ///   scenario, the current solutions to the linear system will be
    ///   placed in the linear problem and return ::Unconverged.
    /// - global convergence has been met. In this case, the current
    ///   solutions to the linear system will be placed in the linear
    ///   problem and the solver manager will return ::Converged
    ///
    /// \return ::ReturnType specifying:
    ///   - ::Converged: the linear problem was solved to the
    ///     specification required by the solver manager.
    ///   - ::Unconverged: the linear problem was not solved to the
    ///     specification desired by the solver manager.
    ReturnType solve();
    
    //@}
    
    /** \name Overridden from Teuchos::Describable */
    //@{
    
    //! Description of the MINRES solver manager
    std::string description() const;
    
    //@}
    
  private:
    /// Validate the given linear problem, by raising
    /// std::invalid_argument (with an informative message) if the
    /// problem is null or its essential components are null.
    static void
    validateProblem (const Teuchos::RCP< LinearProblem< ScalarType, MV, OP > >& problem);

    /// Read the parameters' values from the given parameter list, and
    /// set up solver components as necessary.
    void 
    readParameters (const Teuchos::RCP< const Teuchos::ParameterList >& params);

    //! Linear problem to solve
    Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > problem_;
    
    //! Output manager.
    Teuchos::RCP<OutputManager<ScalarType> > printer_;
    Teuchos::RCP<std::ostream> outputStream_;

    //! The full status test
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > sTest_;

    //! One of sTest_'s components
    Teuchos::RCP<StatusTestMaxIters<ScalarType,MV,OP> > maxIterTest_;

    //! One of sTest_'s components
    Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP> > convTest_;

    //! Status test output
    Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputTest_;

    //! List of current parameters
    Teuchos::RCP< Teuchos::ParameterList > params_;

    //! Current relative residual 2-norm convergence tolerance
    MagnitudeType convtol_;

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
  };


  template<class ScalarType, class MV, class OP>
  Teuchos::RCP< const Teuchos::ParameterList > 
  MinresSolMgr< ScalarType, MV, OP >::defaultParameters()
  {
    using Teuchos::is_null;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::EnhancedNumberValidator;
    typedef Teuchos::ScalarTraits< MagnitudeType > MST;
    typedef Teuchos::ScalarTraits< int > IST;

    // List of parameters accepted by MINRES, and their default values.
    static RCP< const ParameterList > validPL;
  
    // If we haven't already initialized the list of accepted parameters
    // and their default values, do so.
    if (is_null (validPL)) 
      {
	// validPL is a pointer to a const ParameterList, so we need to
	// make a non-const ParameterList first in order to fill it with
	// parameters and their default values.
	RCP< ParameterList > pl (new ParameterList);

	pl->set ("Convergence Tolerance", static_cast< MagnitudeType >(1e-8),
		 "Relative residual tolerance that needs to be achieved by "
		 "the iterative solver, in order for the linear system to be "
		 "declared converged.",
		 rcp (new EnhancedNumberValidator< MagnitudeType > (MST::zero(), MST::rmax())));
	pl->set ("Maximum Iterations", static_cast<int>(1000),
		 "Maximum number of iterations allowed for each right-hand "
		 "side solved.",
		 rcp (new EnhancedNumberValidator< int > (0, INT_MAX)));
	pl->set ("Block Size", static_cast<int>(1),
		 "Number of vectors in each block.  WARNING: The current "
		 "implementation of MINRES only accepts a block size of 1, "
		 "since it can only solve for 1 right-hand side at a time.",
		 rcp (new EnhancedNumberValidator< int > (1, 1)));  // [1,1] is inclusive range
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
	pl->set ("Output Stream", rcp(&std::cout, false),
		 "A reference-counted pointer to the output stream where all "
		 "solver output is sent.  The output stream defaults to stdout.");
	pl->set ("Timer Label", std::string("Belos"),
		 "The string to use as a prefix for the timer labels.");
	validPL = pl;
      }
    return validPL;
  }

  //
  // Empty Constructor
  //
  template<class ScalarType, class MV, class OP>
  MinresSolMgr<ScalarType,MV,OP>::MinresSolMgr () :
    numIters_ (0),
    parametersSet_ (false)
  {
    // Pass in a null parameter list so setParameters grabs the default parameter list.
    Teuchos::RCP< Teuchos::ParameterList > nullParams = Teuchos::null;
    setParameters ( nullParams );
  }

  //
  // Primary constructor (use this one)
  //
  template<class ScalarType, class MV, class OP>
  MinresSolMgr< ScalarType, MV, OP >::
  MinresSolMgr (const Teuchos::RCP< LinearProblem< ScalarType, MV, OP > > &problem,
		const Teuchos::RCP< Teuchos::ParameterList >& params) :
    problem_ (problem),
    numIters_ (0),
    parametersSet_ (false)
  {
    TEST_FOR_EXCEPTION(problem_ == Teuchos::null, std::invalid_argument, "Problem not given to solver manager.");

    // If the parameter list pointer is null, then set the current parameters to the default parameter list.
    if ( !is_null(params) ) {
      setParameters( params );  
    }
  }

  template<class ScalarType, class MV, class OP>
  void
  MinresSolMgr< ScalarType, MV, OP >::
  validateProblem (const Teuchos::RCP< LinearProblem< ScalarType, MV, OP > > &problem) 
  {
    TEST_FOR_EXCEPTION(Teuchos::is_null(problem), std::invalid_argument, 
		       "MINRES requires a non-null LinearProblem object,"
		       "which represents the linear problem to solve.");
    TEST_FOR_EXCEPTION(Teuchos::is_null(problem->getOperator()), 
		       std::invalid_argument, 
		       "MINRES requires a LinearProblem object with a non-null "
		       "operator.");
    TEST_FOR_EXCEPTION(Teuchos::is_null(problem->getRHS()),
		       std::invalid_argument, 
		       "MINRES requires a LinearProblem object with a non-null "
		       "right-hand side.");
  }

  template< class ScalarType, class MV, class OP >
  void
  MinresSolMgr< ScalarType, MV, OP>::
  validateParameters (const Teuchos::RCP< Teuchos::ParameterList >& params)
  {
    using Teuchos::Exceptions::InvalidParameterName;
    using Teuchos::is_null;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

    if (is_null (params)) 
      throw std::logic_error ("MinresSolMgr::validateParameters() requires a "
			      "non-null input, but was given a null input.  Th"
			      "is is likely the result of a bug in the impleme"
			      "ntation of MinresSolMgr.");
    // List of default, known-valid parameters.
    Teuchos::RCP< const Teuchos::ParameterList > defaults = defaultParameters();

    // Validate parameters' values in params, and add default values
    // for parameters that are not in params.
    //
    // NOTE (mfh 06 Dec 2010) This throws
    // Teuchos::Exceptions::InvalidParameterName if params has a
    // parameter that defaults doesn't have (not vice versa).  (That
    // would mean the user provided "extra" parameters.)  We should
    // really just let those go through, for the sake of forwards
    // compatibility.  We could do this by implementing a variant of
    // validateParametersAndSetDefaults() that ignores parameters in
    // params that are not in defaults.
    params->validateParametersAndSetDefaults(*defaults);
  }

  template<class ScalarType, class MV, class OP>
  void 
  MinresSolMgr< ScalarType, MV, OP>::
  readParameters (const Teuchos::RCP< const Teuchos::ParameterList >& params)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::null;
    using Teuchos::is_null;
    using std::string;
    using std::ostream;

    // Set the block size.  The only valid value is 1, for this
    // particular solver.
    blockSize_ = params->get<int> ("Block Size");

    // Change the timer label.
    {
      label_ = params->get< string > ("Timer Label");
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      const string solveLabel = label_ + ": MinresSolMgr total solve time";
      timerSolve_ = Teuchos::TimeMonitor::getNewTimer (solveLabel);
#endif
    }

    // Set verbosity level
    verbosity_ = params->get< int > ("Verbosity");
    // We may not have created the output manager yet.  If we haven't,
    // we will create it below and set the verbosity level
    // appropriately.
    if (! is_null (printer_))
      printer_->setVerbosity (verbosity_);

    // Set output style
    outputStyle_ = params->get<int> ("Output Style");

    // Set output stream
    outputStream_ = params->get< RCP< std::ostream > > ("Output Stream");
    // Tell the output manager about the new output stream.  We
    // may not yet have created the output manager yet, so check
    // here if it is not null.  We will create the output manager
    // if necessary below.
    if (! Teuchos::is_null (printer_))
      printer_->setOStream (outputStream_);

    // Set output frequency level
    if (Belos::StatusTestDetails) 
      {
	outputFreq_ = params->get<int>("Output Frequency");
	// We may not have created the status test output object yet.
	if (! is_null (outputTest_))
	  outputTest_->setOutputFrequency (outputFreq_);
      }

    // Create output manager
    printer_ = rcp (new OutputManager< ScalarType >(verbosity_, outputStream_));
  
    //
    // Set up the convergence tests
    //
    typedef Belos::StatusTestCombo< ScalarType, MV, OP > StatusTestCombo_t;
    typedef Belos::StatusTestGenResNorm< ScalarType, MV, OP > StatusTestResNorm_t;

    // Set convergence tolerance
    convtol_ = params->get< MagnitudeType >("Convergence Tolerance");

    // Residual status test.  It uses the native residual to determine
    // if convergence was achieved.  Initialize it if we haven't yet
    // done so, otherwise tell it the new convergence tolerance.
    if (is_null (convTest_))
      convTest_ = rcp (new StatusTestResNorm_t (convtol_, -1));
    else
      convTest_->setTolerance (convtol_);

    // Set maximum number of iterations.
    maxIters_ = params->get<int>("Maximum Iterations");

    // Maximum number of iterations status test.  It tells the solver to
    // stop iteration, if the maximum number of iterations has been
    // exceeded.  Initialize it if we haven't yet done so, otherwise
    // tell it the new maximum number of iterations.
    if (is_null (maxIterTest_))
      maxIterTest_ = rcp (new StatusTestMaxIters<ScalarType,MV,OP> (maxIters_));
    else
      maxIterTest_->setMaxIters (maxIters_);

    // Create the full status test if we need to.  
    //
    // The full status test: maximum number of iterations reached, OR
    // residual has converged.
    //
    // "If we need to" means if the status test was never created
    // before.  The full status test has pointers to maxIterTest_ and
    // convTest_, so if we changed either the convergence tolerance
    // and/or the maximum number of iterations, the full status test
    // will get the results of the updates.
    if (is_null (sTest_))
      sTest_ = rcp (new StatusTestCombo_t (StatusTestCombo_t::OR, maxIterTest_, convTest_));
  
    // If necessary, create the status test output class.  This class
    // manages and formats the output from the status test.
    if (is_null (outputTest_))
      {
	StatusTestOutputFactory<ScalarType,MV,OP> stoFactory (outputStyle_);
	outputTest_ = stoFactory.create (printer_, sTest_, outputFreq_,
					 Passed+Failed+Undefined);
	// Set the solver string for the output test
	outputTest_->setSolverDesc (std::string(" MINRES "));
      }

    // Create the timer if we need to.
    if (is_null (timerSolve_))
      {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
	const std::string solveLabel = label_ + ": MinresSolMgr total solve time";
	timerSolve_ = Teuchos::TimeMonitor::getNewTimer(solveLabel);
#endif
      }
    // Inform the solver manager that the current parameters were set.
    parametersSet_ = true;
  }


  template<class ScalarType, class MV, class OP>
  void 
  MinresSolMgr< ScalarType, MV, OP>::
  setParameters (const Teuchos::RCP< Teuchos::ParameterList >& params)
  {
    if (Teuchos::is_null (params))
      // We don't need to validate the default parameter values.
      // However, we do need to make a deep copy, since params_ is
      // const and the input params is non-const.
      params_ = Teuchos::rcp (new Teuchos::ParameterList (*defaultParameters()));
    else
      {
	// Validate params.  All the entries that should be there, are
	// there, with default values if the entries were not supplied in
	// the input parameter list, and with valid values otherwise.
	validateParameters (params);
	params_ = params; // both are non-const
      }
    readParameters (params_);
  }


  template<class ScalarType, class MV, class OP>
  ReturnType MinresSolMgr<ScalarType,MV,OP>::solve() 
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor slvtimer(*timerSolve_);
#endif

    // We need a problem to solve, else we can't solve it.
    TEST_FOR_EXCEPTION( !problem_->isProblemSet(),
			MinresSolMgrLinearProblemFailure,
			"Belos::MinresSolMgr::solve(): The linear problem to "
			"solve has not been set, so it is not valid to call the "
			"solve() method.  You should either provide the linear "
			"problem to the solver manager's constructor, or call "
			"setProblem() with a non-null linear problem, before "
			"calling solve()." );

    // Reset the status test for this solve.  
    outputTest_->reset();

    // The linear problem has this many right-hand sides to solve.
    // MINRES can solve one at a time, so we solve for each right-hand
    // side in succession.
    const int numRHS2Solve = MVT::GetNumberVecs (*(problem_->getRHS()));

    // Create MINRES iteration object.  Pass along the solver
    // manager's parameters, which have already been validated.
    RCP< MinresIteration< ScalarType, MV, OP > > minres_iter =
      rcp (new MinresIter< ScalarType, MV, OP > (problem_, printer_, outputTest_, *params_));

    // The index/indices of the right-hand sides for which MINRES did
    // _not_ converge.  Hopefully this is empty after the for loop below!
    // If it is not empty, at least one right-hand side did not converge.
    std::vector<int> notConverged;
    std::vector<int> currentIndices(1);

    numIters_ = 0;

    // Solve for each right-hand side in turn.
    for (int currentRHS = 0; currentRHS < numRHS2Solve; ++currentRHS)
      {
	// Inform the linear problem of the right-hand side(s) currently
	// being solved.  MINRES only knows how to solve linear problems
	// with one right-hand side, so we only include one index, which
	// is the index of the current right-hand side.
	currentIndices[0] = currentRHS;
	problem_->setLSIndex (currentIndices);

	// Reset the number of iterations.
	minres_iter->resetNumIters();
	// Reset the number of calls that the status test output knows about.
	outputTest_->resetNumCalls();
	// Set the new state and initialize the solver.
	MinresIterationState< ScalarType, MV > newstate;

	// Get the residual vector for the current linear system
	// (i.e., for the current right-hand side).
	
	newstate.Y = MVT::CloneViewNonConst (*(Teuchos::rcp_const_cast< MV > (problem_->getInitResVec())), currentIndices);
	minres_iter->initializeMinres (newstate);

	// Attempt to solve for the solution corresponding to the
	// current right-hand side.
	while(1) {
          try {
	    minres_iter->iterate();

            // First check for convergence
	    if (convTest_->getStatus() == Passed) {
	      break;
	    }
            // Now check for max # of iterations
            else if (maxIterTest_->getStatus() == Passed) {
	      // This right-hand side didn't converge!
	      notConverged.push_back (currentRHS);
              break;
            } else {
	      // If we get here, we returned from iterate(), but none of
	      // our status tests Passed.  Something is wrong, and it is
	      // probably our fault.
	      TEST_FOR_EXCEPTION(true, std::logic_error,
                               "Belos::MinresSolMgr::solve(): iterations neither"
                               " converged, nor reached the maximum number of "
                               "iterations " << maxIters_ << ".  That means someth"
                               "ing went wrong.");
	      }
	  } catch (const std::exception &e) {
	    printer_->stream(Errors) 
	      << "Error! Caught std::exception in "
	      << "MinresIter::iterate() at iteration "
	      << minres_iter->getNumIters() << std::endl
	      << e.what() << std::endl;
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
    sTest_->print( printer_->stream(FinalSummary) );

    // Print timing information, if the corresponding compile-time and
    // run-time options are enabled.
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    // Calling summarize() can be expensive, so don't call unless the
    // user wants to print out timing details.  summarize() will do all
    // the work even if it's passed a "black hole" output stream.
    if (verbosity_ & TimingDetails)
      Teuchos::TimeMonitor::summarize (printer_->stream (TimingDetails));
#endif // BELOS_TEUCHOS_TIME_MONITOR
 
    if (notConverged.size() > 0)
      return Unconverged;
    else
      return Converged;
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
