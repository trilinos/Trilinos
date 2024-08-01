// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_PSEUDO_BLOCK_STOCHASTIC_CG_SOLMGR_HPP
#define BELOS_PSEUDO_BLOCK_STOCHASTIC_CG_SOLMGR_HPP

/*! \file BelosPseudoBlockStochasticCGSolMgr.hpp
 *  \brief The Belos::PseudoBlockStochasticCGSolMgr provides a solver manager for the stochastic BlockCG linear solver.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "BelosPseudoBlockStochasticCGIter.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestGenResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosStatusTestOutputFactory.hpp"
#include "BelosOutputManager.hpp"
#ifdef BELOS_TEUCHOS_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

/*! \class Belos::PseudoBlockStochasticCGSolMgr
 *
 *  \brief The Belos::PseudoBlockStochasticCGSolMgr provides a powerful and fully-featured solver manager over the pseudo-block CG iteration.

 \ingroup belos_solver_framework

 \author Chris Siefert, Heidi Thornquist, Chris Baker, and Teri Barth
 */

namespace Belos {

  //! @name PseudoBlockStochasticCGSolMgr Exceptions
  //@{

  /** \brief PseudoBlockStochasticCGSolMgrLinearProblemFailure is thrown when the linear problem is
   * not setup (i.e. setProblem() was not called) when solve() is called.
   *
   * This std::exception is thrown from the PseudoBlockStochasticCGSolMgr::solve() method.
   *
   */
  class PseudoBlockStochasticCGSolMgrLinearProblemFailure : public BelosError {public:
    PseudoBlockStochasticCGSolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  template<class ScalarType, class MV, class OP>
  class PseudoBlockStochasticCGSolMgr : public SolverManager<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;

  public:

    //! @name Constructors/Destructor
    //@{

    /*! \brief Empty constructor for BlockStochasticCGSolMgr.
     * This constructor takes no arguments and sets the default values for the solver.
     * The linear problem must be passed in using setProblem() before solve() is called on this object.
     * The solver values can be changed using setParameters().
     */
    PseudoBlockStochasticCGSolMgr();

    /*! \brief Basic constructor for PseudoBlockStochasticCGSolMgr.
     *
     * This constructor accepts the LinearProblem to be solved in addition
     * to a parameter list of options for the solver manager. These options include the following:
     *   - "Maximum Iterations" - a \c int specifying the maximum number of iterations the underlying solver is allowed to perform.
     *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Belos::Errors
     *   - "Output Style" - a OutputType specifying the style of output. Default: Belos::General
     *   - "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide convergence.
     */
    PseudoBlockStochasticCGSolMgr( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                         const Teuchos::RCP<Teuchos::ParameterList> &pl );

    //! Destructor.
    virtual ~PseudoBlockStochasticCGSolMgr() {};

    //! clone for Inverted Injection (DII)
    Teuchos::RCP<SolverManager<ScalarType, MV, OP> > clone () const override {
      return Teuchos::rcp(new PseudoBlockStochasticCGSolMgr<ScalarType,MV,OP>);
    }
    //@}

    //! @name Accessor methods
    //@{

    const LinearProblem<ScalarType,MV,OP>& getProblem() const override {
      return *problem_;
    }

    /*! \brief Get a parameter list containing the valid parameters for this object.
     */
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

    /*! \brief Get a parameter list containing the current parameters for this object.
     */
    Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters() const override { return params_; }

    /*! \brief Return the timers for this object.
     *
     * The timers are ordered as follows:
     *   - time spent in solve() routine
     */
    Teuchos::Array<Teuchos::RCP<Teuchos::Time> > getTimers() const {
      return Teuchos::tuple(timerSolve_);
    }

    //! Get the iteration count for the most recent call to \c solve().
    int getNumIters() const override {
      return numIters_;
    }

    /*! \brief Return whether a loss of accuracy was detected by this solver during the most current solve.
        \note This flag will be reset the next time solve() is called.
     */
    bool isLOADetected() const override { return false; }

    //@}

    //! @name Set methods
    //@{

    //! Set the linear problem that needs to be solved.
    void setProblem( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem ) override { problem_ = problem; }

    //! Set the parameters the solver manager should use to solve the linear problem.
    void setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params ) override;

    //@}

    //! @name Reset methods
    //@{
    /*! \brief Performs a reset of the solver manager specified by the \c ResetType.  This informs the
     *  solver manager that the solver should prepare for the next call to solve by resetting certain elements
     *  of the iterative solver strategy.
     */
    void reset( const ResetType type ) override { if ((type & Belos::Problem) && !Teuchos::is_null(problem_)) problem_->setProblem(); }
    //@}

    //! @name Solver application methods
    //@{

    /*! \brief This method performs possibly repeated calls to the underlying linear solver's iterate() routine
     * until the problem has been solved (as decided by the solver manager) or the solver manager decides to
     * quit.
     *
     * This method calls PseudoBlockStochasticCGIter::iterate(), which will return either because a specially constructed status test evaluates to
     * ::Passed or an std::exception is thrown.
     *
     * A return from PseudoBlockStochasticCGIter::iterate() signifies one of the following scenarios:
     *    - the maximum number of restarts has been exceeded. In this scenario, the current solutions to the linear system
     *      will be placed in the linear problem and return ::Unconverged.
     *    - global convergence has been met. In this case, the current solutions to the linear system will be placed in the linear
     *      problem and the solver manager will return ::Converged
     *
     * \returns ::ReturnType specifying:
     *     - ::Converged: the linear problem was solved to the specification required by the solver manager.
     *     - ::Unconverged: the linear problem was not solved to the specification desired by the solver manager.
     */
    ReturnType solve() override;

    //@}

    //! Get a copy of the final stochastic vector
    Teuchos::RCP<MV> getStochasticVector() { return Y_;}

    /** \name Overridden from Teuchos::Describable */
    //@{

    /** \brief Method to return description of the block CG solver manager */
    std::string description() const override;

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
    Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP> > convTest_;
    Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputTest_;

    // Current parameter list.
    Teuchos::RCP<Teuchos::ParameterList> params_;

    /// \brief List of valid parameters and their default values.
    ///
    /// This is declared "mutable" because the SolverManager interface
    /// requires that getValidParameters() be declared const, yet we
    /// want to create the valid parameter list only on demand.
    mutable Teuchos::RCP<const Teuchos::ParameterList> validParams_;

    // Default solver values.
    static constexpr int maxIters_default_ = 1000;
    static constexpr bool assertPositiveDefiniteness_default_ = true;
    static constexpr bool showMaxResNormOnly_default_ = false;
    static constexpr int verbosity_default_ = Belos::Errors;
    static constexpr int outputStyle_default_ = Belos::General;
    static constexpr int outputFreq_default_ = -1;
    static constexpr int defQuorum_default_ = 1;
    static constexpr const char * resScale_default_ = "Norm of Initial Residual";
    static constexpr const char * label_default_ = "Belos";

    // Current solver values.
    MagnitudeType convtol_;
    int maxIters_, numIters_;
    int verbosity_, outputStyle_, outputFreq_, defQuorum_;
    bool assertPositiveDefiniteness_, showMaxResNormOnly_;
    std::string resScale_;

    // Timers.
    std::string label_;
    Teuchos::RCP<Teuchos::Time> timerSolve_;

    // Internal state variables.
    bool isSet_;

    // Stashed copy of the stochastic vector
    Teuchos::RCP<MV> Y_;

  };


// Empty Constructor
template<class ScalarType, class MV, class OP>
PseudoBlockStochasticCGSolMgr<ScalarType,MV,OP>::PseudoBlockStochasticCGSolMgr() :
  outputStream_(Teuchos::rcpFromRef(std::cout)),
  convtol_(DefaultSolverParameters::convTol),
  maxIters_(maxIters_default_),
  numIters_(0),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  defQuorum_(defQuorum_default_),
  assertPositiveDefiniteness_(assertPositiveDefiniteness_default_),
  showMaxResNormOnly_(showMaxResNormOnly_default_),
  resScale_(resScale_default_),
  label_(label_default_),
  isSet_(false)
{}

// Basic Constructor
template<class ScalarType, class MV, class OP>
PseudoBlockStochasticCGSolMgr<ScalarType,MV,OP>::
PseudoBlockStochasticCGSolMgr (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                               const Teuchos::RCP<Teuchos::ParameterList> &pl ) :
  problem_(problem),
  outputStream_(Teuchos::rcpFromRef(std::cout)),
  convtol_(DefaultSolverParameters::convTol),
  maxIters_(maxIters_default_),
  numIters_(0),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  defQuorum_(defQuorum_default_),
  assertPositiveDefiniteness_(assertPositiveDefiniteness_default_),
  showMaxResNormOnly_(showMaxResNormOnly_default_),
  resScale_(resScale_default_),
  label_(label_default_),
  isSet_(false)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    problem_.is_null (), std::invalid_argument,
    "Belos::PseudoBlockStochasticCGSolMgr two-argument constructor: "
    "'problem' is null.  You must supply a non-null Belos::LinearProblem "
    "instance when calling this constructor.");

  if (! pl.is_null ()) {
    // Set the parameters using the list that was passed in.
    setParameters (pl);
  }
}

template<class ScalarType, class MV, class OP>
void PseudoBlockStochasticCGSolMgr<ScalarType,MV,OP>::setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params )
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;

  RCP<const ParameterList> defaultParams = getValidParameters();

  // Create the internal parameter list if one doesn't already exist.
  if (params_.is_null()) {
    params_ = parameterList (*defaultParams);
  } else {
    params->validateParameters (*defaultParams);
  }

  // Check for maximum number of iterations
  if (params->isParameter("Maximum Iterations")) {
    maxIters_ = params->get("Maximum Iterations",maxIters_default_);

    // Update parameter in our list and in status test.
    params_->set("Maximum Iterations", maxIters_);
    if (maxIterTest_!=Teuchos::null)
      maxIterTest_->setMaxIters( maxIters_ );
  }

  // Check if positive definiteness assertions are to be performed
  if (params->isParameter("Assert Positive Definiteness")) {
    assertPositiveDefiniteness_ = params->get("Assert Positive Definiteness",assertPositiveDefiniteness_default_);

    // Update parameter in our list.
    params_->set("Assert Positive Definiteness", assertPositiveDefiniteness_);
  }

  // Check to see if the timer label changed.
  if (params->isParameter("Timer Label")) {
    std::string tempLabel = params->get("Timer Label", label_default_);

    // Update parameter in our list and solver timer
    if (tempLabel != label_) {
      label_ = tempLabel;
      params_->set("Timer Label", label_);
      std::string solveLabel = label_ + ": PseudoBlockStochasticCGSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerSolve_ = Teuchos::TimeMonitor::getNewCounter(solveLabel);
#endif
    }
  }

  // Check for a change in verbosity level
  if (params->isParameter("Verbosity")) {
    if (Teuchos::isParameterType<int>(*params,"Verbosity")) {
      verbosity_ = params->get("Verbosity", verbosity_default_);
    } else {
      verbosity_ = (int)Teuchos::getParameter<Belos::MsgType>(*params,"Verbosity");
    }

    // Update parameter in our list.
    params_->set("Verbosity", verbosity_);
    if (printer_ != Teuchos::null)
      printer_->setVerbosity(verbosity_);
  }

  // Check for a change in output style
  if (params->isParameter("Output Style")) {
    if (Teuchos::isParameterType<int>(*params,"Output Style")) {
      outputStyle_ = params->get("Output Style", outputStyle_default_);
    } else {
      outputStyle_ = (int)Teuchos::getParameter<Belos::OutputType>(*params,"Output Style");
    }

    // Reconstruct the convergence test if the explicit residual test is not being used.
    params_->set("Output Style", outputStyle_);
    outputTest_ = Teuchos::null;
  }

  // output stream
  if (params->isParameter("Output Stream")) {
    outputStream_ = Teuchos::getParameter<Teuchos::RCP<std::ostream> >(*params,"Output Stream");

    // Update parameter in our list.
    params_->set("Output Stream", outputStream_);
    if (printer_ != Teuchos::null)
      printer_->setOStream( outputStream_ );
  }

  // frequency level
  if (verbosity_ & Belos::StatusTestDetails) {
    if (params->isParameter("Output Frequency")) {
      outputFreq_ = params->get("Output Frequency", outputFreq_default_);
    }

    // Update parameter in out list and output status test.
    params_->set("Output Frequency", outputFreq_);
    if (outputTest_ != Teuchos::null)
      outputTest_->setOutputFrequency( outputFreq_ );
  }

  // Create output manager if we need to.
  if (printer_ == Teuchos::null) {
    printer_ = Teuchos::rcp( new OutputManager<ScalarType>(verbosity_, outputStream_) );
  }

  // Convergence
  typedef Belos::StatusTestCombo<ScalarType,MV,OP>  StatusTestCombo_t;
  typedef Belos::StatusTestGenResNorm<ScalarType,MV,OP>  StatusTestResNorm_t;

  // Check for convergence tolerance
  if (params->isParameter("Convergence Tolerance")) {
    if (params->isType<MagnitudeType> ("Convergence Tolerance")) {
      convtol_ = params->get ("Convergence Tolerance",
                              static_cast<MagnitudeType> (DefaultSolverParameters::convTol));
    }
    else {
      convtol_ = params->get ("Convergence Tolerance", DefaultSolverParameters::convTol);
    }

    // Update parameter in our list and residual tests.
    params_->set("Convergence Tolerance", convtol_);
    if (convTest_ != Teuchos::null)
      convTest_->setTolerance( convtol_ );
  }

  if (params->isParameter("Show Maximum Residual Norm Only")) {
    showMaxResNormOnly_ = Teuchos::getParameter<bool>(*params,"Show Maximum Residual Norm Only");

    // Update parameter in our list and residual tests
    params_->set("Show Maximum Residual Norm Only", showMaxResNormOnly_);
    if (convTest_ != Teuchos::null)
      convTest_->setShowMaxResNormOnly( showMaxResNormOnly_ );
  }

  // Check for a change in scaling, if so we need to build new residual tests.
  bool newResTest = false;
  {
    // "Residual Scaling" is the old parameter name; "Implicit
    // Residual Scaling" is the new name.  We support both options for
    // backwards compatibility.
    std::string tempResScale = resScale_;
    bool implicitResidualScalingName = false;
    if (params->isParameter ("Residual Scaling")) {
      tempResScale = params->get<std::string> ("Residual Scaling");
    }
    else if (params->isParameter ("Implicit Residual Scaling")) {
      tempResScale = params->get<std::string> ("Implicit Residual Scaling");
      implicitResidualScalingName = true;
    }

    // Only update the scaling if it's different.
    if (resScale_ != tempResScale) {
      Belos::ScaleType resScaleType = convertStringToScaleType( tempResScale );
      resScale_ = tempResScale;

      // Update parameter in our list and residual tests, using the
      // given parameter name.
      if (implicitResidualScalingName) {
        params_->set ("Implicit Residual Scaling", resScale_);
      }
      else {
        params_->set ("Residual Scaling", resScale_);
      }

      if (! convTest_.is_null()) {
        try {
          convTest_->defineScaleForm( resScaleType, Belos::TwoNorm );
        }
        catch (std::exception& e) {
          // Make sure the convergence test gets constructed again.
          newResTest = true;
        }
      }
    }
  }

  // Get the deflation quorum, or number of converged systems before deflation is allowed
  if (params->isParameter("Deflation Quorum")) {
    defQuorum_ = params->get("Deflation Quorum", defQuorum_);
    params_->set("Deflation Quorum", defQuorum_);
    if (convTest_ != Teuchos::null)
      convTest_->setQuorum( defQuorum_ );
  }

  // Create status tests if we need to.

  // Basic test checks maximum iterations and native residual.
  if (maxIterTest_ == Teuchos::null)
    maxIterTest_ = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>( maxIters_ ) );

  // Implicit residual test, using the native residual to determine if convergence was achieved.
  if (convTest_ == Teuchos::null || newResTest) {
    convTest_ = Teuchos::rcp( new StatusTestResNorm_t( convtol_, defQuorum_, showMaxResNormOnly_ ) );
    convTest_->defineScaleForm( convertStringToScaleType( resScale_ ), Belos::TwoNorm );
  }

  if (sTest_ == Teuchos::null || newResTest)
    sTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::OR, maxIterTest_, convTest_ ) );

  if (outputTest_ == Teuchos::null || newResTest) {

    // Create the status test output class.
    // This class manages and formats the output from the status test.
    StatusTestOutputFactory<ScalarType,MV,OP> stoFactory( outputStyle_ );
    outputTest_ = stoFactory.create( printer_, sTest_, outputFreq_, Passed+Failed+Undefined );

    // Set the solver string for the output test
    std::string solverDesc = " Pseudo Block CG ";
    outputTest_->setSolverDesc( solverDesc );

  }

  // Create the timer if we need to.
  if (timerSolve_ == Teuchos::null) {
    std::string solveLabel = label_ + ": PseudoBlockStochasticCGSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    timerSolve_ = Teuchos::TimeMonitor::getNewCounter(solveLabel);
#endif
  }

  // Inform the solver manager that the current parameters were set.
  isSet_ = true;
}


template<class ScalarType, class MV, class OP>
Teuchos::RCP<const Teuchos::ParameterList>
PseudoBlockStochasticCGSolMgr<ScalarType,MV,OP>::getValidParameters() const
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;

  if (validParams_.is_null()) {
    // Set all the valid parameters and their default values.

    // The static_cast is to resolve an issue with older clang versions which
    // would cause the constexpr to link fail. With c++17 the problem is resolved.
    RCP<ParameterList> pl = parameterList ();
    pl->set("Convergence Tolerance", static_cast<MagnitudeType>(DefaultSolverParameters::convTol),
      "The relative residual tolerance that needs to be achieved by the\n"
      "iterative solver in order for the linera system to be declared converged.");
    pl->set("Maximum Iterations", static_cast<int>(maxIters_default_),
      "The maximum number of block iterations allowed for each\n"
      "set of RHS solved.");
    pl->set("Assert Positive Definiteness", static_cast<bool>(assertPositiveDefiniteness_default_),
      "Whether or not to assert that the linear operator\n"
      "and the preconditioner are indeed positive definite.");
    pl->set("Verbosity", static_cast<int>(verbosity_default_),
      "What type(s) of solver information should be outputted\n"
      "to the output stream.");
    pl->set("Output Style", static_cast<int>(outputStyle_default_),
      "What style is used for the solver information outputted\n"
      "to the output stream.");
    pl->set("Output Frequency", static_cast<int>(outputFreq_default_),
      "How often convergence information should be outputted\n"
      "to the output stream.");
    pl->set("Deflation Quorum", static_cast<int>(defQuorum_default_),
      "The number of linear systems that need to converge before\n"
      "they are deflated.  This number should be <= block size.");
    pl->set("Output Stream", Teuchos::rcpFromRef(std::cout),
      "A reference-counted pointer to the output stream where all\n"
      "solver output is sent.");
    pl->set("Show Maximum Residual Norm Only", static_cast<bool>(showMaxResNormOnly_default_),
      "When convergence information is printed, only show the maximum\n"
      "relative residual norm when the block size is greater than one.");
    pl->set("Implicit Residual Scaling", resScale_default_,
      "The type of scaling used in the residual convergence test.");
    // We leave the old name as a valid parameter for backwards
    // compatibility (so that validateParametersAndSetDefaults()
    // doesn't raise an exception if it encounters "Residual
    // Scaling").  The new name was added for compatibility with other
    // solvers, none of which use "Residual Scaling".
    pl->set("Residual Scaling", resScale_default_,
            "The type of scaling used in the residual convergence test.  This "
            "name is deprecated; the new name is \"Implicit Residual Scaling\".");
    pl->set("Timer Label", static_cast<const char *>(label_default_),
      "The string to use as a prefix for the timer labels.");
    validParams_ = pl;
  }
  return validParams_;
}


// solve()
template<class ScalarType, class MV, class OP>
ReturnType PseudoBlockStochasticCGSolMgr<ScalarType,MV,OP>::solve() {

  // Set the current parameters if they were not set before.
  // NOTE:  This may occur if the user generated the solver manager with the default constructor and
  // then didn't set any parameters using setParameters().
  if (!isSet_) { setParameters( params_ ); }

  TEUCHOS_TEST_FOR_EXCEPTION(!problem_->isProblemSet(),PseudoBlockStochasticCGSolMgrLinearProblemFailure,
                     "Belos::PseudoBlockStochasticCGSolMgr::solve(): Linear problem is not ready, setProblem() has not been called.");

  // Create indices for the linear systems to be solved.
  int startPtr = 0;
  int numRHS2Solve = MVT::GetNumberVecs( *(problem_->getRHS()) );
  int numCurrRHS = numRHS2Solve;

  std::vector<int> currIdx( numRHS2Solve ), currIdx2( numRHS2Solve );
  for (int i=0; i<numRHS2Solve; ++i) {
    currIdx[i] = startPtr+i;
    currIdx2[i]=i;
  }

  // Inform the linear problem of the current linear system to solve.
  problem_->setLSIndex( currIdx );

  //////////////////////////////////////////////////////////////////////////////////////
  // Parameter list
  Teuchos::ParameterList plist;

  plist.set("Assert Positive Definiteness",assertPositiveDefiniteness_);

  // Reset the status test.
  outputTest_->reset();

  // Assume convergence is achieved, then let any failed convergence set this to false.
  bool isConverged = true;

  //////////////////////////////////////////////////////////////////////////////////////
  // Pseudo-Block CG solver

  Teuchos::RCP<PseudoBlockStochasticCGIter<ScalarType,MV,OP> > block_cg_iter
    = Teuchos::rcp( new PseudoBlockStochasticCGIter<ScalarType,MV,OP>(problem_,printer_,outputTest_,plist) );

  // Enter solve() iterations
  {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor slvtimer(*timerSolve_);
#endif

    while ( numRHS2Solve > 0 ) {

      // Reset the active / converged vectors from this block
      std::vector<int> convRHSIdx;
      std::vector<int> currRHSIdx( currIdx );
      currRHSIdx.resize(numCurrRHS);

      // Reset the number of iterations.
      block_cg_iter->resetNumIters();

      // Reset the number of calls that the status test output knows about.
      outputTest_->resetNumCalls();

      // Get the current residual for this block of linear systems.
      Teuchos::RCP<MV> R_0 = MVT::CloneViewNonConst( *(Teuchos::rcp_const_cast<MV>(problem_->getInitResVec())), currIdx );

      // Get a new state struct and initialize the solver.
      StochasticCGIterationState<ScalarType,MV> newState;
      newState.R = R_0;
      block_cg_iter->initializeCG(newState);

      while(1) {

        // tell block_gmres_iter to iterate
        try {
          block_cg_iter->iterate();

          ////////////////////////////////////////////////////////////////////////////////////
          //
          // check convergence first
          //
          ////////////////////////////////////////////////////////////////////////////////////
          if ( convTest_->getStatus() == Passed ) {

            // Figure out which linear systems converged.
            std::vector<int> convIdx = Teuchos::rcp_dynamic_cast<StatusTestGenResNorm<ScalarType,MV,OP> >(convTest_)->convIndices();

            // If the number of converged linear systems is equal to the
            // number of current linear systems, then we are done with this block.
            if (convIdx.size() == currRHSIdx.size())
              break;  // break from while(1){block_cg_iter->iterate()}

            // Inform the linear problem that we are finished with this current linear system.
            problem_->setCurrLS();

            // Reset currRHSIdx to have the right-hand sides that are left to converge for this block.
            int have = 0;
            std::vector<int> unconvIdx(currRHSIdx.size());
            for (unsigned int i=0; i<currRHSIdx.size(); ++i) {
              bool found = false;
              for (unsigned int j=0; j<convIdx.size(); ++j) {
                if (currRHSIdx[i] == convIdx[j]) {
                  found = true;
                  break;
                }
              }
              if (!found) {
                currIdx2[have] = currIdx2[i];
                currRHSIdx[have++] = currRHSIdx[i];
              }
            }
            currRHSIdx.resize(have);
            currIdx2.resize(have);

            // Set the remaining indices after deflation.
            problem_->setLSIndex( currRHSIdx );

            // Get the current residual vector.
            std::vector<MagnitudeType> norms;
            R_0 = MVT::CloneCopy( *(block_cg_iter->getNativeResiduals(&norms)),currIdx2 );
            for (int i=0; i<have; ++i) { currIdx2[i] = i; }

            // Set the new state and initialize the solver.
            StochasticCGIterationState<ScalarType,MV> defstate;
            defstate.R = R_0;
            block_cg_iter->initializeCG(defstate);
          }

          ////////////////////////////////////////////////////////////////////////////////////
          //
          // check for maximum iterations
          //
          ////////////////////////////////////////////////////////////////////////////////////
          else if ( maxIterTest_->getStatus() == Passed ) {
            // we don't have convergence
            isConverged = false;
            break;  // break from while(1){block_cg_iter->iterate()}
          }

          ////////////////////////////////////////////////////////////////////////////////////
          //
          // we returned from iterate(), but none of our status tests Passed.
          // something is wrong, and it is probably our fault.
          //
          ////////////////////////////////////////////////////////////////////////////////////

          else {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                               "Belos::PseudoBlockStochasticCGSolMgr::solve(): Invalid return from PseudoBlockStochasticCGIter::iterate().");
          }
        }
        catch (const std::exception &e) {
          printer_->stream(Errors) << "Error! Caught std::exception in PseudoBlockStochasticCGIter::iterate() at iteration "
                                   << block_cg_iter->getNumIters() << std::endl
                                   << e.what() << std::endl;
          throw;
        }
      }

      // Inform the linear problem that we are finished with this block linear system.
      problem_->setCurrLS();

      // Update indices for the linear systems to be solved.
      startPtr += numCurrRHS;
      numRHS2Solve -= numCurrRHS;

      if ( numRHS2Solve > 0 ) {

        numCurrRHS = numRHS2Solve;
        currIdx.resize( numCurrRHS );
        currIdx2.resize( numCurrRHS );
        for (int i=0; i<numCurrRHS; ++i)
          { currIdx[i] = startPtr+i; currIdx2[i] = i; }

        // Set the next indices.
        problem_->setLSIndex( currIdx );
      }
      else {
        currIdx.resize( numRHS2Solve );
      }

    }// while ( numRHS2Solve > 0 )

  }

  // get the final stochastic vector
  Y_=block_cg_iter->getStochasticVector();


  // print final summary
  sTest_->print( printer_->stream(FinalSummary) );

  // print timing information
#ifdef BELOS_TEUCHOS_TIME_MONITOR
  // Calling summarize() can be expensive, so don't call unless the
  // user wants to print out timing details.  summarize() will do all
  // the work even if it's passed a "black hole" output stream.
  if (verbosity_ & TimingDetails)
    Teuchos::TimeMonitor::summarize( printer_->stream(TimingDetails) );
#endif

  // get iteration information for this solve
  numIters_ = maxIterTest_->getNumIters();

  if (!isConverged ) {
    return Unconverged; // return from PseudoBlockStochasticCGSolMgr::solve()
  }
  return Converged; // return from PseudoBlockStochasticCGSolMgr::solve()
}

//  This method requires the solver manager to return a std::string that describes itself.
template<class ScalarType, class MV, class OP>
std::string PseudoBlockStochasticCGSolMgr<ScalarType,MV,OP>::description() const
{
  std::ostringstream oss;
  oss << "Belos::PseudoBlockStochasticCGSolMgr<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
  oss << "{";
  oss << "}";
  return oss.str();
}

} // end Belos namespace

#endif /* BELOS_PSEUDO_BLOCK_CG_SOLMGR_HPP */
