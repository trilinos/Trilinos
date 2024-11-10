// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_PSEUDO_BLOCK_CG_SOLMGR_HPP
#define BELOS_PSEUDO_BLOCK_CG_SOLMGR_HPP

/*! \file BelosPseudoBlockCGSolMgr.hpp
 *  \brief The Belos::PseudoBlockCGSolMgr provides a solver manager for the BlockCG linear solver.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "BelosPseudoBlockCGIter.hpp"
#include "BelosCGSingleRedIter.hpp"
#include "BelosCGIter.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestGenResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosStatusTestOutputFactory.hpp"
#include "BelosOutputManager.hpp"
#include "Teuchos_LAPACK.hpp"
#ifdef BELOS_TEUCHOS_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

/** \example epetra/example/BlockCG/PseudoBlockCGEpetraExFile.cpp
    This is an example of how to use the Belos::PseudoBlockCGSolMgr solver manager using Epetra.
*/
/** \example tpetra/example/BlockCG/PseudoBlockCGTpetraExFile.cpp
    This is an example of how to use the Belos::PseudoBlockCGSolMgr solver manager using Tpetra.
*/
/** \example epetra/example/BlockCG/PseudoBlockPrecCGEpetraExFile.cpp
    This is an example of how to use the Belos::PseudoBlockCGSolMgr solver manager with an Ifpack preconditioner.
*/

/*! \class Belos::PseudoBlockCGSolMgr
 *
 *  \brief The Belos::PseudoBlockCGSolMgr provides a powerful and fully-featured solver manager over the pseudo-block CG iteration.

 \ingroup belos_solver_framework

 \author Heidi Thornquist, Chris Baker, and Teri Barth
 */

namespace Belos {

  //! @name PseudoBlockCGSolMgr Exceptions
  //@{

  /** \brief PseudoBlockCGSolMgrLinearProblemFailure is thrown when the linear problem is
   * not setup (i.e. setProblem() was not called) when solve() is called.
   *
   * This std::exception is thrown from the PseudoBlockCGSolMgr::solve() method.
   *
   */
  class PseudoBlockCGSolMgrLinearProblemFailure : public BelosError {public:
    PseudoBlockCGSolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};


  // Partial specialization for unsupported ScalarType types.
  // This contains a stub implementation.
  template<class ScalarType, class MV, class OP,
           const bool supportsScalarType =
             Belos::Details::LapackSupportsScalar<ScalarType>::value>
  class PseudoBlockCGSolMgr :
    public Details::SolverManagerRequiresLapack<ScalarType, MV, OP,
                                                Belos::Details::LapackSupportsScalar<ScalarType>::value>
  {
    static const bool scalarTypeIsSupported =
      Belos::Details::LapackSupportsScalar<ScalarType>::value;
    typedef Details::SolverManagerRequiresLapack<ScalarType, MV, OP,
                                                 scalarTypeIsSupported> base_type;

  public:
    PseudoBlockCGSolMgr () :
      base_type ()
    {}
    PseudoBlockCGSolMgr (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                         const Teuchos::RCP<Teuchos::ParameterList> &pl) :
      base_type ()
    {}
    virtual ~PseudoBlockCGSolMgr () {}

    Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP> >
    getResidualStatusTest() const { return Teuchos::null; }
  };


  template<class ScalarType, class MV, class OP>
  class PseudoBlockCGSolMgr<ScalarType, MV, OP, true> :
    public Details::SolverManagerRequiresLapack<ScalarType, MV, OP, true>
  {
  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;

  public:

    //! @name Constructors/Destructor
    //@{

    /*! \brief Empty constructor for BlockCGSolMgr.
     * This constructor takes no arguments and sets the default values for the solver.
     * The linear problem must be passed in using setProblem() before solve() is called on this object.
     * The solver values can be changed using setParameters().
     */
    PseudoBlockCGSolMgr();

    /*! \brief Basic constructor for PseudoBlockCGSolMgr.
     *
     * This constructor accepts the LinearProblem to be solved in addition
     * to a parameter list of options for the solver manager. These options include the following:
     *   - "Maximum Iterations" - a \c int specifying the maximum number of iterations the underlying solver is allowed to perform.
     *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Belos::Errors
     *   - "Output Style" - a OutputType specifying the style of output. Default: Belos::General
     *   - "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide convergence.
     *          \param pl [in] ParameterList with construction information
     *                  \htmlonly
     *                  <iframe src="belos_PseudoBlockCG.xml" width=100% scrolling="no" frameborder="0">
     *                  </iframe>
     *                  <hr />
     *                  \endhtmlonly
     */
    PseudoBlockCGSolMgr( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                         const Teuchos::RCP<Teuchos::ParameterList> &pl );

    //! Destructor.
    virtual ~PseudoBlockCGSolMgr() {};

    //! clone for Inverted Injection (DII)
    Teuchos::RCP<SolverManager<ScalarType, MV, OP> > clone () const override {
      return Teuchos::rcp(new PseudoBlockCGSolMgr<ScalarType,MV,OP>);
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


    /// \brief Tolerance achieved by the last \c solve() invocation.
    ///
    /// This is the maximum over all right-hand sides' achieved
    /// convergence tolerances, and is set whether or not the solve
    /// actually managed to achieve the desired convergence tolerance.
    ///
    /// \warning This result may not be meaningful if there was a loss
    ///   of accuracy during the solve.  You should first call \c
    ///   isLOADetected() to check for a loss of accuracy during the
    ///   last solve.
    MagnitudeType achievedTol() const override {
      return achievedTol_;
    }

    //! Get the iteration count for the most recent call to \c solve().
    int getNumIters() const override {
      return numIters_;
    }

    /*! \brief Return whether a loss of accuracy was detected by this solver during the most current solve.
        \note This flag will be reset the next time solve() is called.
     */
    bool isLOADetected() const override { return false; }

    /*! \brief Gets the estimated condition number.
      \note Only works if "Estimate Condition Number" is set on parameterlist
    */
    ScalarType getConditionEstimate() const {return condEstimate_;}
    Teuchos::ArrayRCP<MagnitudeType> getEigenEstimates() const {return eigenEstimates_;}

    //! Return the residual status test
    Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP> >
    getResidualStatusTest() const { return convTest_; }

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
     * This method calls PseudoBlockCGIter::iterate(), which will return either because a specially constructed status test evaluates to
     * ::Passed or an std::exception is thrown.
     *
     * A return from PseudoBlockCGIter::iterate() signifies one of the following scenarios:
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

    /** \name Overridden from Teuchos::Describable */
    //@{

    /** \brief Method to return description of the block CG solver manager */
    std::string description() const override;

    //@}
  private:
    // Compute the condition number estimate
    void compute_condnum_tridiag_sym(Teuchos::ArrayView<MagnitudeType> diag,
                                     Teuchos::ArrayView<MagnitudeType> offdiag,
                                     Teuchos::ArrayRCP<MagnitudeType>& lambdas,
                                     ScalarType & lambda_min,
                                     ScalarType & lambda_max,
                                     ScalarType & ConditionNumber );

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
    static constexpr bool foldConvergenceDetectionIntoAllreduce_default_ = false;
    static constexpr const char * resScale_default_ = "Norm of Initial Residual";
    static constexpr const char * label_default_ = "Belos";
    static constexpr bool genCondEst_default_ = false;

    // Current solver values.
    MagnitudeType convtol_,achievedTol_;
    int maxIters_, numIters_;
    int verbosity_, outputStyle_, outputFreq_, defQuorum_;
    bool assertPositiveDefiniteness_, showMaxResNormOnly_;
    bool foldConvergenceDetectionIntoAllreduce_;
    std::string resScale_;
    bool genCondEst_;
    ScalarType condEstimate_;
    Teuchos::ArrayRCP<MagnitudeType> eigenEstimates_;

    // Timers.
    std::string label_;
    Teuchos::RCP<Teuchos::Time> timerSolve_;

    // Internal state variables.
    bool isSet_;
  };


// Empty Constructor
template<class ScalarType, class MV, class OP>
PseudoBlockCGSolMgr<ScalarType,MV,OP,true>::PseudoBlockCGSolMgr() :
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
  foldConvergenceDetectionIntoAllreduce_(foldConvergenceDetectionIntoAllreduce_default_),
  resScale_(resScale_default_),
  genCondEst_(genCondEst_default_),
  condEstimate_(-Teuchos::ScalarTraits<ScalarType>::one()),
  label_(label_default_),
  isSet_(false)
{}

// Basic Constructor
template<class ScalarType, class MV, class OP>
PseudoBlockCGSolMgr<ScalarType,MV,OP,true>::
PseudoBlockCGSolMgr (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
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
  foldConvergenceDetectionIntoAllreduce_(foldConvergenceDetectionIntoAllreduce_default_),
  resScale_(resScale_default_),
  genCondEst_(genCondEst_default_),
  condEstimate_(-Teuchos::ScalarTraits<ScalarType>::one()),
  label_(label_default_),
  isSet_(false)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    problem_.is_null (), std::invalid_argument,
    "Belos::PseudoBlockCGSolMgr two-argument constructor: "
    "'problem' is null.  You must supply a non-null Belos::LinearProblem "
    "instance when calling this constructor.");

  if (! pl.is_null ()) {
    // Set the parameters using the list that was passed in.
    setParameters (pl);
  }
}

template<class ScalarType, class MV, class OP>
void
PseudoBlockCGSolMgr<ScalarType,MV,OP,true>::
setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<const ParameterList> defaultParams = this->getValidParameters ();

  // Create the internal parameter list if one doesn't already exist.
  // Belos' solvers treat the input ParameterList to setParameters as
  // a "delta" -- that is, a change from the current state -- so the
  // default parameter list (if the input is null) should be empty.
  // This explains also why Belos' solvers copy parameters one by one
  // from the input list to the current list.
  //
  // Belos obfuscates the latter, because it takes the input parameter
  // list by RCP, rather than by (nonconst) reference.  The latter
  // would make more sense, given that it doesn't actually keep the
  // input parameter list.
  //
  // Note, however, that Belos still correctly triggers the "used"
  // field of each parameter in the input list.  While isParameter()
  // doesn't (apparently) trigger the "used" flag, get() certainly
  // does.

  if (params_.is_null ()) {
    // Create an empty list with the same name as the default list.
    params_ = parameterList (defaultParams->name ());
  } else {
    params->validateParameters (*defaultParams);
  }

  // Check for maximum number of iterations
  if (params->isParameter ("Maximum Iterations")) {
    maxIters_ = params->get ("Maximum Iterations", maxIters_default_);

    // Update parameter in our list and in status test.
    params_->set ("Maximum Iterations", maxIters_);
    if (! maxIterTest_.is_null ()) {
      maxIterTest_->setMaxIters (maxIters_);
    }
  }

  // Check if positive definiteness assertions are to be performed
  if (params->isParameter ("Assert Positive Definiteness")) {
    assertPositiveDefiniteness_ =
      params->get ("Assert Positive Definiteness",
                   assertPositiveDefiniteness_default_);

    // Update parameter in our list.
    params_->set ("Assert Positive Definiteness", assertPositiveDefiniteness_);
  }

  if (params->isParameter("Fold Convergence Detection Into Allreduce")) {
    foldConvergenceDetectionIntoAllreduce_ = params->get("Fold Convergence Detection Into Allreduce",
                                                         foldConvergenceDetectionIntoAllreduce_default_);
  }

  // Check to see if the timer label changed.
  if (params->isParameter ("Timer Label")) {
    const std::string tempLabel = params->get ("Timer Label", label_default_);

    // Update parameter in our list and solver timer
    if (tempLabel != label_) {
      label_ = tempLabel;
      params_->set ("Timer Label", label_);
      const std::string solveLabel =
        label_ + ": PseudoBlockCGSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerSolve_ = Teuchos::TimeMonitor::getNewCounter (solveLabel);
#endif
    }
  }

  // Check for a change in verbosity level
  if (params->isParameter ("Verbosity")) {
    if (Teuchos::isParameterType<int> (*params, "Verbosity")) {
      verbosity_ = params->get ("Verbosity", verbosity_default_);
    } else {
      verbosity_ = (int) Teuchos::getParameter<Belos::MsgType> (*params, "Verbosity");
    }

    // Update parameter in our list.
    params_->set ("Verbosity", verbosity_);
    if (! printer_.is_null ()) {
      printer_->setVerbosity (verbosity_);
    }
  }

  // Check for a change in output style
  if (params->isParameter ("Output Style")) {
    if (Teuchos::isParameterType<int> (*params, "Output Style")) {
      outputStyle_ = params->get ("Output Style", outputStyle_default_);
    } else {
      // FIXME (mfh 29 Jul 2015) What if the type is wrong?
      outputStyle_ = (int) Teuchos::getParameter<Belos::OutputType> (*params, "Output Style");
    }

    // Reconstruct the convergence test if the explicit residual test
    // is not being used.
    params_->set ("Output Style", outputStyle_);
    outputTest_ = Teuchos::null;
  }

  // output stream
  if (params->isParameter ("Output Stream")) {
    outputStream_ = params->get<RCP<std::ostream> > ("Output Stream");

    // Update parameter in our list.
    params_->set ("Output Stream", outputStream_);
    if (! printer_.is_null ()) {
      printer_->setOStream (outputStream_);
    }
  }

  // frequency level
  if (verbosity_ & Belos::StatusTestDetails) {
    if (params->isParameter ("Output Frequency")) {
      outputFreq_ = params->get ("Output Frequency", outputFreq_default_);
    }

    // Update parameter in out list and output status test.
    params_->set ("Output Frequency", outputFreq_);
    if (! outputTest_.is_null ()) {
      outputTest_->setOutputFrequency (outputFreq_);
    }
  }

  // Condition estimate
  if (params->isParameter ("Estimate Condition Number")) {
    genCondEst_ = params->get ("Estimate Condition Number", genCondEst_default_);
  }

  // Create output manager if we need to.
  if (printer_.is_null ()) {
    printer_ = rcp (new OutputManager<ScalarType> (verbosity_, outputStream_));
  }

  // Convergence
  typedef Belos::StatusTestCombo<ScalarType,MV,OP> StatusTestCombo_t;
  typedef Belos::StatusTestGenResNorm<ScalarType,MV,OP> StatusTestResNorm_t;

  // Check for convergence tolerance
  if (params->isParameter ("Convergence Tolerance")) {
    if (params->isType<MagnitudeType> ("Convergence Tolerance")) {
      convtol_ = params->get ("Convergence Tolerance",
                              static_cast<MagnitudeType> (DefaultSolverParameters::convTol));
    }
    else {
      convtol_ = params->get ("Convergence Tolerance", DefaultSolverParameters::convTol);
    }

    // Update parameter in our list and residual tests.
    params_->set ("Convergence Tolerance", convtol_);
    if (! convTest_.is_null ()) {
      convTest_->setTolerance (convtol_);
    }
  }

  if (params->isParameter ("Show Maximum Residual Norm Only")) {
    showMaxResNormOnly_ = params->get<bool> ("Show Maximum Residual Norm Only");

    // Update parameter in our list and residual tests
    params_->set ("Show Maximum Residual Norm Only", showMaxResNormOnly_);
    if (! convTest_.is_null ()) {
      convTest_->setShowMaxResNormOnly (showMaxResNormOnly_);
    }
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
      const Belos::ScaleType resScaleType =
        convertStringToScaleType (tempResScale);
      resScale_ = tempResScale;

      // Update parameter in our list and residual tests, using the
      // given parameter name.
      if (implicitResidualScalingName) {
        params_->set ("Implicit Residual Scaling", resScale_);
      }
      else {
        params_->set ("Residual Scaling", resScale_);
      }

      if (! convTest_.is_null ()) {
        try {
          convTest_->defineScaleForm (resScaleType, Belos::TwoNorm);
        }
        catch (std::exception& e) {
          // Make sure the convergence test gets constructed again.
          newResTest = true;
        }
      }
    }
  }

  // Get the deflation quorum, or number of converged systems before deflation is allowed
  if (params->isParameter ("Deflation Quorum")) {
    defQuorum_ = params->get ("Deflation Quorum", defQuorum_);
    params_->set ("Deflation Quorum", defQuorum_);
    if (! convTest_.is_null ()) {
      convTest_->setQuorum( defQuorum_ );
    }
  }

  // Create status tests if we need to.

  // Basic test checks maximum iterations and native residual.
  if (maxIterTest_.is_null ()) {
    maxIterTest_ = rcp (new StatusTestMaxIters<ScalarType,MV,OP> (maxIters_));
  }

  // Implicit residual test, using the native residual to determine if convergence was achieved.
  if (convTest_.is_null () || newResTest) {
    convTest_ = rcp (new StatusTestResNorm_t (convtol_, defQuorum_, showMaxResNormOnly_));
    convTest_->defineScaleForm (convertStringToScaleType (resScale_), Belos::TwoNorm);
  }

  if (sTest_.is_null () || newResTest) {
    sTest_ = rcp (new StatusTestCombo_t (StatusTestCombo_t::OR, maxIterTest_, convTest_));
  }

  if (outputTest_.is_null () || newResTest) {
    // Create the status test output class.
    // This class manages and formats the output from the status test.
    StatusTestOutputFactory<ScalarType,MV,OP> stoFactory (outputStyle_);
    outputTest_ = stoFactory.create (printer_, sTest_, outputFreq_,
                                     Passed+Failed+Undefined);

    // Set the solver string for the output test
    const std::string solverDesc = " Pseudo Block CG ";
    outputTest_->setSolverDesc (solverDesc);
  }

  // Create the timer if we need to.
  if (timerSolve_.is_null ()) {
    const std::string solveLabel =
      label_ + ": PseudoBlockCGSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    timerSolve_ = Teuchos::TimeMonitor::getNewCounter (solveLabel);
#endif
  }

  // Inform the solver manager that the current parameters were set.
  isSet_ = true;
}


template<class ScalarType, class MV, class OP>
Teuchos::RCP<const Teuchos::ParameterList>
PseudoBlockCGSolMgr<ScalarType,MV,OP,true>::getValidParameters() const
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;

  if (validParams_.is_null()) {
    // Set all the valid parameters and their default values.
    RCP<ParameterList> pl = parameterList ();
    pl->set("Convergence Tolerance", static_cast<MagnitudeType>(DefaultSolverParameters::convTol),
      "The relative residual tolerance that needs to be achieved by the\n"
      "iterative solver in order for the linear system to be declared converged.");
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
    pl->set("Estimate Condition Number", static_cast<bool>(genCondEst_default_),
      "Whether or not to estimate the condition number of the preconditioned system.");
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
    pl->set("Fold Convergence Detection Into Allreduce",static_cast<bool>(foldConvergenceDetectionIntoAllreduce_default_),
      "Merge the allreduce for convergence detection with the one for CG.\n"
      "This saves one all-reduce, but incurs more computation.");
    validParams_ = pl;
  }
  return validParams_;
}


// solve()
template<class ScalarType, class MV, class OP>
ReturnType PseudoBlockCGSolMgr<ScalarType,MV,OP,true>::solve ()
{
  const char prefix[] = "Belos::PseudoBlockCGSolMgr::solve: ";

  // Set the current parameters if they were not set before.
  // NOTE:  This may occur if the user generated the solver manager with the default constructor and
  // then didn't set any parameters using setParameters().
  if (!isSet_) { setParameters( params_ ); }

  TEUCHOS_TEST_FOR_EXCEPTION
    (! problem_->isProblemSet (), PseudoBlockCGSolMgrLinearProblemFailure,
     prefix << "The linear problem to solve is not ready.  You must call "
     "setProblem() on the Belos::LinearProblem instance before telling the "
     "Belos solver to solve it.");

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
  // Parameter list (iteration)
  Teuchos::ParameterList plist;

  plist.set("Assert Positive Definiteness",assertPositiveDefiniteness_);
  if(genCondEst_) plist.set("Max Size For Condest",maxIters_);

  // Reset the status test.
  outputTest_->reset();

  // Assume convergence is achieved, then let any failed convergence set this to false.
  bool isConverged = true;

  //////////////////////////////////////////////////////////////////////////////////////
  // Pseudo-Block CG solver
  Teuchos::RCP<CGIteration<ScalarType,MV,OP> > block_cg_iter;
  if (numRHS2Solve == 1) {
    plist.set("Fold Convergence Detection Into Allreduce",
              foldConvergenceDetectionIntoAllreduce_);
    block_cg_iter =
      Teuchos::rcp (new CGIter<ScalarType,MV,OP> (problem_, printer_, outputTest_, convTest_, plist));
  } else {
    block_cg_iter =
      Teuchos::rcp (new PseudoBlockCGIter<ScalarType,MV,OP> (problem_, printer_, outputTest_, plist));
  }

  // Setup condition estimate
  block_cg_iter->setDoCondEst(genCondEst_);
  bool condEstPerf = false;

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
      CGIterationState<ScalarType,MV> newState;
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

            // Compute condition estimate if the very first linear system in the block has converged.
            if (currRHSIdx[0] != 0 && genCondEst_ && !condEstPerf)
            {
              // Compute the estimate.
              ScalarType l_min, l_max;
              Teuchos::ArrayView<MagnitudeType> diag    = block_cg_iter->getDiag();
              Teuchos::ArrayView<MagnitudeType> offdiag = block_cg_iter->getOffDiag();
              compute_condnum_tridiag_sym(diag,offdiag,eigenEstimates_,l_min,l_max,condEstimate_);

              // Make sure not to do more condition estimate computations for this solve.
              block_cg_iter->setDoCondEst(false); 
              condEstPerf = true;
            }

            // Set the remaining indices after deflation.
            problem_->setLSIndex( currRHSIdx );

            // Get the current residual vector.
            std::vector<MagnitudeType> norms;
            R_0 = MVT::CloneCopy( *(block_cg_iter->getNativeResiduals(&norms)),currIdx2 );
            for (int i=0; i<have; ++i) { currIdx2[i] = i; }

            // Set the new state and initialize the solver.
            CGIterationState<ScalarType,MV> defstate;
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
                               "Belos::PseudoBlockCGSolMgr::solve(): Invalid return from PseudoBlockCGIter::iterate().");
          }
        }
        catch (const StatusTestNaNError& e) {
          // A NaN was detected in the solver.  Set the solution to zero and return unconverged.
          achievedTol_ = MT::one();
          Teuchos::RCP<MV> X = problem_->getLHS();
          MVT::MvInit( *X, SCT::zero() );
          printer_->stream(Warnings) << "Belos::PseudoBlockCGSolMgr::solve(): Warning! NaN has been detected!" 
                                     << std::endl;
          return Unconverged;
        }
        catch (const std::exception &e) {
          printer_->stream(Errors) << "Error! Caught std::exception in PseudoBlockCGIter::iterate() at iteration "
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

  // Save the convergence test value ("achieved tolerance") for this
  // solve.
  const std::vector<MagnitudeType>* pTestValues = convTest_->getTestValue();
  if (pTestValues != NULL && pTestValues->size () > 0) {
    achievedTol_ = *std::max_element (pTestValues->begin(), pTestValues->end());
  }

  // Do condition estimate, if needed
  if (genCondEst_ && !condEstPerf) {
    ScalarType l_min, l_max;
    Teuchos::ArrayView<MagnitudeType> diag    = block_cg_iter->getDiag();
    Teuchos::ArrayView<MagnitudeType> offdiag = block_cg_iter->getOffDiag();
    compute_condnum_tridiag_sym(diag,offdiag,eigenEstimates_,l_min,l_max,condEstimate_);
    condEstPerf = true;
  }

  if (! isConverged) {
    return Unconverged; // return from PseudoBlockCGSolMgr::solve()
  }
  return Converged; // return from PseudoBlockCGSolMgr::solve()
}

//  This method requires the solver manager to return a std::string that describes itself.
template<class ScalarType, class MV, class OP>
std::string PseudoBlockCGSolMgr<ScalarType,MV,OP,true>::description() const
{
  std::ostringstream oss;
  oss << "Belos::PseudoBlockCGSolMgr<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
  oss << "{";
  oss << "}";
  return oss.str();
}


template<class ScalarType, class MV, class OP>
void
PseudoBlockCGSolMgr<ScalarType,MV,OP,true>::
compute_condnum_tridiag_sym (Teuchos::ArrayView<MagnitudeType> diag,
                             Teuchos::ArrayView<MagnitudeType> offdiag,
                             Teuchos::ArrayRCP<MagnitudeType>& lambdas,
                             ScalarType & lambda_min,
                             ScalarType & lambda_max,
                             ScalarType & ConditionNumber )
{
  typedef Teuchos::ScalarTraits<ScalarType> STS;

  /* Copied from az_cg.c: compute_condnum_tridiag_sym */
  /* diag ==      ScalarType vector of size N, containing the diagonal
     elements of A
     offdiag ==   ScalarType vector of size N-1, containing the offdiagonal
     elements of A. Note that A is supposed to be symmatric
  */
  int info = 0;
  const int N = diag.size ();
  ScalarType scalar_dummy;
  std::vector<MagnitudeType> mag_dummy(4*N);
  char char_N = 'N';
  Teuchos::LAPACK<int,ScalarType> lapack;

  lambdas.resize(N, 0.0);
  lambda_min = STS::one ();
  lambda_max = STS::one ();
  if( N > 2 ) {
    lapack.PTEQR (char_N, N, diag.getRawPtr (), offdiag.getRawPtr (),
                  &scalar_dummy, 1, &mag_dummy[0], &info);
    TEUCHOS_TEST_FOR_EXCEPTION
      (info < 0, std::logic_error, "Belos::PseudoBlockCGSolMgr::"
       "compute_condnum_tridiag_sym: LAPACK's _PTEQR failed with info = "
       << info << " < 0.  This suggests there might be a bug in the way Belos "
       "is calling LAPACK.  Please report this to the Belos developers.");
    for (int k = 0; k < N; k++) {
      lambdas[k] = diag[N - 1 - k];
    }
    lambda_min = Teuchos::as<ScalarType> (diag[N-1]);
    lambda_max = Teuchos::as<ScalarType> (diag[0]);
  }

  // info > 0 means that LAPACK's eigensolver didn't converge.  This
  // is unlikely but may be possible.  In that case, the best we can
  // do is use the eigenvalues that it computes, as long as lambda_max
  // >= lambda_min.
  if (STS::real (lambda_max) < STS::real (lambda_min)) {
    ConditionNumber = STS::one ();
  }
  else {
    // It's OK for the condition number to be Inf.
    ConditionNumber = lambda_max / lambda_min;
  }

} /* compute_condnum_tridiag_sym */





} // end Belos namespace

#endif /* BELOS_PSEUDO_BLOCK_CG_SOLMGR_HPP */
