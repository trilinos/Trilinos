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

#ifndef BELOS_TFQMR_SOLMGR_HPP
#define BELOS_TFQMR_SOLMGR_HPP

/*! \file BelosTFQMRSolMgr.hpp
 *  \brief The Belos::TFQMRSolMgr provides a solver manager for the TFQMR linear solver.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "BelosTFQMRIter.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestGenResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosStatusTestOutputFactory.hpp"
#include "BelosOutputManager.hpp"
#ifdef BELOS_TEUCHOS_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

/** \example TFQMR/TFQMREpetraExFile.cpp
    This is an example of how to use the Belos::TFQMRSolMgr solver manager.
*/

/*! \class Belos::TFQMRSolMgr
 *
 *  \brief The Belos::TFQMRSolMgr provides a powerful and fully-featured solver manager over the TFQMR linear solver.

 \ingroup belos_solver_framework

 \author Heidi Thornquist
*/

namespace Belos {

  //! @name TFQMRSolMgr Exceptions
  //@{

  /** \brief TFQMRSolMgrLinearProblemFailure is thrown when the linear problem is
   * not setup (i.e. setProblem() was not called) when solve() is called.
   *
   * This std::exception is thrown from the TFQMRSolMgr::solve() method.
   *
   */
  class TFQMRSolMgrLinearProblemFailure : public BelosError {public:
    TFQMRSolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief TFQMRSolMgrOrthoFailure is thrown when the orthogonalization manager is
   * unable to generate orthonormal columns from the initial basis vectors.
   *
   * This std::exception is thrown from the TFQMRSolMgr::solve() method.
   *
   */
  class TFQMRSolMgrOrthoFailure : public BelosError {public:
    TFQMRSolMgrOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  template<class ScalarType, class MV, class OP>
  class TFQMRSolMgr : public SolverManager<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;

  public:

    //! @name Constructors/Destructor
    //@{

    /*! \brief Empty constructor for TFQMRSolMgr.
     * This constructor takes no arguments and sets the default values for the solver.
     * The linear problem must be passed in using setProblem() before solve() is called on this object.
     * The solver values can be changed using setParameters().
     */
     TFQMRSolMgr();

    /*! \brief Basic constructor for TFQMRSolMgr.
     *
     * This constructor accepts the LinearProblem to be solved in addition
     * to a parameter list of options for the solver manager. These options include the following:
     *   - "Maximum Iterations" - an \c int specifying the maximum number of iterations the
     *                            underlying solver is allowed to perform. Default: 1000
     *   - "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms
     *                               must reach to decide convergence. Default: 1e-8.
     *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Belos::Errors
     *   - "Output Style" - a OutputType specifying the style of output. Default: Belos::General
     *   - "Output Stream" - a reference-counted pointer to the output stream where all
     *                       solver output is sent.  Default: Teuchos::rcp(&std::cout,false)
     *   - "Output Frequency" - an \c int specifying how often convergence information should be
     *                          outputted.  Default: -1 (never)
     *   - "Timer Label" - a \c std::string to use as a prefix for the timer labels.  Default: "Belos"
     */
    TFQMRSolMgr( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                 const Teuchos::RCP<Teuchos::ParameterList> &pl );

    //! Destructor.
    virtual ~TFQMRSolMgr() {};

    //! clone for Inverted Injection (DII)
    Teuchos::RCP<SolverManager<ScalarType, MV, OP> > clone () const override {
      return Teuchos::rcp(new TFQMRSolMgr<ScalarType,MV,OP>);
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
    MagnitudeType achievedTol() const override {
      return achievedTol_;
    }

    //! Get the iteration count for the most recent call to \c solve().
    int getNumIters() const override {
      return numIters_;
    }

    /// \brief Whether loss of accuracy was detected during the last \c solve() invocation.
    ///
    /// In solvers that can detect a loss of accuracy, this method
    /// would say whether the solver detected it in the most recent \c
    /// solve() invocation.  However, our TFQMR implementation does
    /// not currently detect a loss of accuracy, so this method always
    /// returns false.
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

    /*! \brief This method performs possibly repeated calls to the underlying linear solver's
     *         iterate() routine until the problem has been solved (as decided by the solver manager)
     *         or the solver manager decides to quit.
     *
     * This method calls TFQMRIter::iterate(), which will return either because a
     * specially constructed status test evaluates to ::Passed or an std::exception is thrown.
     *
     * A return from TFQMRIter::iterate() signifies one of the following scenarios:
     *    - the maximum number of iterations has been exceeded. In this scenario, the current solutions
     *      to the linear system will be placed in the linear problem and return ::Unconverged.
     *    - global convergence has been met. In this case, the current solutions to the linear system
     *      will be placed in the linear problem and the solver manager will return ::Converged
     *
     * \returns ::ReturnType specifying:
     *     - ::Converged: the linear problem was solved to the specification required by the solver manager.
     *     - ::Unconverged: the linear problem was not solved to the specification desired by the solver manager.
     */
    ReturnType solve() override;

    //@}
    /** \name Overridden from Teuchos::Describable */
    //@{

    /** \brief Method to return description of the TFQMR solver manager */
    std::string description() const override;
    //@}

  private:

    // Method for checking current status test against defined linear problem.
    bool checkStatusTest();

    // Linear problem.
    Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > problem_;

    // Output manager.
    Teuchos::RCP<OutputManager<ScalarType> > printer_;
    Teuchos::RCP<std::ostream> outputStream_;

    // Status test.
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > sTest_;
    Teuchos::RCP<StatusTestMaxIters<ScalarType,MV,OP> > maxIterTest_;
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > convTest_;
    Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP> > expConvTest_, impConvTest_;
    Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputTest_;

    // Current parameter list.
    Teuchos::RCP<Teuchos::ParameterList> params_;

    // Default solver values.
    static constexpr int maxIters_default_ = 1000;
    static constexpr bool expResTest_default_ = false;
    static constexpr int verbosity_default_ = Belos::Errors;
    static constexpr int outputStyle_default_ = Belos::General;
    static constexpr int outputFreq_default_ = -1;
    static constexpr const char * impResScale_default_ = "Norm of Preconditioned Initial Residual";
    static constexpr const char * expResScale_default_ = "Norm of Initial Residual";
    static constexpr const char * label_default_ = "Belos";
    static constexpr std::ostream * outputStream_default_ = &std::cout;

    // Current solver values.
    MagnitudeType convtol_, impTolScale_, achievedTol_;
    int maxIters_, numIters_;
    int verbosity_, outputStyle_, outputFreq_;
    int blockSize_;
    bool expResTest_;
    std::string impResScale_, expResScale_;

    // Timers.
    std::string label_;
    Teuchos::RCP<Teuchos::Time> timerSolve_;

    // Internal state variables.
    bool isSet_, isSTSet_;
  };


// Empty Constructor
template<class ScalarType, class MV, class OP>
TFQMRSolMgr<ScalarType,MV,OP>::TFQMRSolMgr() :
  outputStream_(Teuchos::rcp(outputStream_default_,false)),
  convtol_(DefaultSolverParameters::convTol),
  impTolScale_(DefaultSolverParameters::impTolScale),
  achievedTol_(Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>::zero()),
  maxIters_(maxIters_default_),
  numIters_(0),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  blockSize_(1),
  expResTest_(expResTest_default_),
  impResScale_(impResScale_default_),
  expResScale_(expResScale_default_),
  label_(label_default_),
  isSet_(false),
  isSTSet_(false)
{}


// Basic Constructor
template<class ScalarType, class MV, class OP>
TFQMRSolMgr<ScalarType,MV,OP>::TFQMRSolMgr(
                                             const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                                             const Teuchos::RCP<Teuchos::ParameterList> &pl ) :
  problem_(problem),
  outputStream_(Teuchos::rcp(outputStream_default_,false)),
  convtol_(DefaultSolverParameters::convTol),
  impTolScale_(DefaultSolverParameters::impTolScale),
  achievedTol_(Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>::zero()),
  maxIters_(maxIters_default_),
  numIters_(0),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  blockSize_(1),
  expResTest_(expResTest_default_),
  impResScale_(impResScale_default_),
  expResScale_(expResScale_default_),
  label_(label_default_),
  isSet_(false),
  isSTSet_(false)
{
  TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null, std::invalid_argument, "Problem not given to solver manager.");

  // If the parameter list pointer is null, then set the current parameters to the default parameter list.
  if ( !is_null(pl) ) {
    setParameters( pl );
  }
}

template<class ScalarType, class MV, class OP>
void TFQMRSolMgr<ScalarType,MV,OP>::setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params )
{
  // Create the internal parameter list if ones doesn't already exist.
  if (params_ == Teuchos::null) {
    params_ = Teuchos::rcp( new Teuchos::ParameterList(*getValidParameters()) );
  }
  else {
    params->validateParameters(*getValidParameters());
  }

  // Check for maximum number of iterations
  if (params->isParameter("Maximum Iterations")) {
    maxIters_ = params->get("Maximum Iterations",maxIters_default_);

    // Update parameter in our list and in status test.
    params_->set("Maximum Iterations", maxIters_);
    if (maxIterTest_!=Teuchos::null)
      maxIterTest_->setMaxIters( maxIters_ );
  }

  // Check for blocksize
  if (params->isParameter("Block Size")) {
    blockSize_ = params->get("Block Size",1);
    TEUCHOS_TEST_FOR_EXCEPTION(blockSize_ != 1, std::invalid_argument,
                       "Belos::TFQMRSolMgr: \"Block Size\" must be 1.");

    // Update parameter in our list.
    params_->set("Block Size", blockSize_);
  }

  // Check to see if the timer label changed.
  if (params->isParameter("Timer Label")) {
    std::string tempLabel = params->get("Timer Label", label_default_);

    // Update parameter in our list and solver timer
    if (tempLabel != label_) {
      label_ = tempLabel;
      params_->set("Timer Label", label_);
      std::string solveLabel = label_ + ": TFQMRSolMgr total solve time";
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
    isSTSet_ = false;
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

  // Check for convergence tolerance
  if (params->isParameter("Convergence Tolerance")) {
    if (params->isType<MagnitudeType> ("Convergence Tolerance")) {
      convtol_ = params->get ("Convergence Tolerance",
                              static_cast<MagnitudeType> (DefaultSolverParameters::convTol));
    }
    else {
      convtol_ = params->get ("Convergence Tolerance", DefaultSolverParameters::convTol);
    }

    // Update parameter in our list.
    params_->set("Convergence Tolerance", convtol_);
    isSTSet_ = false;
  }

  // Check for implicit residual scaling
  if (params->isParameter("Implicit Tolerance Scale Factor")) {
    if (params->isType<MagnitudeType> ("Implicit Tolerance Scale Factor")) {
      impTolScale_ = params->get ("Implicit Tolerance Scale Factor",
                                  static_cast<MagnitudeType> (DefaultSolverParameters::impTolScale));

    }
    else {
      impTolScale_ = params->get ("Implicit Tolerance Scale Factor",
                                  DefaultSolverParameters::impTolScale);
    }

    // Update parameter in our list.
    params_->set("Implicit Tolerance Scale Factor", impTolScale_);
    isSTSet_ = false;
  }

  // Check for a change in scaling, if so we need to build new residual tests.
  if (params->isParameter("Implicit Residual Scaling")) {
    std::string tempImpResScale = Teuchos::getParameter<std::string>( *params, "Implicit Residual Scaling" );

    // Only update the scaling if it's different.
    if (impResScale_ != tempImpResScale) {
      impResScale_ = tempImpResScale;

      // Update parameter in our list and residual tests
      params_->set("Implicit Residual Scaling", impResScale_);

      // Make sure the convergence test gets constructed again.
      isSTSet_ = false;
    }
  }

  if (params->isParameter("Explicit Residual Scaling")) {
    std::string tempExpResScale = Teuchos::getParameter<std::string>( *params, "Explicit Residual Scaling" );

    // Only update the scaling if it's different.
    if (expResScale_ != tempExpResScale) {
      expResScale_ = tempExpResScale;

      // Update parameter in our list and residual tests
      params_->set("Explicit Residual Scaling", expResScale_);

      // Make sure the convergence test gets constructed again.
      isSTSet_ = false;
    }
  }

  if (params->isParameter("Explicit Residual Test")) {
    expResTest_ = Teuchos::getParameter<bool>( *params,"Explicit Residual Test" );

    // Reconstruct the convergence test if the explicit residual test is not being used.
    params_->set("Explicit Residual Test", expResTest_);
    if (expConvTest_ == Teuchos::null) {
      isSTSet_ = false;
    }
  }

  // Create the timer if we need to.
  if (timerSolve_ == Teuchos::null) {
    std::string solveLabel = label_ + ": TFQMRSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    timerSolve_ = Teuchos::TimeMonitor::getNewCounter(solveLabel);
#endif
  }

  // Inform the solver manager that the current parameters were set.
  isSet_ = true;
}


// Check the status test versus the defined linear problem
template<class ScalarType, class MV, class OP>
bool TFQMRSolMgr<ScalarType,MV,OP>::checkStatusTest() {

  typedef Belos::StatusTestCombo<ScalarType,MV,OP>  StatusTestCombo_t;
  typedef Belos::StatusTestGenResNorm<ScalarType,MV,OP>  StatusTestGenResNorm_t;

  // Basic test checks maximum iterations and native residual.
  maxIterTest_ = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>( maxIters_ ) );

  if (expResTest_) {

    // Implicit residual test, using the native residual to determine if convergence was achieved.
    Teuchos::RCP<StatusTestGenResNorm_t> tmpImpConvTest =
      Teuchos::rcp( new StatusTestGenResNorm_t( impTolScale_*convtol_ ) );
    tmpImpConvTest->defineScaleForm( convertStringToScaleType(impResScale_), Belos::TwoNorm );
    impConvTest_ = tmpImpConvTest;

    // Explicit residual test once the native residual is below the tolerance
    Teuchos::RCP<StatusTestGenResNorm_t> tmpExpConvTest =
      Teuchos::rcp( new StatusTestGenResNorm_t( convtol_ ) );
    tmpExpConvTest->defineResForm( StatusTestGenResNorm_t::Explicit, Belos::TwoNorm );
    tmpExpConvTest->defineScaleForm( convertStringToScaleType(expResScale_), Belos::TwoNorm );
    expConvTest_ = tmpExpConvTest;

    // The convergence test is a combination of the "cheap" implicit test and explicit test.
    convTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::SEQ, impConvTest_, expConvTest_ ) );
  }
  else {

    // Implicit residual test, using the native residual to determine if convergence was achieved.
    Teuchos::RCP<StatusTestGenResNorm_t> tmpImpConvTest =
      Teuchos::rcp( new StatusTestGenResNorm_t( convtol_ ) );
    tmpImpConvTest->defineScaleForm( convertStringToScaleType(impResScale_), Belos::TwoNorm );
    impConvTest_ = tmpImpConvTest;

    // Set the explicit and total convergence test to this implicit test that checks for accuracy loss.
    expConvTest_ = impConvTest_;
    convTest_ = impConvTest_;
  }
  sTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::OR, maxIterTest_, convTest_ ) );

  // Create the status test output class.
  // This class manages and formats the output from the status test.
  StatusTestOutputFactory<ScalarType,MV,OP> stoFactory( outputStyle_ );
  outputTest_ = stoFactory.create( printer_, sTest_, outputFreq_, Passed+Failed+Undefined );

  // Set the solver string for the output test
  std::string solverDesc = " TFQMR ";
  outputTest_->setSolverDesc( solverDesc );


  // The status test is now set.
  isSTSet_ = true;

  return false;
}


template<class ScalarType, class MV, class OP>
Teuchos::RCP<const Teuchos::ParameterList>
TFQMRSolMgr<ScalarType,MV,OP>::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;

  // Set all the valid parameters and their default values.
  if(is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

    // The static_cast is to resolve an issue with older clang versions which
    // would cause the constexpr to link fail. With c++17 the problem is resolved.
    pl->set("Convergence Tolerance", static_cast<MagnitudeType>(DefaultSolverParameters::convTol),
      "The relative residual tolerance that needs to be achieved by the\n"
      "iterative solver in order for the linear system to be declared converged.");
    pl->set("Implicit Tolerance Scale Factor", static_cast<MagnitudeType>(DefaultSolverParameters::impTolScale),
      "The scale factor used by the implicit residual test when explicit residual\n"
      "testing is used.  May enable faster convergence when TFQMR bound is too loose.");
    pl->set("Maximum Iterations", static_cast<int>(maxIters_default_),
      "The maximum number of block iterations allowed for each\n"
      "set of RHS solved.");
    pl->set("Verbosity", static_cast<int>(verbosity_default_),
      "What type(s) of solver information should be outputted\n"
      "to the output stream.");
    pl->set("Output Style", static_cast<int>(outputStyle_default_),
      "What style is used for the solver information outputted\n"
      "to the output stream.");
    pl->set("Output Frequency", static_cast<int>(outputFreq_default_),
      "How often convergence information should be outputted\n"
      "to the output stream.");
    pl->set("Output Stream", Teuchos::rcp(outputStream_default_,false),
      "A reference-counted pointer to the output stream where all\n"
      "solver output is sent.");
    pl->set("Explicit Residual Test", static_cast<bool>(expResTest_default_),
      "Whether the explicitly computed residual should be used in the convergence test.");
    pl->set("Implicit Residual Scaling", static_cast<const char *>(impResScale_default_),
      "The type of scaling used in the implicit residual convergence test.");
    pl->set("Explicit Residual Scaling", static_cast<const char *>(expResScale_default_),
      "The type of scaling used in the explicit residual convergence test.");
    pl->set("Timer Label", static_cast<const char *>(label_default_),
      "The string to use as a prefix for the timer labels.");
    validPL = pl;
  }
  return validPL;
}


// solve()
template<class ScalarType, class MV, class OP>
ReturnType TFQMRSolMgr<ScalarType,MV,OP>::solve() {

  // Set the current parameters if they were not set before.
  // NOTE:  This may occur if the user generated the solver manager with the default constructor and
  // then didn't set any parameters using setParameters().
  if (!isSet_) {
    setParameters(Teuchos::parameterList(*getValidParameters()));
  }

  TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null,TFQMRSolMgrLinearProblemFailure,
                     "Belos::TFQMRSolMgr::solve(): Linear problem is not a valid object.");

  TEUCHOS_TEST_FOR_EXCEPTION(!problem_->isProblemSet(),TFQMRSolMgrLinearProblemFailure,
                     "Belos::TFQMRSolMgr::solve(): Linear problem is not ready, setProblem() has not been called.");

  if (!isSTSet_) {
    TEUCHOS_TEST_FOR_EXCEPTION( checkStatusTest(),TFQMRSolMgrLinearProblemFailure,
                        "Belos::TFQMRSolMgr::solve(): Linear problem and requested status tests are incompatible.");
  }

  // Create indices for the linear systems to be solved.
  int startPtr = 0;
  int numRHS2Solve = MVT::GetNumberVecs( *(problem_->getRHS()) );
  int numCurrRHS = blockSize_;

  std::vector<int> currIdx, currIdx2;

  //  The index set is generated that informs the linear problem that some linear systems are augmented.
  currIdx.resize( blockSize_ );
  currIdx2.resize( blockSize_ );
  for (int i=0; i<numCurrRHS; ++i)
    { currIdx[i] = startPtr+i; currIdx2[i]=i; }

  // Inform the linear problem of the current linear system to solve.
  problem_->setLSIndex( currIdx );

  //////////////////////////////////////////////////////////////////////////////////////
  // Parameter list
  Teuchos::ParameterList plist;
  plist.set("Block Size",blockSize_);

  // Reset the status test.
  outputTest_->reset();

  // Assume convergence is achieved, then let any failed convergence set this to false.
  bool isConverged = true;

  //////////////////////////////////////////////////////////////////////////////////////
  // TFQMR solver

  Teuchos::RCP<TFQMRIter<ScalarType,MV,OP> > tfqmr_iter =
    Teuchos::rcp( new TFQMRIter<ScalarType,MV,OP>(problem_,printer_,outputTest_,plist) );

  // Enter solve() iterations
  {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor slvtimer(*timerSolve_);
#endif

    while ( numRHS2Solve > 0 ) {
      //
      // Reset the active / converged vectors from this block
      std::vector<int> convRHSIdx;
      std::vector<int> currRHSIdx( currIdx );
      currRHSIdx.resize(numCurrRHS);

      // Reset the number of iterations.
      tfqmr_iter->resetNumIters();

      // Reset the number of calls that the status test output knows about.
      outputTest_->resetNumCalls();

      // Get the current residual for this block of linear systems.
      Teuchos::RCP<MV> R_0 = MVT::CloneViewNonConst( *(Teuchos::rcp_const_cast<MV>(problem_->getInitPrecResVec())), currIdx );

      // Set the new state and initialize the solver.
      TFQMRIterState<ScalarType,MV> newstate;
      newstate.R = R_0;
      tfqmr_iter->initializeTFQMR(newstate);

      while(1) {

        // tell tfqmr_iter to iterate
        try {
          tfqmr_iter->iterate();

          ////////////////////////////////////////////////////////////////////////////////////
          //
          // check convergence first
          //
          ////////////////////////////////////////////////////////////////////////////////////
          if ( convTest_->getStatus() == Passed ) {
            // We have convergence of the linear system.
            break;  // break from while(1){tfqmr_iter->iterate()}
          }
          ////////////////////////////////////////////////////////////////////////////////////
          //
          // check for maximum iterations
          //
          ////////////////////////////////////////////////////////////////////////////////////
          else if ( maxIterTest_->getStatus() == Passed ) {
            // we don't have convergence
            isConverged = false;
            break;  // break from while(1){tfqmr_iter->iterate()}
          }

          ////////////////////////////////////////////////////////////////////////////////////
          //
          // we returned from iterate(), but none of our status tests Passed.
          // something is wrong, and it is probably our fault.
          //
          ////////////////////////////////////////////////////////////////////////////////////

          else {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                               "Belos::TFQMRSolMgr::solve(): Invalid return from TFQMRIter::iterate().");
          }
        }
        catch (const std::exception &e) {
          printer_->stream(Errors) << "Error! Caught std::exception in TFQMRIter::iterate() at iteration "
                                   << tfqmr_iter->getNumIters() << std::endl
                                   << e.what() << std::endl;
          throw;
        }
      }

      // Update the current solution with the update computed by the iteration object.
      problem_->updateSolution( tfqmr_iter->getCurrentUpdate(), true );

      // Inform the linear problem that we are finished with this block linear system.
      problem_->setCurrLS();

      // Update indices for the linear systems to be solved.
      startPtr += numCurrRHS;
      numRHS2Solve -= numCurrRHS;
      if ( numRHS2Solve > 0 ) {
        numCurrRHS = blockSize_;

        currIdx.resize( blockSize_ );
        currIdx2.resize( blockSize_ );
        for (int i=0; i<numCurrRHS; ++i)
          { currIdx[i] = startPtr+i; currIdx2[i] = i; }
        // Set the next indices.
        problem_->setLSIndex( currIdx );

        // Set the new blocksize for the solver.
        tfqmr_iter->setBlockSize( blockSize_ );
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
  // solve.  For this solver, convTest_ may either be a single
  // (implicit) residual norm test, or a combination of two residual
  // norm tests.  In the latter case, the master convergence test
  // convTest_ is a SEQ combo of the implicit resp. explicit tests.
  // If the implicit test never passes, then the explicit test won't
  // ever be executed.  This manifests as
  // expConvTest_->getTestValue()->size() < 1.  We deal with this case
  // by using the values returned by impConvTest_->getTestValue().
  {
    // We'll fetch the vector of residual norms one way or the other.
    const std::vector<MagnitudeType>* pTestValues = NULL;
    if (expResTest_) {
      pTestValues = expConvTest_->getTestValue();
      if (pTestValues == NULL || pTestValues->size() < 1) {
        pTestValues = impConvTest_->getTestValue();
      }
    }
    else {
      // Only the implicit residual norm test is being used.
      pTestValues = impConvTest_->getTestValue();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(pTestValues == NULL, std::logic_error,
      "Belos::TFQMRSolMgr::solve(): The implicit convergence test's "
      "getTestValue() method returned NULL.  Please report this bug to the "
      "Belos developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(pTestValues->size() < 1, std::logic_error,
      "Belos::TMQMRSolMgr::solve(): The implicit convergence test's "
      "getTestValue() method returned a vector of length zero.  Please report "
      "this bug to the Belos developers.");

    // FIXME (mfh 12 Dec 2011) Does pTestValues really contain the
    // achieved tolerances for all vectors in the current solve(), or
    // just for the vectors from the last deflation?
    achievedTol_ = *std::max_element (pTestValues->begin(), pTestValues->end());
  }

  if (!isConverged) {
    return Unconverged; // return from TFQMRSolMgr::solve()
  }
  return Converged; // return from TFQMRSolMgr::solve()
}

//  This method requires the solver manager to return a std::string that describes itself.
template<class ScalarType, class MV, class OP>
std::string TFQMRSolMgr<ScalarType,MV,OP>::description() const
{
  std::ostringstream oss;
  oss << "Belos::TFQMRSolMgr<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
  oss << "{}";
  return oss.str();
}

} // end Belos namespace

#endif /* BELOS_TFQMR_SOLMGR_HPP */
