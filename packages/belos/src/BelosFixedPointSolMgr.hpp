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

#ifndef BELOS_FIXEDPOINT_SOLMGR_HPP
#define BELOS_FIXEDPOINT_SOLMGR_HPP

/*! \file BelosFixedPointSolMgr.hpp
 *  \brief The Belos::FixedPointSolMgr provides a solver manager for the FixedPoint linear solver.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "BelosCGIter.hpp"
#include "BelosFixedPointIter.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestGenResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosStatusTestOutputFactory.hpp"
#include "BelosOutputManager.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#ifdef BELOS_TEUCHOS_TIME_MONITOR
#  include "Teuchos_TimeMonitor.hpp"
#endif
#include <algorithm>

/*! \class Belos::FixedPointSolMgr
 *
 *  \brief The Belos::FixedPointSolMgr provides a powerful and fully-featured solver manager over the FixedPoint linear solver.

 \ingroup belos_solver_framework

 */

namespace Belos {

  //! @name FixedPointSolMgr Exceptions
  //@{

  /** \brief FixedPointSolMgrLinearProblemFailure is thrown when the linear problem is
   * not setup (i.e. setProblem() was not called) when solve() is called.
   *
   * This std::exception is thrown from the FixedPointSolMgr::solve() method.
   *
   */
  class FixedPointSolMgrLinearProblemFailure : public BelosError {public:
    FixedPointSolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  template<class ScalarType, class MV, class OP>
  class FixedPointSolMgr : public SolverManager<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;

  public:

    //! @name Constructors/Destructor
    //@{

    /*! \brief Empty constructor for FixedPointSolMgr.
     * This constructor takes no arguments and sets the default values for the solver.
     * The linear problem must be passed in using setProblem() before solve() is called on this object.
     * The solver values can be changed using setParameters().
     */
     FixedPointSolMgr();

    /*! \brief Basic constructor for FixedPointSolMgr.
     *
     * This constructor accepts the LinearProblem to be solved in addition
     * to a parameter list of options for the solver manager. These options include the following:
     *   - "Block Size" - an \c int specifying the block size to be used by the underlying block
     *                    conjugate-gradient solver. Default: 1
     *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Belos::Errors
     *   - "Output Style" - a OutputType specifying the style of output. Default: Belos::General
     *   - "Output Stream" - a reference-counted pointer to the output stream where all
     *                       solver output is sent.  Default: Teuchos::rcp(&std::cout,false)
     *   - "Output Frequency" - an \c int specifying how often convergence information should be
     *                          outputted.  Default: -1 (never)
     *   - "Show Maximum Residual Norm Only" - a \c bool specifying whether that only the maximum
     *                                         relative residual norm is printed if convergence
     *                                         information is printed. Default: false
     *   - "Timer Label" - a \c std::string to use as a prefix for the timer labels.  Default: "Belos"
     */
    FixedPointSolMgr( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                   const Teuchos::RCP<Teuchos::ParameterList> &pl );

    //! Destructor.
    virtual ~FixedPointSolMgr() {};

    //! clone for Inverted Injection (DII)
    Teuchos::RCP<SolverManager<ScalarType, MV, OP> > clone () const override {
      return Teuchos::rcp(new FixedPointSolMgr<ScalarType,MV,OP>);
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

    /*! \brief Return whether a loss of accuracy was detected by this solver during the most current solve.
     */
    bool isLOADetected() const override { return false; }
    //@}

    //! @name Set methods
    //@{
   
    //! Set the linear problem that needs to be solved. 
    void setProblem( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem ) override { problem_ = problem; }
   
    //! Set the parameters the solver manager should use to solve the linear problem. 
    void setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params ) override;
    
    //! Set user-defined convergence status test.
    void replaceUserConvStatusTest( const Teuchos::RCP<StatusTestResNorm<ScalarType,MV,OP> > &userConvStatusTest )
    {

      convTest_ = userConvStatusTest;

      typedef Belos::StatusTestCombo<ScalarType,MV,OP>  StatusTestCombo_t;
      sTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::OR, maxIterTest_, convTest_ ) );

      StatusTestOutputFactory<ScalarType,MV,OP> stoFactory( outputStyle_ );
      outputTest_ = stoFactory.create( printer_, sTest_, outputFreq_, Passed+Failed+Undefined );

      std::string solverDesc = " Fixed Point ";
      outputTest_->setSolverDesc( solverDesc );
    }

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
     * This method calls FixedPointIter::iterate() or CGIter::iterate(), which will return either because a
     * specially constructed status test evaluates to ::Passed or an std::exception is thrown.
     *
     * A return from FixedPointIter::iterate() signifies one of the following scenarios:
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

    /** \brief Method to return description of the block CG solver manager */
    std::string description() const override;
    //@}

  private:

    //! The linear problem to solve.
    Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > problem_;

    //! Output manager, that handles printing of different kinds of messages.
    Teuchos::RCP<OutputManager<ScalarType> > printer_;
    //! Output stream to which the output manager prints.
    Teuchos::RCP<std::ostream> outputStream_;

    /// \brief Aggregate stopping criterion.
    ///
    /// This is an OR combination of the maximum iteration count test
    /// (\c maxIterTest_) and convergence test (\c convTest_).
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > sTest_;

    //! Maximum iteration count stopping criterion.
    Teuchos::RCP<StatusTestMaxIters<ScalarType,MV,OP> > maxIterTest_;

    //! Convergence stopping criterion.
    Teuchos::RCP<StatusTestResNorm<ScalarType,MV,OP> > convTest_;

    //! Output "status test" that controls all the other status tests.
    Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputTest_;

    //! Current parameter list.
    Teuchos::RCP<Teuchos::ParameterList> params_;

    //
    // Default solver parameters.
    //
    static constexpr int maxIters_default_ = 1000;
    static constexpr bool showMaxResNormOnly_default_ = false;
    static constexpr int blockSize_default_ = 1;
    static constexpr int verbosity_default_ = Belos::Errors;
    static constexpr int outputStyle_default_ = Belos::General;
    static constexpr int outputFreq_default_ = -1;
    static constexpr const char * label_default_ = "Belos";
    static constexpr std::ostream * outputStream_default_ = &std::cout;

    //
    // Current solver parameters and other values.
    //

    //! Convergence tolerance (read from parameter list).
    MagnitudeType convtol_;

    /// \brief Tolerance achieved by the last \c solve() invocation.
    ///
    /// This is the maximum over all right-hand sides' achieved
    /// convergence tolerances, and is set whether or not the solve
    /// actually managed to achieve the desired convergence tolerance.
    MagnitudeType achievedTol_;

    //! Maximum iteration count (read from parameter list).
    int maxIters_;

    //! Number of iterations taken by the last \c solve() invocation.
    int numIters_;

    int blockSize_, verbosity_, outputStyle_, outputFreq_;
    bool showMaxResNormOnly_;

    //! Prefix label for all the timers.
    std::string label_;

    //! Solve timer.
    Teuchos::RCP<Teuchos::Time> timerSolve_;

    //! Whether or not the parameters have been set (via \c setParameters()).
    bool isSet_;
  };


// Empty Constructor
template<class ScalarType, class MV, class OP>
FixedPointSolMgr<ScalarType,MV,OP>::FixedPointSolMgr() :
  outputStream_(Teuchos::rcp(outputStream_default_,false)),
  convtol_(DefaultSolverParameters::convTol),
  achievedTol_(Teuchos::ScalarTraits<MagnitudeType>::zero()),
  maxIters_(maxIters_default_),
  numIters_(0),
  blockSize_(blockSize_default_),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  showMaxResNormOnly_(showMaxResNormOnly_default_),
  label_(label_default_),
  isSet_(false)
{}


// Basic Constructor
template<class ScalarType, class MV, class OP>
FixedPointSolMgr<ScalarType,MV,OP>::
FixedPointSolMgr(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
              const Teuchos::RCP<Teuchos::ParameterList> &pl) :
  problem_(problem),
  outputStream_(Teuchos::rcp(outputStream_default_,false)),
  convtol_(DefaultSolverParameters::convTol),
  achievedTol_(Teuchos::ScalarTraits<MagnitudeType>::zero()),
  maxIters_(maxIters_default_),
  numIters_(0),
  blockSize_(blockSize_default_),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  showMaxResNormOnly_(showMaxResNormOnly_default_),
  label_(label_default_),
  isSet_(false)
{
  TEUCHOS_TEST_FOR_EXCEPTION(problem_.is_null(), std::invalid_argument,
    "FixedPointSolMgr's constructor requires a nonnull LinearProblem instance.");

  // If the user passed in a nonnull parameter list, set parameters.
  // Otherwise, the next solve() call will use default parameters,
  // unless the user calls setParameters() first.
  if (! pl.is_null()) {
    setParameters (pl);
  }
}

template<class ScalarType, class MV, class OP>
void
FixedPointSolMgr<ScalarType,MV,OP>::
setParameters (const Teuchos::RCP<Teuchos::ParameterList> &params)
{
  // Create the internal parameter list if one doesn't already exist.
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
    blockSize_ = params->get("Block Size",blockSize_default_);
    TEUCHOS_TEST_FOR_EXCEPTION(blockSize_ <= 0, std::invalid_argument,
                       "Belos::FixedPointSolMgr: \"Block Size\" must be strictly positive.");

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
      std::string solveLabel = label_ + ": FixedPointSolMgr total solve time";
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

    // Update parameter in our list.
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

  // Create status tests if we need to.

  // Basic test checks maximum iterations and native residual.
  if (maxIterTest_ == Teuchos::null)
    maxIterTest_ = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>( maxIters_ ) );

  // Implicit residual test, using the native residual to determine if convergence was achieved.
  if (convTest_ == Teuchos::null)
    convTest_ = Teuchos::rcp( new StatusTestResNorm_t( convtol_, 1 ) );

  if (sTest_ == Teuchos::null)
    sTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::OR, maxIterTest_, convTest_ ) );

  if (outputTest_ == Teuchos::null) {

    // Create the status test output class.
    // This class manages and formats the output from the status test.
    StatusTestOutputFactory<ScalarType,MV,OP> stoFactory( outputStyle_ );
    outputTest_ = stoFactory.create( printer_, sTest_, outputFreq_, Passed+Failed+Undefined );

    // Set the solver string for the output test
    std::string solverDesc = " Fixed Point ";
    outputTest_->setSolverDesc( solverDesc );

  }

  // Create the timer if we need to.
  if (timerSolve_ == Teuchos::null) {
    std::string solveLabel = label_ + ": FixedPointSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    timerSolve_ = Teuchos::TimeMonitor::getNewCounter(solveLabel);
#endif
  }

  // Inform the solver manager that the current parameters were set.
  isSet_ = true;
}


template<class ScalarType, class MV, class OP>
Teuchos::RCP<const Teuchos::ParameterList>
FixedPointSolMgr<ScalarType,MV,OP>::getValidParameters() const
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
    pl->set("Maximum Iterations", static_cast<int>(maxIters_default_),
      "The maximum number of block iterations allowed for each\n"
      "set of RHS solved.");
    pl->set("Block Size", static_cast<int>(blockSize_default_),
      "The number of vectors in each block.");
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
    pl->set("Show Maximum Residual Norm Only", static_cast<bool>(showMaxResNormOnly_default_),
      "When convergence information is printed, only show the maximum\n"
      "relative residual norm when the block size is greater than one.");
    pl->set("Timer Label", static_cast<const char *>(label_default_),
      "The string to use as a prefix for the timer labels.");
    validPL = pl;
  }
  return validPL;
}


// solve()
template<class ScalarType, class MV, class OP>
ReturnType FixedPointSolMgr<ScalarType,MV,OP>::solve() {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;

  // Set the current parameters if they were not set before.  NOTE:
  // This may occur if the user generated the solver manager with the
  // default constructor and then didn't set any parameters using
  // setParameters().
  if (!isSet_) {
    setParameters(Teuchos::parameterList(*getValidParameters()));
  }

  Teuchos::BLAS<int,ScalarType> blas;
  Teuchos::LAPACK<int,ScalarType> lapack;

  TEUCHOS_TEST_FOR_EXCEPTION( !problem_->isProblemSet(),
    FixedPointSolMgrLinearProblemFailure,
    "Belos::FixedPointSolMgr::solve(): Linear problem is not ready, setProblem() "
    "has not been called.");

  // Create indices for the linear systems to be solved.
  int startPtr = 0;
  int numRHS2Solve = MVT::GetNumberVecs( *(problem_->getRHS()) );
  int numCurrRHS = ( numRHS2Solve < blockSize_) ? numRHS2Solve : blockSize_;

  std::vector<int> currIdx, currIdx2;
  currIdx.resize( blockSize_ );
  currIdx2.resize( blockSize_ );
  for (int i=0; i<numCurrRHS; ++i)
    { currIdx[i] = startPtr+i; currIdx2[i]=i; }
  for (int i=numCurrRHS; i<blockSize_; ++i)
    { currIdx[i] = -1; currIdx2[i] = i; }

  // Inform the linear problem of the current linear system to solve.
  problem_->setLSIndex( currIdx );

  ////////////////////////////////////////////////////////////////////////////
  // Set up the parameter list for the Iteration subclass.
  Teuchos::ParameterList plist;
  plist.set("Block Size",blockSize_);

  // Reset the output status test (controls all the other status tests).
  outputTest_->reset();

  // Assume convergence is achieved, then let any failed convergence
  // set this to false.  "Innocent until proven guilty."
  bool isConverged = true;

  ////////////////////////////////////////////////////////////////////////////
  // Set up the FixedPoint Iteration subclass.

  RCP<FixedPointIteration<ScalarType,MV,OP> > block_fp_iter;
  block_fp_iter = rcp (new FixedPointIter<ScalarType,MV,OP> (problem_, printer_, outputTest_, plist));

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
      block_fp_iter->resetNumIters();

      // Reset the number of calls that the status test output knows about.
      outputTest_->resetNumCalls();

      // Get the current residual for this block of linear systems.
      RCP<MV> R_0 = MVT::CloneViewNonConst( *(rcp_const_cast<MV>(problem_->getInitResVec())), currIdx );

      // Set the new state and initialize the solver.
      FixedPointIterationState<ScalarType,MV> newstate;
      newstate.R = R_0;
      block_fp_iter->initializeFixedPoint(newstate);

      while(1) {

        // tell block_fp_iter to iterate
        try {
          block_fp_iter->iterate();
          //
          // Check whether any of the linear systems converged.
          //
          if (convTest_->getStatus() == Passed) {
            // At least one of the linear system(s) converged.
            //
            // Get the column indices of the linear systems that converged.
            std::vector<int> convIdx = convTest_->convIndices();

            // If the number of converged linear systems equals the
            // number of linear systems currently being solved, then
            // we are done with this block.
            if (convIdx.size() == currRHSIdx.size())
              break;  // break from while(1){block_fp_iter->iterate()}

            // Inform the linear problem that we are finished with
            // this current linear system.
            problem_->setCurrLS();

            // Reset currRHSIdx to contain the right-hand sides that
            // are left to converge for this block.
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
              else {
              }
            }
            currRHSIdx.resize(have);
            currIdx2.resize(have);

            // Set the remaining indices after deflation.
            problem_->setLSIndex( currRHSIdx );

            // Get the current residual vector.
            std::vector<MagnitudeType> norms;
            R_0 = MVT::CloneCopy( *(block_fp_iter->getNativeResiduals(&norms)),currIdx2 );
            for (int i=0; i<have; ++i) { currIdx2[i] = i; }

            // Set the new blocksize for the solver.
            block_fp_iter->setBlockSize( have );

            // Set the new state and initialize the solver.
            FixedPointIterationState<ScalarType,MV> defstate;
            defstate.R = R_0;
            block_fp_iter->initializeFixedPoint(defstate);
          }
          //
          // None of the linear systems converged.  Check whether the
          // maximum iteration count was reached.
          //
          else if (maxIterTest_->getStatus() == Passed) {
            isConverged = false; // None of the linear systems converged.
            break;  // break from while(1){block_fp_iter->iterate()}
          }
          //
          // iterate() returned, but none of our status tests Passed.
          // This indicates a bug.
          //
          else {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
              "Belos::FixedPointSolMgr::solve(): Neither the convergence test nor "
              "the maximum iteration count test passed.  Please report this bug "
              "to the Belos developers.");
          }
        }
        catch (const std::exception &e) {
          std::ostream& err = printer_->stream (Errors);
          err << "Error! Caught std::exception in FixedPointIteration::iterate() at "
              << "iteration " << block_fp_iter->getNumIters() << std::endl
              << e.what() << std::endl;
          throw;
        }
      }

      // Inform the linear problem that we are finished with this
      // block linear system.
      problem_->setCurrLS();

      // Update indices for the linear systems to be solved.
      startPtr += numCurrRHS;
      numRHS2Solve -= numCurrRHS;
      if ( numRHS2Solve > 0 ) {
        numCurrRHS = ( numRHS2Solve < blockSize_) ? numRHS2Solve : blockSize_;


        currIdx.resize( blockSize_ );
        currIdx2.resize( blockSize_ );
        for (int i=0; i<numCurrRHS; ++i)
          { currIdx[i] = startPtr+i; currIdx2[i] = i; }
        for (int i=numCurrRHS; i<blockSize_; ++i)
          { currIdx[i] = -1; currIdx2[i] = i; }

        // Set the next indices.
        problem_->setLSIndex( currIdx );

        // Set the new blocksize for the solver.
        block_fp_iter->setBlockSize( blockSize_ );
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
  // Calling summarize() requires communication in general, so don't
  // call it unless the user wants to print out timing details.
  // summarize() will do all the work even if it's passed a "black
  // hole" output stream.
  if (verbosity_ & TimingDetails) {
    Teuchos::TimeMonitor::summarize( printer_->stream(TimingDetails) );
  }
#endif

  // Save the iteration count for this solve.
  numIters_ = maxIterTest_->getNumIters();

  // Save the convergence test value ("achieved tolerance") for this solve.
  {
    // testValues is nonnull and not persistent.
    const std::vector<MagnitudeType>* pTestValues = convTest_->getTestValue();

    TEUCHOS_TEST_FOR_EXCEPTION(pTestValues == NULL, std::logic_error,
      "Belos::FixedPointSolMgr::solve(): The convergence test's getTestValue() "
      "method returned NULL.  Please report this bug to the Belos developers.");

    TEUCHOS_TEST_FOR_EXCEPTION(pTestValues->size() < 1, std::logic_error,
      "Belos::FixedPointSolMgr::solve(): The convergence test's getTestValue() "
      "method returned a vector of length zero.  Please report this bug to the "
      "Belos developers.");

    // FIXME (mfh 12 Dec 2011) Does pTestValues really contain the
    // achieved tolerances for all vectors in the current solve(), or
    // just for the vectors from the last deflation?
    achievedTol_ = *std::max_element (pTestValues->begin(), pTestValues->end());
  }

  if (!isConverged) {
    return Unconverged; // return from FixedPointSolMgr::solve()
  }
  return Converged; // return from FixedPointSolMgr::solve()
}

//  This method requires the solver manager to return a std::string that describes itself.
template<class ScalarType, class MV, class OP>
std::string FixedPointSolMgr<ScalarType,MV,OP>::description() const
{
  std::ostringstream oss;
  oss << "Belos::FixedPointSolMgr<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
  return oss.str();
}

} // end Belos namespace

#endif /* BELOS_FIXEDPOINT_SOLMGR_HPP */
