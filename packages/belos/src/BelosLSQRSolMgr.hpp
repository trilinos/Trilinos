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

/*! \file BelosLSQRSolMgr.hpp
 *  \brief The LSQRSolMgr provides a solver manager for the LSQR linear solver.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "BelosLSQRIteration.hpp"
#include "BelosLSQRIter.hpp"
#include "BelosDGKSOrthoManager.hpp"
#include "BelosICGSOrthoManager.hpp"
#include "BelosIMGSOrthoManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosLSQRStatusTest.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosStatusTestOutputFactory.hpp"
#include "BelosOutputManager.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"

#ifdef BELOS_TEUCHOS_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

/*! \class Belos::LSQRSolMgr
 *
 *  \brief The LSQRSolMgr managers the LSQR linear least squares solver.

 \ingroup belos_solver_framework

 \author David Day

 */

namespace Belos {

  
//! @name LSQRSolMgr Exceptions
//@{

/** \brief Belos::LSQRSolMgrLinearProblemFailure is thrown when the linear problem is
 * not setup (i.e. setProblem() was not called) when solve() is called.
 *
 * This std::exception is thrown from the LSQRSolMgr::solve() method.
 *
 */
class LSQRSolMgrLinearProblemFailure : public BelosError {public:
    LSQRSolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg)
  {}};

/** \brief LSQRSolMgrOrthoFailure is thrown when the orthogonalization manager is
 * unable to generate orthonormal columns from the initial basis vectors.
 *
 * This std::exception is thrown from the LSQRSolMgr::solve() method.
 *
 */
class LSQRSolMgrOrthoFailure : public BelosError {public:
    LSQRSolMgrOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
  {}};

/** \brief LSQRSolMgrBlockSizeFailure is thrown when the linear problem has
 * more than one RHS.  This is unique to single vector methods.
 *
 * This std::exception is thrown from the LSQRSolMgr::solve() method.
 *
 */
class LSQRSolMgrBlockSizeFailure : public BelosError {public:
    LSQRSolMgrBlockSizeFailure(const std::string& what_arg) : BelosError(what_arg)
  {}};

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
   *   - "Maximum Iterations" - an \c int specifying the maximum number of iterations the 
   *                            underlying solver is allowed to perform. Default: 1000
   *   - "Condition Limit" - a \c MagnitudeType specifying the upper limit of the estimate of
   *                         the norm of Abar to decide convergence. Default: 0.
   *   - "Term Iter Max" - the number of consecutive successful iterations required before
   *                       convergence is declared.
   *   - "Rel RHS Err" - an estimate of the error in the data defining the RHS.
   *   - "Rel Mat Err" - an estimate of the error in the data defining the matrix.
   *   - "Orthogonalization" - a \c std::string specifying the desired orthogonalization:  
   *                           DGKS ,ICGS, and IMGS. Default: "DGKS"
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
    return tuple(timerSolve_);
  }
  
  //! Get the iteration count for the most recent call to \c solve().
  int getNumIters() const {
    return numIters_;
  }

  //! Get the for the most recent call to \c solve().
  MagnitudeType getMatCondNum () const {
    return matCondNum_; 
  }

  //! Get the for the most recent call to \c solve().
  MagnitudeType getMatNorm () const {
    return matNorm_;
  }

  //! Get the for the most recent call to \c solve().
  MagnitudeType getResNorm () const {
    return resNorm_;
  }

  //! Get the for the most recent call to \c solve().
  MagnitudeType getMatResNorm () const {
    return matResNorm_;
  }


  
  /*! \brief Return whether a loss of accuracy was detected by this solver during the most current solve.
   */
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



  // Method to convert std::string to enumerated type for residual.
  Belos::ScaleType convertStringToScaleType( std::string& scaleType ) {
    if (scaleType == "Norm of Initial Residual") {
      return Belos::NormOfInitRes;
    } else if (scaleType == "Norm of Preconditioned Initial Residual") {
      return Belos::NormOfPrecInitRes;
    } else if (scaleType == "Norm of RHS") {
      return Belos::NormOfRHS;
    } else if (scaleType == "None") {
      return Belos::None;
    } else
      TEST_FOR_EXCEPTION( true ,std::logic_error,
        "Belos::LSQRSolMgr(): Invalid residual scaling type.");
  }


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

    
  // Default solver values.
  static const MagnitudeType lambda_default_;
  static const MagnitudeType relRhsErr_default_;
  static const MagnitudeType relMatErr_default_;
  static const MagnitudeType condMax_default_;
  static const int maxIters_default_; 
  static const int termIterMax_default_;
  static const std::string orthoType_default_;
  static const MagnitudeType orthoKappa_default_;
  static const int verbosity_default_;
  static const int outputStyle_default_;
  static const int outputFreq_default_;
  static const std::string label_default_;
  static const Teuchos::RCP<std::ostream> outputStream_default_;
  
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
  bool isSTSet_;
  bool loaDetected_;

};


// Default solver values.
template<class ScalarType, class MV, class OP>
const typename LSQRSolMgr<ScalarType,MV,OP>::MagnitudeType LSQRSolMgr<ScalarType,MV,OP>::lambda_default_ = 
  Teuchos::ScalarTraits<MagnitudeType>::zero();

template<class ScalarType, class MV, class OP>
const typename LSQRSolMgr<ScalarType,MV,OP>::MagnitudeType LSQRSolMgr<ScalarType,MV,OP>::relRhsErr_default_ = 
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType(10) * 
  Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>::squareroot (Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>::eps());

template<class ScalarType, class MV, class OP>
const typename LSQRSolMgr<ScalarType,MV,OP>::MagnitudeType LSQRSolMgr<ScalarType,MV,OP>::relMatErr_default_ = 
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType(10) * 
  Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>::squareroot (Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>::eps());

template<class ScalarType, class MV, class OP>
const typename LSQRSolMgr<ScalarType,MV,OP>::MagnitudeType LSQRSolMgr<ScalarType,MV,OP>::condMax_default_ = 
 Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>::one() /
 Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>::eps() ;

template<class ScalarType, class MV, class OP>
const int LSQRSolMgr<ScalarType,MV,OP>::maxIters_default_ = 1000;

template<class ScalarType, class MV, class OP>
const int LSQRSolMgr<ScalarType,MV,OP>::termIterMax_default_ = 1;

template<class ScalarType, class MV, class OP>
const std::string LSQRSolMgr<ScalarType,MV,OP>::orthoType_default_ = "DGKS";

template<class ScalarType, class MV, class OP>
const typename LSQRSolMgr<ScalarType,MV,OP>::MagnitudeType LSQRSolMgr<ScalarType,MV,OP>::orthoKappa_default_ = -1.0;

template<class ScalarType, class MV, class OP>
const int LSQRSolMgr<ScalarType,MV,OP>::verbosity_default_ = Belos::Errors;

template<class ScalarType, class MV, class OP>
const int LSQRSolMgr<ScalarType,MV,OP>::outputStyle_default_ = Belos::General;

template<class ScalarType, class MV, class OP>
const int LSQRSolMgr<ScalarType,MV,OP>::outputFreq_default_ = -1;

template<class ScalarType, class MV, class OP>
const std::string LSQRSolMgr<ScalarType,MV,OP>::label_default_ = "Belos";

template<class ScalarType, class MV, class OP>
const Teuchos::RCP<std::ostream> LSQRSolMgr<ScalarType,MV,OP>::outputStream_default_ = Teuchos::rcp(&std::cout,false);


// Empty Constructor
template<class ScalarType, class MV, class OP>
LSQRSolMgr<ScalarType,MV,OP>::LSQRSolMgr() :
  outputStream_(outputStream_default_),
  lambda_(lambda_default_),
  relRhsErr_(relRhsErr_default_),
  relMatErr_(relMatErr_default_),
  condMax_(condMax_default_),
  maxIters_(maxIters_default_),
  termIterMax_(termIterMax_default_),
  orthoType_(orthoType_default_),
  orthoKappa_(orthoKappa_default_),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  label_(label_default_),
  isSet_(false),
  isSTSet_(false),
  loaDetected_(false)
{}


// Basic Constructor
template<class ScalarType, class MV, class OP>
LSQRSolMgr<ScalarType,MV,OP>::LSQRSolMgr( 
  const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
  const Teuchos::RCP<Teuchos::ParameterList> &pl ) : 
  problem_(problem),
  outputStream_(outputStream_default_),
  lambda_(lambda_default_),
  relRhsErr_(relRhsErr_default_),
  relMatErr_(relMatErr_default_),
  condMax_(condMax_default_),
  maxIters_(maxIters_default_),
  termIterMax_(termIterMax_default_),
  orthoType_(orthoType_default_),
  orthoKappa_(orthoKappa_default_),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  label_(label_default_),
  isSet_(false),
  isSTSet_(false),
  loaDetected_(false)
{
  TEST_FOR_EXCEPTION(problem_ == Teuchos::null, std::invalid_argument, "Problem not given to solver manager.");

  // If the parameter list pointer is null, 
  // then set the current parameters to the default parameter list.
  if ( !is_null(pl) ) {
    setParameters( pl );  
  }
}


template<class ScalarType, class MV, class OP>
Teuchos::RCP<const Teuchos::ParameterList> 
LSQRSolMgr<ScalarType,MV,OP>::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;
  
  // Set all the valid parameters and their default values.
  if(is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("Output Stream", outputStream_default_,
      "is a reference-counted pointer to the output stream receiving\n"
      "all solver output.");
    pl->set("Lambda", lambda_default_, "is the damping parameter.");
    pl->set("Rel RHS Err", relRhsErr_default_,
	    "estimates the error in the data defining the right-\n"
	    "hand side.");
    pl->set("Rel Mat Err", relMatErr_default_,
	    "estimates the error in the data defining the matrix.");
    pl->set("Condition limit", condMax_default_,
      "bounds the estimated condition number of Abar.");
    pl->set("Maximum Iterations", maxIters_default_,
      "allows at most the maximum number of iterations.");
    pl->set("Term Iter Max", termIterMax_default_,
      "consecutive iterations meeting thresholds are necessary for\n"
       "for convergence.");
    pl->set("Orthogonalization", orthoType_default_,
      "uses orthogonalization of either DGKS, ICGS or IMGS.");
    pl->set("Orthogonalization Constant",orthoKappa_default_,
      "is the threshold used by DGKS orthogonalization to determine\n"
      "whether or not to repeat classical Gram-Schmidt.");
    pl->set("Verbosity", verbosity_default_,
      "type(s) of solver information are outputted to the output\n"
      "stream.");
    pl->set("Output Style", outputStyle_default_,
      "the style used for the solver information outputted to the\n"
      "output stream.");
    pl->set("Output Frequency", outputFreq_default_,
      "is the frequency at which information is written to the\n"
      "output stream.");  
    pl->set("Timer Label", label_default_,
      "is the string to use as a prefix for the timer labels.");
    //  pl->set("Restart Timers", restartTimers_);
    validPL = pl;
  }
  return validPL;
}


template<class ScalarType, class MV, class OP>
void LSQRSolMgr<ScalarType,MV,OP>::setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params )
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
  using Teuchos::Exceptions::InvalidParameter;
  using Teuchos::Exceptions::InvalidParameterName;
  using Teuchos::Exceptions::InvalidParameterType;

  // Create the internal parameter list if ones doesn't already exist.
  if (params_ == Teuchos::null) {
    params_ = Teuchos::rcp( new Teuchos::ParameterList(*getValidParameters()) );
  }
  else {
    params->validateParameters(*getValidParameters());
  }

  // Check for damping value lambda
  if (params->isParameter("Lambda")) {
    lambda_ = params->get("Lambda",lambda_default_);

    // Update parameter in our list and in status test.
    params_->set("Lambda",lambda_);
  }

  // Check for maximum number of iterations
  if (params->isParameter("Maximum Iterations")) {
    maxIters_ = params->get("Maximum Iterations",maxIters_default_);

    // Update parameter in our list and in status test.
    params_->set("Maximum Iterations", maxIters_);
    if (maxIterTest_!=Teuchos::null)
      maxIterTest_->setMaxIters( maxIters_ );
  }

  // Check to see if the timer label changed.
  if (params->isParameter("Timer Label")) {
    std::string tempLabel = params->get("Timer Label", label_default_);

    // Update parameter in our list and solver timer
    if (tempLabel != label_) {
      label_ = tempLabel;
      params_->set("Timer Label", label_);
      std::string solveLabel = label_ + ": LSQRSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerSolve_ = Teuchos::TimeMonitor::getNewTimer(solveLabel);
#endif
    }
  }
  // Check for a change in verbosity level
  if (params->isParameter ("Verbosity")) {
    if (isParameterType<int> (*params, "Verbosity")) {   //'isParameterType' was not declared in this scope
      verbosity_ = params->get ("Verbosity", verbosity_default_);
    } else {
      verbosity_ = (int) getParameter<Belos::MsgType> (*params, "Verbosity"); 
    } //'getParameter' was not declared in this scope
    // Update parameter in our list.
    params_->set ("Verbosity", verbosity_);
    // If the output manager (printer_) is null, then we will
    // instantiate it later with the correct verbosity.
    if (! printer_.is_null())
      printer_->setVerbosity (verbosity_);
  }


  // Check for a change in output style
  if (params->isParameter ("Output Style")) {
    if (isParameterType<int> (*params, "Output Style")) {
      outputStyle_ = params->get ("Output Style", outputStyle_default_);
    } else {
      outputStyle_ = (int) getParameter<OutputType> (*params, "Output Style");
    }

    // Update parameter in our list.
    params_->set ("Output Style", outputStyle_);
    // We will (re)instantiate the output status test afresh below.
    outputTest_ = null;
  }

  // Get the output stream for the output manager.
  //
  // FIXME (mfh 28 Feb 2011) While storing the output stream in the
  // parameter list (either as an RCP or as a nonconst reference) is
  // convenient and safe for programming, it makes it nearly
  // impossible to serialize the parameter list, read it back in from
  // the serialized representation, and get the same output stream as
  // before.  However, a general solution is likely impossible,
  // because output streams may be arbitrary constructed objects.  
  //
  // In case the output stream can't be read back in, we default to
  // stdout (std::cout), just to ensure reasonable behavior.
  if (params->isParameter ("Output Stream")) {
    try {
      outputStream_ = getParameter<RCP<std::ostream> > (*params, "Output Stream");
    } catch (InvalidParameter&) {
      outputStream_ = rcpFromRef (std::cout);
    }
    // We assume that a null output stream indicates that the user
    // doesn't want to print anything, so we replace it with a "black
    // hole" stream that prints nothing sent to it.  (We can't use a
    // null output stream, since the output manager always sends
    // things it wants to print to the output stream.)
    if (outputStream_.is_null())
      outputStream_ = rcp (new Teuchos::oblackholestream);

    // Update parameter in our list.
    params_->set ("Output Stream", outputStream_);
    // If the output manager (printer_) is null, then we will
    // instantiate it later with the correct output stream.
    if (! printer_.is_null())
      printer_->setOStream (outputStream_);
  }

  // frequency level
  if (verbosity_ & Belos::StatusTestDetails) {
    if (params->isParameter ("Output Frequency")) {
      outputFreq_ = params->get ("Output Frequency", outputFreq_default_);
    }

    // Update parameter in out list and output status test.
    params_->set("Output Frequency", outputFreq_);
    if (! outputTest_.is_null())
      outputTest_->setOutputFrequency (outputFreq_);
  }

  // Create output manager if we need to, using the verbosity level
  // and output stream that we fetched above.  We do this here because
  // instantiating an OrthoManager using OrthoManagerFactory requires
  // a valid OutputManager.
  if (printer_.is_null()) {
    printer_ = rcp (new OutputManager<ScalarType> (verbosity_, outputStream_));
  }

  // Check if the orthogonalization changed.
  if (params->isParameter("Orthogonalization")) {
    std::string tempOrthoType = params->get("Orthogonalization",orthoType_default_);
    TEST_FOR_EXCEPTION( tempOrthoType != "DGKS" && tempOrthoType != "ICGS" && tempOrthoType != "IMGS", 
                        std::invalid_argument,
			"LSQRSolMgr: \"Orthogonalization\" must be either \"DGKS\", \"ICGS\", or \"IMGS\".");
    if (tempOrthoType != orthoType_) {
      orthoType_ = tempOrthoType;
      // Create orthogonalization manager
      if (orthoType_=="DGKS") {
	if (orthoKappa_ <= 0) {
	  ortho_ = Teuchos::rcp( new Belos::DGKSOrthoManager<ScalarType,MV,OP>( label_ ) );
	}
	else {
	  ortho_ = Teuchos::rcp( new Belos::DGKSOrthoManager<ScalarType,MV,OP>( label_ ) );
	  Teuchos::rcp_dynamic_cast<Belos::DGKSOrthoManager<ScalarType,MV,OP> >(ortho_)->setDepTol( orthoKappa_ );
	}
      }
      else if (orthoType_=="ICGS") {
	ortho_ = Teuchos::rcp( new Belos::ICGSOrthoManager<ScalarType,MV,OP>( label_ ) );
      } 
      else if (orthoType_=="IMGS") {
	ortho_ = Teuchos::rcp( new Belos::IMGSOrthoManager<ScalarType,MV,OP>( label_ ) );
      } 
    }  
  }

  // Check which orthogonalization constant to use.
  if (params->isParameter("Orthogonalization Constant")) {
    orthoKappa_ = params->get("Orthogonalization Constant",orthoKappa_default_);

    // Update parameter in our list.
    params_->set("Orthogonalization Constant",orthoKappa_);
    if (orthoType_=="DGKS") {
      if (orthoKappa_ > 0 && ortho_ != Teuchos::null) {
	Teuchos::rcp_dynamic_cast<Belos::DGKSOrthoManager<ScalarType,MV,OP> >(ortho_)->setDepTol( orthoKappa_ );
      }
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
    printer_ = Teuchos::rcp( new Belos::OutputManager<ScalarType>(verbosity_, outputStream_) );
  }  
  
  // Check for convergence tolerance
  if (params->isParameter("Condition Limit")) {
    condMax_ = params->get("Condition Limit",condMax_default_);

    // Update parameter in our list and residual tests.
    params_->set("Condition Limit", condMax_);
    if (convTest_ != Teuchos::null)
      convTest_->setCondLim( condMax_ );
  }

  // Check for number of consecutive passed iterations
  if (params->isParameter("Term Iter Max")) {
    termIterMax_ = params->get("Term Iter Max", termIterMax_default_);

    // Update parameter in our list and residual tests
    params_->set("Term Iter Max", termIterMax_);
    if (convTest_ != Teuchos::null)
      convTest_->setTermIterMax( termIterMax_ );
  }

  // Check for relative RHS error
  if (params->isParameter("Rel RHS Err")) {
    relRhsErr_ = params->get("Rel RHS Err",relRhsErr_default_);

    // Update parameter in our list and residual tests
    params_->set("Rel RHS Err", relRhsErr_);
    if (convTest_ != Teuchos::null)
      convTest_->setRelRhsErr( relRhsErr_ );
  }

  // Check for relative matrix error
  if (params->isParameter("Rel Mat Err")) {
    relMatErr_ = params->get("Rel Mat Err",relMatErr_default_);

    // Update parameter in our list and residual tests
    params_->set("Rel Mat Err", relMatErr_);
    if (convTest_ != Teuchos::null)
      convTest_->setRelMatErr( relMatErr_ );
  }

  // Maximum Iterations and native residual test
  if (maxIterTest_ == Teuchos::null)
    maxIterTest_ = Teuchos::rcp( new Belos::StatusTestMaxIters<ScalarType,MV,OP>( maxIters_ ) );

  if (convTest_ == Teuchos::null)
    convTest_ = Teuchos::rcp( new LSQRStatusTest<ScalarType,MV,OP>(condMax_, termIterMax_, relRhsErr_, relMatErr_) );
  
  typedef Belos::StatusTestCombo<ScalarType,MV,OP>  StatusTestCombo_t;

  if (sTest_ == Teuchos::null)
    sTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::OR, maxIterTest_, convTest_ ) );
  
  if (outputTest_ == Teuchos::null) {

    // Create the status test output class.
    // This class manages and formats the output from the status test.
    Belos::StatusTestOutputFactory<ScalarType,MV,OP> stoFactory( outputStyle_ );
    outputTest_ = stoFactory.create( printer_, sTest_, outputFreq_, Belos::Passed+Belos::Failed+Belos::Undefined );

    // Set the solver string for the output test
    std::string solverDesc = " LSQR ";
    outputTest_->setSolverDesc( solverDesc );

  }

  // Create orthogonalization manager if we need to.
  if (ortho_ == Teuchos::null) {
    if (orthoType_=="DGKS") {
      if (orthoKappa_ <= 0) {
	ortho_ = Teuchos::rcp( new Belos::DGKSOrthoManager<ScalarType,MV,OP>( label_ ) );
      }
      else {
	ortho_ = Teuchos::rcp( new Belos::DGKSOrthoManager<ScalarType,MV,OP>( label_ ) );
	Teuchos::rcp_dynamic_cast<Belos::DGKSOrthoManager<ScalarType,MV,OP> >(ortho_)->setDepTol( orthoKappa_ );
      }
    }
    else if (orthoType_=="ICGS") {
      ortho_ = Teuchos::rcp( new Belos::ICGSOrthoManager<ScalarType,MV,OP>( label_ ) );
    } 
    else if (orthoType_=="IMGS") {
      ortho_ = Teuchos::rcp( new Belos::IMGSOrthoManager<ScalarType,MV,OP>( label_ ) );
    } 
    else {
      TEST_FOR_EXCEPTION(orthoType_!="ICGS"&&orthoType_!="DGKS"&&orthoType_!="IMGS",std::logic_error,
			 "LSQRSolMgr(): Invalid orthogonalization type.");
    }  
  }

  // Create the timer if we need to.
  if (timerSolve_ == Teuchos::null) {
    std::string solveLabel = label_ + ": LSQRSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    timerSolve_ = Teuchos::TimeMonitor::getNewTimer(solveLabel);
#endif
  }

  // Inform the solver manager that the current parameters were set.
  isSet_ = true;
}


template<class ScalarType, class MV, class OP>
Belos::ReturnType LSQRSolMgr<ScalarType,MV,OP>::solve() {

  // Set the current parameters if they were not set before.
  // NOTE:  This may occur if the user generated the solver manager with the default constructor
  // but did not set any parameters using setParameters().
  if (!isSet_) {
    setParameters(Teuchos::parameterList(*getValidParameters()));
  }

  //Teuchos::BLAS<int,ScalarType> blas;
  //Teuchos::LAPACK<int,ScalarType> lapack;
  
  TEST_FOR_EXCEPTION(!problem_->isProblemSet(),LSQRSolMgrLinearProblemFailure,
                     "LSQRSolMgr::solve(): Linear problem is not ready, setProblem() has not been called.");
  TEST_FOR_EXCEPTION(MVT::GetNumberVecs( *(problem_->getRHS()) ) != 1, LSQRSolMgrBlockSizeFailure, "LSQRSolMgr::solve(): Incorrect number of RHS vectors, should be exactly 1.");

  // test isFlexible might go here.
  // isSTSet abbreviates is Status Test set.

  // Next the right-hand sides to solve are identified.  Among other things,
  // this enables getCurrLHSVec() to get the current initial guess vector,
  // and getCurrRHSVec() to get the current right-hand side (in Iter).
  std::vector<int> currRHSIdx(1, 0); problem_->setLSIndex(currRHSIdx);

  // Reset the status test.  
  outputTest_->reset();
  // Assume convergence is achieved, then let any failed convergence set this
  // to false.
  bool isConverged = true;	

  // FIXME: Currently we are setting the initial guess to zero, since the
  // solver doesn't yet know how to handle a nonzero initial guess.  This 
  // could be fixed by rewriting the solver to work with the residual and
  // a delta.  

  // In a least squares problem with a nonzero initial guess, the minimzation
  // problem involves the distance (in a norm depending on the preconditioner)
  // between the solution and the the initial guess.

  //////////////////////////////////////////////////////////////////////////////////////
  // LSQR solver
  Teuchos::ParameterList plist; // Parameter list
  Teuchos::RCP<LSQRIter<ScalarType,MV,OP> > lsqr_iter = Teuchos::rcp( new LSQRIter<ScalarType,MV,OP>(problem_, printer_, outputTest_, plist) );
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

    ////////////////////////////////////////////////////////////////////////////////////
    //
    // check convergence first
    //
    ////////////////////////////////////////////////////////////////////////////////////
    if ( convTest_->getStatus() == Belos::Passed ) {
    }
    else if ( maxIterTest_->getStatus() == Belos::Passed ) {
      isConverged = false;
    }
    ////////////////////////////////////////////////////////////////////////////////////
    //
    // return from iterate(), no status tests Passed
    // Something is wrong, and it is probably our fault.
    //
    ////////////////////////////////////////////////////////////////////////////////////
    else {
      TEST_FOR_EXCEPTION(true,std::logic_error,
			 "LSQRSolMgr::solve(): Invalid return from LSQRIteration::iterate().");
    }
  }
  catch (const std::exception &e) {
    printer_->stream(Belos::Errors) << "Error! Caught std::exception in LSQRIter::iterate() at iteration " 
				    << lsqr_iter->getNumIters() << std::endl 
				    << e.what() << std::endl;
    throw;
  }
      
  // identify current linear system as solved LinearProblem 
  problem_->setCurrLS();
  // print final summary
  sTest_->print( printer_->stream(Belos::FinalSummary) );

  // Print timing information, if the corresponding compile-time and
  // run-time options are enabled.
#ifdef BELOS_TEUCHOS_TIME_MONITOR
  // Calling summarize() can be expensive, so don't call unless the
  // user wants to print out timing details.  summarize() will do all
  // the work even if it's passed a "black hole" output stream.
  if (verbosity_ & TimingDetails)
    Teuchos::TimeMonitor::summarize( printer_->stream(Belos::TimingDetails) );
#endif // BELOS_TEUCHOS_TIME_MONITOR

  // A posteriori solve information
  numIters_ = maxIterTest_->getNumIters();
  matCondNum_ = convTest_->getMatCondNum();
  matNorm_ = convTest_->getMatNorm();
  resNorm_ = convTest_->getResidNorm();
  matResNorm_ = convTest_->getLSResidNorm();

 
  if (!isConverged) {
    return Belos::Unconverged; // return from LSQRSolMgr::solve() 
  }
  return Belos::Converged; // return from LSQRSolMgr::solve() 
}

//  LSQRSolMgr requires the solver manager to return an eponyous std::string.
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
